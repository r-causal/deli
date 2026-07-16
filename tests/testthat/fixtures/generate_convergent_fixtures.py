# /// script
# requires-python = ">=3.9"
# dependencies = [
#     "numpy",
#     "scipy",
# ]
# ///
"""Regenerate convergent fixtures for ee_aft_weibull and ee_tobit,
replacing previous fixtures that had non-convergent theta values.

Run from the Delicatessen root directory:
    uv run r-pkg/deli/tests/testthat/fixtures/generate_convergent_fixtures.py

Outputs JSON fixture files to r-pkg/deli/tests/testthat/fixtures/
"""

import json
import os
import sys

import numpy as np

# Ensure the local delicatessen package is importable
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..', '..', '..', '..'))

from delicatessen import MEstimator
from delicatessen.estimating_equations import (
    ee_aft,
    ee_robust_regression,
    ee_tobit,
)

FIXTURES_DIR = os.path.dirname(__file__)


def to_serializable(obj):
    """Convert numpy types to JSON-serializable Python types."""
    if isinstance(obj, np.ndarray):
        return obj.tolist()
    if isinstance(obj, (np.float64, np.float32)):
        return float(obj)
    if isinstance(obj, (np.int64, np.int32)):
        return int(obj)
    return obj


def extract_results(mestimator, alpha=0.05):
    """Extract standard results from a fitted MEstimator."""
    results = {
        'theta': to_serializable(mestimator.theta),
        'variance': to_serializable(mestimator.variance),
        'asymptotic_variance': to_serializable(mestimator.asymptotic_variance),
        'bread': to_serializable(mestimator.bread),
    }
    ci = mestimator.confidence_intervals(alpha=alpha)
    results['ci_lower'] = to_serializable(ci[:, 0])
    results['ci_upper'] = to_serializable(ci[:, 1])
    return results


def save_fixture(name, data):
    """Save fixture data as JSON."""
    filepath = os.path.join(FIXTURES_DIR, f'{name}.json')
    with open(filepath, 'w') as f:
        json.dump(data, f, indent=2)
    print(f'  Saved: {filepath}')


def generate_aft_weibull():
    """Generate ee_aft_weibull fixture with Weibull survival data.

    Uses the approach from the Delicatessen docs/examples.
    """
    print('Generating ee_aft_weibull fixture...')

    np.random.seed(101)
    n = 200

    # Binary treatment covariate
    X1 = np.random.binomial(1, 0.5, n)

    # Generate Weibull survival times
    # Using numpy weibull parameterization
    shape_param = 1.5
    scale_factor = np.exp(1.0 + 0.5 * X1)
    T_event = scale_factor * np.random.weibull(shape_param, n)

    # Random censoring
    C = np.random.exponential(3, n)
    delta = (T_event <= C).astype(int)
    t = np.minimum(T_event, C)
    t = np.maximum(t, 1e-8)  # floor small values

    censoring_rate = 1 - np.mean(delta)
    print(f'  Censoring rate: {censoring_rate:.2%}')
    print(f'  Time range: [{np.min(t):.3f}, {np.max(t):.3f}]')

    # Design matrix: intercept + treatment
    X = np.column_stack([np.ones(n), X1])

    # theta = [beta_0, beta_1, log(sigma)]
    init = [1.0, 0.5, 0.0]

    def psi(theta):
        return ee_aft(theta, X=X, t=t, delta=delta, distribution='weibull')

    # Try multiple solvers
    for solver in ['lm', 'hybr']:
        try:
            estr = MEstimator(psi, init=init)
            estr.estimate(solver=solver)
            if all(abs(th) < 100 for th in estr.theta):
                print(f'  Converged with solver={solver}')
                print(f'  theta: {estr.theta}')
                save_fixture('ee_aft_weibull', {
                    't': t.tolist(),
                    'delta': delta.astype(int).tolist(),
                    'X': X.tolist(),
                    'distribution': 'weibull',
                    'init': init,
                    **extract_results(estr),
                })
                return True
        except Exception as e:
            print(f'  solver={solver} failed: {e}')

    print('  WARNING: ee_aft_weibull did not converge with any solver')
    return False


def generate_tobit():
    """Generate ee_tobit fixture with left-censored normal data."""
    print('Generating ee_tobit fixture...')

    np.random.seed(42)
    n = 300

    # Single continuous covariate
    X1 = np.random.normal(0, 1, n)

    # True parameters
    true_beta0 = 2.0
    true_beta1 = 1.0
    true_sigma = 1.0

    # Generate latent outcome
    Y_latent = true_beta0 + true_beta1 * X1 + np.random.normal(0, true_sigma, n)

    # Left-censor at 0
    lower = 0.0
    Y_censored = np.maximum(Y_latent, lower)

    censored_pct = np.mean(Y_censored == lower)
    print(f'  Censoring rate: {censored_pct:.2%}')

    # Design matrix with intercept
    X = np.column_stack([np.ones(n), X1])

    # theta = [beta_0, beta_1, log(sigma)]
    init = [true_beta0, true_beta1, np.log(true_sigma)]

    def psi(theta):
        return ee_tobit(theta, X=X, y=Y_censored, lower=lower)

    # Try multiple solvers
    for solver in ['lm', 'hybr']:
        try:
            estr = MEstimator(psi, init=init)
            estr.estimate(solver=solver)
            if all(abs(th) < 100 for th in estr.theta):
                print(f'  Converged with solver={solver}')
                print(f'  theta: {estr.theta}')
                save_fixture('ee_tobit', {
                    'y': Y_censored.tolist(),
                    'X': X.tolist(),
                    'lower': lower,
                    'init': init,
                    **extract_results(estr),
                })
                return True
        except Exception as e:
            print(f'  solver={solver} failed: {e}')

    print('  WARNING: ee_tobit did not converge with any solver')
    return False


def generate_robust_regression_huber():
    """Generate ee_robust_regression fixture with a solvable Huber problem.

    The design places most residuals inside the Huber threshold and adds a
    handful of outliers so the robust fit departs from ordinary least squares.
    Starting from the origin the Jacobian is well conditioned, so the root
    finder converges to non-degenerate estimates.
    """
    print('Generating ee_robust_regression_huber fixture...')

    np.random.seed(20240607)
    n = 100

    # Single continuous covariate
    X1 = np.random.normal(0, 1, n)

    # True parameters
    true_beta0 = 1.0
    true_beta1 = 2.0

    # Outcome with a few positive outliers to exercise the robust loss
    y = true_beta0 + true_beta1 * X1 + np.random.normal(0, 1, n)
    y[:5] += 12.0

    # Design matrix with intercept
    X = np.column_stack([np.ones(n), X1])

    # Huber tuning constant (classic 95% efficiency value)
    k = 1.345
    init = [0.0, 0.0]

    def psi(theta):
        return ee_robust_regression(theta, X=X, y=y, model='linear',
                                    k=k, loss='huber')

    estr = MEstimator(psi, init=init)
    estr.estimate(solver='lm')
    print(f'  theta: {estr.theta}')

    save_fixture('ee_robust_regression_huber', {
        'X': X.tolist(),
        'y': y.tolist(),
        'model': 'linear',
        'k': k,
        'loss': 'huber',
        'init': init,
        **extract_results(estr),
    })
    return True


def main():
    aft_ok = generate_aft_weibull()
    tobit_ok = generate_tobit()
    huber_ok = generate_robust_regression_huber()

    if aft_ok and tobit_ok and huber_ok:
        print('\nDone! All fixtures regenerated successfully.')
    else:
        if not aft_ok:
            print('\nWARNING: ee_aft_weibull fixture still non-convergent')
        if not tobit_ok:
            print('\nWARNING: ee_tobit fixture still non-convergent')
        if not huber_ok:
            print('\nWARNING: ee_robust_regression_huber fixture still non-convergent')


if __name__ == '__main__':
    main()
