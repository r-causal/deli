# /// script
# requires-python = ">=3.9"
# dependencies = [
#     "numpy",
#     "scipy",
# ]
# ///
"""Reference fixtures for the Levenberg-Marquardt solver (solver='lm').

Two problems are recorded so the R port can pin parity of the LM solver
against scipy.optimize.root(method='lm'):

  * ee_solver_lm_linear: a well-conditioned linear regression where every
    root-finder should agree to tight tolerance.
  * ee_solver_lm_logistic: a logistic regression, a nonlinear system that
    exercises the least-squares behaviour of Levenberg-Marquardt.

Each fixture stores the data, initial values, and the fitted theta,
variance, asymptotic variance, bread, and confidence intervals produced by
solver='lm'.

Run from the Delicatessen root directory:
    uv run --with-editable . python \
        r-pkg/deli/tests/testthat/fixtures/generate_solver_lm_fixtures.py
"""

import json
import os

import numpy as np

from delicatessen import MEstimator

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


def extract_results(estr, alpha=0.05):
    """Extract the fitted results from a solved MEstimator."""
    results = {
        'theta': to_serializable(estr.theta),
        'variance': to_serializable(estr.variance),
        'asymptotic_variance': to_serializable(estr.asymptotic_variance),
        'bread': to_serializable(estr.bread),
    }
    ci = estr.confidence_intervals(alpha=alpha)
    results['ci_lower'] = to_serializable(ci[:, 0])
    results['ci_upper'] = to_serializable(ci[:, 1])
    return results


def save_fixture(name, data):
    """Write fixture data to disk as JSON."""
    filepath = os.path.join(FIXTURES_DIR, f'{name}.json')
    with open(filepath, 'w') as f:
        json.dump(data, f, indent=2)
    print(f'  Saved: {filepath}')


def generate_linear():
    """Well-conditioned linear regression solved with solver='lm'."""
    print("Generating ee_solver_lm_linear fixture...")

    rng = np.random.default_rng(20240607)
    n = 200
    x1 = rng.normal(0, 1, n)
    x2 = rng.binomial(1, 0.5, n)
    X = np.column_stack([np.ones(n), x1, x2])
    y = 1.0 + 0.5 * x1 - 0.75 * x2 + rng.normal(0, 1, n)

    init = [0.0, 0.0, 0.0]

    def psi(theta):
        residuals = y - X @ np.asarray(theta)
        return (X * residuals[:, None]).T

    estr = MEstimator(psi, init=init)
    estr.estimate(solver='lm')
    print(f'  theta: {estr.theta}')

    save_fixture('ee_solver_lm_linear', {
        'X': X.tolist(),
        'y': y.tolist(),
        'init': init,
        **extract_results(estr),
    })


def generate_logistic():
    """Logistic regression solved with solver='lm'."""
    print("Generating ee_solver_lm_logistic fixture...")

    rng = np.random.default_rng(20240608)
    n = 250
    x1 = rng.normal(0, 1, n)
    x2 = rng.binomial(1, 0.5, n)
    X = np.column_stack([np.ones(n), x1, x2])
    eta = -0.5 + 1.0 * x1 + 0.8 * x2
    prob = 1.0 / (1.0 + np.exp(-eta))
    y = rng.binomial(1, prob, n).astype(float)

    init = [0.0, 0.0, 0.0]

    def psi(theta):
        pred = 1.0 / (1.0 + np.exp(-(X @ np.asarray(theta))))
        residuals = y - pred
        return (X * residuals[:, None]).T

    estr = MEstimator(psi, init=init)
    estr.estimate(solver='lm')
    print(f'  theta: {estr.theta}')

    save_fixture('ee_solver_lm_logistic', {
        'X': X.tolist(),
        'y': y.tolist(),
        'init': init,
        **extract_results(estr),
    })


def main():
    generate_linear()
    generate_logistic()
    print('\nDone.')


if __name__ == '__main__':
    main()
