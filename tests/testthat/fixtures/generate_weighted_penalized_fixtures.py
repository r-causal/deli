# /// script
# requires-python = ">=3.9"
# dependencies = [
#     "numpy",
#     "scipy",
# ]
# ///
"""Generate reference results from Python Delicatessen for weighted penalized
regression, pinning the placement of the penalty term inside the weight.

Python weights the entire penalized estimating equation, including the penalty
term. Every penalized regression EE returns

    w * (((y - pred_y) * X).T - penalty_terms[:, None])

so the penalty is multiplied by the per-observation weight, not subtracted after
weighting the score. With unit weights the placement is invisible, but with
non-unit weights it changes both the estimating-function value at any theta and
the solved coefficients, variance, and bread.

These fixtures use non-unit integer weights that vary substantially (1 to 5), so
sum(weights) differs markedly from n. Each fixture records:

  * the solved theta, variance, asymptotic_variance, bread, meat, and CIs, and
  * ee_sum_at_init: the row sums of the raw estimating functions evaluated at a
    deliberately non-zero init.

The init is non-zero so that the penalty term is non-zero there. This makes
ee_sum_at_init a direct, solve-free pin on the weighted-penalty form: the
weighted-penalty and unweighted-penalty implementations differ at init by
penalty_term_b * (sum(weights) - n) in each row.

Python's approximate LASSO and elastic-net emit a non-differentiability
UserWarning; that is expected and the warnings are silenced while generating.

Run from the Delicatessen root directory:
    uv run r-pkg/deli/tests/testthat/fixtures/generate_weighted_penalized_fixtures.py

Outputs JSON fixture files to r-pkg/deli/tests/testthat/fixtures/
"""

import json
import os
import sys
import warnings

import numpy as np

# Ensure the local delicatessen package is importable
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..', '..', '..', '..'))

from delicatessen import MEstimator
from delicatessen.estimating_equations import (
    ee_ridge_regression,
    ee_bridge_regression,
    ee_lasso_regression,
    ee_dlasso_regression,
    ee_elasticnet_regression,
    ee_additive_regression,
)
from delicatessen.utilities import additive_design_matrix

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
        'meat': to_serializable(mestimator.meat),
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


def generate_weighted_linear_data():
    """Simulate a linear model with deterministic non-unit integer weights.

    The weights are drawn from 1 to 5 (inclusive), so they vary substantially
    and sum(weights) is well away from n. That gap is exactly what discriminates
    the weighted-penalty form from an unweighted-penalty one.
    """
    rng = np.random.default_rng(20240706)
    n = 100

    x1 = rng.normal(0, 1, n)
    x2 = rng.normal(0, 1, n)
    y = 1.0 + 1.0 * x1 + 1.0 * x2 + rng.normal(0, 1, n)

    weights = rng.integers(1, 6, size=n).astype(float)   # integer weights in {1,...,5}
    X = np.column_stack([np.ones(n), x1, x2])
    return n, X, y, weights


# Non-zero init so the penalty term is non-zero at init (the intercept stays
# unpenalized via penalty[0] = 0). This makes ee_sum_at_init discriminate.
INIT = [0.5, 0.5, -0.3]
# Penalise the two slopes but not the intercept, with a strong penalty so the
# weighted-versus-unweighted divergence is large relative to any tolerance.
PENALTY = [0., 10., 10.]


def ee_row_sums_at_init(psi):
    """Row sums of the raw estimating functions evaluated at INIT."""
    ee_init = np.asarray(psi(INIT))          # b-by-n estimating functions
    return ee_init.shape[0], np.sum(ee_init, axis=1)


def generate_ridge_fixture():
    print('Generating ee_ridge_regression (weighted) fixture...')
    n, X, y, weights = generate_weighted_linear_data()

    def psi(theta):
        return ee_ridge_regression(theta, X=X, y=y, model='linear',
                                   penalty=PENALTY, weights=weights)

    ee_nrow, ee_sum_at_init = ee_row_sums_at_init(psi)
    estr = MEstimator(psi, init=INIT)
    estr.estimate(solver='lm', maxiter=20000)

    save_fixture('ee_ridge_regression_weighted', {
        'y': y.tolist(),
        'X': X.tolist(),
        'weights': weights.tolist(),
        'model': 'linear',
        'penalty': PENALTY,
        'init': INIT,
        'ee_nrow': int(ee_nrow),
        'ee_sum_at_init': to_serializable(ee_sum_at_init),
        **extract_results(estr),
    })


def generate_bridge_fixture():
    print('Generating ee_bridge_regression (weighted) fixture...')
    n, X, y, weights = generate_weighted_linear_data()
    gamma = 2.3

    def psi(theta):
        return ee_bridge_regression(theta, X=X, y=y, model='linear',
                                    penalty=PENALTY, gamma=gamma,
                                    weights=weights)

    ee_nrow, ee_sum_at_init = ee_row_sums_at_init(psi)
    estr = MEstimator(psi, init=INIT)
    estr.estimate(solver='lm', maxiter=20000)

    save_fixture('ee_bridge_regression_weighted', {
        'y': y.tolist(),
        'X': X.tolist(),
        'weights': weights.tolist(),
        'model': 'linear',
        'penalty': PENALTY,
        'gamma': gamma,
        'init': INIT,
        'ee_nrow': int(ee_nrow),
        'ee_sum_at_init': to_serializable(ee_sum_at_init),
        **extract_results(estr),
    })


def generate_lasso_fixture():
    print('Generating ee_lasso_regression (weighted) fixture...')
    n, X, y, weights = generate_weighted_linear_data()

    def psi(theta):
        return ee_lasso_regression(theta, X=X, y=y, model='linear',
                                   penalty=PENALTY, weights=weights)

    # The approximate LASSO is not everywhere differentiable; silence the
    # expected UserWarning while evaluating and solving.
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        ee_nrow, ee_sum_at_init = ee_row_sums_at_init(psi)
        estr = MEstimator(psi, init=INIT)
        estr.estimate(solver='lm', maxiter=20000)

    save_fixture('ee_lasso_regression_weighted', {
        'y': y.tolist(),
        'X': X.tolist(),
        'weights': weights.tolist(),
        'model': 'linear',
        'penalty': PENALTY,
        'init': INIT,
        'ee_nrow': int(ee_nrow),
        'ee_sum_at_init': to_serializable(ee_sum_at_init),
        **extract_results(estr),
    })


def generate_dlasso_fixture():
    print('Generating ee_dlasso_regression (weighted) fixture...')
    n, X, y, weights = generate_weighted_linear_data()

    def psi(theta):
        return ee_dlasso_regression(theta, X=X, y=y, model='linear',
                                    penalty=PENALTY, weights=weights)

    ee_nrow, ee_sum_at_init = ee_row_sums_at_init(psi)
    estr = MEstimator(psi, init=INIT)
    estr.estimate(solver='lm', maxiter=20000)

    save_fixture('ee_dlasso_regression_weighted', {
        'y': y.tolist(),
        'X': X.tolist(),
        'weights': weights.tolist(),
        'model': 'linear',
        'penalty': PENALTY,
        'init': INIT,
        'ee_nrow': int(ee_nrow),
        'ee_sum_at_init': to_serializable(ee_sum_at_init),
        **extract_results(estr),
    })


def generate_elasticnet_fixture():
    print('Generating ee_elasticnet_regression (weighted) fixture...')
    n, X, y, weights = generate_weighted_linear_data()
    ratio = 0.5

    def psi(theta):
        return ee_elasticnet_regression(theta, X=X, y=y, model='linear',
                                        penalty=PENALTY, ratio=ratio,
                                        weights=weights)

    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        ee_nrow, ee_sum_at_init = ee_row_sums_at_init(psi)
        estr = MEstimator(psi, init=INIT)
        estr.estimate(solver='lm', maxiter=20000)

    save_fixture('ee_elasticnet_regression_weighted', {
        'y': y.tolist(),
        'X': X.tolist(),
        'weights': weights.tolist(),
        'model': 'linear',
        'penalty': PENALTY,
        'ratio': ratio,
        'init': INIT,
        'ee_nrow': int(ee_nrow),
        'ee_sum_at_init': to_serializable(ee_sum_at_init),
        **extract_results(estr),
    })


def generate_additive_fixture():
    """Weighted additive regression, which delegates to ee_bridge_regression.

    This pins the weighted-penalty form through the additive delegation path.
    The spline specification matches the pattern in the unweighted additive
    fixture so R can reconstruct the same design.
    """
    print('Generating ee_additive_regression (weighted) fixture...')

    rng = np.random.default_rng(20240706)
    n = 100
    x1 = rng.uniform(-3, 3, n)
    y = np.sin(x1) + 0.5 * x1 + rng.normal(0, 0.5, n)
    weights = rng.integers(1, 6, size=n).astype(float)
    X = np.column_stack([np.ones(n), x1])

    specs = [None, {"knots": [-2, -1, 0, 1, 2], "penalty": 5,
                    "natural": True, "power": 3}]

    Xa = additive_design_matrix(X=X, specifications=specs)
    n_params = Xa.shape[1]

    # Non-zero init so the spline penalty terms are non-zero at init.
    init = [0.3] * n_params

    def psi(theta):
        return ee_additive_regression(theta, X=X, y=y, model='linear',
                                      specifications=specs, weights=weights)

    ee_init = np.asarray(psi(init))
    ee_nrow = ee_init.shape[0]
    ee_sum_at_init = np.sum(ee_init, axis=1)

    estr = MEstimator(psi, init=init)
    estr.estimate(solver='lm', maxiter=20000)

    save_fixture('ee_additive_regression_weighted', {
        'y': y.tolist(),
        'X': X.tolist(),
        'weights': weights.tolist(),
        'model': 'linear',
        'specifications': specs,
        'n_params': int(n_params),
        'init': init,
        'ee_nrow': int(ee_nrow),
        'ee_sum_at_init': to_serializable(ee_sum_at_init),
        **extract_results(estr),
    })


def main():
    print('Generating weighted penalized regression fixtures...')
    generate_ridge_fixture()
    generate_bridge_fixture()
    generate_lasso_fixture()
    generate_dlasso_fixture()
    generate_elasticnet_fixture()
    generate_additive_fixture()
    print('Done.')


if __name__ == '__main__':
    main()
