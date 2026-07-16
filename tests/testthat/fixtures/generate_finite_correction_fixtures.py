# /// script
# requires-python = ">=3.9"
# dependencies = [
#     "numpy",
#     "scipy",
# ]
# ///
"""Generate fixtures for t-distribution inference under finite_correction.

When ``finite_correction`` is set, delicatessen switches confidence intervals
and P-values from the standard normal distribution to the t-distribution with
``df = n_obs - n_params``. These fixtures capture small-sample scenarios where
the difference between the normal and t distributions is large enough to be
detected well above the 1e-6 comparison tolerance.

Run from the Delicatessen root directory:
    uv run --with-editable . python r-pkg/deli/tests/testthat/fixtures/generate_finite_correction_fixtures.py

Outputs JSON fixture files to r-pkg/deli/tests/testthat/fixtures/
"""

import json
import os
import sys

import numpy as np

# Ensure the local delicatessen package is importable
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..', '..', '..', '..'))

from delicatessen import GMMEstimator, MEstimator
from delicatessen.estimating_equations import ee_regression

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


def extract_inference(estr, alpha=0.05, null=0):
    """Extract point estimates and t-based inference from a fitted estimator."""
    ci = estr.confidence_intervals(alpha=alpha)
    return {
        'theta': to_serializable(estr.theta),
        'variance': to_serializable(estr.variance),
        'n_obs': to_serializable(estr.n_obs),
        'n_params': to_serializable(estr.n_params),
        'df': to_serializable(estr.n_obs - estr.n_params),
        'alpha': alpha,
        'null': null,
        'ci_lower': to_serializable(ci[:, 0]),
        'ci_upper': to_serializable(ci[:, 1]),
        'z_scores': to_serializable(estr.z_scores(null=null)),
        'p_values': to_serializable(estr.p_values(null=null)),
        's_values': to_serializable(estr.s_values(null=null)),
    }


def save_fixture(name, data):
    """Save fixture data as JSON."""
    filepath = os.path.join(FIXTURES_DIR, f'{name}.json')
    with open(filepath, 'w') as f:
        json.dump(data, f, indent=2)
    print(f'  Saved: {filepath}')


def generate_mean():
    """Small-sample mean with HC1 correction (df = n - 1)."""
    print('Generating ee_mean_finite_correction fixture...')

    np.random.seed(20240706)
    n = 12
    y = np.round(np.random.normal(5.0, 2.0, n), 4)

    def psi(theta):
        return (y - theta[0]).reshape(1, -1)

    estr = MEstimator(psi, init=[0.0], finite_correction='HC1')
    estr.estimate(solver='lm')
    print(f'  n_obs={estr.n_obs}, n_params={estr.n_params}, df={estr.n_obs - estr.n_params}')
    print(f'  theta: {estr.theta}')

    # A null close to theta keeps the z-score moderate (near 2.4) so that the
    # t-based and normal-based P-values are separated well above the comparison
    # tolerance. A null of 0 would push both P-values below the tolerance floor
    # and make the P-value test pass against either distribution.
    save_fixture('ee_mean_finite_correction', {
        'y': y.tolist(),
        'init': [0.0],
        'finite_correction': 'HC1',
        **extract_inference(estr, null=4.0),
    })


def generate_regression():
    """Small-sample linear regression with HC1 correction (df = n - p)."""
    print('Generating ee_regression_finite_correction fixture...')

    np.random.seed(1234)
    n = 20
    x1 = np.random.normal(0, 1, n)
    x2 = np.random.normal(0, 1, n)
    y = 1.0 + 0.5 * x1 - 0.75 * x2 + np.random.normal(0, 1, n)
    X = np.column_stack([np.ones(n), x1, x2])

    def psi(theta):
        return ee_regression(theta, X=X, y=y, model='linear')

    estr = MEstimator(psi, init=[0.0, 0.0, 0.0], finite_correction='HC1')
    estr.estimate(solver='lm')
    print(f'  n_obs={estr.n_obs}, n_params={estr.n_params}, df={estr.n_obs - estr.n_params}')
    print(f'  theta: {estr.theta}')

    save_fixture('ee_regression_finite_correction', {
        'X': X.tolist(),
        'y': y.tolist(),
        'model': 'linear',
        'init': [0.0, 0.0, 0.0],
        'finite_correction': 'HC1',
        **extract_inference(estr),
    })


def generate_gmm_overid():
    """Over-identified GMM (instrumental-variables) fit with HC1 correction.

    The moment system stacks three instrument-weighted residual conditions for
    a two-parameter linear model, so there are more moment conditions than
    parameters. Under ``finite_correction`` the GMMEstimator switches confidence
    intervals and P-values to the t-distribution with ``df = n_obs - n_params``,
    where ``n_params`` is the number of parameters (the length of ``init``, here
    two), not the number of moment conditions.

    A small endogenous slope leaves both P-values in the interior of (0, 1), so
    the comparison is above the absolute tolerance floor and separates the t and
    normal reference distributions.
    """
    print('Generating gmm_overid_finite_correction fixture...')

    np.random.seed(20240706)
    n = 20
    z1 = np.round(np.random.normal(0, 1, n), 4)
    z2 = np.round(np.random.normal(0, 1, n), 4)
    u = np.random.normal(0, 1, n)
    x = np.round(0.8 * z1 + 0.6 * z2 + u, 4)
    y = np.round(1.0 + 0.7 * u + np.random.normal(0, 1, n), 4)

    # Design matrix (intercept, endogenous regressor) has two columns; the
    # instrument matrix (intercept, two instruments) has three, so the moment
    # system is over-identified by one condition.
    X = np.column_stack([np.ones(n), x])
    Z = np.column_stack([np.ones(n), z1, z2])

    def psi(theta):
        resid = y - X @ theta
        return (Z * resid[:, None]).T

    estr = GMMEstimator(psi, init=[0.0, 0.0], finite_correction='HC1')
    estr.estimate(solver='bfgs')
    print(f'  n_obs={estr.n_obs}, n_params={estr.n_params}, df={estr.n_obs - estr.n_params}')
    print(f'  theta: {estr.theta}')

    save_fixture('gmm_overid_finite_correction', {
        'x': x.tolist(),
        'z1': z1.tolist(),
        'z2': z2.tolist(),
        'y': y.tolist(),
        'init': [0.0, 0.0],
        'finite_correction': 'HC1',
        **extract_inference(estr),
    })


def main():
    generate_mean()
    generate_regression()
    generate_gmm_overid()
    print('\nDone! Finite-correction inference fixtures generated.')


if __name__ == '__main__':
    main()
