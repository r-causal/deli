# /// script
# requires-python = ">=3.9"
# dependencies = [
#     "numpy",
#     "scipy",
# ]
# ///
"""Generate fixtures for the top-level compute_sandwich function.

These fixtures call ``delicatessen.compute_sandwich`` directly, rather than
through ``MEstimator``, so that the R port can pin the user-facing contract of
its own ``compute_sandwich`` export. The recorded sandwich is the asymptotic
variance (the bread and meat are each scaled by ``n`` inside compute_sandwich,
so the returned matrix is not divided by ``n``). Dividing the recorded sandwich
by ``n`` gives the variance on the standard-error scale.

For each problem the parameter vector is solved with MEstimator first, and then
the sandwich is computed at that root with each derivative method the R port
supports: exact automatic differentiation and the by-hand central, forward, and
backward finite differences. The SciPy forward-difference default ('approx') is
not recorded because the R port does not replicate that scheme.

Run from the Delicatessen root directory:
    uv run --with-editable . python r-pkg/deli/tests/testthat/fixtures/generate_sandwich_fixtures.py

Outputs JSON fixture files to r-pkg/deli/tests/testthat/fixtures/
"""

import json
import os
import sys

import numpy as np

# Ensure the local delicatessen package is importable
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..', '..', '..', '..'))

from delicatessen import MEstimator, compute_sandwich
from delicatessen.estimating_equations import ee_mean_variance, ee_regression

FIXTURES_DIR = os.path.dirname(__file__)

# Derivative methods shared by the Python and R implementations
DERIV_METHODS = ['exact', 'capprox', 'fapprox', 'bapprox']


def to_serializable(obj):
    """Convert numpy types to JSON-serializable Python types."""
    if isinstance(obj, np.ndarray):
        return obj.tolist()
    if isinstance(obj, (np.float64, np.float32)):
        return float(obj)
    if isinstance(obj, (np.int64, np.int32)):
        return int(obj)
    return obj


def sandwich_by_method(psi, theta):
    """Compute the sandwich at theta for each shared derivative method."""
    return {
        method: to_serializable(compute_sandwich(psi, theta=theta, deriv_method=method))
        for method in DERIV_METHODS
    }


def save_fixture(name, data):
    """Save fixture data as JSON."""
    filepath = os.path.join(FIXTURES_DIR, f'{name}.json')
    with open(filepath, 'w') as f:
        json.dump(data, f, indent=2)
    print(f'  Saved: {filepath}')


def generate_mean_variance():
    """Mean and variance estimating equation with two parameters."""
    print('Generating sandwich_mean_variance fixture...')

    np.random.seed(20240706)
    n = 15
    y = np.round(np.random.normal(4.0, 2.0, n), 4)

    def psi(theta):
        return ee_mean_variance(theta=theta, y=y)

    estr = MEstimator(psi, init=[0.0, 1.0])
    estr.estimate(solver='lm')
    theta = estr.theta
    print(f'  n={n}, theta={theta}')

    save_fixture('sandwich_mean_variance', {
        'y': y.tolist(),
        'init': [0.0, 1.0],
        'theta': to_serializable(theta),
        'n_obs': n,
        'sandwich': sandwich_by_method(psi, theta),
    })


def generate_regression():
    """Linear regression estimating equation with three parameters."""
    print('Generating sandwich_regression fixture...')

    np.random.seed(1234)
    n = 30
    x1 = np.round(np.random.normal(0, 1, n), 4)
    x2 = np.round(np.random.normal(0, 1, n), 4)
    y = np.round(2.0 + 1.5 * x1 - 0.8 * x2 + np.random.normal(0, 1, n), 4)
    X = np.column_stack([np.ones(n), x1, x2])

    def psi(theta):
        return ee_regression(theta, X=X, y=y, model='linear')

    estr = MEstimator(psi, init=[0.0, 0.0, 0.0])
    estr.estimate(solver='lm')
    theta = estr.theta
    print(f'  n={n}, theta={theta}')

    save_fixture('sandwich_regression', {
        'X': X.tolist(),
        'y': y.tolist(),
        'model': 'linear',
        'init': [0.0, 0.0, 0.0],
        'theta': to_serializable(theta),
        'n_obs': n,
        'sandwich': sandwich_by_method(psi, theta),
    })


def main():
    generate_mean_variance()
    generate_regression()
    print('\nDone! compute_sandwich fixtures generated.')


if __name__ == '__main__':
    main()
