# /// script
# requires-python = ">=3.9"
# dependencies = [
#     "numpy",
#     "scipy",
# ]
# ///
"""Generate fixtures for delta_method in exact (autodiff) mode.

Python's ``delta_method`` defaults to ``deriv_method='exact'`` (forward-mode
automatic differentiation). These fixtures pin the exact-mode Jacobian and the
resulting delta-method covariance for a handful of nontrivial transforms, so the
R port can assert that its own exact mode agrees with Python near machine
precision.

The parameter vector and its covariance come from an actual logistic-regression
M-estimation fit, so the numbers are realistic rather than handcrafted. The
transforms are chosen to exercise the autodiff surface that a zero Jacobian
would silently break: element-wise math (``exp``, ``log``), products of
parameters, a matrix-vector product ``X %*% theta``, and the inverse-logit
transform ``1 / (1 + exp(-eta))``.

Run from the Delicatessen root directory:
    uv run --with-editable . python r-pkg/deli/tests/testthat/fixtures/generate_delta_method_fixtures.py

Outputs JSON fixture files to r-pkg/deli/tests/testthat/fixtures/
"""

import json
import os
import sys

import numpy as np

# Ensure the local delicatessen package is importable
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..', '..', '..', '..'))

from delicatessen import MEstimator, delta_method
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


# A fixed design matrix used by the matrix-product transform. Two output rows so
# the transform returns a length-2 vector and the Jacobian is 2-by-3.
X_TRANSFORM = np.array([
    [1.0, 0.5, -0.25],
    [1.0, -1.0, 2.0],
])


def transform_scalar_mix(theta):
    """Vector-valued transform mixing exp, a product, and a log of a square."""
    return [
        np.exp(theta[0]),
        theta[0] * theta[1],
        np.log(theta[2] ** 2 + 1.0),
    ]


def transform_matrix_invlogit(theta):
    """Inverse-logit of a matrix-vector product X %*% theta (length-2 output)."""
    eta = np.dot(X_TRANSFORM, theta)
    return 1.0 / (1.0 + np.exp(-eta))


def main():
    # Fit a logistic-regression M-estimator to obtain a realistic theta and its
    # covariance matrix (intercept plus two coefficients).
    rng = np.random.default_rng(20260706)
    n = 400
    x1 = rng.normal(size=n)
    x2 = rng.normal(size=n)
    design = np.column_stack([np.ones(n), x1, x2])
    linpred = 0.4 - 0.8 * x1 + 0.6 * x2
    prob = 1.0 / (1.0 + np.exp(-linpred))
    y = rng.binomial(1, prob)

    def psi(theta):
        return ee_regression(theta, X=design, y=y, model='logistic')

    estr = MEstimator(psi, init=[0.0, 0.0, 0.0])
    estr.estimate(solver='lm')

    theta = np.asarray(estr.theta)
    covariance = np.asarray(estr.variance)

    fixtures = {}

    # Transform 1: element-wise mixed transform.
    g1 = np.asarray(transform_scalar_mix(theta))
    var1_exact = delta_method(theta, transform_scalar_mix, covariance,
                              deriv_method='exact')
    var1_capprox = delta_method(theta, transform_scalar_mix, covariance,
                                deriv_method='capprox')
    fixtures['scalar_mix'] = {
        'estimate': to_serializable(g1),
        'variance': to_serializable(var1_exact),
        'variance_capprox': to_serializable(var1_capprox),
    }

    # Transform 2: inverse-logit of a matrix-vector product.
    g2 = np.asarray(transform_matrix_invlogit(theta))
    var2_exact = delta_method(theta, transform_matrix_invlogit, covariance,
                              deriv_method='exact')
    var2_capprox = delta_method(theta, transform_matrix_invlogit, covariance,
                                deriv_method='capprox')
    fixtures['matrix_invlogit'] = {
        'X': to_serializable(X_TRANSFORM),
        'estimate': to_serializable(g2),
        'variance': to_serializable(var2_exact),
        'variance_capprox': to_serializable(var2_capprox),
    }

    out = {
        'theta': to_serializable(theta),
        'covariance': to_serializable(covariance),
        'transforms': fixtures,
    }

    path = os.path.join(FIXTURES_DIR, 'delta_method_exact.json')
    with open(path, 'w') as f:
        json.dump(out, f, indent=2)
    print('Wrote', path)


if __name__ == '__main__':
    main()
