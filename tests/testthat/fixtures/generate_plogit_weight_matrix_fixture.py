# /// script
# requires-python = ">=3.9"
# dependencies = [
#     "numpy",
#     "scipy",
# ]
# ///
"""Generate a reference result for ee_plogit with a time-varying (2D) weight
matrix.

Python's ee_plogit accepts an n-by-K weight matrix in addition to the usual
1D weight vector. Each column weights one time interval of the discrete-time
design matrix S. Internally, generate_weights returns the matrix unchanged and
ee_plogit validates the shape (rows == number of observations, columns ==
number of time intervals), transposes it to K-by-n, and multiplies it into the
K-by-n residual matrix elementwise before summing across intervals. See
delicatessen/estimating_equations/survival.py (the `if weights.ndim == 2`
block) and processing.py generate_weights.

This fixture fits the breast cancer data used by the other plogit fixtures with
a Gompertz-like linear time specification (S = intercept + linear time, so
K = max(t) = 225 intervals). The weight matrix is drawn from a fixed RNG stream
so it differs across both rows and columns: a fit that mishandles the axis
mapping (for example, transposing the matrix or collapsing it to a single axis)
cannot reproduce the recorded theta. The unweighted (unit-weight) theta is
recorded on the same data so the R test can assert the weighted fit is
genuinely different (discrimination guard). The smallest per-parameter gap
between the weighted and unweighted theta is about 0.00296 (the statin
coefficient); the intercept differs by about 0.219.

Run from the Delicatessen root directory:
    uv run --with-editable . python \
        r-pkg/deli/tests/testthat/fixtures/generate_plogit_weight_matrix_fixture.py

Outputs ee_plogit_weight_matrix.json to r-pkg/deli/tests/testthat/fixtures/
"""

import json
import os
import sys

import numpy as np

# Ensure the local delicatessen package is importable
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..', '..', '..', '..'))

from delicatessen import MEstimator
from delicatessen.estimating_equations import ee_plogit
from delicatessen.data import load_breast_cancer

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


def generate_plogit_weight_matrix_fixture():
    """Generate fixture for ee_plogit with a 2D (time-varying) weight matrix."""
    print('Generating ee_plogit weight-matrix fixture...')

    # Breast cancer data, the same source used by the other plogit fixtures
    d = load_breast_cancer()
    delta = d[:, 0].astype(int)
    t = d[:, 1].astype(int)
    stain = d[:, 2]
    X = stain[:, None]                       # n-by-1 covariate design matrix

    # Gompertz-like linear time specification: intercept + linear term. With S
    # provided, the number of intervals K equals max(t).
    t_steps = np.arange(1, np.max(t) + 1)
    s_matrix = np.concatenate(
        [np.ones(t_steps.shape)[:, None], t_steps[:, None]], axis=1
    )
    n = d.shape[0]                           # number of observations
    k = t_steps.shape[0]                     # number of time intervals (== max(t))

    # Time-varying weight matrix, n-by-K, differing across both rows and columns.
    # Drawn from a dedicated deterministic RNG stream and bounded away from zero.
    weights = np.random.default_rng(2024).uniform(0.5, 2.0, size=(n, k))

    init = [0., -3., 0.]                     # statin coeff, time intercept, time slope

    def psi_weighted(theta):
        return ee_plogit(theta, X=X, t=t, delta=delta, S=s_matrix, weights=weights)

    def psi_unweighted(theta):
        return ee_plogit(theta, X=X, t=t, delta=delta, S=s_matrix)

    estr = MEstimator(psi_weighted, init=init)
    estr.estimate(solver='lm')

    estr_unweighted = MEstimator(psi_unweighted, init=init)
    estr_unweighted.estimate(solver='lm')

    save_fixture('ee_plogit_weight_matrix', {
        'delta': delta.tolist(),
        't': t.tolist(),
        'X': X.tolist(),
        'S': s_matrix.tolist(),
        'weights': weights.tolist(),        # n-by-K nested list (loads as an n-by-K matrix in R)
        'n': int(n),
        'k': int(k),
        'init': init,
        'theta_unweighted': to_serializable(estr_unweighted.theta),
        **extract_results(estr),
    })


def main():
    generate_plogit_weight_matrix_fixture()
    print('Done.')


if __name__ == '__main__':
    main()
