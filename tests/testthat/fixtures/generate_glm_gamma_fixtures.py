# /// script
# requires-python = ">=3.9"
# dependencies = [
#     "numpy",
#     "scipy",
# ]
# ///
"""Generate reference results from Python Delicatessen for the ee_glm gamma
distribution, including its digamma-based shape nuisance estimating equation.

The gamma GLM estimates one extra parameter (the log of the shape alpha, whose
reciprocal is the GLM dispersion phi = 1 / alpha), so theta has length
X.shape[1] + 1. Python appends the nuisance estimating equation and treats the
shape honestly in the sandwich variance. The gamma variance function is
v(mu) = mu^2, so the beta score equations coincide with the tweedie distribution
at power p = 2.

The fixtures cover the log link, a non-canonical inverse link, and a weighted
log-link fit. The weighted fixture pins Python's deliberate asymmetry: the gamma
nuisance row is weighted (regression.py:334), unlike the negative binomial one.

Run from the Delicatessen root directory:
    uv run --with-editable . python \
        r-pkg/deli/tests/testthat/fixtures/generate_glm_gamma_fixtures.py

Outputs JSON fixture files to r-pkg/deli/tests/testthat/fixtures/
"""

import json
import os
import sys

import numpy as np

# Ensure the local delicatessen package is importable
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..', '..', '..', '..'))

from delicatessen import MEstimator
from delicatessen.estimating_equations import ee_glm

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
    # Confidence intervals
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


def generate_positive_continuous_data():
    """Simulate strictly positive continuous data (no zeros).

    The gamma nuisance equation involves log(y), so positive support is
    required. The mean follows a log link, mu = exp(0.5 + 0.5*X1 - 0.3*X2), with
    gamma-distributed responses (mean mu, all values strictly > 0).
    """
    np.random.seed(7)
    n = 200

    X1 = np.random.normal(0, 1, n)
    X2 = np.random.normal(0, 1, n)

    mu = np.exp(0.5 + 0.5 * X1 - 0.3 * X2)

    shape = 2.0
    y = np.random.gamma(shape, mu / shape)   # mean mu, all values strictly > 0

    return n, X1, X2, y


def generate_glm_gamma_fixture():
    """Generate fixture for ee_glm with gamma distribution, log link."""
    print('Generating ee_glm (gamma/log) fixture...')

    n, X1, X2, y = generate_positive_continuous_data()

    # Design matrix with intercept
    X = np.column_stack([np.ones(n), X1, X2])

    # theta = [beta0, beta1, beta2, log(alpha)]; length is X.shape[1] + 1
    init = [0., 0., 0., 0.]

    def psi(theta):
        return ee_glm(theta, X=X, y=y, distribution='gamma', link='log')

    # Raw estimating-function row sums at the starting values. The final row is
    # the shape nuisance equation, so this pins the exact form of the
    # estimating equation, not just the solution.
    ee_init = psi(init)                       # (b+1)-by-n estimating functions
    ee_sum_at_init = np.sum(ee_init, axis=1)  # sum over observations per row

    estr = MEstimator(psi, init=init)
    estr.estimate(solver='lm')

    save_fixture('ee_glm_gamma_log', {
        'y': y.tolist(),
        'X': X.tolist(),
        'distribution': 'gamma',
        'link': 'log',
        'init': init,
        'ee_nrow': int(ee_init.shape[0]),
        'ee_sum_at_init': to_serializable(ee_sum_at_init),
        **extract_results(estr),
    })


def generate_glm_gamma_weighted_fixture():
    """Generate fixture for a weighted ee_glm gamma (log link).

    Uses non-uniform integer weights so that the weighting genuinely changes the
    solution. Python weights the shape nuisance row as well as the beta scores,
    so this fixture fails if the R port drops the weights from the nuisance row.
    Integer weights also let the R test cross-check against the equivalent fit on
    the row-expanded (repeated) data.
    """
    print('Generating ee_glm (gamma/log, weighted) fixture...')

    n, X1, X2, y = generate_positive_continuous_data()

    # Deterministic non-uniform integer weights, drawn from a dedicated stream.
    weights = np.random.default_rng(11).integers(1, 5, size=n).astype(float)

    X = np.column_stack([np.ones(n), X1, X2])
    init = [0., 0., 0., 0.]

    def psi(theta):
        return ee_glm(theta, X=X, y=y, distribution='gamma', link='log',
                      weights=weights)

    estr = MEstimator(psi, init=init)
    estr.estimate(solver='lm', maxiter=5000)

    save_fixture('ee_glm_gamma_weighted_log', {
        'y': y.tolist(),
        'X': X.tolist(),
        'weights': weights.tolist(),
        'distribution': 'gamma',
        'link': 'log',
        'init': init,
        **extract_results(estr),
    })


def generate_inverse_link_data():
    """Simulate gamma data with a strictly positive linear predictor.

    The inverse link sets mu = 1 / (X beta), so the linear predictor must stay
    strictly positive for a positive mean. X1 is bounded so eta = 0.5 + 0.3*X1
    stays well away from zero. Responses are gamma-distributed (mean mu > 0).
    """
    np.random.seed(19)
    n = 200

    X1 = np.random.uniform(-1, 1, n)
    eta = 0.5 + 0.3 * X1                      # strictly positive linear predictor
    mu = 1.0 / eta

    shape = 2.0
    y = np.random.gamma(shape, mu / shape)   # mean mu, all values strictly > 0

    return n, X1, y


def generate_glm_gamma_inverse_fixture():
    """Generate fixture for ee_glm gamma with the non-canonical inverse link."""
    print('Generating ee_glm (gamma/inverse) fixture...')

    n, X1, y = generate_inverse_link_data()

    X = np.column_stack([np.ones(n), X1])
    # theta = [beta0, beta1, log(alpha)]; length is X.shape[1] + 1. Start near
    # the data-generating linear predictor so the inverse link stays positive.
    init = [0.5, 0.3, 0.]

    def psi(theta):
        return ee_glm(theta, X=X, y=y, distribution='gamma', link='inverse')

    estr = MEstimator(psi, init=init)
    estr.estimate(solver='lm', maxiter=5000)

    save_fixture('ee_glm_gamma_inverse', {
        'y': y.tolist(),
        'X': X.tolist(),
        'distribution': 'gamma',
        'link': 'inverse',
        'init': init,
        **extract_results(estr),
    })


def main():
    print('Generating gamma GLM fixtures...')
    generate_glm_gamma_fixture()
    generate_glm_gamma_weighted_fixture()
    generate_glm_gamma_inverse_fixture()
    print('Done.')


if __name__ == '__main__':
    main()
