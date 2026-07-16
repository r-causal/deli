# /// script
# requires-python = ">=3.9"
# dependencies = [
#     "numpy",
#     "scipy",
# ]
# ///
"""Generate reference results from Python Delicatessen for the ee_glm tweedie
distribution across meaningful hyperparameter (power) values.

The tweedie distribution uses the ``hyperparameter`` argument as the power p in
the variance function v(mu) = mu**p. Unlike the negative binomial and gamma
distributions, tweedie does NOT estimate an additional parameter, so theta has
length X.shape[1] and no nuisance estimating equation is appended. The power is
a fixed (not estimated) hyperparameter bounded to be non-negative.

Notable equivalences encoded by these fixtures:
  * p = 1 gives v(mu) = mu, identical to the Poisson distribution.
  * p = 2 gives v(mu) = mu**2, matching the gamma distribution variance (the
    beta score equations coincide; gamma additionally estimates a dispersion
    nuisance, tweedie does not).
  * 1 < p < 2 is the compound Poisson-gamma regime (a point mass at zero plus a
    continuous positive part).

Run from the Delicatessen root directory:
    uv run --with-editable . python \
        r-pkg/deli/tests/testthat/fixtures/generate_glm_tweedie_fixtures.py

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


def generate_compound_poisson_gamma_data():
    """Simulate proper tweedie (compound Poisson-gamma) data with power p=1.5.

    Uses the standard compound Poisson-gamma parameterization of the tweedie
    family. Each observation is a sum of a Poisson number of gamma draws, so the
    response has an exact point mass at zero plus a continuous positive part.
    The mean follows a log link, mu = exp(0.5 + 0.5*X1 - 0.3*X2).
    """
    np.random.seed(42)
    n = 200

    X1 = np.random.normal(0, 1, n)
    X2 = np.random.normal(0, 1, n)

    mu = np.exp(0.5 + 0.5 * X1 - 0.3 * X2)

    p = 1.5                             # power in the compound regime (1, 2)
    phi = 1.0                           # dispersion
    # Poisson rate, gamma shape, and gamma scale for the compound representation
    lam = mu ** (2 - p) / (phi * (2 - p))
    shape = (2 - p) / (p - 1)
    scale = phi * (p - 1) * mu ** (p - 1)

    y = np.zeros(n)
    for i in range(n):
        n_i = np.random.poisson(lam[i])          # number of gamma events
        if n_i > 0:
            y[i] = np.random.gamma(n_i * shape, scale[i])  # sum of gamma draws

    return n, X1, X2, y


def generate_positive_continuous_data():
    """Simulate strictly positive continuous data (no zeros).

    Positive support is required for the p=2 tweedie/gamma comparison, since the
    gamma nuisance equation involves log(y). The mean follows a log link,
    mu = exp(0.5 + 0.5*X1 - 0.3*X2), with gamma-distributed responses.
    """
    np.random.seed(7)
    n = 200

    X1 = np.random.normal(0, 1, n)
    X2 = np.random.normal(0, 1, n)

    mu = np.exp(0.5 + 0.5 * X1 - 0.3 * X2)

    shape = 2.0
    y = np.random.gamma(shape, mu / shape)   # mean mu, all values strictly > 0

    return n, X1, X2, y


def generate_tweedie_fixture(name, X, y, hyperparameter, with_shape=False):
    """Fit tweedie ee_glm on the given data and save a fixture.

    When ``with_shape`` is True, the raw estimating-function row count and the
    row sums at the starting values are recorded to pin the exact form of the
    estimating equation (tweedie must return exactly X.shape[1] rows, with no
    nuisance row appended).
    """
    print(f'Generating {name} (tweedie, hyperparameter={hyperparameter})...')

    init = [0., 0., 0.]                  # length X.shape[1]; no extra parameter

    def psi(theta):
        return ee_glm(theta, X=X, y=y,
                      distribution='tweedie', link='log',
                      hyperparameter=hyperparameter)

    payload = {
        'y': to_serializable(y),
        'X': X.tolist(),
        'distribution': 'tweedie',
        'link': 'log',
        'hyperparameter': hyperparameter,
        'init': init,
    }

    if with_shape:
        ee_init = psi(init)                       # b-by-n estimating functions
        ee_sum_at_init = np.sum(ee_init, axis=1)  # sum over observations per row
        payload['ee_nrow'] = int(ee_init.shape[0])
        payload['ee_sum_at_init'] = to_serializable(ee_sum_at_init)

    estr = MEstimator(psi, init=init)
    estr.estimate(solver='lm')
    payload.update(extract_results(estr))

    save_fixture(name, payload)


def main():
    print('Generating tweedie GLM fixtures...')

    # Compound Poisson-gamma data (has zeros) for the main p=1.5 case and for
    # the p=1 / Poisson equivalence. The shape assertions pin the EE form here.
    n, X1, X2, y_cpg = generate_compound_poisson_gamma_data()
    X_cpg = np.column_stack([np.ones(n), X1, X2])
    generate_tweedie_fixture('ee_glm_tweedie_p15_log', X_cpg, y_cpg, 1.5,
                             with_shape=True)
    generate_tweedie_fixture('ee_glm_tweedie_p1_log', X_cpg, y_cpg, 1.0)

    # Strictly positive data (no zeros) for the p=2 / gamma equivalence.
    n, X1, X2, y_pos = generate_positive_continuous_data()
    X_pos = np.column_stack([np.ones(n), X1, X2])
    generate_tweedie_fixture('ee_glm_tweedie_p2_log', X_pos, y_pos, 2.0)

    print('Done.')


if __name__ == '__main__':
    main()
