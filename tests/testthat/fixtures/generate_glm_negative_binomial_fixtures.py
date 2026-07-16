# /// script
# requires-python = ">=3.9"
# dependencies = [
#     "numpy",
#     "scipy",
# ]
# ///
"""Generate reference results from Python Delicatessen for the ee_glm negative
binomial distribution, including its polygamma nuisance estimating equation.

The negative binomial GLM estimates one extra parameter (the log of the
dispersion alpha), so theta has length X.shape[1] + 1. Python appends the
polygamma-based nuisance estimating equation and treats the dispersion
honestly in the sandwich variance.

Run from the Delicatessen root directory:
    uv run --with-editable . python \
        r-pkg/deli/tests/testthat/fixtures/generate_glm_negative_binomial_fixtures.py

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


def generate_negative_binomial_data():
    """Generate overdispersed count data for the negative binomial GLM.

    The mean follows a log link, mu = exp(0.5 + 0.5*X1 - 0.3*X2), and the
    counts are drawn with dispersion so that Var(Y) = mu + alpha*mu^2 with
    alpha = 1 / r. numpy parameterizes negative_binomial(r, p) with
    p = r / (r + mu), giving mean mu and dispersion 1 / r.
    """
    np.random.seed(42)
    n = 200

    X1 = np.random.normal(0, 1, n)
    X2 = np.random.normal(0, 1, n)

    mu = np.exp(0.5 + 0.5 * X1 - 0.3 * X2)

    r = 2.0                             # NB shape; dispersion alpha = 1 / r = 0.5
    p = r / (r + mu)                    # numpy success probability per trial
    Y_nb = np.random.negative_binomial(r, p)

    return n, X1, X2, Y_nb


def generate_glm_negative_binomial_fixture():
    """Generate fixture for ee_glm with negative binomial distribution, log link."""
    print('Generating ee_glm (negative_binomial/log) fixture...')

    n, X1, X2, Y_nb = generate_negative_binomial_data()

    # Design matrix with intercept
    X = np.column_stack([np.ones(n), X1, X2])

    # theta = [beta0, beta1, beta2, log(alpha)]; length is X.shape[1] + 1
    init = [0., 0., 0., 0.]

    def psi(theta):
        return ee_glm(theta, X=X, y=Y_nb,
                      distribution='negative_binomial', link='log')

    # Raw estimating-function row sums at the starting values. The final row is
    # the polygamma nuisance equation for the dispersion parameter, so this
    # pins the exact form of the estimating equation, not just the solution.
    ee_init = psi(init)                       # (b+1)-by-n estimating functions
    ee_sum_at_init = np.sum(ee_init, axis=1)  # sum over observations per row

    estr = MEstimator(psi, init=init)
    estr.estimate(solver='lm')

    save_fixture('ee_glm_negative_binomial_log', {
        'y': Y_nb.tolist(),
        'X': X.tolist(),
        'distribution': 'negative_binomial',
        'link': 'log',
        'init': init,
        'ee_nrow': int(ee_init.shape[0]),
        'ee_sum_at_init': to_serializable(ee_sum_at_init),
        **extract_results(estr),
    })


def generate_glm_nb_alias_fixture():
    """Generate fixture for ee_glm using the 'nb' distribution alias, log link.

    This exercises the alias handling. The estimated theta and variance must
    match the full 'negative_binomial' name exactly.
    """
    print("Generating ee_glm (nb alias/log) fixture...")

    n, X1, X2, Y_nb = generate_negative_binomial_data()

    X = np.column_stack([np.ones(n), X1, X2])
    init = [0., 0., 0., 0.]

    def psi(theta):
        return ee_glm(theta, X=X, y=Y_nb,
                      distribution='nb', link='log')

    estr = MEstimator(psi, init=init)
    estr.estimate(solver='lm')

    save_fixture('ee_glm_nb_alias_log', {
        'y': Y_nb.tolist(),
        'X': X.tolist(),
        'distribution': 'nb',
        'link': 'log',
        'init': init,
        **extract_results(estr),
    })


def generate_glm_negative_binomial_weighted_fixture():
    """Generate fixture for a weighted ee_glm negative binomial (log link).

    Uses non-uniform integer weights so that the weighting genuinely changes the
    solution. Python weights the beta score rows but deliberately leaves the
    dispersion nuisance row unweighted (regression.py:341), so this fixture fails
    if the R port ever applies the weights to the nuisance row.
    """
    print('Generating ee_glm (negative_binomial/log, weighted) fixture...')

    n, X1, X2, Y_nb = generate_negative_binomial_data()

    # Deterministic non-uniform integer weights, drawn from a dedicated stream.
    weights = np.random.default_rng(11).integers(1, 5, size=n).astype(float)

    X = np.column_stack([np.ones(n), X1, X2])
    init = [0., 0., 0., 0.]

    def psi(theta):
        return ee_glm(theta, X=X, y=Y_nb,
                      distribution='negative_binomial', link='log',
                      weights=weights)

    estr = MEstimator(psi, init=init)
    estr.estimate(solver='lm', maxiter=5000)

    save_fixture('ee_glm_negative_binomial_weighted_log', {
        'y': Y_nb.tolist(),
        'X': X.tolist(),
        'weights': weights.tolist(),
        'distribution': 'negative_binomial',
        'link': 'log',
        'init': init,
        **extract_results(estr),
    })


def main():
    print('Generating negative binomial GLM fixtures...')
    generate_glm_negative_binomial_fixture()
    generate_glm_nb_alias_fixture()
    generate_glm_negative_binomial_weighted_fixture()
    print('Done.')


if __name__ == '__main__':
    main()
