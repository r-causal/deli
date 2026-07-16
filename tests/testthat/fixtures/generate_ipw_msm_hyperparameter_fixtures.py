# /// script
# requires-python = ">=3.9"
# dependencies = [
#     "numpy",
#     "scipy",
# ]
# ///
"""Generate reference results from Python Delicatessen for ee_ipw_msm with a
tweedie marginal structural model, across meaningful hyperparameter (power)
values.

Marginal structural model outcome families and the theta partition
--------------------------------------------------------------------
ee_ipw_msm partitions theta as ``alpha = theta[:V.shape[1]]`` (the marginal
structural model coefficients) and ``beta = theta[V.shape[1]:]`` (the logistic
propensity score coefficients). The MSM slot is fixed at exactly ``V.shape[1]``
elements. It fits the MSM by handing ``alpha`` to ee_glm as the full theta of a
weighted GLM, weighted by the inverse probability weights.

This partition reserves no slot for an estimated nuisance parameter. The gamma
and negative binomial GLM families read the LAST element of their theta as the
log of a nuisance parameter (gamma shape / negative binomial dispersion), which
would consume one of the MSM coefficient slots and leave ee_glm's design matrix
``V`` (``V.shape[1]`` columns) multiplying only ``V.shape[1] - 1`` coefficients.
Python therefore CANNOT fit a gamma or negative binomial MSM through
ee_ipw_msm: it raises a numpy shape-alignment error regardless of the length of
``init`` (a longer init merely misassigns the extra element to the propensity
model, which then mismatches ``W``). These fixtures deliberately do not attempt
those families; the R tests pin the same non-support.

The tweedie family is different: its power is a FIXED (not estimated)
hyperparameter, so tweedie appends no nuisance estimating equation and needs
exactly ``V.shape[1]`` MSM coefficients. The hyperparameter flows straight
through ee_ipw_msm to ee_glm's variance function v(mu) = mu**p. Tweedie is thus
the outcome family through which the hyperparameter argument is exercised here.

Design of the discrimination
----------------------------
When the MSM is saturated in the exposure (e.g. V = [intercept, binary A]), the
tweedie power cancels out of the estimating equations: the solution reduces to
inverse-probability-weighted group means, which do not depend on p, so theta and
its sandwich variance are identical for every hyperparameter. To make the
hyperparameter genuinely matter, the MSM here includes a continuous
effect-modifier column, so the variance weighting v(mu) = mu**p actually
reweights observations and the fit changes with p. Two fixtures on the SAME data
(p = 1.5 and p = 1.8) let the R tests assert a quantified, non-vacuous change.

Run from the Delicatessen root directory:
    uv run --with-editable . python \
        r-pkg/deli/tests/testthat/fixtures/generate_ipw_msm_hyperparameter_fixtures.py

Outputs JSON fixture files to r-pkg/deli/tests/testthat/fixtures/
"""

import json
import os
import sys

import numpy as np

# Ensure the local delicatessen package is importable
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..', '..', '..', '..'))

from delicatessen import MEstimator
from delicatessen.estimating_equations import ee_ipw_msm

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


def generate_confounded_tweedie_data():
    """Simulate confounded data with a compound Poisson-gamma (tweedie) outcome.

    W1 is a continuous confounder of both the binary exposure A and the outcome.
    Vc is a continuous effect-modifier that enters the marginal structural model
    (so the MSM design is NOT saturated in A, which is what makes the tweedie
    power p matter). The outcome is drawn from the compound Poisson-gamma
    representation of the tweedie family with power p = 1.5, giving a response
    with an exact point mass at zero plus a continuous positive part.
    """
    np.random.seed(20240706)
    n = 600

    W1 = np.random.normal(0, 1, n)                 # continuous confounder
    Vc = np.random.normal(0, 1, n)                 # continuous MSM effect-modifier

    # Confounded binary exposure: Pr(A = 1 | W1) via a logistic model
    pa = 1 / (1 + np.exp(-(-0.3 + 0.9 * W1)))
    A = np.random.binomial(1, pa, n)

    # Conditional mean with a log link; confounding enters through W1
    mu = np.exp(0.5 + 0.4 * A + 0.3 * Vc + 0.35 * W1)

    p = 1.5                                         # power in the compound regime
    phi = 1.0                                       # dispersion
    lam = mu ** (2 - p) / (phi * (2 - p))
    shape = (2 - p) / (p - 1)
    scale = phi * (p - 1) * mu ** (p - 1)

    Y = np.zeros(n)
    for i in range(n):
        k = np.random.poisson(lam[i])              # number of gamma events
        if k > 0:
            Y[i] = np.random.gamma(k * shape, scale[i])

    C = np.ones(n)
    W = np.column_stack([C, W1])                    # propensity score design
    V = np.column_stack([C, A, Vc])                 # marginal structural model design
    return n, Y, A, W, V


def generate_ipw_msm_tweedie_fixture(name, Y, A, W, V, hyperparameter,
                                     with_shape=False):
    """Fit a tweedie MSM via ee_ipw_msm on the given data and save a fixture.

    ``init`` has length V.shape[1] + W.shape[1]: the MSM coefficients followed by
    the propensity score coefficients, with no slot for a nuisance parameter
    (tweedie estimates none). When ``with_shape`` is True, the raw stacked
    estimating-function row count and the row sums at the starting values are
    recorded to pin the exact form of the stacked estimating equation.
    """
    print(f'Generating {name} (tweedie MSM, hyperparameter={hyperparameter})...')

    init = [0.] * (V.shape[1] + W.shape[1])

    def psi(theta):
        return ee_ipw_msm(theta, y=Y, A=A, W=W, V=V,
                          distribution='tweedie', link='log',
                          hyperparameter=hyperparameter)

    payload = {
        'y': to_serializable(Y),
        'a': to_serializable(A),
        'W': W.tolist(),
        'V': V.tolist(),
        'distribution': 'tweedie',
        'link': 'log',
        'hyperparameter': hyperparameter,
        'init': init,
    }

    if with_shape:
        ee_init = np.asarray(psi(init))            # (c+b)-by-n stacked equations
        ee_sum_at_init = np.sum(ee_init, axis=1)   # sum over observations per row
        payload['ee_nrow'] = int(ee_init.shape[0])
        payload['ee_sum_at_init'] = to_serializable(ee_sum_at_init)

    estr = MEstimator(psi, init=init)
    estr.estimate(solver='lm')
    payload.update(extract_results(estr))

    save_fixture(name, payload)


def main():
    print('Generating ee_ipw_msm tweedie hyperparameter fixtures...')

    n, Y, A, W, V = generate_confounded_tweedie_data()

    # Two fixtures on the SAME data. p = 1.5 is the primary case (in the compound
    # Poisson-gamma regime and records the stacked-EE shape); p = 1.8 supports a
    # quantified discrimination that the hyperparameter changes the fit.
    generate_ipw_msm_tweedie_fixture('ee_ipw_msm_tweedie_p15', Y, A, W, V, 1.5,
                                     with_shape=True)
    generate_ipw_msm_tweedie_fixture('ee_ipw_msm_tweedie_p18', Y, A, W, V, 1.8)

    print('Done.')


if __name__ == '__main__':
    main()
