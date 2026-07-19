# /// script
# requires-python = ">=3.9"
# dependencies = [
#     "numpy",
#     "scipy",
# ]
# ///
"""Generate a weighted reference result for ee_rogan_gladen_extended.

The unweighted extended fixture is produced by
generate_all_remaining_fixtures.py. This companion adds observation weights so
the R port can verify that the weights reach the sensitivity and specificity
nuisance models, not only the corrected-mean row.

Run against the local Delicatessen source:
    uv run generate_rogan_gladen_extended_weighted_fixture.py
"""

import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..', '..', '..', '..'))

from delicatessen import MEstimator
from delicatessen.estimating_equations import ee_rogan_gladen_extended

FIXTURES_DIR = os.path.dirname(os.path.abspath(__file__))


def to_serializable(obj):
    if isinstance(obj, np.ndarray):
        return obj.tolist()
    if isinstance(obj, (np.float64, np.float32)):
        return float(obj)
    if isinstance(obj, (np.int64, np.int32)):
        return int(obj)
    return obj


def extract_results(mestimator, alpha=0.05):
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
    filepath = os.path.join(FIXTURES_DIR, f'{name}.json')
    with open(filepath, 'w') as f:
        json.dump(data, f, indent=2)
    print(f'  Saved: {filepath}')


def main():
    # Same Cole et al. example as the unweighted extended fixture, with an
    # intercept-only design matrix.
    Y_star = np.array([0, 1] + [0, 1, 0, 1])
    Y = np.array([0, 0] + [0, 0, 1, 1])
    S = np.array([1, 1] + [0, 0, 0, 0])
    counts = np.array([270, 680] + [71, 18, 38, 203])

    Y_star_exp = np.repeat(Y_star, counts)
    Y_exp = np.repeat(Y, counts)
    S_exp = np.repeat(S, counts)

    X = np.ones((len(Y_exp), 1))

    # Non-constant weights so the weighted nuisance fits differ from the
    # unweighted fits. Mismeasured positives carry twice the weight.
    weights = np.where(Y_star_exp == 1, 2.0, 1.0)

    init = [0.5, 1., 1.]

    def psi(theta):
        return ee_rogan_gladen_extended(theta, y=Y_exp, y_star=Y_star_exp,
                                        r=S_exp, X=X, weights=weights)

    estr = MEstimator(psi, init=init)
    estr.estimate(solver='lm')

    save_fixture('ee_rogan_gladen_extended_weighted', {
        'y': Y_exp.tolist(),
        'y_star': Y_star_exp.tolist(),
        'r': S_exp.tolist(),
        'X': X.tolist(),
        'weights': weights.tolist(),
        'init': init,
        **extract_results(estr),
    })


if __name__ == '__main__':
    main()
