# /// script
# requires-python = ">=3.9"
# dependencies = [
#     "numpy",
#     "scipy",
# ]
# ///
"""Generate reference results for measurement error estimating equations in the
presence of NaN values in the positions Python masks before multiplying by the
sample indicator.

Python replaces the gold-standard measurement with a placeholder (-999) in the
main study (R=1) rows before it is multiplied by the sample indicator, so NaN
values there do not propagate. These fixtures capture that tolerance: the
gold-standard measurement (``y`` / ``a``) carries literal NaN in the main-study
rows, and the naive outcome carries NaN in the external rows for the regression
calibration example.

Run from the Delicatessen root directory:
    uv run r-pkg/deli/tests/testthat/fixtures/generate_measurement_na_fixtures.py

Outputs JSON fixture files to r-pkg/deli/tests/testthat/fixtures/
"""

import json
import os
import sys

import numpy as np

# Ensure the local delicatessen package is importable
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..', '..', '..', '..'))

from delicatessen import MEstimator
from delicatessen.estimating_equations import (
    ee_rogan_gladen,
    ee_rogan_gladen_extended,
    ee_regression_calibration,
    ee_regression,
)

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


def nan_to_null(values):
    """Serialize an input vector so masked NaN positions become JSON null.

    Standard JSON has no NaN literal, so the masked positions are written as
    null. When read back in R (via jsonlite) they become NA, which is the exact
    input the R estimating function receives in the corresponding test.
    """
    return [None if isinstance(v, float) and np.isnan(v) else v
            for v in np.asarray(values).tolist()]


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
    """Save fixture data as JSON (allowing NaN literals)."""
    filepath = os.path.join(FIXTURES_DIR, f'{name}.json')
    with open(filepath, 'w') as f:
        json.dump(data, f, indent=2)
    print(f'  Saved: {filepath}')


def generate_rogan_gladen_na_fixture():
    """ee_rogan_gladen with NaN for the gold-standard Y in the main study."""
    print('Generating ee_rogan_gladen_na fixture...')

    # Cole et al. (2023) example, but the main study gold-standard Y is NaN
    Y_star = np.array([0, 1] + [0, 1, 0, 1])
    Y = np.array([np.nan, np.nan] + [0, 0, 1, 1])  # NaN in the main study (R=1)
    S = np.array([1, 1] + [0, 0, 0, 0])            # R indicator (1=main, 0=external)
    counts = np.array([270, 680] + [71, 18, 38, 203])

    Y_star_exp = np.repeat(Y_star, counts).astype(float)
    Y_exp = np.repeat(Y, counts).astype(float)
    S_exp = np.repeat(S, counts)

    init = [0.5, 0.5, 0.75, 0.75]

    def psi(theta):
        return ee_rogan_gladen(theta, y=Y_exp, y_star=Y_star_exp, r=S_exp)

    estr = MEstimator(psi, init=init)
    estr.estimate(solver='lm')

    save_fixture('ee_rogan_gladen_na', {
        'y': nan_to_null(Y_exp),
        'y_star': Y_star_exp.tolist(),
        'r': S_exp.tolist(),
        'init': init,
        **extract_results(estr),
    })


def generate_rogan_gladen_extended_na_fixture():
    """ee_rogan_gladen_extended with NaN for the gold-standard Y in the main study."""
    print('Generating ee_rogan_gladen_extended_na fixture...')

    Y_star = np.array([0, 1] + [0, 1, 0, 1])
    Y = np.array([np.nan, np.nan] + [0, 0, 1, 1])
    S = np.array([1, 1] + [0, 0, 0, 0])
    counts = np.array([270, 680] + [71, 18, 38, 203])

    Y_star_exp = np.repeat(Y_star, counts).astype(float)
    Y_exp = np.repeat(Y, counts).astype(float)
    S_exp = np.repeat(S, counts)

    X = np.ones((len(Y_exp), 1))

    init = [0.5, 1., 1.]

    def psi(theta):
        return ee_rogan_gladen_extended(theta, y=Y_exp, y_star=Y_star_exp,
                                        r=S_exp, X=X)

    estr = MEstimator(psi, init=init)
    estr.estimate(solver='lm')

    save_fixture('ee_rogan_gladen_extended_na', {
        'y': nan_to_null(Y_exp),
        'y_star': Y_star_exp.tolist(),
        'r': S_exp.tolist(),
        'X': X.tolist(),
        'init': init,
        **extract_results(estr),
    })


def generate_regression_calibration_na_fixture():
    """ee_regression_calibration with NaN for the gold-standard A in the main study."""
    print('Generating ee_regression_calibration_na fixture...')

    A_star = np.array([0, 1, 0, 1] + [1, 0, 1, 0])
    A = np.array([np.nan, np.nan, np.nan, np.nan] + [1, 1, 0, 0])  # NaN in main (R=1)
    Y_vals = np.array([0, 0, 1, 1] + [0, 1, np.nan, np.nan])       # NaN in external (R=0)
    S = np.array([1, 1, 1, 1] + [0, 0, 0, 0])
    counts = np.array([266, 67, 400, 267] + [90, 10, 20, 80])

    A_star_exp = np.repeat(A_star, counts).astype(float)
    A_exp = np.repeat(A, counts).astype(float)
    Y_exp = np.repeat(Y_vals, counts).astype(float)
    S_exp = np.repeat(S, counts)

    X_main = np.column_stack([np.ones(len(A_star_exp)), A_star_exp])

    init = [0., 0.75, 0., 0., 0.]

    def psi(theta):
        theta_calib = theta[:3]
        theta_main = theta[3:]

        # Mask the naive outcome in the external data before it is used
        y_safe = np.where(S_exp == 0, -999, Y_exp).astype(float)

        ee_logit = ee_regression(theta=theta_main, y=y_safe,
                                 X=X_main, model='logistic')
        ee_logit = ee_logit * S_exp

        ee_calib = ee_regression_calibration(theta=theta_calib,
                                             beta=theta_main[1],
                                             a=A_exp,
                                             a_star=A_star_exp,
                                             r=S_exp)

        return np.vstack([ee_calib, ee_logit])

    estr = MEstimator(psi, init=init)
    estr.estimate(solver='lm')

    save_fixture('ee_regression_calibration_na', {
        'a': nan_to_null(A_exp),
        'a_star': A_star_exp.tolist(),
        'y': nan_to_null(Y_exp),
        'r': S_exp.tolist(),
        'init': init,
        **extract_results(estr),
    })


def main():
    generators = [
        generate_rogan_gladen_na_fixture,
        generate_rogan_gladen_extended_na_fixture,
        generate_regression_calibration_na_fixture,
    ]
    for gen in generators:
        gen()
    print('Done.')


if __name__ == '__main__':
    main()
