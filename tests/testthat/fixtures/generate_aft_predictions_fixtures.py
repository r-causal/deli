# /// script
# requires-python = ">=3.9"
# dependencies = [
#     "numpy",
#     "scipy",
# ]
# ///
"""Generate fixtures for aft_predictions_function.

Fits accelerated failure time models with ee_aft on the deterministic
Collett (2015) breast cancer data (no RNG), then records function-level
predictions with delta-method point-wise confidence intervals across
distributions, measures, covariate patterns, and alpha levels.

Python's aft_predictions_function computes the delta-method variance with
exact forward-mode automatic differentiation (delta_method defaults to
deriv_method='exact'). The recorded theta and covariance are also stored so
the R test can feed identical inputs to aft_predictions_function, isolating
the prediction and delta-method logic from AFT fitting parity.

Run from the Delicatessen root directory:
    uv run --with-editable . python \
        r-pkg/deli/tests/testthat/fixtures/generate_aft_predictions_fixtures.py

Outputs aft_predictions_function.json to the fixtures directory.
"""

import json
import os
import sys

import numpy as np

# Ensure the local delicatessen package is importable
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..', '..', '..', '..'))

from delicatessen import MEstimator
from delicatessen.estimating_equations import ee_aft
from delicatessen.utilities import aft_predictions_function
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


def save_fixture(name, data):
    """Save fixture data as JSON."""
    filepath = os.path.join(FIXTURES_DIR, f'{name}.json')
    with open(filepath, 'w') as f:
        json.dump(data, f, indent=2)
    print(f'  Saved: {filepath}')


def fit_aft(t, delta, X, distribution, init):
    """Fit an AFT model, returning theta and covariance."""
    def psi(theta):
        return ee_aft(theta, X=X, t=t, delta=delta, distribution=distribution)

    for solver in ['lm', 'hybr']:
        estr = MEstimator(psi, init=init)
        estr.estimate(solver=solver)
        if all(abs(th) < 100 for th in estr.theta):
            print(f'  {distribution}: converged with solver={solver}')
            return estr.theta, estr.variance
    raise RuntimeError(f'{distribution} did not converge')


def main():
    dat = load_breast_cancer()
    delta = dat[:, 0]
    t = dat[:, 1]
    stain = dat[:, 2]
    X = np.column_stack([np.ones(dat.shape[0]), stain])

    # Fit models for each distribution.
    # Exponential has no scale parameter, so init has 2 elements.
    model_specs = {
        'weibull': [5.0, 0.0, 0.0],
        'log-normal': [5.0, 0.0, 0.0],
        'log-logistic': [5.0, 0.0, 0.0],
        'exponential': [5.0, 0.0],
    }

    models = {}
    for dist, init in model_specs.items():
        theta, covariance = fit_aft(t, delta, X, dist, init)
        models[dist] = {
            'theta': to_serializable(theta),
            'covariance': to_serializable(covariance),
        }

    times = [50.0, 100.0, 150.0, 200.0]
    measures = ['survival', 'risk', 'hazard', 'cumulative_hazard', 'density']
    patterns = {
        'neg': [[1.0, 0.0]],   # negative stain
        'pos': [[1.0, 1.0]],   # positive stain
    }

    scenarios = []
    for dist in model_specs:
        theta = np.asarray(models[dist]['theta'])
        covariance = np.asarray(models[dist]['covariance'])
        for measure in measures:
            for pat_name, Xpat in patterns.items():
                for alpha in [0.05, 0.10]:
                    preds = aft_predictions_function(
                        X=Xpat, times=times, theta=theta,
                        covariance=covariance, distribution=dist,
                        measure=measure, alpha=alpha,
                    )
                    scenarios.append({
                        'distribution': dist,
                        'measure': measure,
                        'pattern': pat_name,
                        'X': Xpat,
                        'times': times,
                        'alpha': alpha,
                        'predicted': to_serializable(preds[:, 0]),
                        'variance': to_serializable(preds[:, 1]),
                        'lower': to_serializable(preds[:, 2]),
                        'upper': to_serializable(preds[:, 3]),
                    })

    save_fixture('aft_predictions_function', {
        'models': models,
        'times': times,
        'scenarios': scenarios,
    })
    print(f'  {len(scenarios)} scenarios recorded')


if __name__ == '__main__':
    main()
