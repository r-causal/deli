# /// script
# requires-python = ">=3.9"
# dependencies = [
#     "numpy",
#     "scipy",
# ]
# ///
"""Generate fixtures for survival_predictions (bd-1v9o.9 / bd-1v9o.10).

Fits parametric survival models with ee_survival_model on the deterministic
Collett (2015) breast cancer data (no RNG), then records the function-level
predictions with delta-method point-wise confidence intervals across
distributions, measures, times, and alpha levels.

Python's survival_predictions computes the delta-method variance with exact
forward-mode automatic differentiation: it calls delta_method, which defaults
to deriv_method='exact'. The recorded theta and covariance are stored so the R
test can feed identical inputs to survival_predictions, isolating the
prediction and delta-method logic from survival-model fitting parity. Because
Python uses exact autodiff internally, the recorded variances and confidence
limits pin the exact-mode values.

Run from the Delicatessen root directory:
    uv run --with-editable . python \
        r-pkg/deli/tests/testthat/fixtures/generate_survival_predictions_fixture.py

Outputs survival_predictions.json to the fixtures directory.
"""

import json
import os
import sys

import numpy as np

# Ensure the local delicatessen package is importable
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..', '..', '..', '..'))

from delicatessen import MEstimator
from delicatessen.estimating_equations import ee_survival_model
from delicatessen.utilities import survival_predictions
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


def fit_survival(t, delta, distribution, init):
    """Fit a parametric survival model, returning theta and covariance."""
    def psi(theta):
        return ee_survival_model(theta, t=t, delta=delta,
                                 distribution=distribution)

    for solver in ['lm', 'hybr']:
        estr = MEstimator(psi, init=init)
        estr.estimate(solver=solver)
        print(f'  {distribution}: solver={solver}, theta={estr.theta}')
        return estr.theta, estr.variance
    raise RuntimeError(f'{distribution} did not converge')


def main():
    dat = load_breast_cancer()
    delta = dat[:, 0]
    t = dat[:, 1]

    # Fit models for each distribution. The Weibull carries a shape parameter
    # (theta = [lambda, gamma]); the exponential fixes gamma = 1
    # (theta = [lambda]).
    model_specs = {
        'weibull': [0.001, 1.0],
        'exponential': [0.01],
    }

    models = {}
    for dist, init in model_specs.items():
        theta, covariance = fit_survival(t, delta, dist, init)
        models[dist] = {
            'theta': to_serializable(theta),
            'covariance': to_serializable(covariance),
        }

    times = [25.0, 50.0, 100.0, 150.0, 200.0]
    measures = ['survival', 'risk', 'hazard', 'cumulative_hazard', 'density']

    scenarios = []
    for dist in model_specs:
        theta = np.asarray(models[dist]['theta'])
        covariance = np.asarray(models[dist]['covariance'])
        for measure in measures:
            for alpha in [0.05, 0.10]:
                preds = survival_predictions(
                    times=times, theta=theta, covariance=covariance,
                    distribution=dist, measure=measure, alpha=alpha,
                )
                scenarios.append({
                    'distribution': dist,
                    'measure': measure,
                    'times': times,
                    'alpha': alpha,
                    'predicted': to_serializable(preds[:, 0]),
                    'variance': to_serializable(preds[:, 1]),
                    'lower': to_serializable(preds[:, 2]),
                    'upper': to_serializable(preds[:, 3]),
                })

    save_fixture('survival_predictions', {
        'models': models,
        'times': times,
        'scenarios': scenarios,
    })
    print(f'  {len(scenarios)} scenarios recorded')


if __name__ == '__main__':
    main()
