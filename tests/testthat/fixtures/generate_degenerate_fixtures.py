# /// script
# requires-python = ">=3.9"
# dependencies = [
#     "numpy",
#     "scipy",
# ]
# ///
"""Regenerate degenerate fixtures (ee_mean_geometric, ee_percentile_50,
ee_loglogistic) with better data and init values.

Run from the Delicatessen root directory:
    uv run r-pkg/deli/tests/testthat/fixtures/generate_degenerate_fixtures.py

Outputs JSON fixture files to r-pkg/deli/tests/testthat/fixtures/
"""

import json
import os
import sys

import numpy as np

# Ensure the local delicatessen package is importable
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..', '..', '..', '..'))

from delicatessen import MEstimator
from delicatessen.estimating_equations import ee_mean_geometric, ee_percentile, ee_loglogistic

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


def generate_mean_geometric():
    """Generate ee_mean_geometric fixture with positive data and good init."""
    print('Generating ee_mean_geometric fixture...')

    np.random.seed(42)
    # Positive data (geometric mean requires positive values)
    y = np.exp(np.random.normal(1.5, 0.5, 100))

    init = [float(np.mean(y))]  # reasonable starting point

    def psi(theta):
        return ee_mean_geometric(theta=theta, y=y)

    estr = MEstimator(psi, init=init)
    estr.estimate(solver='lm')

    print(f'  theta: {estr.theta}')

    save_fixture('ee_mean_geometric', {
        'y': y.tolist(),
        'init': init,
        **extract_results(estr),
    })


def generate_percentile_50():
    """Generate ee_percentile fixture with proper solver settings."""
    print('Generating ee_percentile_50 fixture...')

    np.random.seed(42)
    y = np.random.normal(5, 2, 200)

    init = [np.median(y)]  # start near the true median

    def psi(theta):
        return ee_percentile(theta=theta, y=y, q=0.5)

    estr = MEstimator(psi, init=init)
    # Percentile EE is non-smooth; need special solver settings
    estr.estimate(solver='hybr', tolerance=1e-3, dx=1)

    print(f'  theta: {estr.theta}')

    save_fixture('ee_percentile_50', {
        'y': y.tolist(),
        'q': 0.5,
        'init': [float(init[0])],
        **extract_results(estr),
    })


def generate_loglogistic():
    """Generate ee_loglogistic fixture with clear dose-response curve."""
    print('Generating ee_loglogistic fixture...')

    np.random.seed(42)
    n_per_dose = 15
    doses = np.array([0, 1, 5, 10, 25, 50, 100])

    # True parameters: lower=2, upper=10, ed50=20, steepness=2
    true_lower = 2.0
    true_upper = 10.0
    true_ed50 = 20.0
    true_steep = 2.0

    dose_all = np.repeat(doses, n_per_dose)
    # Generate response from true model + noise
    dose_calc = np.where(dose_all == 0, 1e-9, dose_all)
    rho = (dose_calc / true_ed50) ** true_steep
    true_response = true_lower + (true_upper - true_lower) / (1 + rho)
    response_all = true_response + np.random.normal(0, 0.3, len(dose_all))

    init = [true_lower, true_upper, true_ed50, true_steep]

    def psi(theta):
        return ee_loglogistic(theta=theta, dose=dose_all, response=response_all)

    estr = MEstimator(psi, init=init)
    estr.estimate(solver='lm')

    print(f'  theta: {estr.theta}')

    save_fixture('ee_loglogistic', {
        'dose': dose_all.tolist(),
        'response': response_all.tolist(),
        'init': init,
        **extract_results(estr),
    })


def main():
    generate_mean_geometric()
    generate_percentile_50()
    generate_loglogistic()
    print('\nDone! Degenerate fixtures regenerated.')


if __name__ == '__main__':
    main()
