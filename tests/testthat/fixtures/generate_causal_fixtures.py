# /// script
# requires-python = ">=3.9"
# dependencies = [
#     "numpy",
#     "scipy",
# ]
# ///
"""Generate reference results from Python Delicatessen for causal inference
estimating equations (ee_gformula, ee_ipw, ee_aipw).

Run from the Delicatessen root directory:
    uv run r-pkg/deli/tests/testthat/fixtures/generate_causal_fixtures.py

Outputs JSON fixture files to r-pkg/deli/tests/testthat/fixtures/
"""

import json
import os
import sys

import numpy as np

# Ensure the local delicatessen package is importable
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..', '..', '..', '..'))

from delicatessen import MEstimator
from delicatessen.estimating_equations import ee_gformula, ee_ipw, ee_aipw

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


def generate_causal_data():
    """Generate a simple binary treatment dataset for causal inference fixtures.

    Data generating process:
    - W1, W2: continuous confounders (standard normal)
    - A: binary treatment, depends on W1 and W2
    - Y: binary outcome, depends on A, W1, W2
    """
    np.random.seed(42)
    n = 200

    # Confounders
    W1 = np.random.normal(0, 1, n)
    W2 = np.random.normal(0, 1, n)

    # Treatment model: logit(P(A=1)) = -0.5 + 0.5*W1 + 0.3*W2
    prob_a = 1 / (1 + np.exp(-(-0.5 + 0.5 * W1 + 0.3 * W2)))
    A = np.random.binomial(1, prob_a)

    # Outcome model: logit(P(Y=1)) = -0.2 + 0.8*A + 0.5*W1 - 0.3*W2
    prob_y = 1 / (1 + np.exp(-(-0.2 + 0.8 * A + 0.5 * W1 - 0.3 * W2)))
    Y = np.random.binomial(1, prob_y)

    return n, W1, W2, A, Y


def generate_gformula_fixture():
    """Generate fixture for ee_gformula (average causal effect)."""
    print('Generating ee_gformula fixture...')

    n, W1, W2, A, Y = generate_causal_data()

    # Design matrices for outcome model: intercept, A, W1, W2
    X = np.column_stack([np.ones(n), A, W1, W2])

    # Counterfactual data: set A=1 for all
    X1 = X.copy()
    X1[:, 1] = 1

    # Counterfactual data: set A=0 for all
    X0 = X.copy()
    X0[:, 1] = 0

    # When X0 is provided: theta = [ATE, mu1, mu0, beta0, beta1, beta2, beta3]
    init = [0., 0.5, 0.5, 0., 0., 0., 0.]

    def psi(theta):
        return ee_gformula(theta, y=Y, X=X, X1=X1, X0=X0)

    estr = MEstimator(psi, init=init)
    estr.estimate(solver='lm')

    save_fixture('ee_gformula', {
        'y': Y.tolist(),
        'X': X.tolist(),
        'X1': X1.tolist(),
        'X0': X0.tolist(),
        'a': A.tolist(),
        'init': init,
        **extract_results(estr),
    })


def generate_ipw_fixture():
    """Generate fixture for ee_ipw (inverse probability weighting)."""
    print('Generating ee_ipw fixture...')

    n, W1, W2, A, Y = generate_causal_data()

    # Propensity score model design matrix: intercept, W1, W2
    W = np.column_stack([np.ones(n), W1, W2])

    # theta = [ATE, mu1, mu0, alpha0, alpha1, alpha2]
    init = [0., 0.5, 0.5, 0., 0., 0.]

    def psi(theta):
        return ee_ipw(theta, y=Y, A=A, W=W)

    estr = MEstimator(psi, init=init)
    estr.estimate(solver='lm')

    save_fixture('ee_ipw', {
        'y': Y.tolist(),
        'a': A.tolist(),
        'W': W.tolist(),
        'init': init,
        **extract_results(estr),
    })


def generate_aipw_fixture():
    """Generate fixture for ee_aipw (augmented inverse probability weighting)."""
    print('Generating ee_aipw fixture...')

    n, W1, W2, A, Y = generate_causal_data()

    # Propensity score model design matrix: intercept, W1, W2
    W = np.column_stack([np.ones(n), W1, W2])

    # Outcome model design matrix: intercept, A, W1, W2
    X = np.column_stack([np.ones(n), A, W1, W2])

    # Counterfactual data: set A=1 for all
    X1 = X.copy()
    X1[:, 1] = 1

    # Counterfactual data: set A=0 for all
    X0 = X.copy()
    X0[:, 1] = 0

    # theta = [ATE, mu1, mu0, alpha0, alpha1, alpha2, beta0, beta1, beta2, beta3]
    # 3 causal params + 3 propensity score params + 4 outcome model params = 10
    init = [0., 0.5, 0.5, 0., 0., 0., 0., 0., 0., 0.]

    def psi(theta):
        return ee_aipw(theta, y=Y, A=A, W=W, X=X, X1=X1, X0=X0)

    estr = MEstimator(psi, init=init)
    estr.estimate(solver='lm')

    save_fixture('ee_aipw', {
        'y': Y.tolist(),
        'a': A.tolist(),
        'W': W.tolist(),
        'X': X.tolist(),
        'X1': X1.tolist(),
        'X0': X0.tolist(),
        'init': init,
        **extract_results(estr),
    })


def main():
    os.makedirs(FIXTURES_DIR, exist_ok=True)
    generate_gformula_fixture()
    generate_ipw_fixture()
    generate_aipw_fixture()
    print('\nDone! All causal fixtures written.')


if __name__ == '__main__':
    main()
