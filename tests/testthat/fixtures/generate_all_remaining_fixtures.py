# /// script
# requires-python = ">=3.9"
# dependencies = [
#     "numpy",
#     "scipy",
# ]
# ///
"""Generate reference results from Python Delicatessen for all remaining
estimating equations that do not yet have JSON fixtures.

Run from the Delicatessen root directory:
    uv run r-pkg/deli/tests/testthat/fixtures/generate_all_remaining_fixtures.py

Outputs JSON fixture files to r-pkg/deli/tests/testthat/fixtures/
"""

import json
import os
import sys
import traceback

import numpy as np

# Ensure the local delicatessen package is importable
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..', '..', '..', '..'))

from delicatessen import MEstimator
from delicatessen.estimating_equations import (
    ee_ipw_msm,
    ee_gestimation_snmm,
    ee_iv_causal,
    ee_2sls,
    ee_mean_sensitivity_analysis,
    ee_glm,
    ee_mlogit,
    ee_beta_regression,
    ee_tobit,
    ee_lasso_regression,
    ee_dlasso_regression,
    ee_elasticnet_regression,
    ee_bridge_regression,
    ee_additive_regression,
    ee_survival_model,
    ee_aft,
    ee_plogit,
    ee_rogan_gladen,
    ee_rogan_gladen_extended,
    ee_regression_calibration,
    ee_regression,
    ee_emax,
    ee_emax_ed,
    ee_loglogistic,
    ee_loglogistic_ed,
)
from delicatessen.data import load_inderjit, load_breast_cancer
from delicatessen.utilities import inverse_logit

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


###############################################################################
# Shared data generators
###############################################################################

def generate_causal_data():
    """Generate a simple binary treatment dataset for causal inference fixtures.

    Data generating process:
    - W1, W2: continuous confounders (standard normal)
    - A: binary treatment, depends on W1 and W2
    - Y_continuous: continuous outcome, depends on A, W1, W2
    - Y_binary: binary outcome, depends on A, W1, W2
    """
    np.random.seed(42)
    n = 200

    # Confounders
    W1 = np.random.normal(0, 1, n)
    W2 = np.random.normal(0, 1, n)

    # Treatment model: logit(P(A=1)) = -0.5 + 0.5*W1 + 0.3*W2
    prob_a = 1 / (1 + np.exp(-(-0.5 + 0.5 * W1 + 0.3 * W2)))
    A = np.random.binomial(1, prob_a)

    # Binary outcome: logit(P(Y=1)) = -0.2 + 0.8*A + 0.5*W1 - 0.3*W2
    prob_y = 1 / (1 + np.exp(-(-0.2 + 0.8 * A + 0.5 * W1 - 0.3 * W2)))
    Y_binary = np.random.binomial(1, prob_y)

    # Continuous outcome: Y = 5 + 2*A + 1.5*W1 - 0.8*W2 + noise
    Y_continuous = 5 + 2 * A + 1.5 * W1 - 0.8 * W2 + np.random.normal(0, 1, n)

    return n, W1, W2, A, Y_binary, Y_continuous


def generate_regression_data():
    """Generate regression data for penalized and general regression fixtures."""
    np.random.seed(42)
    n = 200

    X1 = np.random.normal(0, 1, n)
    X2 = np.random.normal(0, 1, n)

    # Linear outcome
    Y_linear = 0.5 + 2 * X1 - 1 * X2 + np.random.normal(0, 1, n)

    # Binary outcome
    prob_y = 1 / (1 + np.exp(-(0.5 + 2 * X1 - 1 * X2)))
    Y_binary = np.random.binomial(1, prob_y)

    # Poisson outcome
    Y_poisson = np.random.poisson(lam=np.exp(0.5 + 0.5 * X1 - 0.3 * X2))

    return n, X1, X2, Y_linear, Y_binary, Y_poisson


###############################################################################
# Causal inference fixtures (remaining 5)
###############################################################################

def generate_ipw_msm_fixture():
    """Generate fixture for ee_ipw_msm (IPW marginal structural model)."""
    print('Generating ee_ipw_msm fixture...')

    n, W1, W2, A, Y_binary, Y_continuous = generate_causal_data()

    # Propensity score design matrix: intercept, W1, W2
    W = np.column_stack([np.ones(n), W1, W2])

    # MSM design matrix: intercept, A
    V = np.column_stack([np.ones(n), A])

    # theta = [MSM params (2), PS params (3)]
    # MSM: normal distribution, identity link for continuous outcome
    init = [0., 0., 0., 0., 0.]

    def psi(theta):
        return ee_ipw_msm(theta, y=Y_continuous, A=A, W=W, V=V,
                          distribution='normal', link='identity')

    estr = MEstimator(psi, init=init)
    estr.estimate(solver='lm')

    save_fixture('ee_ipw_msm', {
        'y': Y_continuous.tolist(),
        'a': A.tolist(),
        'W': W.tolist(),
        'V': V.tolist(),
        'distribution': 'normal',
        'link': 'identity',
        'init': init,
        **extract_results(estr),
    })


def generate_gestimation_snmm_fixture():
    """Generate fixture for ee_gestimation_snmm (g-estimation of SNMM)."""
    print('Generating ee_gestimation_snmm fixture...')

    n, W1, W2, A, Y_binary, Y_continuous = generate_causal_data()

    # PS model design matrix: intercept, W1, W2
    W = np.column_stack([np.ones(n), W1, W2])

    # Effect modifier matrix for SMM: intercept only (constant treatment effect)
    V = np.ones((n, 1))

    # theta = [SMM params (1), PS params (3)] = 4 total (inefficient g-estimator)
    init = [0., 0., 0., 0.]

    def psi(theta):
        return ee_gestimation_snmm(theta, y=Y_continuous, A=A,
                                   W=W, V=V, model='linear')

    estr = MEstimator(psi, init=init)
    estr.estimate(solver='lm')

    save_fixture('ee_gestimation_snmm', {
        'y': Y_continuous.tolist(),
        'a': A.tolist(),
        'W': W.tolist(),
        'V': V.tolist(),
        'model': 'linear',
        'init': init,
        **extract_results(estr),
    })


def generate_iv_causal_fixture():
    """Generate fixture for ee_iv_causal (instrumental variable)."""
    print('Generating ee_iv_causal fixture...')

    np.random.seed(42)
    n = 500

    # Instrument Z: binary, randomized
    Z = np.random.binomial(1, 0.5, n)

    # Unmeasured confounder U
    U = np.random.normal(0, 1, n)

    # Treatment A depends on Z and U (relevance + confounding)
    prob_a = 1 / (1 + np.exp(-(U + Z)))
    A = np.random.binomial(1, prob_a)

    # Outcome Y depends on A and U (causal effect = 2)
    Y = 2 * A - U + np.random.normal(0, 1, n)

    # theta = [causal effect, mean of Z]
    init = [0., 0.5]

    def psi(theta):
        return ee_iv_causal(theta, y=Y, A=A, Z=Z)

    estr = MEstimator(psi, init=init)
    estr.estimate(solver='lm')

    save_fixture('ee_iv_causal', {
        'y': Y.tolist(),
        'a': A.tolist(),
        'z': Z.tolist(),
        'init': init,
        **extract_results(estr),
    })


def generate_2sls_fixture():
    """Generate fixture for ee_2sls (two-stage least squares)."""
    print('Generating ee_2sls fixture...')

    np.random.seed(42)
    n = 500

    # Instrument Z: binary
    Z = np.random.binomial(1, 0.5, n)

    # Unmeasured confounder U
    U = np.random.normal(0, 1, n)

    # Treatment A depends on Z and U
    prob_a = 1 / (1 + np.exp(-(U + Z)))
    A = np.random.binomial(1, prob_a)

    # Outcome Y depends on A and U (causal effect = 2)
    Y = 2 * A - U + np.random.normal(0, 1, n)

    # Without exogenous covariates (W=None)
    # Z must be 2d for ee_2sls
    Z_matrix = Z[:, None]

    # theta = [second-stage (1 param), first-stage (1 param)] = 2 total
    init = [0., 0.]

    def psi(theta):
        return ee_2sls(theta, y=Y, A=A, Z=Z_matrix)

    estr = MEstimator(psi, init=init)
    estr.estimate(solver='lm')

    save_fixture('ee_2sls', {
        'y': Y.tolist(),
        'a': A.tolist(),
        'z': Z_matrix.tolist(),
        'init': init,
        **extract_results(estr),
    })


def generate_mean_sensitivity_analysis_fixture():
    """Generate fixture for ee_mean_sensitivity_analysis.

    SKIPPED: This function requires a callable H_function argument which
    cannot be serialized to JSON for use in R fixtures. The R implementation
    would need to independently define this function, making a direct
    comparison fixture less useful. Additionally, this EE involves missing
    data patterns that add complexity.
    """
    print('Generating ee_mean_sensitivity_analysis fixture...')

    np.random.seed(42)
    n = 200

    # Generate data with missing values
    W = np.random.binomial(1, p=0.5, size=n)
    Y_full = 200. - 35 * W + np.random.normal(scale=5, size=n)

    # Missing data mechanism
    prob_obs = 1 / (1 + np.exp(-(2 + W - 0.01 * Y_full)))
    M = np.random.binomial(1, p=prob_obs)

    # Set missing Y to 0 (will be masked by delta)
    Y = np.where(M == 0, 0, Y_full)

    # Design matrix for selection model
    X = np.column_stack([np.ones(n), W])

    # Sensitivity parameter alpha = 0 (MAR case)
    alpha_sens = 0.0
    q_eval = alpha_sens * Y  # q(Y; alpha) = alpha * Y, all zeros for MAR

    # theta = [mean, gamma_0, gamma_1]
    init = [200., 0., 0.]

    def psi(theta):
        return ee_mean_sensitivity_analysis(theta, y=Y, delta=M, X=X,
                                            q_eval=q_eval,
                                            H_function=inverse_logit)

    estr = MEstimator(psi, init=init)
    estr.estimate(solver='lm')

    save_fixture('ee_mean_sensitivity_analysis', {
        'y': Y.tolist(),
        'y_full': Y_full.tolist(),
        'delta': M.tolist(),
        'X': X.tolist(),
        'alpha_sens': alpha_sens,
        'q_eval': q_eval.tolist(),
        'H_function': 'inverse_logit',
        'init': init,
        **extract_results(estr),
    })


###############################################################################
# Penalized regression fixtures (4)
###############################################################################

def generate_lasso_regression_fixture():
    """Generate fixture for ee_lasso_regression (LASSO regression)."""
    print('Generating ee_lasso_regression fixture...')

    n, X1, X2, Y_linear, Y_binary, Y_poisson = generate_regression_data()

    # Design matrix with intercept
    X = np.column_stack([np.ones(n), X1, X2])

    # Don't penalize the intercept
    penalty = [0., 5., 5.]

    init = [0.01, 0.01, 0.01]

    def psi(theta):
        return ee_lasso_regression(theta, X=X, y=Y_linear,
                                   model='linear', penalty=penalty)

    estr = MEstimator(psi, init=init)
    estr.estimate(solver='lm', maxiter=20000)

    save_fixture('ee_lasso_regression', {
        'y': Y_linear.tolist(),
        'X': X.tolist(),
        'model': 'linear',
        'penalty': penalty,
        'init': init,
        **extract_results(estr),
    })


def generate_dlasso_regression_fixture():
    """Generate fixture for ee_dlasso_regression (differentiable LASSO)."""
    print('Generating ee_dlasso_regression fixture...')

    n, X1, X2, Y_linear, Y_binary, Y_poisson = generate_regression_data()

    # Design matrix with intercept
    X = np.column_stack([np.ones(n), X1, X2])

    # Don't penalize the intercept
    penalty = [0., 5., 5.]

    init = [0.01, 0.01, 0.01]

    def psi(theta):
        return ee_dlasso_regression(theta, X=X, y=Y_linear,
                                    model='linear', penalty=penalty)

    estr = MEstimator(psi, init=init)
    estr.estimate(solver='lm', maxiter=20000)

    save_fixture('ee_dlasso_regression', {
        'y': Y_linear.tolist(),
        'X': X.tolist(),
        'model': 'linear',
        'penalty': penalty,
        'init': init,
        **extract_results(estr),
    })


def generate_elasticnet_regression_fixture():
    """Generate fixture for ee_elasticnet_regression (elastic net)."""
    print('Generating ee_elasticnet_regression fixture...')

    n, X1, X2, Y_linear, Y_binary, Y_poisson = generate_regression_data()

    # Design matrix with intercept
    X = np.column_stack([np.ones(n), X1, X2])

    # Don't penalize the intercept
    penalty = [0., 5., 5.]
    ratio = 0.5

    init = [0.01, 0.01, 0.01]

    def psi(theta):
        return ee_elasticnet_regression(theta, X=X, y=Y_linear,
                                        model='linear', penalty=penalty,
                                        ratio=ratio)

    estr = MEstimator(psi, init=init)
    estr.estimate(solver='lm', maxiter=20000)

    save_fixture('ee_elasticnet_regression', {
        'y': Y_linear.tolist(),
        'X': X.tolist(),
        'model': 'linear',
        'penalty': penalty,
        'ratio': ratio,
        'init': init,
        **extract_results(estr),
    })


def generate_bridge_regression_fixture():
    """Generate fixture for ee_bridge_regression (bridge penalty)."""
    print('Generating ee_bridge_regression fixture...')

    n, X1, X2, Y_linear, Y_binary, Y_poisson = generate_regression_data()

    # Design matrix with intercept
    X = np.column_stack([np.ones(n), X1, X2])

    # Don't penalize the intercept
    penalty = [0., 5., 5.]
    gamma = 2.3

    init = [0., 0., 0.]

    def psi(theta):
        return ee_bridge_regression(theta, X=X, y=Y_linear,
                                    model='linear', penalty=penalty,
                                    gamma=gamma)

    estr = MEstimator(psi, init=init)
    estr.estimate(solver='lm')

    save_fixture('ee_bridge_regression', {
        'y': Y_linear.tolist(),
        'X': X.tolist(),
        'model': 'linear',
        'penalty': penalty,
        'gamma': gamma,
        'init': init,
        **extract_results(estr),
    })


###############################################################################
# Remaining regression fixtures (4)
###############################################################################

def generate_glm_normal_fixture():
    """Generate fixture for ee_glm with normal distribution, identity link."""
    print('Generating ee_glm (normal/identity) fixture...')

    n, X1, X2, Y_linear, Y_binary, Y_poisson = generate_regression_data()

    # Design matrix with intercept
    X = np.column_stack([np.ones(n), X1, X2])

    init = [0., 0., 0.]

    def psi(theta):
        return ee_glm(theta, X=X, y=Y_linear,
                      distribution='normal', link='identity')

    estr = MEstimator(psi, init=init)
    estr.estimate(solver='lm')

    save_fixture('ee_glm_normal_identity', {
        'y': Y_linear.tolist(),
        'X': X.tolist(),
        'distribution': 'normal',
        'link': 'identity',
        'init': init,
        **extract_results(estr),
    })


def generate_glm_poisson_fixture():
    """Generate fixture for ee_glm with Poisson distribution, log link."""
    print('Generating ee_glm (poisson/log) fixture...')

    n, X1, X2, Y_linear, Y_binary, Y_poisson = generate_regression_data()

    # Design matrix with intercept
    X = np.column_stack([np.ones(n), X1, X2])

    init = [0., 0., 0.]

    def psi(theta):
        return ee_glm(theta, X=X, y=Y_poisson,
                      distribution='poisson', link='log')

    estr = MEstimator(psi, init=init)
    estr.estimate(solver='lm')

    save_fixture('ee_glm_poisson_log', {
        'y': Y_poisson.tolist(),
        'X': X.tolist(),
        'distribution': 'poisson',
        'link': 'log',
        'init': init,
        **extract_results(estr),
    })


def generate_mlogit_fixture():
    """Generate fixture for ee_mlogit (multinomial logistic regression)."""
    print('Generating ee_mlogit fixture...')

    np.random.seed(42)
    n = 200

    # Covariate
    W = np.random.normal(0, 1, n)

    # Generate 3-category outcome
    # Category probabilities depend on W
    log_odds_2 = -0.5 + 0.8 * W  # log-odds of cat 2 vs cat 1
    log_odds_3 = 0.3 - 0.5 * W   # log-odds of cat 3 vs cat 1
    denom = 1 + np.exp(log_odds_2) + np.exp(log_odds_3)
    p1 = 1 / denom
    p2 = np.exp(log_odds_2) / denom
    p3 = np.exp(log_odds_3) / denom

    # Sample from multinomial
    Y_cat = np.zeros(n, dtype=int)
    for i in range(n):
        Y_cat[i] = np.random.choice([1, 2, 3], p=[p1[i], p2[i], p3[i]])

    # Create indicator matrix (reference = category 1)
    Y1 = (Y_cat == 1).astype(float)
    Y2 = (Y_cat == 2).astype(float)
    Y3 = (Y_cat == 3).astype(float)
    Y_matrix = np.column_stack([Y1, Y2, Y3])

    # Design matrix: intercept, W
    X = np.column_stack([np.ones(n), W])

    # theta: 2 params per non-reference category = 2 * 2 = 4
    init = [0., 0., 0., 0.]

    def psi(theta):
        return ee_mlogit(theta, X=X, y=Y_matrix)

    estr = MEstimator(psi, init=init)
    estr.estimate(solver='lm')

    save_fixture('ee_mlogit', {
        'y': Y_matrix.tolist(),
        'y_cat': Y_cat.tolist(),
        'X': X.tolist(),
        'init': init,
        **extract_results(estr),
    })


def generate_beta_regression_fixture():
    """Generate fixture for ee_beta_regression (beta regression)."""
    print('Generating ee_beta_regression fixture...')

    np.random.seed(42)
    n = 200

    # Covariate
    W = np.random.normal(0, 1, n)

    # Beta-distributed outcome with mean depending on W
    # logit(mu) = 0.5 + 0.8*W
    mu = 1 / (1 + np.exp(-(0.5 + 0.8 * W)))
    phi = 10  # precision
    alpha_param = mu * phi
    beta_param = (1 - mu) * phi
    Y = np.random.beta(alpha_param, beta_param)

    # Clip to avoid exact 0 or 1
    Y = np.clip(Y, 0.001, 0.999)

    # Design matrix: intercept, W
    X = np.column_stack([np.ones(n), W])

    # theta: 2 regression params + 1 log(precision) = 3
    init = [0., 0., 0.]

    def psi(theta):
        return ee_beta_regression(theta, X=X, y=Y)

    estr = MEstimator(psi, init=init)
    estr.estimate(solver='lm')

    save_fixture('ee_beta_regression', {
        'y': Y.tolist(),
        'X': X.tolist(),
        'init': init,
        **extract_results(estr),
    })


def generate_tobit_fixture():
    """Generate fixture for ee_tobit (Tobit regression)."""
    print('Generating ee_tobit fixture...')

    np.random.seed(42)
    n = 300

    # Covariates
    X1 = np.random.normal(0, 1, n)

    # Latent outcome
    Y_latent = 0.5 + 2 * X1 + np.random.normal(0, 1, n)

    # Censoring: left-censor at 0
    lower = 0.0
    Y_censored = np.clip(Y_latent, a_min=lower, a_max=None)

    # Design matrix with intercept
    X = np.column_stack([np.ones(n), X1])

    # theta: 2 regression params + 1 log(sigma) = 3
    init = [np.mean(Y_censored), 0., 0.]

    def psi(theta):
        return ee_tobit(theta, X=X, y=Y_censored, lower=lower)

    estr = MEstimator(psi, init=init)
    estr.estimate(solver='lm')

    save_fixture('ee_tobit', {
        'y': Y_censored.tolist(),
        'X': X.tolist(),
        'lower': lower,
        'init': init,
        **extract_results(estr),
    })


###############################################################################
# Additive regression fixture (1)
###############################################################################

def generate_additive_regression_fixture():
    """Generate fixture for ee_additive_regression (GAM).

    NOTE: ee_additive_regression requires a 'specifications' argument that
    is a list of dicts controlling spline knots, penalty, etc. This is complex
    to replicate in R. We generate a simple fixture with a small number of
    knots for testing purposes.
    """
    print('Generating ee_additive_regression fixture...')

    np.random.seed(42)
    n = 300

    # Generate data with non-linear relationship
    X1 = np.random.uniform(-3, 3, n)
    Y = np.sin(X1) + 0.5 * X1 + np.random.normal(0, 0.5, n)

    # Design matrix: intercept, X1
    X = np.column_stack([np.ones(n), X1])

    # Specifications: no spline for intercept, spline for X1
    specs = [None, {"knots": [-2, -1, 0, 1, 2], "penalty": 3, "natural": True, "power": 3}]

    # Determine number of parameters using the additive_design_matrix utility
    from delicatessen.utilities import additive_design_matrix
    Xa = additive_design_matrix(X=X, specifications=specs)
    n_params = Xa.shape[1]

    init = [0.] * n_params

    def psi(theta):
        return ee_additive_regression(theta, X=X, y=Y,
                                      model='linear',
                                      specifications=specs)

    estr = MEstimator(psi, init=init)
    estr.estimate(solver='lm')

    save_fixture('ee_additive_regression', {
        'y': Y.tolist(),
        'X': X.tolist(),
        'model': 'linear',
        'specifications': specs,
        'n_params': n_params,
        'init': init,
        **extract_results(estr),
    })


###############################################################################
# Survival fixtures (3)
###############################################################################

def generate_survival_model_fixture():
    """Generate fixture for ee_survival_model (parametric survival)."""
    print('Generating ee_survival_model fixture...')

    np.random.seed(42)
    n = 200

    # Generate Weibull survival times
    T_event = 0.8 * np.random.weibull(a=0.8, size=n)

    # Generate censoring times
    C = np.random.weibull(a=1, size=n)
    C = np.where(C > 5, 5, C)

    # Observed times and event indicators
    delta = np.where(T_event < C, 1, 0)
    t = np.where(delta == 1, T_event, C)

    # Weibull model: theta = [lambda, gamma]
    init = [1., 1.]

    def psi(theta):
        return ee_survival_model(theta, t=t, delta=delta,
                                 distribution='weibull')

    estr = MEstimator(psi, init=init)
    estr.estimate(solver='lm')

    save_fixture('ee_survival_model_weibull', {
        't': t.tolist(),
        'delta': delta.tolist(),
        'distribution': 'weibull',
        'init': init,
        **extract_results(estr),
    })


def generate_aft_fixture():
    """Generate fixture for ee_aft (accelerated failure time model)."""
    print('Generating ee_aft fixture...')

    np.random.seed(42)
    n = 200

    # Covariate: binary treatment
    X1 = np.random.binomial(1, 0.5, n)

    # Generate Weibull survival times with AFT structure
    T_event = (1/1.25 + 1/np.exp(0.5) * X1) * np.random.weibull(a=0.75, size=n)

    # Censoring
    C = np.random.weibull(a=1, size=n)
    C = np.where(C > 10, 10, C)

    delta = np.where(T_event < C, 1, 0)
    t = np.where(delta == 1, T_event, C)

    # Design matrix: intercept, X1
    X = np.column_stack([np.ones(n), X1])

    # theta = [beta_0, beta_X, log(sigma)] = 3 params
    init = [2., 0., 0.]

    def psi(theta):
        return ee_aft(theta, X=X, t=t, delta=delta, distribution='weibull')

    estr = MEstimator(psi, init=init)
    estr.estimate(solver='lm')

    save_fixture('ee_aft_weibull', {
        't': t.tolist(),
        'delta': delta.tolist(),
        'X': X.tolist(),
        'distribution': 'weibull',
        'init': init,
        **extract_results(estr),
    })


def generate_plogit_fixture():
    """Generate fixture for ee_plogit (pooled logistic regression)."""
    print('Generating ee_plogit fixture...')

    # Use the built-in breast cancer data
    d = load_breast_cancer()
    delta = d[:, 0].astype(int)
    t = d[:, 1].astype(int)
    stain = d[:, 2]

    # X is the covariate design matrix (no intercept for time)
    X = stain[:, None]

    # Use a simple linear time specification for S
    # This corresponds to a Gompertz-like model
    t_steps = np.arange(1, np.max(t) + 1)
    intercept_s = np.ones(t_steps.shape)[:, None]
    s_matrix = np.concatenate([intercept_s, t_steps[:, None]], axis=1)

    # theta = [stain coeff, intercept_time, slope_time] = 3
    init = [0., -3., 0.]

    def psi(theta):
        return ee_plogit(theta, X=X, t=t, delta=delta, S=s_matrix)

    estr = MEstimator(psi, init=init)
    estr.estimate(solver='lm')

    save_fixture('ee_plogit', {
        'delta': delta.tolist(),
        't': t.tolist(),
        'X': X.tolist(),
        'S': s_matrix.tolist(),
        'init': init,
        **extract_results(estr),
    })


def generate_plogit_predict_fixture():
    """Generate fixture for plogit_predict (predictions from pooled logistic).

    Fits ee_plogit end to end on the breast cancer data with a Gompertz-like
    time specification, then records plogit_predict output at a handful of
    times. Used to pin the ee_plogit -> plogit_predict pipeline with the
    consistent time/event semantics (Python: t = observed times, delta =
    event indicator).
    """
    print('Generating plogit_predict fixture...')
    from delicatessen.utilities import plogit_predict

    # Same data and time design as generate_plogit_fixture
    d = load_breast_cancer()
    delta = d[:, 0].astype(int)
    t = d[:, 1].astype(int)
    stain = d[:, 2]
    X = stain[:, None]

    t_steps = np.arange(1, np.max(t) + 1)
    intercept_s = np.ones(t_steps.shape)[:, None]
    s_matrix = np.concatenate([intercept_s, t_steps[:, None]], axis=1)

    init = [0., -3., 0.]

    def psi(theta):
        return ee_plogit(theta, X=X, t=t, delta=delta, S=s_matrix)

    estr = MEstimator(psi, init=init)
    estr.estimate(solver='lm')

    # Predicted risk at several interior times (all in [1, max(t)])
    times_to_predict = [25, 50, 100, 150, 200]
    measure = 'risk'
    predictions = plogit_predict(theta=estr.theta, t=t, delta=delta, X=X,
                                 S=s_matrix, times_to_predict=times_to_predict,
                                 measure=measure)

    save_fixture('plogit_predict', {
        'delta': delta.tolist(),
        't': t.tolist(),
        'X': X.tolist(),
        'S': s_matrix.tolist(),
        'init': init,
        'theta': to_serializable(estr.theta),
        'times_to_predict': times_to_predict,
        'measure': measure,
        'predictions': to_serializable(predictions),
    })


###############################################################################
# Measurement error fixtures (3)
###############################################################################

def generate_rogan_gladen_fixture():
    """Generate fixture for ee_rogan_gladen (misclassification correction)."""
    print('Generating ee_rogan_gladen fixture...')

    # Replicate the published Cole et al. (2023) example
    # Main study data (R=1): only Y_star observed
    # External data (R=0): both Y and Y_star observed
    Y_star = np.array([0, 1] + [0, 1, 0, 1])
    Y = np.array([0, 0] + [0, 0, 1, 1])  # placeholder for main study
    S = np.array([1, 1] + [0, 0, 0, 0])  # R indicator (1=main, 0=external)
    counts = np.array([270, 680] + [71, 18, 38, 203])

    # Expand data
    Y_star_exp = np.repeat(Y_star, counts)
    Y_exp = np.repeat(Y, counts)
    S_exp = np.repeat(S, counts)

    # theta = [corrected mean, naive mean, sensitivity, specificity]
    init = [0.5, 0.5, 0.75, 0.75]

    def psi(theta):
        return ee_rogan_gladen(theta, y=Y_exp, y_star=Y_star_exp, r=S_exp)

    estr = MEstimator(psi, init=init)
    estr.estimate(solver='lm')

    save_fixture('ee_rogan_gladen', {
        'y': Y_exp.tolist(),
        'y_star': Y_star_exp.tolist(),
        'r': S_exp.tolist(),
        'init': init,
        **extract_results(estr),
    })


def generate_rogan_gladen_extended_fixture():
    """Generate fixture for ee_rogan_gladen_extended."""
    print('Generating ee_rogan_gladen_extended fixture...')

    # Same Cole et al. example but with intercept-only design matrix
    Y_star = np.array([0, 1] + [0, 1, 0, 1])
    Y = np.array([0, 0] + [0, 0, 1, 1])
    S = np.array([1, 1] + [0, 0, 0, 0])
    counts = np.array([270, 680] + [71, 18, 38, 203])

    Y_star_exp = np.repeat(Y_star, counts)
    Y_exp = np.repeat(Y, counts)
    S_exp = np.repeat(S, counts)

    # Design matrix for Se/Sp models: intercept only
    X = np.ones((len(Y_exp), 1))

    # theta = [corrected mean, sens (logit), spec (logit)] = 3
    init = [0.5, 1., 1.]

    def psi(theta):
        return ee_rogan_gladen_extended(theta, y=Y_exp, y_star=Y_star_exp,
                                        r=S_exp, X=X)

    estr = MEstimator(psi, init=init)
    estr.estimate(solver='lm')

    save_fixture('ee_rogan_gladen_extended', {
        'y': Y_exp.tolist(),
        'y_star': Y_star_exp.tolist(),
        'r': S_exp.tolist(),
        'X': X.tolist(),
        'init': init,
        **extract_results(estr),
    })


def generate_regression_calibration_fixture():
    """Generate fixture for ee_regression_calibration."""
    print('Generating ee_regression_calibration fixture...')

    # Following the docstring example from the Python source
    A_star = np.array([0, 1, 0, 1] + [1, 0, 1, 0])
    A = np.array([0, 0, 0, 0] + [1, 1, 0, 0])  # placeholder for main
    Y_vals = np.array([0, 0, 1, 1] + [0, 1, 0, 0])  # placeholder for external
    S = np.array([1, 1, 1, 1] + [0, 0, 0, 0])
    counts = np.array([266, 67, 400, 267] + [90, 10, 20, 80])

    A_star_exp = np.repeat(A_star, counts)
    A_exp = np.repeat(A, counts)
    Y_exp = np.repeat(Y_vals, counts)
    S_exp = np.repeat(S, counts)

    # First fit the naive model in the main study only
    X_main = np.column_stack([np.ones(len(A_star_exp)), A_star_exp])
    main_mask = S_exp == 1

    # Stack: calibration EE (3 params) + naive logistic (2 params) = 5
    init = [0., 0.75, 0., 0., 0.]

    def psi(theta):
        theta_calib = theta[:3]
        theta_main = theta[3:]

        # Replace NaN-like values for Y in external data
        y_safe = np.where(S_exp == 0, -999, Y_exp).astype(float)

        # Naive logistic regression (main study only)
        ee_logit = ee_regression(theta=theta_main, y=y_safe,
                                 X=X_main, model='logistic')
        ee_logit = ee_logit * S_exp  # Only main study contributes

        # Regression calibration
        ee_calib = ee_regression_calibration(theta=theta_calib,
                                             beta=theta_main[1],
                                             a=A_exp,
                                             a_star=A_star_exp,
                                             r=S_exp)

        return np.vstack([ee_calib, ee_logit])

    estr = MEstimator(psi, init=init)
    estr.estimate(solver='lm')

    save_fixture('ee_regression_calibration', {
        'a': A_exp.tolist(),
        'a_star': A_star_exp.tolist(),
        'y': Y_exp.tolist(),
        'r': S_exp.tolist(),
        'init': init,
        **extract_results(estr),
    })


###############################################################################
# Pharmacokinetics fixtures (2)
###############################################################################

def generate_emax_ed_fixture():
    """Generate fixture for ee_emax + ee_emax_ed (stacked for ED25)."""
    print('Generating ee_emax_ed fixture...')

    d = load_inderjit()
    dose = d[:, 1]
    # Invert response for Emax (increasing dose -> increasing effect)
    response = np.max(d[:, 0]) - d[:, 0]

    # theta = [e_zero, e_max, ed50, ed25]
    init = [np.min(response), np.max(response),
            np.median(dose), np.median(dose)]

    def psi(theta):
        emax_ee = ee_emax(theta=theta[:3], dose=dose, response=response)
        ed_25 = ee_emax_ed(theta[3], dose=dose, delta=0.25,
                           ed50=theta[2])
        return np.vstack((emax_ee, ed_25))

    estr = MEstimator(psi, init=init)
    estr.estimate(solver='lm')

    save_fixture('ee_emax_ed', {
        'dose': dose.tolist(),
        'response': response.tolist(),
        'delta': 0.25,
        'init': init,
        **extract_results(estr),
    })


def generate_loglogistic_ed_fixture():
    """Generate fixture for ee_loglogistic + ee_loglogistic_ed (stacked for ED25)."""
    print('Generating ee_loglogistic_ed fixture...')

    d = load_inderjit()
    dose = d[:, 1]
    response = d[:, 0]

    # 3-parameter log-logistic (lower=0) + ED25
    # theta = [upper, ed50, steepness, ed25]
    lower = 0
    init = [np.max(response),
            (np.max(response) + np.min(response)) / 2,
            (np.max(response) + np.min(response)) / 2,
            (np.max(response) + np.min(response)) / 2]

    def psi(theta):
        # 3PL: fix lower=0, estimate upper, ed50, steepness
        ll_ee = ee_loglogistic(theta=[lower] + list(theta[:3]),
                               dose=dose, response=response)
        ll_ee = ll_ee[1:, :]  # Drop lower limit EE

        ed_25 = ee_loglogistic_ed(theta[3], dose=dose, delta=0.25,
                                  lower=lower, upper=theta[0],
                                  ed50=theta[1], steepness=theta[2])
        return np.vstack((ll_ee, ed_25))

    estr = MEstimator(psi, init=init)
    estr.estimate(solver='lm')

    save_fixture('ee_loglogistic_ed', {
        'dose': dose.tolist(),
        'response': response.tolist(),
        'lower': lower,
        'delta': 0.25,
        'init': init,
        **extract_results(estr),
    })


###############################################################################
# Main
###############################################################################

def main():
    os.makedirs(FIXTURES_DIR, exist_ok=True)

    generators = [
        # Causal (remaining 5)
        generate_ipw_msm_fixture,
        generate_gestimation_snmm_fixture,
        generate_iv_causal_fixture,
        generate_2sls_fixture,
        generate_mean_sensitivity_analysis_fixture,

        # Penalized regression (4)
        generate_lasso_regression_fixture,
        generate_dlasso_regression_fixture,
        generate_elasticnet_regression_fixture,
        generate_bridge_regression_fixture,

        # Remaining regression (4 + 1 Poisson GLM)
        generate_glm_normal_fixture,
        generate_glm_poisson_fixture,
        generate_mlogit_fixture,
        generate_beta_regression_fixture,
        generate_tobit_fixture,

        # Additive regression (1)
        generate_additive_regression_fixture,

        # Survival (3)
        generate_survival_model_fixture,
        generate_aft_fixture,
        generate_plogit_fixture,
        generate_plogit_predict_fixture,

        # Measurement (3)
        generate_rogan_gladen_fixture,
        generate_rogan_gladen_extended_fixture,
        generate_regression_calibration_fixture,

        # Pharmacokinetics (2)
        generate_emax_ed_fixture,
        generate_loglogistic_ed_fixture,
    ]

    succeeded = 0
    failed = 0

    for gen in generators:
        try:
            gen()
            succeeded += 1
        except Exception as e:
            print(f'  FAILED: {gen.__name__}: {e}')
            traceback.print_exc()
            failed += 1

    print(f'\nDone! {succeeded} fixtures generated, {failed} failed.')


if __name__ == '__main__':
    main()
