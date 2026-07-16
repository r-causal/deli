# /// script
# requires-python = ">=3.9"
# dependencies = [
#     "numpy",
#     "scipy",
# ]
# ///
"""Reference fixtures for exact-mode (forward-mode autodiff) sandwich variance.

Python Delicatessen computes the bread with deriv_method='exact' here, so the
bread, variance, and confidence intervals reflect exact derivatives rather than
the finite-difference default. The R port compares its own exact-mode fits
against these values.

Run from the Delicatessen root directory:
    uv run --with-editable . python r-pkg/deli/tests/testthat/fixtures/generate_exact_mode_fixtures.py

Outputs JSON fixture files to r-pkg/deli/tests/testthat/fixtures/
"""

import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..', '..', '..', '..'))

from delicatessen import MEstimator
from delicatessen.estimating_equations import (
    ee_regression,
    ee_glm,
    ee_mlogit,
    ee_beta_regression,
    ee_ridge_regression,
    ee_dlasso_regression,
    ee_robust_regression,
    ee_tobit,
    ee_aft,
    ee_plogit,
    ee_ipw,
    ee_gformula,
    ee_aipw,
    ee_gestimation_snmm,
    ee_2sls,
    ee_rogan_gladen_extended,
    ee_regression_calibration,
)

FIXTURES_DIR = os.path.dirname(__file__)


def to_serializable(obj):
    if isinstance(obj, np.ndarray):
        return obj.tolist()
    if isinstance(obj, (np.float64, np.float32)):
        return float(obj)
    if isinstance(obj, (np.int64, np.int32)):
        return int(obj)
    return obj


def extract_results(estr, alpha=0.05):
    results = {
        'theta': to_serializable(estr.theta),
        'variance': to_serializable(estr.variance),
        'asymptotic_variance': to_serializable(estr.asymptotic_variance),
        'bread': to_serializable(estr.bread),
    }
    ci = estr.confidence_intervals(alpha=alpha)
    results['ci_lower'] = to_serializable(ci[:, 0])
    results['ci_upper'] = to_serializable(ci[:, 1])
    return results


def save_fixture(name, data):
    filepath = os.path.join(FIXTURES_DIR, f'{name}.json')
    with open(filepath, 'w') as f:
        json.dump(data, f, indent=2)
    print(f'  Saved: {filepath}')


def design(seed, n=200):
    rng = np.random.default_rng(seed)
    x1 = rng.normal(size=n)
    x2 = rng.normal(size=n)
    X = np.column_stack([np.ones(n), x1, x2])
    return rng, X, x1, x2


def fit_exact(psi, init, maxiter=5000):
    estr = MEstimator(psi, init=init)
    estr.estimate(solver='lm', deriv_method='exact', maxiter=maxiter)
    return estr


def generate_regression():
    for model, seed in (('linear', 11), ('logistic', 12), ('poisson', 13)):
        rng, X, x1, x2 = design(seed)
        eta = 0.7 + 1.2 * x1 - 0.6 * x2
        if model == 'linear':
            y = eta + rng.normal(size=X.shape[0])
        elif model == 'logistic':
            y = rng.binomial(1, 1 / (1 + np.exp(-eta)))
        else:
            y = rng.poisson(np.exp(eta - 0.5))

        def psi(theta, X=X, y=y, model=model):
            return ee_regression(theta, X=X, y=y, model=model)

        estr = fit_exact(psi, [0., 0., 0.])
        print(f'Generating ee_regression_exact_{model}: theta={estr.theta}')
        save_fixture(f'ee_regression_exact_{model}', {
            'X': X.tolist(), 'y': to_serializable(y), 'model': model,
            'init': [0., 0., 0.], **extract_results(estr),
        })


def generate_glm_gamma():
    rng, X, x1, x2 = design(21)
    mu = np.exp(0.6 + 0.4 * x1 - 0.3 * x2)
    y = rng.gamma(shape=2.0, scale=mu / 2.0)

    def psi(theta):
        return ee_glm(theta, X=X, y=y, distribution='gamma', link='log')

    estr = fit_exact(psi, [0., 0., 0., 0.])
    print(f'Generating ee_glm_exact_gamma_log: theta={estr.theta}')
    save_fixture('ee_glm_exact_gamma_log', {
        'X': X.tolist(), 'y': to_serializable(y),
        'distribution': 'gamma', 'link': 'log',
        'init': [0., 0., 0., 0.], **extract_results(estr),
    })


def generate_ridge():
    rng, X, x1, x2 = design(31)
    y = 0.7 + 1.2 * x1 - 0.6 * x2 + rng.normal(size=X.shape[0])
    penalty = [0., 5., 5.]

    def psi(theta):
        return ee_ridge_regression(theta, X=X, y=y, model='linear', penalty=penalty)

    estr = fit_exact(psi, [0., 0., 0.])
    print(f'Generating ee_ridge_exact: theta={estr.theta}')
    save_fixture('ee_ridge_exact', {
        'X': X.tolist(), 'y': to_serializable(y),
        'model': 'linear', 'penalty': penalty,
        'init': [0., 0., 0.], **extract_results(estr),
    })


def generate_dlasso():
    rng, X, x1, x2 = design(41)
    y = 0.7 + 1.2 * x1 - 0.6 * x2 + rng.normal(size=X.shape[0])
    penalty = [0., 5., 5.]

    def psi(theta):
        return ee_dlasso_regression(theta, X=X, y=y, model='linear', penalty=penalty)

    estr = fit_exact(psi, [0.1, 0.1, 0.1])
    print(f'Generating ee_dlasso_exact: theta={estr.theta}')
    save_fixture('ee_dlasso_exact', {
        'X': X.tolist(), 'y': to_serializable(y),
        'model': 'linear', 'penalty': penalty,
        'init': [0.1, 0.1, 0.1], **extract_results(estr),
    })


def generate_robust_cauchy():
    rng, X, x1, x2 = design(51)
    y = 0.7 + 1.2 * x1 - 0.6 * x2 + rng.normal(size=X.shape[0])
    k = 1.345

    def psi(theta):
        return ee_robust_regression(theta, X=X, y=y, model='linear', k=k, loss='cauchy')

    estr = fit_exact(psi, [0.7, 1.2, -0.6], maxiter=20000)
    print(f'Generating ee_robust_exact_cauchy: theta={estr.theta}')
    save_fixture('ee_robust_exact_cauchy', {
        'X': X.tolist(), 'y': to_serializable(y),
        'model': 'linear', 'k': k, 'loss': 'cauchy',
        'init': [0.7, 1.2, -0.6], **extract_results(estr),
    })


def generate_robust_clipping_and_branching():
    # Losses whose score clips (huber) or branches (tukey, andrew, hampel).
    # Wider residual spread (scale=2) drives a meaningful fraction of the
    # residuals into the clipped or redescended regions, so the exact-mode
    # bread reflects the branch selection rather than a single smooth arm.
    # Hampel's a, b, and c parameters map to the R port's length-3 k vector
    # c(a, b, c); in Python they are the lower, upper, and k arguments.
    scale = 2.0
    losses = (
        ('huber', 61, dict(k=1.345), 'ee_robust_exact_huber'),
        ('tukey', 62, dict(k=4.685), 'ee_robust_exact_tukey'),
        ('andrew', 63, dict(k=1.339), 'ee_robust_exact_andrew'),
        ('hampel', 64, dict(k=6.0, lower=1.5, upper=3.0), 'ee_robust_exact_hampel'),
    )
    for loss, seed, kw, name in losses:
        rng, X, x1, x2 = design(seed)
        y = 0.7 + 1.2 * x1 - 0.6 * x2 + rng.normal(scale=scale, size=X.shape[0])

        def psi(theta, X=X, y=y, kw=kw, loss=loss):
            return ee_robust_regression(theta, X=X, y=y, model='linear', loss=loss, **kw)

        estr = fit_exact(psi, [0.7, 1.2, -0.6], maxiter=20000)
        print(f'Generating {name}: theta={estr.theta}')
        # The R port packs Hampel's (a, b, c) into a single length-3 k vector.
        if loss == 'hampel':
            k_field = [kw['lower'], kw['upper'], kw['k']]
        else:
            k_field = kw['k']
        save_fixture(name, {
            'X': X.tolist(), 'y': to_serializable(y),
            'model': 'linear', 'k': k_field, 'loss': loss,
            'init': [0.7, 1.2, -0.6], **extract_results(estr),
        })


def generate_glm_probit():
    # The probit inverse link routes the linear predictor through the standard
    # normal CDF and PDF, both of which must carry tangents under exact mode.
    rng, X, x1, x2 = design(71)
    eta = 0.5 + 1.0 * x1 - 0.7 * x2
    from scipy.stats import norm
    y = rng.binomial(1, norm.cdf(eta))

    def psi(theta):
        return ee_glm(theta, X=X, y=y, distribution='binomial', link='probit')

    estr = fit_exact(psi, [0., 0., 0.])
    print(f'Generating ee_glm_exact_probit: theta={estr.theta}')
    save_fixture('ee_glm_exact_probit', {
        'X': X.tolist(), 'y': to_serializable(y),
        'distribution': 'binomial', 'link': 'probit',
        'init': [0., 0., 0.], **extract_results(estr),
    })


def generate_tobit():
    # Left-censored outcome. The Tobit score divides the standard normal PDF by
    # the CDF of the scaled distance to the censoring limit, so both the PDF and
    # the CDF must carry tangents under exact mode, including through the CDF
    # clip.
    rng, X, x1, x2 = design(81)
    latent = 0.7 + 1.2 * x1 - 0.6 * x2 + rng.normal(size=X.shape[0])
    lower = -0.5
    y = np.where(latent <= lower, lower, latent)
    # The R port solves the root with rootSolve, which diverges from an all-zero
    # start on this score, so begin near the data-generating values. The
    # derivative is what the exact-mode comparison exercises, not the solver.
    init = [0.5, 1.0, -0.5, 0.2]

    def psi(theta):
        return ee_tobit(theta, X=X, y=y, lower=lower)

    estr = fit_exact(psi, init)
    print(f'Generating ee_tobit_exact: theta={estr.theta}')
    save_fixture('ee_tobit_exact', {
        'X': X.tolist(), 'y': to_serializable(y), 'lower': lower,
        'init': init, **extract_results(estr),
    })


def generate_aft_lognormal():
    # Accelerated failure time with a log-normal error. The survival term of the
    # score divides the standard normal PDF by one minus the CDF of the scaled
    # log-time, so both must carry tangents under exact mode.
    rng, X, x1, x2 = design(91)
    sigma_true = 0.7
    log_t = 0.7 + 1.2 * x1 - 0.6 * x2 + sigma_true * rng.normal(size=X.shape[0])
    t_event = np.exp(log_t)
    c_time = np.exp(rng.normal(loc=2.5, scale=1.0, size=X.shape[0]))
    t = np.minimum(t_event, c_time)
    delta = (t_event <= c_time).astype(float)
    # The log-normal AFT score is highly non-linear in theta and the solver
    # diverges from an all-zero start, so begin near the data-generating values
    # (theta[-1] = -log(sigma_true)). The fit then isolates the derivative rather
    # than the solver path.
    init = [0.7, 1.2, -0.6, -float(np.log(sigma_true))]

    def psi(theta):
        return ee_aft(theta, X=X, t=t, delta=delta, distribution='log-normal')

    estr = fit_exact(psi, init)
    print(f'Generating ee_aft_exact_lognormal: theta={estr.theta}')
    save_fixture('ee_aft_exact_lognormal', {
        'X': X.tolist(), 't': to_serializable(t), 'delta': to_serializable(delta),
        'distribution': 'log-normal',
        'init': init, **extract_results(estr),
    })


def _binary_causal(seed, n=300):
    # Confounded binary treatment and binary outcome with continuous confounders.
    rng = np.random.default_rng(seed)
    w1 = rng.normal(size=n)
    w2 = rng.normal(size=n)
    a = rng.binomial(1, 1 / (1 + np.exp(-(-0.5 + 0.5 * w1 + 0.3 * w2))))
    y = rng.binomial(1, 1 / (1 + np.exp(-(-0.2 + 0.8 * a + 0.5 * w1 - 0.3 * w2))))
    return n, w1, w2, a, y


def generate_ipw():
    n, w1, w2, a, y = _binary_causal(201)
    W = np.column_stack([np.ones(n), w1, w2])
    init = [0., 0.5, 0.5, 0., 0., 0.]

    def psi(theta):
        return ee_ipw(theta, y=y, A=a, W=W)

    estr = fit_exact(psi, init)
    print(f'Generating ee_ipw_exact: theta={estr.theta}')
    save_fixture('ee_ipw_exact', {
        'y': to_serializable(y), 'a': to_serializable(a), 'W': W.tolist(),
        'init': init, **extract_results(estr),
    })


def generate_gformula():
    n, w1, w2, a, y = _binary_causal(202)
    X = np.column_stack([np.ones(n), a, w1, w2])
    X1 = X.copy(); X1[:, 1] = 1
    X0 = X.copy(); X0[:, 1] = 0
    init = [0., 0.5, 0.5, 0., 0., 0., 0.]

    def psi(theta):
        return ee_gformula(theta, y=y, X=X, X1=X1, X0=X0)

    estr = fit_exact(psi, init)
    print(f'Generating ee_gformula_exact: theta={estr.theta}')
    save_fixture('ee_gformula_exact', {
        'y': to_serializable(y), 'X': X.tolist(), 'X1': X1.tolist(), 'X0': X0.tolist(),
        'init': init, **extract_results(estr),
    })


def generate_aipw():
    n, w1, w2, a, y = _binary_causal(203)
    W = np.column_stack([np.ones(n), w1, w2])
    X = np.column_stack([np.ones(n), a, w1, w2])
    X1 = X.copy(); X1[:, 1] = 1
    X0 = X.copy(); X0[:, 1] = 0
    init = [0., 0.5, 0.5, 0., 0., 0., 0., 0., 0., 0.]

    def psi(theta):
        return ee_aipw(theta, y=y, A=a, W=W, X=X, X1=X1, X0=X0)

    estr = fit_exact(psi, init)
    print(f'Generating ee_aipw_exact: theta={estr.theta}')
    save_fixture('ee_aipw_exact', {
        'y': to_serializable(y), 'a': to_serializable(a), 'W': W.tolist(),
        'X': X.tolist(), 'X1': X1.tolist(), 'X0': X0.tolist(),
        'init': init, **extract_results(estr),
    })


def generate_gestimation():
    # Inefficient (no outcome model) linear g-estimation with a two-column
    # structural mean model, so the SMM parameter subset is a tangent vector
    # rather than a lone scalar.
    rng = np.random.default_rng(204)
    n = 300
    w1 = rng.normal(size=n)
    v1 = rng.binomial(1, 0.5, n)
    a = rng.binomial(1, 1 / (1 + np.exp(-(0.25 + 0.5 * v1 + w1))))
    y = 5 - 2 * a - 0.5 * a * v1 + 3 * v1 + w1 + rng.normal(size=n)
    W = np.column_stack([np.ones(n), v1, w1])
    V = np.column_stack([np.ones(n), v1])
    init = [0.] * (V.shape[1] + W.shape[1])

    def psi(theta):
        return ee_gestimation_snmm(theta, y=y, A=a, W=W, V=V, model='linear')

    estr = fit_exact(psi, init)
    print(f'Generating ee_gestimation_exact: theta={estr.theta}')
    save_fixture('ee_gestimation_exact', {
        'y': to_serializable(y), 'a': to_serializable(a), 'W': W.tolist(), 'V': V.tolist(),
        'model': 'linear', 'init': init, **extract_results(estr),
    })


def generate_gestimation_single_v():
    # Inefficient (no outcome model) linear g-estimation with a single-column
    # structural mean model (V is the intercept alone), so the SMM parameter
    # subset theta[0:1] indexes down to a lone scalar rather than a tangent
    # vector. This is the case where the scalar branch of the matrix product
    # M %*% phi must carry tangents under exact mode.
    rng = np.random.default_rng(207)
    n = 300
    w1 = rng.normal(size=n)
    a = rng.binomial(1, 1 / (1 + np.exp(-(0.25 + w1))))
    y = 5 - 2 * a + w1 + rng.normal(size=n)
    W = np.column_stack([np.ones(n), w1])
    V = np.ones((n, 1))
    init = [0.] * (V.shape[1] + W.shape[1])

    def psi(theta):
        return ee_gestimation_snmm(theta, y=y, A=a, W=W, V=V, model='linear')

    estr = fit_exact(psi, init)
    print(f'Generating ee_gestimation_single_v_exact: theta={estr.theta}')
    save_fixture('ee_gestimation_single_v_exact', {
        'y': to_serializable(y), 'a': to_serializable(a), 'W': W.tolist(), 'V': V.tolist(),
        'model': 'linear', 'init': init, **extract_results(estr),
    })


def generate_gestimation_efficient():
    # Efficient (outcome model provided) linear g-estimation. The outcome model
    # regresses h(phi) = Y - (V*A) phi on X, and because h(phi) is theta-derived
    # it carries tangents under exact mode. The outcome design X routes that
    # tangent-carrying vector through ee_regression as the y argument.
    rng = np.random.default_rng(208)
    n = 300
    w1 = rng.normal(size=n)
    v1 = rng.binomial(1, 0.5, n)
    a = rng.binomial(1, 1 / (1 + np.exp(-(0.25 + 0.5 * v1 + w1))))
    y = 5 - 2 * a - 0.5 * a * v1 + 3 * v1 + w1 + rng.normal(size=n)
    W = np.column_stack([np.ones(n), v1, w1])
    V = np.column_stack([np.ones(n), v1])
    X = np.column_stack([np.ones(n), v1, w1])
    init = [0.] * (V.shape[1] + W.shape[1] + X.shape[1])

    def psi(theta):
        return ee_gestimation_snmm(theta, y=y, A=a, W=W, V=V, X=X, model='linear')

    estr = fit_exact(psi, init)
    print(f'Generating ee_gestimation_efficient_exact: theta={estr.theta}')
    save_fixture('ee_gestimation_efficient_exact', {
        'y': to_serializable(y), 'a': to_serializable(a),
        'W': W.tolist(), 'V': V.tolist(), 'X': X.tolist(),
        'model': 'linear', 'init': init, **extract_results(estr),
    })


def generate_2sls_with_w():
    # Two-stage least squares with exogenous covariates W. The second-stage
    # design stacks predicted A (theta-derived, tangent-carrying) with W and
    # routes it through ee_regression as the X argument.
    rng = np.random.default_rng(211)
    n = 500
    x_cov = rng.normal(size=n)
    z_raw = rng.binomial(1, 0.5, n)
    u = rng.normal(size=n)
    a = 0.5 * z_raw + 0.3 * x_cov + u + rng.normal(scale=0.5, size=n)
    y = 2 * a - u + 0.5 * x_cov + rng.normal(size=n)
    Z = np.column_stack([z_raw])
    W = np.column_stack([np.ones(n), x_cov])

    # Ordinary least squares starting values so the solver converges reliably.
    dmatrix1 = np.column_stack([Z, W])
    alpha_init = np.linalg.lstsq(dmatrix1, a, rcond=None)[0]
    a_hat_init = dmatrix1 @ alpha_init
    dmatrix2 = np.column_stack([a_hat_init, W])
    beta_init = np.linalg.lstsq(dmatrix2, y, rcond=None)[0]
    init = list(beta_init) + list(alpha_init)

    def psi(theta):
        return ee_2sls(theta, y=y, A=a, Z=Z, W=W)

    estr = fit_exact(psi, init)
    print(f'Generating ee_2sls_exact_w: theta={estr.theta}')
    save_fixture('ee_2sls_exact_w', {
        'y': to_serializable(y), 'a': to_serializable(a),
        'Z': Z.tolist(), 'W': W.tolist(),
        'init': [float(v) for v in init], **extract_results(estr),
    })


def generate_2sls_no_w():
    # Two-stage least squares without exogenous covariates. The intercept lives
    # in Z, so the second stage regresses Y on predicted A alone: the design is a
    # single tangent-carrying column and the second-stage coefficient is a lone
    # scalar, exercising the scalar branch of the tangent-aware matrix product.
    rng = np.random.default_rng(212)
    n = 500
    z_raw = rng.binomial(1, 0.5, n)
    u = rng.normal(size=n)
    a = 0.5 * z_raw + u + rng.normal(scale=0.5, size=n)
    y = 2 * a - u + rng.normal(size=n)
    Z = np.column_stack([np.ones(n), z_raw])

    dmatrix1 = Z
    alpha_init = np.linalg.lstsq(dmatrix1, a, rcond=None)[0]
    a_hat_init = dmatrix1 @ alpha_init
    dmatrix2 = a_hat_init.reshape(-1, 1)
    beta_init = np.linalg.lstsq(dmatrix2, y, rcond=None)[0]
    init = list(beta_init) + list(alpha_init)

    def psi(theta):
        return ee_2sls(theta, y=y, A=a, Z=Z)

    estr = fit_exact(psi, init)
    print(f'Generating ee_2sls_exact_now: theta={estr.theta}')
    save_fixture('ee_2sls_exact_now', {
        'y': to_serializable(y), 'a': to_serializable(a),
        'Z': Z.tolist(),
        'init': [float(v) for v in init], **extract_results(estr),
    })


def generate_rogan_gladen_extended():
    # Well-conditioned validation/main design: sensitivity and specificity both
    # depend on a covariate, so the extended correction is identified.
    rng = np.random.default_rng(205)
    n = 600
    r = rng.binomial(1, 0.55, n)
    x = rng.normal(size=n)
    X = np.column_stack([np.ones(n), x])
    y_true = rng.binomial(1, 0.4, n)
    sens = 1 / (1 + np.exp(-(1.1 + 0.5 * x)))
    spec = 1 / (1 + np.exp(-(1.3 - 0.4 * x)))
    y_star = rng.binomial(1, np.where(y_true == 1, sens, 1 - spec)).astype(float)
    y = y_true.astype(float)
    init = [0.4, 0.5, 0.0, 0.5, 0.0]

    def psi(theta):
        return ee_rogan_gladen_extended(theta, y=y, y_star=y_star, r=r, X=X)

    estr = fit_exact(psi, init)
    print(f'Generating ee_rogan_gladen_extended_exact: theta={estr.theta}')
    save_fixture('ee_rogan_gladen_extended_exact', {
        'y': to_serializable(y), 'y_star': to_serializable(y_star),
        'r': to_serializable(r), 'X': X.tolist(),
        'init': init, **extract_results(estr),
    })


def generate_regression_calibration():
    # Binary mismeasured action with a gold standard observed in the validation
    # subset. The calibration factor gamma_0 is well away from zero.
    rng = np.random.default_rng(206)
    n = 600
    r = rng.binomial(1, 0.5, n)
    a_star = rng.binomial(1, 0.5, n).astype(float)
    a = rng.binomial(1, np.where(a_star == 1, 0.8, 0.2)).astype(float)
    beta = 1.5
    init = [1.5, 0.6, 0.2]

    def psi(theta):
        return ee_regression_calibration(theta, beta=beta, a=a, a_star=a_star, r=r)

    estr = fit_exact(psi, init)
    print(f'Generating ee_regression_calibration_exact: theta={estr.theta}')
    save_fixture('ee_regression_calibration_exact', {
        'beta': beta, 'a': to_serializable(a), 'a_star': to_serializable(a_star),
        'r': to_serializable(r), 'init': init, **extract_results(estr),
    })


def generate_mlogit():
    # Multinomial logistic regression with three categories. The stacked score
    # multiplies each category residual by the design, and the shared denominator
    # couples the categories, so every parameter carries tangents through the
    # exponentials under exact mode.
    rng, X, x1, x2 = design(301)
    lp1 = 0.3 + 0.5 * x1 - 0.4 * x2
    lp2 = -0.2 + 0.8 * x1 + 0.3 * x2
    denom = 1 + np.exp(lp1) + np.exp(lp2)
    probs = np.column_stack([1 / denom, np.exp(lp1) / denom, np.exp(lp2) / denom])
    cats = np.array([rng.choice(3, p=p) for p in probs])
    y = np.column_stack([(cats == k).astype(float) for k in range(3)])
    init = [0.] * (X.shape[1] * (y.shape[1] - 1))

    def psi(theta):
        return ee_mlogit(theta, X=X, y=y)

    estr = fit_exact(psi, init)
    print(f'Generating ee_mlogit_exact: theta={estr.theta}')
    save_fixture('ee_mlogit_exact', {
        'X': X.tolist(), 'y': y.tolist(),
        'init': init, **extract_results(estr),
    })


def generate_beta_regression():
    # Beta regression with the mean-precision parameterization. The residual runs
    # the predicted mean through the digamma function twice, and the precision row
    # reshapes a tangent-carrying vector into a single-row matrix, so both the
    # regression coefficients and the log-precision parameter carry tangents.
    rng, X, x1, x2 = design(302)
    eta = 0.3 + 0.8 * x1 - 0.5 * x2
    mu = 1 / (1 + np.exp(-eta))
    phi = 5.0
    y = rng.beta(mu * phi, (1 - mu) * phi)
    y = np.clip(y, 1e-4, 1 - 1e-4)
    init = [0., 0., 0., 0.]

    def psi(theta):
        return ee_beta_regression(theta, X=X, y=y)

    estr = fit_exact(psi, init)
    print(f'Generating ee_beta_regression_exact: theta={estr.theta}')
    save_fixture('ee_beta_regression_exact', {
        'X': X.tolist(), 'y': to_serializable(y),
        'init': init, **extract_results(estr),
    })


def generate_plogit():
    # Pooled logistic regression for discrete-time survival. The design has no
    # intercept because the disjoint time indicators supply one. The covariate
    # log-odds are tiled across the unit-time intervals with a replication, and
    # the interval residuals are summed down the columns; both operations carry
    # tangents under exact mode. The intercept-free covariate keeps the model
    # identified. X is stored as a two-dimensional column so the R port reads it
    # as a single-column matrix rather than a vector.
    rng = np.random.default_rng(303)
    n = 300
    x1 = rng.normal(size=n)
    X = x1[:, None]
    t = rng.integers(1, 5, size=n).astype(float)
    delta = rng.binomial(1, 0.6, n).astype(float)
    n_time_steps = np.unique(t[delta == 1]).shape[0]
    init = [0.] * (X.shape[1] + n_time_steps)

    def psi(theta):
        return ee_plogit(theta, X=X, t=t, delta=delta)

    estr = fit_exact(psi, init)
    print(f'Generating ee_plogit_exact: theta={estr.theta}')
    save_fixture('ee_plogit_exact', {
        'X': X.tolist(), 'time': to_serializable(t), 'event': to_serializable(delta),
        'init': init, **extract_results(estr),
    })


if __name__ == '__main__':
    generate_regression()
    generate_mlogit()
    generate_beta_regression()
    generate_plogit()
    generate_glm_gamma()
    generate_ridge()
    generate_dlasso()
    generate_robust_cauchy()
    generate_robust_clipping_and_branching()
    generate_glm_probit()
    generate_tobit()
    generate_aft_lognormal()
    generate_ipw()
    generate_gformula()
    generate_aipw()
    generate_gestimation()
    generate_gestimation_single_v()
    generate_gestimation_efficient()
    generate_2sls_with_w()
    generate_2sls_no_w()
    generate_rogan_gladen_extended()
    generate_regression_calibration()
    print('Done.')
