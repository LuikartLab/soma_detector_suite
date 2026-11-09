"""
mixed_stats.py — the Brain-paper statistical model, dependency-free (numpy/scipy).

Prina, Goyette et al., Brain 2026:  Stata
    xtmixed y i.genotype || animalid:, variance robust covariance(exchangeable)
    pwcompare genotype, mcompare(scheffe) effects
i.e. a linear mixed model with a random intercept per animal fitted by maximum
likelihood, fixed-effect standard errors from the cluster-robust (sandwich)
estimator clustered on animal, and Scheffé-adjusted pairwise contrasts
(for z statistics: p = P[chi2(k-1) >= z^2], k = number of groups).

fit_random_intercept(y, group, x)      -> dict (beta, cov_robust, sigma2_u, sigma2_e, loglik, ...)
genotype_model(values, animals, conds, reference) -> dict with the contrasts table

Validated against the Stata output stored with the published datasheets
(Supplemental Datasheets.xlsx: Figures 1c, 1e, 1g, 1i, 2): contrasts, robust
SEs, z and Scheffé p reproduce to the printed precision.
"""
import numpy as np
from scipy import optimize, stats


def _prepare(y, group, X):
    y = np.asarray(y, float); X = np.asarray(X, float); group = np.asarray(group)
    keys, inv = np.unique(group, return_inverse=True)
    idx = [np.where(inv == g)[0] for g in range(len(keys))]
    return y, X, idx, keys


def _profile(params, y, X, idx):
    """Negative log-likelihood (ML) at log-variances params = (log s2u, log s2e),
    plus the GLS beta and per-cluster pieces at that point."""
    s2u, s2e = np.exp(params)
    A = np.zeros((X.shape[1], X.shape[1])); b = np.zeros(X.shape[1])
    pieces = []
    for ix in idx:
        n = len(ix); Xi = X[ix]; yi = y[ix]
        # V^-1 = (I - c J)/s2e with c = s2u / (s2e + n s2u)
        c = s2u / (s2e + n * s2u)
        XtVi = (Xi.T - c * np.outer(Xi.sum(0), np.ones(n))) / s2e      # X' V^-1
        A += XtVi @ Xi; b += XtVi @ yi
        pieces.append((ix, c, XtVi))
    try:
        beta = np.linalg.solve(A, b)
    except np.linalg.LinAlgError:
        beta = np.linalg.pinv(A) @ b
    ll = 0.0
    for ix, c, XtVi in pieces:
        n = len(ix); r = y[ix] - X[ix] @ beta
        quad = (r @ r - c * r.sum() ** 2) / s2e
        logdet = (n - 1) * np.log(s2e) + np.log(s2e + n * s2u)
        ll += -0.5 * (logdet + quad + n * np.log(2 * np.pi))
    return -ll, beta, A, pieces


def fit_random_intercept(y, group, X, robust=True):
    """Maximum-likelihood random-intercept model y = X beta + u[group] + e."""
    y, X, idx, _keys = _prepare(y, group, X)
    # starting values from OLS residual variance
    beta0 = np.linalg.lstsq(X, y, rcond=None)[0]; v = np.var(y - X @ beta0)
    start = np.log([max(v * 0.3, 1e-8), max(v * 0.7, 1e-8)])
    lv = np.log(max(v, 1e-8))
    bounds = [(lv - 25, lv + 12), (lv - 25, lv + 12)]         # keep the simplex away from degenerate variances
    f = lambda p: _profile(p, y, X, idx)[0]
    res = optimize.minimize(f, start, method='Nelder-Mead', bounds=bounds,
                            options={'xatol': 1e-10, 'fatol': 1e-12, 'maxiter': 20000})
    res2 = optimize.minimize(f, res.x, method='L-BFGS-B', bounds=bounds, options={'ftol': 1e-14, 'gtol': 1e-10})
    if res2.fun <= res.fun:
        res = res2
    nll, beta, A, pieces = _profile(res.x, y, X, idx)
    s2u, s2e = np.exp(res.x)
    Ainv = np.linalg.inv(A)
    out = {'beta': beta, 'cov_model': Ainv, 'sigma2_u': s2u, 'sigma2_e': s2e, 'loglik': -nll,
           'n_obs': len(y), 'n_groups': len(idx), 'converged': bool(res.success)}
    if robust:
        # cluster-robust sandwich, clustered on the random-effect groups, with
        # Stata's small-sample factor G/(G-1)
        B = np.zeros_like(A)
        for ix, _c, XtVi in pieces:
            s = XtVi @ (y[ix] - X[ix] @ beta); B += np.outer(s, s)
        G = len(idx)
        out['cov_robust'] = Ainv @ B @ Ainv * (G / (G - 1.0))
    return out


def scheffe_p(z, k):
    """Scheffé-adjusted p-value of a z statistic among k groups."""
    return float(stats.chi2.sf(z * z, k - 1))


def genotype_model(values, animals, conds, reference, pairwise=False):
    """The paper's model on per-cell values: contrasts of every condition versus
    the reference (and optionally all pairwise contrasts), with robust SE, z,
    unadjusted and Scheffé p.  Returns dict(model=..., contrasts=[...])."""
    values = np.asarray(values, float); animals = np.asarray(animals).astype(str); conds = np.asarray(conds).astype(str)
    reference = str(reference)
    ok = np.isfinite(values); values, animals, conds = values[ok], animals[ok], conds[ok]
    if reference not in set(conds):
        raise ValueError('reference condition %r not present' % reference)
    levels = [reference] + sorted(c for c in np.unique(conds) if c != reference)
    k = len(levels)
    X = np.column_stack([np.ones(len(values))] + [(conds == lv).astype(float) for lv in levels[1:]])
    fit = fit_random_intercept(values, animals, X)
    cov = fit['cov_robust']; beta = fit['beta']
    rows = []

    def contrast(cvec, a, b):
        est = float(cvec @ beta); se = float(np.sqrt(cvec @ cov @ cvec)); z = est / se if se > 0 else np.nan
        rows.append({'contrast': '%s vs %s' % (a, b), 'group': a, 'reference': b, 'estimate': est, 'se_robust': se, 'z': z,
                     'p': float(2 * stats.norm.sf(abs(z))), 'p_scheffe': scheffe_p(z, k),
                     'n_group': int((conds == a).sum()), 'n_ref': int((conds == b).sum()),
                     'animals_group': int(len(np.unique(animals[conds == a]))), 'animals_ref': int(len(np.unique(animals[conds == b])))})
    for j, lv in enumerate(levels[1:], 1):
        c = np.zeros(k); c[j] = 1.0; contrast(c, lv, reference)
    if pairwise:
        for a in range(1, k):
            for b in range(a + 1, k):
                c = np.zeros(k); c[b] = 1.0; c[a] = -1.0; contrast(c, levels[b], levels[a])
    means = {lv: float(values[conds == lv].mean()) for lv in levels}
    n_cells = {lv: int((conds == lv).sum()) for lv in levels}
    n_anim = {lv: int(len(np.unique(animals[conds == lv]))) for lv in levels}
    return {'levels': levels, 'k': k, 'n_obs': fit['n_obs'], 'n_animals': fit['n_groups'], 'loglik': fit['loglik'],
            'sigma2_animal': fit['sigma2_u'], 'sigma2_residual': fit['sigma2_e'], 'intercept': float(beta[0]),
            'group_means': means, 'n_cells': n_cells, 'n_animals_by_group': n_anim,
            'beta': beta, 'cov_robust': cov, 'contrasts': rows, 'converged': fit['converged']}


def contrast(result, a, b):
    """Contrast b − a between two conditions of a fitted genotype_model():
    dict(estimate, se_robust, z, p, p_scheffe)."""
    levels, beta, cov, k = result['levels'], result['beta'], result['cov_robust'], result['k']
    c = np.zeros(k)
    if b != levels[0]:
        c[levels.index(b)] = 1.0
    if a != levels[0]:
        c[levels.index(a)] -= 1.0
    est = float(c @ beta); se = float(np.sqrt(max(c @ cov @ c, 0.0)))
    z = est / se if se > 0 else float('nan')
    return {'estimate': est, 'se_robust': se, 'z': z, 'p': float(2 * stats.norm.sf(abs(z))) if se > 0 else float('nan'),
            'p_scheffe': scheffe_p(z, k) if se > 0 else float('nan')}


def stars(p):
    return '****' if p < 0.0001 else '***' if p < 0.001 else '**' if p < 0.01 else '*' if p < 0.05 else 'ns'
