#!/usr/bin/env python3
"""Replication statistics for a (discovery, replication) pair of DE results.

Every function here takes two DESeq2 results frames and **first intersects their testable
genes**.  That intersection is not cosmetic: `build_cohort` filters for expression per
cohort and DESeq2 applies independent filtering per run, so without it the "replication
rate" would be diluted by genes that the replication study never had the opportunity to
call.  This is the single most common way to get this number wrong.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from scipy import stats

FDR = 0.05
LFC = 1.0


# --------------------------------------------------------------------------------------
# gene selection rules (Shi et al. 2008 show the rule alone moves overlap 20% -> 90%)
# --------------------------------------------------------------------------------------
def select(df: pd.DataFrame, rule: str) -> pd.Index:
    """Return the genes a study would 'report' under a given selection rule."""
    if rule == "joint":            # FDR < 0.05 AND |log2FC| > 1   (the hypothesis's rule)
        m = (df.padj < FDR) & (df.log2FoldChange.abs() > LFC)
        return df.index[m]
    if rule == "fdr":              # FDR < 0.05 only
        return df.index[df.padj < FDR]
    if rule == "lfc":              # |log2FC| > 1 only, ignoring significance
        return df.index[df.log2FoldChange.abs() > LFC]
    if rule.startswith("topp"):    # top-N by p-value
        n = int(rule[4:])
        return df.nsmallest(n, "pvalue").index
    if rule.startswith("topfc"):   # top-N by |log2FC| among FDR<0.05 (Shi's FC-ranking)
        n = int(rule[5:])
        sub = df[df.padj < FDR]
        return sub.log2FoldChange.abs().sort_values(ascending=False).index[:n]
    raise ValueError(rule)


RULES = ["joint", "fdr", "lfc", "topp100", "topp500", "topfc100", "topfc500"]


# --------------------------------------------------------------------------------------
# Storey's pi1
# --------------------------------------------------------------------------------------
def storey_pi0(p: np.ndarray, lambdas: np.ndarray | None = None) -> float:
    """Storey & Tibshirani (2003) bootstrap estimator of pi0 (proportion of true nulls).

    The bootstrap variant is used rather than a fixed lambda because the p-value sets here
    are small (hundreds of discovery hits) and a fixed lambda=0.5 is unstable at that size.
    """
    p = np.asarray(p, dtype=float)
    p = p[np.isfinite(p)]
    if p.size < 20:
        return np.nan
    if lambdas is None:
        lambdas = np.arange(0.05, 0.96, 0.05)
    ps = np.sort(p)
    m = ps.size

    def _pi0(sorted_p, n):
        # count of p > lam via binary search on the sorted vector
        cnt = n - np.searchsorted(sorted_p, lambdas, side="right")
        return np.maximum(cnt / (n * (1 - lambdas)), 1e-8)

    pi0s = _pi0(ps, m)
    min_pi0 = np.quantile(pi0s, 0.1)
    # 100 bootstrap replicates; pick the lambda minimising estimated MSE against min_pi0.
    rng = np.random.default_rng(0)
    B = 100
    boots = np.sort(rng.choice(ps, size=(B, m), replace=True), axis=1)
    cnt = m - np.array([np.searchsorted(b, lambdas, side="right") for b in boots])
    pi0b = np.maximum(cnt / (m * (1 - lambdas)), 1e-8)
    mse = ((pi0b - min_pi0) ** 2).sum(axis=0)
    return float(np.clip(pi0s[np.argmin(mse)], 0.0, 1.0))


def storey_pi1(p) -> float:
    """pi1 = 1 - pi0: the estimated fraction of genuinely non-null tests."""
    pi0 = storey_pi0(np.asarray(p, dtype=float))
    return np.nan if not np.isfinite(pi0) else 1.0 - pi0


# --------------------------------------------------------------------------------------
# pairwise replication statistics
# --------------------------------------------------------------------------------------
def pair_stats(dis: pd.DataFrame, rep: pd.DataFrame, rules=RULES) -> dict:
    """All gene-level replication statistics for one ordered (discovery -> replication) pair.

    Parameters
    ----------
    dis, rep : DESeq2 results frames (index = gene symbol) for the discovery and
               replication studies.

    Returns a flat dict of metrics.  Keys prefixed `rep_<rule>` are replication rates under
    each selection rule; `n_dis_<rule>` are the corresponding discovery-hit counts.
    """
    common = dis.index.intersection(rep.index)
    d, r = dis.loc[common], rep.loc[common]
    out = {"n_common": len(common), "n_dis_tested": len(dis), "n_rep_tested": len(rep)}

    for rule in rules:
        sd, sr = select(d, rule), select(r, rule)
        n = len(sd)
        out[f"n_dis_{rule}"] = n
        out[f"n_rep_{rule}"] = len(sr)
        if n == 0:
            out[f"rep_{rule}"] = np.nan
            out[f"repdir_{rule}"] = np.nan
            continue
        inter = sd.intersection(sr)
        out[f"rep_{rule}"] = len(inter) / n
        # Same-direction replication: the stricter and more meaningful version.
        same = np.sign(d.loc[inter, "log2FoldChange"]) == np.sign(r.loc[inter, "log2FoldChange"])
        out[f"repdir_{rule}"] = float(same.sum()) / n

    # --- matched comparators for the pathway arm -------------------------------------
    sd = select(d, "joint")
    if len(sd) >= 5:
        # (a) sign concordance among discovery hits  <-> pathway NES sign concordance
        out["gene_sign_conc"] = float(
            (np.sign(d.loc[sd, "log2FoldChange"]) == np.sign(r.loc[sd, "log2FoldChange"])).mean())
        # (b) both-significant (FDR only, either direction) <-> pathway both-significant
        out["gene_bothsig_fdr"] = float((r.loc[sd, "padj"] < FDR).mean())
        # (c) both-significant AND same sign
        out["gene_bothsig_sign"] = float(
            ((r.loc[sd, "padj"] < FDR) &
             (np.sign(d.loc[sd, "log2FoldChange"]) == np.sign(r.loc[sd, "log2FoldChange"]))).mean())
        # (d) pi1 in the replication study restricted to discovery hits
        out["pi1"] = storey_pi1(r.loc[sd, "pvalue"].to_numpy())
        out["rep_median_abs_lfc"] = float(r.loc[sd, "log2FoldChange"].abs().median())
        out["dis_median_abs_lfc"] = float(d.loc[sd, "log2FoldChange"].abs().median())
        # Effect-size shrinkage on replication (Errington 2021 reported a median 85% drop)
        out["lfc_shrinkage"] = 1.0 - out["rep_median_abs_lfc"] / max(out["dis_median_abs_lfc"], 1e-9)
    else:
        for k in ("gene_sign_conc", "gene_bothsig_fdr", "gene_bothsig_sign", "pi1",
                  "rep_median_abs_lfc", "dis_median_abs_lfc", "lfc_shrinkage"):
            out[k] = np.nan

    # --- threshold-free agreement over all commonly tested genes ----------------------
    if len(common) > 100:
        out["spearman_lfc"] = float(stats.spearmanr(d.log2FoldChange, r.log2FoldChange).statistic)
        out["spearman_stat"] = float(stats.spearmanr(d.stat, r.stat).statistic)
        out["pearson_lfc"] = float(np.corrcoef(d.log2FoldChange, r.log2FoldChange)[0, 1])
        out["global_sign_conc"] = float((np.sign(d.log2FoldChange) == np.sign(r.log2FoldChange)).mean())
    else:
        for k in ("spearman_lfc", "spearman_stat", "pearson_lfc", "global_sign_conc"):
            out[k] = np.nan
    return out


# --------------------------------------------------------------------------------------
# inference helpers
# --------------------------------------------------------------------------------------
def _pool_weighted(values, type_codes, n_types, mult, agg):
    """Weighted pooled statistic under a bootstrap multiplicity vector.

    Vectorised with bincount rather than by materialising the resampled frame: the
    bootstrap runs thousands of iterations per metric, and the naive groupby version is
    two orders of magnitude slower for no difference in the answer.
    """
    if agg == "type_equal":
        num = np.bincount(type_codes, weights=mult * values, minlength=n_types)
        den = np.bincount(type_codes, weights=mult, minlength=n_types)
        ok = den > 0
        if not ok.any():
            return np.nan
        return float((num[ok] / den[ok]).mean())
    tot = mult.sum()
    return float((mult * values).sum() / tot) if tot > 0 else np.nan


def cluster_bootstrap(df: pd.DataFrame, value: str,
                      cluster_cols=("discovery", "replication"),
                      n_boot: int = 10000, seed: int = 42, agg: str = "type_equal",
                      type_col: str = "cancer_type"):
    """Bootstrap distribution of a pooled statistic, resampling **cohorts**, not pairs.

    Pair-level observations are not independent: with 13 LIHC cohorts, a single aberrant
    cohort contributes 24 of the 156 LIHC pairs.  Resampling cohorts with replacement and
    weighting each pair by the multiplicity of its members propagates that dependence.

    agg='type_equal' averages within cancer type before pooling, so that LIHC/BRCA/PRAD
    (30 of 51 cohorts) cannot silently determine the headline number.

    Returns (point_estimate, bootstrap_array).
    """
    rng = np.random.default_rng(seed)
    d = df.dropna(subset=[value])
    if d.empty:
        return np.nan, np.array([])

    values = d[value].to_numpy(dtype=float)
    tcat = pd.Categorical(d[type_col])
    type_codes = tcat.codes.astype(np.intp)
    n_types = len(tcat.categories)

    cols = list(dict.fromkeys(cluster_cols))
    cohort_cat = pd.Categorical(pd.concat([d[c] for c in cols]))
    n_coh = len(cohort_cat.categories)
    # per-pair cohort indices, so a resample maps to multiplicities by integer lookup
    idx = [cohort_cat.codes[i * len(d):(i + 1) * len(d)].astype(np.intp) for i in range(len(cols))]

    point = _pool_weighted(values, type_codes, n_types, np.ones(len(d)), agg)

    boots = np.empty(n_boot)
    for b in range(n_boot):
        draw = rng.integers(0, n_coh, n_coh)
        counts = np.bincount(draw, minlength=n_coh).astype(float)
        mult = counts[idx[0]]
        for j in range(1, len(idx)):
            mult = mult * counts[idx[j]]
        boots[b] = _pool_weighted(values, type_codes, n_types, mult, agg)
    boots = boots[np.isfinite(boots)]
    return point, boots


def summarise(df, value, **kw):
    """Point estimate plus percentile bootstrap 95% CI."""
    point, boots = cluster_bootstrap(df, value, **kw)
    if boots.size == 0:
        return dict(point=point, lo=np.nan, hi=np.nan, n=int(df[value].notna().sum()))
    return dict(point=point, lo=float(np.percentile(boots, 2.5)),
                hi=float(np.percentile(boots, 97.5)),
                n=int(df[value].notna().sum()))


def bootstrap_p(df, value, threshold, direction, **kw):
    """One-sided percentile-bootstrap p-value for pooled(value) < / > threshold."""
    _, boots = cluster_bootstrap(df, value, **kw)
    if boots.size == 0:
        return np.nan
    p = np.mean(boots >= threshold) if direction == "less" else np.mean(boots <= threshold)
    return float(max(p, 1.0 / boots.size))


def paired_test(a, b):
    """Wilcoxon signed-rank on paired observations, with the rank-biserial effect size."""
    a = np.asarray(a, float); b = np.asarray(b, float)
    m = np.isfinite(a) & np.isfinite(b)
    a, b = a[m], b[m]
    d = a - b
    d = d[d != 0]
    if d.size < 6:
        return dict(n=int(m.sum()), median_diff=float(np.median(a - b)) if m.sum() else np.nan,
                    p=np.nan, rbc=np.nan)
    w = stats.wilcoxon(a, b, zero_method="wilcox")
    r = stats.rankdata(np.abs(d))
    rbc = (r[d > 0].sum() - r[d < 0].sum()) / r.sum()
    return dict(n=int(m.sum()), median_diff=float(np.median(a - b)),
                p=float(w.pvalue), rbc=float(rbc))
