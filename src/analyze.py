#!/usr/bin/env python3
"""Phase 3 analysis: hypothesis tests, tables, and figures.

Reads the pairwise metric tables produced by `run_metrics.py` and `run_pathway_metrics.py`
and writes:

  results/tables/*.tsv     every table quoted in REPORT.md
  results/summary.json     machine-readable headline numbers
  figures/*.png            all figures

Inference throughout uses a **cohort-level cluster bootstrap** with cancer types weighted
equally (see `metrics.cluster_bootstrap`).  Wilcoxon signed-rank tests are reported next to
it as a secondary check and are explicitly labelled anticonservative, because pair-level
observations share cohorts.
"""
from __future__ import annotations

import os
import sys
import json

import numpy as np
import pandas as pd
from scipy import stats as sps

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import design
import metrics as M

ROOT = design.ROOT
MET = os.path.join(ROOT, "results", "metrics")
TAB = os.path.join(ROOT, "results", "tables")
FIG = os.path.join(ROOT, "figures")
N_BOOT = 5000

ARM_LABEL = {
    "rsplit": "Random split\n(same cohort)",
    "tsplit": "Site-disjoint split\n(same cohort)",
    "cross_m10": "Independent studies\n(matched 10v10)",
    "cross_full": "Independent studies\n(native n)",
    "shuffle": "Shuffled labels\n(negative control)",
}
ARM_ORDER = ["rsplit", "tsplit", "cross_m10", "cross_full", "shuffle"]
ARM_COLOR = {"rsplit": "#3b7dd8", "tsplit": "#5fbf8f", "cross_m10": "#d9822b",
             "cross_full": "#c0504d", "shuffle": "#9b9b9b"}

# Palette for gene / pathway / random-set comparisons
UNIT_COLOR = {"gene": "#d9822b", "hallmark": "#3b7dd8", "random": "#9b9b9b",
              "kegg": "#5fbf8f", "reactome": "#8e6bbf", "gobp": "#c0504d"}


# ======================================================================================
# loading
# ======================================================================================
def load_all(qc_filter=True):
    qc = pd.read_csv(os.path.join(MET, "cohort_qc.tsv"), sep="\t")
    good = set(qc.cohort[qc.qc_pass])

    genes = pd.concat([
        pd.read_csv(os.path.join(MET, "pairs_full.tsv"), sep="\t"),
        pd.read_csv(os.path.join(MET, "pairs_m10.tsv"), sep="\t"),
        pd.read_csv(os.path.join(MET, "splits.tsv"), sep="\t"),
    ], ignore_index=True)
    paths = pd.read_csv(os.path.join(MET, "pathways.tsv"), sep="\t")

    for df in (genes, paths):
        df["ok"] = df.status.astype(str) == "ok"
    genes = genes[genes.ok].copy()
    paths = paths[paths.ok].copy()

    if qc_filter:
        for name, df in (("genes", genes), ("paths", paths)):
            keep = df.discovery.isin(good) & df.replication.isin(good)
            print(f"  QC filter drops {int((~keep).sum())}/{len(df)} {name} rows")
        genes = genes[genes.discovery.isin(good) & genes.replication.isin(good)].copy()
        paths = paths[paths.discovery.isin(good) & paths.replication.isin(good)].copy()
    return genes, paths, qc


def clusters_for(arm):
    """Split arms compare a cohort with itself, so the bootstrap clusters on one column."""
    return ("discovery",) if arm in ("rsplit", "tsplit", "shuffle") else ("discovery", "replication")


def summarise_arm(df, arm, value, n_boot=N_BOOT):
    sub = df[df.arm == arm]
    if sub.empty or sub[value].notna().sum() == 0:
        return dict(arm=arm, metric=value, point=np.nan, lo=np.nan, hi=np.nan, n=0)
    s = M.summarise(sub, value, cluster_cols=clusters_for(arm), n_boot=n_boot)
    s.update(arm=arm, metric=value)
    return s


# ======================================================================================
# Table 1 - headline replication by arm
# ======================================================================================
def table_arms(genes, out="T1_replication_by_arm.tsv"):
    rows = []
    for arm in ARM_ORDER:
        for value in ["rep_joint", "repdir_joint", "rep_fdr", "gene_sign_conc",
                      "gene_bothsig_fdr", "pi1", "spearman_lfc", "spearman_stat",
                      "lfc_shrinkage", "global_sign_conc"]:
            rows.append(summarise_arm(genes, arm, value))
    t = pd.DataFrame(rows)
    t["n_pairs"] = t["n"]
    t.to_csv(os.path.join(TAB, out), sep="\t", index=False)
    return t


# ======================================================================================
# Table 2 - per cancer type
# ======================================================================================
def table_by_type(genes, out="T2_replication_by_cancer_type.tsv"):
    rows = []
    for arm in ("cross_m10", "cross_full"):
        sub = genes[genes.arm == arm]
        for ct, g in sub.groupby("cancer_type"):
            rows.append(dict(arm=arm, cancer_type=ct, n_pairs=len(g),
                             n_cohorts=g.discovery.nunique(),
                             rep_joint=g.rep_joint.mean(),
                             rep_joint_sd=g.rep_joint.std(),
                             gene_sign_conc=g.gene_sign_conc.mean(),
                             pi1=g.pi1.mean(),
                             spearman_lfc=g.spearman_lfc.mean(),
                             median_n_dis_joint=g.n_dis_joint.median()))
    t = pd.DataFrame(rows).sort_values(["arm", "cancer_type"])
    t.to_csv(os.path.join(TAB, out), sep="\t", index=False)
    return t


# ======================================================================================
# Table 3 - selection-rule sweep
# ======================================================================================
def table_rules(genes, out="T3_selection_rules.tsv"):
    rows = []
    for arm in ARM_ORDER:
        for rule in M.RULES:
            s = summarise_arm(genes, arm, f"rep_{rule}", n_boot=2000)
            s["rule"] = rule
            sub = genes[genes.arm == arm]
            s["median_list_size"] = float(sub[f"n_dis_{rule}"].median()) if len(sub) else np.nan
            rows.append(s)
    t = pd.DataFrame(rows)
    t.to_csv(os.path.join(TAB, out), sep="\t", index=False)
    return t


# ======================================================================================
# Table 4 - gene vs pathway vs random-set, matched comparators (H3)
# ======================================================================================
def table_units(genes, paths, out="T4_gene_vs_pathway.tsv"):
    """Pooled estimates of the three matched comparators for genes, real gene sets, and
    size-matched random gene sets."""
    rows = []
    keys = ["arm", "cancer_type", "discovery", "replication", "seed"]
    for arm in ARM_ORDER:
        g = genes[genes.arm == arm]
        if len(g):
            for gm, name in (("gene_sign_conc", "sign_conc"),
                             ("gene_bothsig_fdr", "bothsig"),
                             ("gene_bothsig_sign", "bothsig_sign")):
                s = summarise_arm(genes, arm, gm)
                rows.append(dict(arm=arm, unit="gene", comparator=name,
                                 point=s["point"], lo=s["lo"], hi=s["hi"], n=s["n"]))
        for coll in ("hallmark", "random", "kegg", "reactome", "gobp"):
            p = paths[(paths.arm == arm) & (paths.collection == coll)]
            if p.empty:
                continue
            for pm in ("sign_conc", "bothsig", "bothsig_sign"):
                s = M.summarise(p, pm, cluster_cols=clusters_for(arm), n_boot=N_BOOT)
                rows.append(dict(arm=arm, unit=coll, comparator=pm,
                                 point=s["point"], lo=s["lo"], hi=s["hi"], n=s["n"]))
    t = pd.DataFrame(rows)
    t.to_csv(os.path.join(TAB, out), sep="\t", index=False)
    return t


def paired_frame(genes, paths, arm, coll="hallmark"):
    """Merge gene- and pathway-level statistics on the pair key for paired testing."""
    keys = ["arm", "cancer_type", "discovery", "replication", "seed"]
    g = genes[genes.arm == arm][keys + ["gene_sign_conc", "gene_bothsig_fdr",
                                        "gene_bothsig_sign", "rep_joint"]]
    p = paths[(paths.arm == arm) & (paths.collection == coll)][
        keys + ["sign_conc", "bothsig", "bothsig_sign", "n_dis_sig"]]
    if arm in ("rsplit", "tsplit", "shuffle"):
        # split arms have two directions per rep; include direction in the key
        g = genes[genes.arm == arm][keys + ["direction", "gene_sign_conc", "gene_bothsig_fdr",
                                            "gene_bothsig_sign", "rep_joint"]]
        p = paths[(paths.arm == arm) & (paths.collection == coll)][
            keys + ["direction", "sign_conc", "bothsig", "bothsig_sign", "n_dis_sig"]]
        keys = keys + ["direction"]
    return g.merge(p, on=keys, how="inner", suffixes=("", "_p"))


def table_h3_tests(genes, paths, out="T5_H3_paired_tests.tsv"):
    """Paired tests for H3: is the pathway advantage real once comparators are matched and
    a size-matched random-set null is in place?"""
    rows = []
    for arm in ("cross_m10", "cross_full"):
        hp = paired_frame(genes, paths, arm, "hallmark")
        rp = paired_frame(genes, paths, arm, "random")
        pairs = [
            ("hallmark_vs_gene", "sign_conc", hp, "sign_conc", "gene_sign_conc"),
            ("hallmark_vs_gene", "bothsig", hp, "bothsig", "gene_bothsig_fdr"),
            ("hallmark_vs_gene", "bothsig_sign", hp, "bothsig_sign", "gene_bothsig_sign"),
            ("random_vs_gene", "sign_conc", rp, "sign_conc", "gene_sign_conc"),
            ("random_vs_gene", "bothsig", rp, "bothsig", "gene_bothsig_fdr"),
        ]
        for name, comp, df, a, b in pairs:
            if df.empty:
                continue
            r = M.paired_test(df[a], df[b])
            r.update(arm=arm, contrast=name, comparator=comp)
            rows.append(r)
        # hallmark vs random gene sets, merged on the pair key
        keyc = [c for c in ("arm", "cancer_type", "discovery", "replication", "seed", "direction")
                if c in hp.columns and c in rp.columns]
        m = hp.merge(rp, on=keyc, suffixes=("_h", "_r"))
        for comp in ("sign_conc", "bothsig", "bothsig_sign"):
            if f"{comp}_h" in m and f"{comp}_r" in m:
                r = M.paired_test(m[f"{comp}_h"], m[f"{comp}_r"])
                r.update(arm=arm, contrast="hallmark_vs_random", comparator=comp)
                rows.append(r)
    t = pd.DataFrame(rows)
    if len(t):
        # Bonferroni over the tests within each arm
        t["p_bonf"] = np.minimum(t.groupby("arm")["p"].transform(lambda s: s * s.notna().sum()), 1.0)
    t.to_csv(os.path.join(TAB, out), sep="\t", index=False)
    return t


# ======================================================================================
# Table 6 - what drives the replication rate?
# ======================================================================================
def table_variance(genes, out="T6_variance_decomposition.tsv"):
    """Attribute variation in the replication rate to the discovery vs the replication study.

    A one-way ANOVA-style R^2 on cohort identity.  If replication failure were a property of
    the *discovery finding* (a false positive), the discovery cohort would dominate.  If it
    is a property of the *replicating study's sensitivity*, the replication cohort dominates.
    """
    import statsmodels.formula.api as smf
    rows = []
    for arm in ("cross_m10", "cross_full"):
        d = genes[genes.arm == arm].dropna(subset=["rep_joint"])
        if len(d) < 20:
            continue
        r2 = {}
        for name, formula in (("discovery_only", "rep_joint ~ C(discovery)"),
                              ("replication_only", "rep_joint ~ C(replication)"),
                              ("both", "rep_joint ~ C(discovery)+C(replication)")):
            r2[name] = float(smf.ols(formula, d).fit().rsquared)
        rows.append(dict(arm=arm, n_pairs=len(d), **r2,
                         rho_rep_with_replication_hits=float(
                             d[["rep_joint", "n_rep_joint"]].corr(method="spearman").iloc[0, 1]),
                         rho_rep_with_discovery_hits=float(
                             d[["rep_joint", "n_dis_joint"]].corr(method="spearman").iloc[0, 1]),
                         mean_pi1=float(d.pi1.mean()),
                         mean_rep=float(d.rep_joint.mean())))
    t = pd.DataFrame(rows)
    t.to_csv(os.path.join(TAB, out), sep="\t", index=False)
    return t


# ======================================================================================
# Hypothesis tests
# ======================================================================================
def hypothesis_tests(genes, paths):
    res = {}
    for arm in ("cross_m10", "cross_full"):
        g = genes[genes.arm == arm]
        h = paths[(paths.arm == arm) & (paths.collection == "hallmark")]
        cc = clusters_for(arm)
        res[f"H1_{arm}"] = dict(
            **M.summarise(g, "rep_joint", cluster_cols=cc, n_boot=N_BOOT),
            threshold=0.40, direction="less",
            p_bootstrap=M.bootstrap_p(g, "rep_joint", 0.40, "less",
                                      cluster_cols=cc, n_boot=N_BOOT),
            p_wilcoxon_secondary=(float(sps.wilcoxon(g.rep_joint.dropna() - 0.40,
                                                     alternative="less").pvalue)
                                  if g.rep_joint.notna().sum() > 10 else np.nan))
        res[f"H2_{arm}"] = dict(
            **M.summarise(h, "sign_conc", cluster_cols=cc, n_boot=N_BOOT),
            threshold=0.70, direction="greater",
            p_bootstrap=M.bootstrap_p(h, "sign_conc", 0.70, "greater",
                                      cluster_cols=cc, n_boot=N_BOOT),
            p_wilcoxon_secondary=(float(sps.wilcoxon(h.sign_conc.dropna() - 0.70,
                                                     alternative="greater").pvalue)
                                  if h.sign_conc.notna().sum() > 10 else np.nan))
    return res


# ======================================================================================
# Figures
# ======================================================================================
def _style(ax):
    ax.spines[["top", "right"]].set_visible(False)
    ax.grid(axis="y", color="#e6e6e6", lw=0.8)
    ax.set_axisbelow(True)


def fig_arms(genes, t1):
    """F1: replication rate and sign concordance across the four designs + control."""
    fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.6))
    for ax, value, title, thr in (
            (axes[0], "rep_joint", "Gene-level replication rate\n(FDR<0.05 & |log2FC|>1)", 0.40),
            (axes[1], "gene_sign_conc", "log2FC sign concordance\namong discovery hits", 0.70)):
        xs, labels, colors = [], [], []
        for i, arm in enumerate(ARM_ORDER):
            v = genes.loc[genes.arm == arm, value].dropna()
            if v.empty:
                continue
            xs.append(v); labels.append(ARM_LABEL[arm]); colors.append(ARM_COLOR[arm])
        bp = ax.boxplot(xs, patch_artist=True, widths=0.6, showfliers=False,
                        medianprops=dict(color="black", lw=1.6))
        for patch, c in zip(bp["boxes"], colors):
            patch.set_facecolor(c); patch.set_alpha(0.55); patch.set_edgecolor(c)
        for i, v in enumerate(xs, 1):
            j = np.random.default_rng(0).normal(0, 0.055, len(v))
            ax.scatter(i + j, v, s=5, color="#333333", alpha=0.25, linewidths=0, zorder=3)
        ax.axhline(thr, ls="--", color="#c0504d", lw=1.2)
        ax.text(0.02, thr + 0.015, f"pre-registered {thr:.0%}", color="#c0504d",
                fontsize=8, transform=ax.get_yaxis_transform())
        ax.set_xticks(range(1, len(labels) + 1)); ax.set_xticklabels(labels, fontsize=8)
        ax.set_ylim(0, 1.02); ax.set_ylabel("proportion"); ax.set_title(title, fontsize=10)
        _style(ax)
    fig.suptitle("Replication of bulk RNA-seq differential expression, by study design",
                 fontsize=12, y=0.99)
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "F1_replication_by_design.png"), dpi=180)
    plt.close(fig)


def fig_by_type(genes):
    """F2: per-cancer-type replication, matched 10v10 arm."""
    sub = genes[genes.arm == "cross_m10"]
    if sub.empty:
        return
    order = sub.groupby("cancer_type").rep_joint.mean().sort_values().index
    fig, ax = plt.subplots(figsize=(9, 4.4))
    data = [sub.loc[sub.cancer_type == ct, "rep_joint"].dropna() for ct in order]
    bp = ax.boxplot(data, patch_artist=True, widths=0.62, showfliers=False,
                    medianprops=dict(color="black", lw=1.5))
    for p in bp["boxes"]:
        p.set_facecolor("#d9822b"); p.set_alpha(0.5); p.set_edgecolor("#d9822b")
    for i, v in enumerate(data, 1):
        j = np.random.default_rng(1).normal(0, 0.06, len(v))
        ax.scatter(i + j, v, s=5, color="#333", alpha=0.3, linewidths=0, zorder=3)
    ax.axhline(0.40, ls="--", color="#c0504d", lw=1.2)
    n = sub.groupby("cancer_type").size()
    ax.set_xticks(range(1, len(order) + 1))
    ax.set_xticklabels([f"{ct}\n({n[ct]})" for ct in order], fontsize=8)
    ax.set_ylabel("replication rate"); ax.set_ylim(0, 1.02)
    ax.set_title("Gene-level replication by cancer type (independent studies, matched 10v10)\n"
                 "numbers in parentheses = ordered cohort pairs x subsampling seeds", fontsize=10)
    _style(ax)
    fig.tight_layout(); fig.savefig(os.path.join(FIG, "F2_replication_by_cancer_type.png"), dpi=180)
    plt.close(fig)


def fig_rules(t3):
    """F3: how the gene-selection rule moves the replication rate (Shi et al. 2008)."""
    fig, ax = plt.subplots(figsize=(9.5, 4.4))
    rules = M.RULES
    labels = {"joint": "FDR<0.05 &\n|log2FC|>1", "fdr": "FDR<0.05", "lfc": "|log2FC|>1",
              "topp100": "top 100\nby p", "topp500": "top 500\nby p",
              "topfc100": "top 100\nby |FC|", "topfc500": "top 500\nby |FC|"}
    w = 0.16
    for k, arm in enumerate(["rsplit", "tsplit", "cross_m10", "cross_full", "shuffle"]):
        s = t3[(t3.arm == arm)].set_index("rule").reindex(rules)
        if s["point"].isna().all():
            continue
        x = np.arange(len(rules)) + (k - 2) * w
        ax.bar(x, s["point"], width=w, color=ARM_COLOR[arm], alpha=0.85,
               label=ARM_LABEL[arm].replace("\n", " "))
        ax.errorbar(x, s["point"], yerr=[s["point"] - s["lo"], s["hi"] - s["point"]],
                    fmt="none", ecolor="#333", elinewidth=0.8, capsize=2)
    ax.axhline(0.40, ls="--", color="#c0504d", lw=1.2)
    ax.set_xticks(np.arange(len(rules))); ax.set_xticklabels([labels[r] for r in rules], fontsize=8)
    ax.set_ylabel("replication rate"); ax.set_ylim(0, 1.02)
    ax.legend(fontsize=8, frameon=False, ncol=3)
    ax.set_title("The gene-selection rule moves the replication rate as much as the study design does",
                 fontsize=10)
    _style(ax)
    fig.tight_layout(); fig.savefig(os.path.join(FIG, "F3_selection_rules.png"), dpi=180)
    plt.close(fig)


def fig_units(t4):
    """F4: the H3 test - matched comparators for genes, hallmark sets, and random sets."""
    fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.4), sharey=True)
    for ax, arm in zip(axes, ("cross_m10", "cross_full")):
        sub = t4[(t4.arm == arm)]
        comps = ["sign_conc", "bothsig", "bothsig_sign"]
        units = ["gene", "hallmark", "random"]
        w = 0.26
        for k, u in enumerate(units):
            s = sub[sub.unit == u].set_index("comparator").reindex(comps)
            x = np.arange(len(comps)) + (k - 1) * w
            ax.bar(x, s["point"], width=w, color=UNIT_COLOR[u], alpha=0.9,
                   label={"gene": "genes", "hallmark": "hallmark gene sets",
                          "random": "size-matched random sets"}[u])
            ax.errorbar(x, s["point"], yerr=[s["point"] - s["lo"], s["hi"] - s["point"]],
                        fmt="none", ecolor="#333", elinewidth=0.8, capsize=2)
        ax.axhline(0.70, ls="--", color="#c0504d", lw=1.2)
        ax.set_xticks(np.arange(len(comps)))
        ax.set_xticklabels(["same direction", "both significant",
                            "both significant\n& same direction"], fontsize=8)
        ax.set_ylim(0, 1.05)
        ax.set_title({"cross_m10": "Independent studies, matched 10v10",
                      "cross_full": "Independent studies, native n"}[arm], fontsize=10)
        _style(ax)
    axes[0].set_ylabel("proportion of discovery-significant units")
    axes[0].legend(fontsize=8, frameon=False, loc="lower left")
    fig.suptitle("Genes vs pathways vs random gene sets, on matched comparators", fontsize=12)
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "F4_gene_vs_pathway_vs_random.png"), dpi=180)
    plt.close(fig)


def fig_decomposition(t1):
    """F5: the sampling-vs-heterogeneity decomposition, all arms at n=10v10."""
    arms = ["shuffle", "cross_m10", "tsplit", "rsplit"]
    fig, ax = plt.subplots(figsize=(8.2, 4.4))
    s = t1[t1.metric == "rep_joint"].set_index("arm").reindex(arms)
    y = np.arange(len(arms))
    ax.barh(y, s["point"], color=[ARM_COLOR[a] for a in arms], alpha=0.85, height=0.6)
    ax.errorbar(s["point"], y, xerr=[s["point"] - s["lo"], s["hi"] - s["point"]],
                fmt="none", ecolor="#333", elinewidth=0.9, capsize=3)
    for yi, (p, n) in enumerate(zip(s["point"], s["n"])):
        if np.isfinite(p):
            ax.text(p + 0.015, yi, f"{p:.2f}  (n={int(n)})", va="center", fontsize=9)
    ax.set_yticks(y); ax.set_yticklabels([ARM_LABEL[a].replace("\n", " ") for a in arms], fontsize=9)
    ax.axvline(0.40, ls="--", color="#c0504d", lw=1.2)
    ax.set_xlim(0, 1.0); ax.set_xlabel("gene-level replication rate (FDR<0.05 & |log2FC|>1)")
    ax.set_title("Decomposition at identical sample size (10 tumour vs 10 normal)\n"
                 "gap between random split and independent studies = between-study heterogeneity",
                 fontsize=10)
    ax.spines[["top", "right"]].set_visible(False)
    ax.grid(axis="x", color="#e6e6e6", lw=0.8); ax.set_axisbelow(True)
    fig.tight_layout(); fig.savefig(os.path.join(FIG, "F5_decomposition.png"), dpi=180)
    plt.close(fig)


def fig_pi1(genes):
    """F6: replication rate vs pi1 - how much of the shortfall is sub-threshold signal."""
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.3))
    sub = genes[genes.arm.isin(["cross_m10", "cross_full"])]
    ax = axes[0]
    for arm in ("cross_m10", "cross_full"):
        s = sub[sub.arm == arm]
        ax.scatter(s.rep_joint, s.pi1, s=8, alpha=0.4, color=ARM_COLOR[arm],
                   label=ARM_LABEL[arm].replace("\n", " "), linewidths=0)
    ax.plot([0, 1], [0, 1], ls=":", color="#888", lw=1)
    ax.set_xlabel("replication rate (thresholded)"); ax.set_ylabel(r"$\pi_1$ in replication study")
    ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.legend(fontsize=8, frameon=False, loc="lower right")
    ax.set_title(r"Most discovery hits are non-null in the replication study" "\n"
                 r"even when they miss the threshold ($\pi_1 \gg$ replication rate)", fontsize=10)
    _style(ax)

    ax = axes[1]
    for arm in ("cross_m10", "cross_full"):
        s = sub[sub.arm == arm]
        ax.scatter(s.spearman_lfc, s.rep_joint, s=8, alpha=0.4, color=ARM_COLOR[arm],
                   linewidths=0, label=ARM_LABEL[arm].replace("\n", " "))
    ax.set_xlabel("Spearman correlation of log2FC (all common genes)")
    ax.set_ylabel("replication rate")
    ax.set_ylim(0, 1); ax.legend(fontsize=8, frameon=False, loc="upper left")
    ax.set_title("Thresholded replication vs threshold-free agreement", fontsize=10)
    _style(ax)
    fig.tight_layout(); fig.savefig(os.path.join(FIG, "F6_pi1_and_correlation.png"), dpi=180)
    plt.close(fig)


def fig_collections(t4):
    """F7: sensitivity of the pathway result to the gene-set collection."""
    sub = t4[(t4.arm == "cross_m10") & (t4.comparator.isin(["sign_conc", "bothsig"]))]
    units = [u for u in ["hallmark", "kegg", "reactome", "gobp", "random"]
             if u in set(sub.unit)]
    if not units:
        return
    fig, ax = plt.subplots(figsize=(8.5, 4.2))
    w = 0.38
    for k, comp in enumerate(["sign_conc", "bothsig"]):
        s = sub[sub.comparator == comp].set_index("unit").reindex(units)
        x = np.arange(len(units)) + (k - 0.5) * w
        ax.bar(x, s["point"], width=w, alpha=0.9,
               color="#3b7dd8" if comp == "sign_conc" else "#8e6bbf",
               label="same NES direction" if comp == "sign_conc" else "both significant")
        ax.errorbar(x, s["point"], yerr=[s["point"] - s["lo"], s["hi"] - s["point"]],
                    fmt="none", ecolor="#333", elinewidth=0.8, capsize=2)
    ax.axhline(0.70, ls="--", color="#c0504d", lw=1.2)
    ax.set_xticks(np.arange(len(units)))
    ax.set_xticklabels([{"hallmark": "hallmark\n(50)", "kegg": "KEGG\n(186)",
                         "reactome": "Reactome\n(1692)", "gobp": "GO:BP\n(7647)",
                         "random": "random sets\n(200, size-matched)"}.get(u, u) for u in units],
                       fontsize=8)
    ax.set_ylabel("proportion"); ax.set_ylim(0, 1.05)
    ax.legend(fontsize=8, frameon=False)
    ax.set_title("Pathway replication by gene-set collection (independent studies, matched 10v10)",
                 fontsize=10)
    _style(ax)
    fig.tight_layout(); fig.savefig(os.path.join(FIG, "F7_collections.png"), dpi=180)
    plt.close(fig)


def fig_size_effect(genes):
    """F8: replication rate against the harmonic-mean sample size of the pair (native n)."""
    sub = genes[genes.arm == "cross_full"].copy()
    if sub.empty:
        return
    cohorts = design.cohort_table().set_index("cohort")
    n = (cohorts.n_tumor + cohorts.n_normal)
    sub["n_dis"] = sub.discovery.map(n); sub["n_rep"] = sub.replication.map(n)
    sub["n_harm"] = 2 / (1 / sub.n_dis + 1 / sub.n_rep)
    fig, ax = plt.subplots(figsize=(7.2, 4.3))
    sc = ax.scatter(sub.n_harm, sub.rep_joint, c=sub.n_dis_joint, s=14, alpha=0.7,
                    cmap="viridis", norm=matplotlib.colors.LogNorm(), linewidths=0)
    plt.colorbar(sc, ax=ax, label="discovery DE genes")
    ax.set_xscale("log"); ax.set_xlabel("harmonic-mean cohort size (samples)")
    ax.set_ylabel("replication rate"); ax.set_ylim(0, 1)
    ax.axhline(0.40, ls="--", color="#c0504d", lw=1.2)
    m = np.isfinite(sub.n_harm) & np.isfinite(sub.rep_joint)
    rho = sps.spearmanr(sub.n_harm[m], sub.rep_joint[m])
    ax.set_title(f"Replication rises with cohort size (native n)\n"
                 f"Spearman rho = {rho.statistic:.2f}", fontsize=10)
    _style(ax)
    fig.tight_layout(); fig.savefig(os.path.join(FIG, "F8_size_effect.png"), dpi=180)
    plt.close(fig)


# ======================================================================================
def main():
    os.makedirs(TAB, exist_ok=True); os.makedirs(FIG, exist_ok=True)
    print("loading metric tables ...")
    genes, paths, qc = load_all()
    print(f"  gene rows {len(genes)}  pathway rows {len(paths)}")
    print(genes.groupby("arm").size())

    print("T1 replication by arm ...");        t1 = table_arms(genes)
    print("T2 by cancer type ...");            t2 = table_by_type(genes)
    print("T3 selection rules ...");           t3 = table_rules(genes)
    print("T4 gene vs pathway vs random ...");  t4 = table_units(genes, paths)
    print("T5 H3 paired tests ...");           t5 = table_h3_tests(genes, paths)
    print("T6 variance decomposition ...");    t6 = table_variance(genes)
    print("hypothesis tests ...");             ht = hypothesis_tests(genes, paths)

    # sensitivity: no QC filtering
    genes_all, paths_all, _ = load_all(qc_filter=False)
    sens = {arm: M.summarise(genes_all[genes_all.arm == arm], "rep_joint",
                             cluster_cols=clusters_for(arm), n_boot=2000)
            for arm in ("cross_m10", "cross_full")}

    print("figures ...")
    fig_arms(genes, t1); fig_by_type(genes); fig_rules(t3); fig_units(t4)
    fig_decomposition(t1); fig_pi1(genes); fig_collections(t4); fig_size_effect(genes)

    summary = dict(
        n_cohorts_total=int(len(qc)), n_cohorts_qc_pass=int(qc.qc_pass.sum()),
        qc_failures=list(qc.cohort[~qc.qc_pass]),
        n_gene_comparisons=int(len(genes)), n_pathway_comparisons=int(len(paths)),
        arm_counts={k: int(v) for k, v in genes.groupby("arm").size().items()},
        hypothesis_tests=ht,
        sensitivity_no_qc_filter=sens,
        headline={r.arm + "|" + r.metric: dict(point=r.point, lo=r.lo, hi=r.hi, n=int(r.n))
                  for r in t1.itertuples()},
        variance_decomposition=t6.to_dict(orient="records"),
    )
    with open(os.path.join(ROOT, "results", "summary.json"), "w") as fh:
        json.dump(summary, fh, indent=2, default=float)

    pd.set_option("display.width", 200)
    print("\n=== T1 (rep_joint / gene_sign_conc / pi1) ===")
    print(t1[t1.metric.isin(["rep_joint", "gene_sign_conc", "pi1", "spearman_lfc"])]
          [["arm", "metric", "point", "lo", "hi", "n"]].to_string(index=False))
    print("\n=== hypothesis tests ===")
    print(json.dumps(ht, indent=2, default=float))
    print("\n=== T4 (matched comparators) ===")
    print(t4[t4.arm.isin(["cross_m10", "cross_full"])].to_string(index=False))
    print("\n=== T5 (paired H3 tests) ===")
    print(t5.to_string(index=False))
    print("\n=== T6 (variance decomposition) ===")
    print(t6.to_string(index=False))
    print("\nwrote results/tables/*.tsv, results/summary.json, figures/*.png")


if __name__ == "__main__":
    main()
