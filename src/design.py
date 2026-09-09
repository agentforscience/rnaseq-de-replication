#!/usr/bin/env python3
"""Study design: which cohorts enter the analysis, which pairs, which control splits.

Kept separate from execution so that the design is auditable in one place and identical
across the full-size arm, the matched-n arm, and the within-cohort control arms.
"""
from __future__ import annotations

import os
import sys

import numpy as np
import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "code"))

# --- pre-registered constants -------------------------------------------------------
MIN_PER_GROUP = 8          # a cohort must have >= this many tumour and normal samples
MATCHED_N = 10             # sample size used for every power-matched comparison
N_SUBSAMPLE_SEEDS = 5      # independent 10v10 draws per cohort in the matched arm
N_SPLIT_REPS = 15          # split-pairs per TCGA cohort for each control type
MIN_NORMALS_FOR_SPLIT = 25 # a cohort needs this many normals to give two disjoint 10-sets
FDR = 0.05
LFC = 1.0
SEED = 42


def cohort_table() -> pd.DataFrame:
    """Cohorts eligible for the analysis: >=2 cohorts of the same cancer type, and at least
    MIN_PER_GROUP tumour and normal samples each."""
    c = pd.read_csv(os.path.join(ROOT, "datasets/manifest/cohorts.tsv"), sep="\t")
    c = c[(c.n_tumor >= MIN_PER_GROUP) & (c.n_normal >= MIN_PER_GROUP)].copy()
    counts = c.cancer_type.value_counts()
    c = c[c.cancer_type.isin(counts[counts >= 2].index)].copy()
    return c.sort_values(["cancer_type", "n_tumor"], ascending=[True, False]).reset_index(drop=True)


def cohort_pairs(cohorts: pd.DataFrame) -> pd.DataFrame:
    """All ordered (discovery, replication) pairs of distinct cohorts within a cancer type."""
    rows = []
    for ct, g in cohorts.groupby("cancer_type"):
        names = list(g.cohort)
        for a in names:
            for b in names:
                if a != b:
                    rows.append(dict(cancer_type=ct, discovery=a, replication=b))
    return pd.DataFrame(rows)


def split_cohorts(cohorts: pd.DataFrame) -> list:
    """TCGA cohorts large enough to be split into two disjoint 10v10 halves."""
    c = cohorts[(cohorts.source == "tcga") &
                (cohorts.n_normal >= MIN_NORMALS_FOR_SPLIT) &
                (cohorts.n_tumor >= 2 * MATCHED_N)]
    return list(c.cohort)


def _samples() -> pd.DataFrame:
    return pd.read_csv(os.path.join(ROOT, "datasets/manifest/samples.tsv"), sep="\t", dtype=str)


def random_split_runs(cohort: str, rep: int, n: int = MATCHED_N):
    """Two disjoint (tumour n, normal n) sample sets drawn at random from one cohort.

    This is the sampling-noise ceiling: the two 'studies' differ only by which patients
    happened to be drawn.  Returns (runs_A, runs_B).
    """
    man = _samples()
    man = man[man.cohort == cohort]
    rng = np.random.default_rng(SEED + 1000 * rep + hash(cohort) % 997)
    a, b = [], []
    for cond in ("tumor", "normal"):
        pool = man.run_id[man.condition == cond].to_numpy()
        pick = rng.choice(pool, 2 * n, replace=False)
        a += list(pick[:n]); b += list(pick[n:])
    return tuple(a), tuple(b)


def tss_split_runs(cohort: str, rep: int, n: int = MATCHED_N):
    """Two disjoint (tumour n, normal n) sets drawn from *disjoint tissue-source sites*.

    TCGA barcodes carry a tissue-source-site code identifying the contributing institution.
    Partitioning sites (rather than patients) adds collection centre, local protocol,
    population and batch differences while holding cancer type, assay and processing
    pipeline exactly fixed.  Returns (runs_A, runs_B) or None if no feasible partition is
    found within the retry budget.
    """
    man = _samples()
    man = man[(man.cohort == cohort) & man.tss.notna()]
    if man.empty:
        return None
    sites = sorted(man.tss.unique())
    rng = np.random.default_rng(SEED + 7919 * rep + hash(cohort) % 997)

    for _ in range(400):
        perm = rng.permutation(sites)
        cut = rng.integers(1, len(perm))
        s1, s2 = set(perm[:cut]), set(perm[cut:])
        g1 = man[man.tss.isin(s1)]
        g2 = man[man.tss.isin(s2)]
        ok = all((g.condition == c).sum() >= n for g in (g1, g2) for c in ("tumor", "normal"))
        if not ok:
            continue
        out = []
        for g in (g1, g2):
            runs = []
            for cond in ("tumor", "normal"):
                pool = g.run_id[g.condition == cond].to_numpy()
                runs += list(rng.choice(pool, n, replace=False))
            out.append(tuple(runs))
        return out[0], out[1]
    return None
