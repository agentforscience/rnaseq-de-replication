#!/usr/bin/env python3
"""Pre-ranked gene-set enrichment for every cached DE result, plus the random-set null.

Method choice.  `gseapy.prerank` implements the fgsea algorithm on a ranked list; the
ranking statistic is the DESeq2 Wald statistic, which encodes both effect size and
precision.  Geistlinger et al. (2021) show that re-running the DE model inside the GSEA
permutation loop (sample-permutation GSEA) costs large amounts of compute without changing
conclusions, so gene-permutation prerank is used throughout — and, importantly, it is a
*competitive* null, which is the family that does not call every set significant.

The random-set null.  For each DE result we additionally score `N_RANDOM` gene sets whose
sizes are drawn to match the hallmark size distribution but whose members are sampled
uniformly from that run's tested genes.  Applying the identical concordance metrics to
these sets answers the question the pathway literature usually leaves open: is pathway
robustness *biological coherence*, or is it just the variance reduction you get for free by
averaging 50-500 noisy genes?  (Venet et al. 2011.)
"""
from __future__ import annotations

import os
import sys
import warnings

import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
GSEA_DIR = os.path.join(ROOT, "results", "gsea")
GS_DIR = os.path.join(ROOT, "datasets", "genesets")

COLLECTIONS = {
    "hallmark": "h.all.v2023.2.Hs.symbols.gmt",      # primary: 50 curated, non-redundant sets
    "kegg": "c2.cp.kegg_legacy.v2023.2.Hs.symbols.gmt",
    "reactome": "c2.cp.reactome.v2023.2.Hs.symbols.gmt",
    "gobp": "c5.go.bp.v2023.2.Hs.symbols.gmt",
}
MIN_SIZE, MAX_SIZE = 15, 500
N_PERM = 1000
N_RANDOM = 200
RANDOM_SEED = 7


def read_gmt(path: str) -> dict:
    sets = {}
    with open(path) as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) > 2:
                sets[parts[0]] = [g for g in parts[2:] if g]
    return sets


def hallmark_sizes(background: set) -> list:
    """Sizes of hallmark sets after restriction to the background, for size matching."""
    h = read_gmt(os.path.join(GS_DIR, COLLECTIONS["hallmark"]))
    sizes = [len(set(v) & background) for v in h.values()]
    return [s for s in sizes if MIN_SIZE <= s <= MAX_SIZE]


def random_sets(background: list, sizes: list, n: int = N_RANDOM, seed: int = RANDOM_SEED) -> dict:
    """`n` random gene sets with sizes resampled from `sizes`, drawn from `background`.

    The seed is fixed and independent of the DE run, so the *same* random sets are scored in
    the discovery and replication studies — exactly as a real gene-set collection would be.
    """
    rng = np.random.default_rng(seed)
    bg = np.asarray(sorted(background))
    out = {}
    for i in range(n):
        k = int(rng.choice(sizes))
        k = min(k, len(bg))
        out[f"RANDOM_{i:04d}"] = list(rng.choice(bg, k, replace=False))
    return out


def rank_vector(de: pd.DataFrame) -> pd.Series:
    """Ranking statistic: DESeq2 Wald statistic, de-duplicated by gene symbol."""
    x = de["stat"].replace([np.inf, -np.inf], np.nan).dropna()
    return x.groupby(level=0).mean().sort_values(ascending=False)


def run_prerank(de: pd.DataFrame, gene_sets, threads: int = 1, seed: int = 1) -> pd.DataFrame:
    import gseapy as gp
    rnk = rank_vector(de)
    pre = gp.prerank(rnk=rnk, gene_sets=gene_sets, permutation_num=N_PERM,
                     min_size=MIN_SIZE, max_size=MAX_SIZE, threads=threads,
                     seed=seed, outdir=None, verbose=False)
    d = pre.res2d.copy()
    d["Term"] = d["Term"].astype(str)
    d = d.set_index("Term")
    for c in ("ES", "NES", "NOM p-val", "FDR q-val"):
        if c in d.columns:
            d[c] = pd.to_numeric(d[c], errors="coerce")
    return d[["ES", "NES", "NOM p-val", "FDR q-val"]].rename(
        columns={"NOM p-val": "pval", "FDR q-val": "fdr"})


def gsea_path(tag: str, collection: str) -> str:
    return os.path.join(GSEA_DIR, f"{tag}__{collection}.parquet")


def _worker(args):
    """Run all requested collections (plus the random-set null) for one DE tag."""
    tag, collections, do_random = args
    from de_core import load_de
    try:
        de = load_de(tag)
        rnk_genes = set(rank_vector(de).index)
        for coll in collections:
            out = gsea_path(tag, coll)
            if os.path.exists(out):
                continue
            res = run_prerank(de, os.path.join(GS_DIR, COLLECTIONS[coll]))
            res.to_parquet(out)
        if do_random:
            out = gsea_path(tag, "random")
            if not os.path.exists(out):
                sizes = hallmark_sizes(rnk_genes)
                rs = random_sets(sorted(rnk_genes), sizes)
                run_prerank(de, rs).to_parquet(out)
        return (tag, "ok", "")
    except Exception as exc:  # noqa: BLE001
        return (tag, "fail", f"{type(exc).__name__}: {exc}")


def run_gsea(tags, collections=("hallmark",), do_random=False, workers: int = 28,
             label: str = "gsea") -> pd.DataFrame:
    import multiprocessing as mp
    import time
    os.makedirs(GSEA_DIR, exist_ok=True)

    def done(tag):
        ok = all(os.path.exists(gsea_path(tag, c)) for c in collections)
        return ok and (not do_random or os.path.exists(gsea_path(tag, "random")))

    todo = [t for t in tags if not done(t)]
    print(f"[{label}] {len(tags)} tags, {len(tags)-len(todo)} cached, {len(todo)} to run",
          flush=True)
    if not todo:
        return pd.DataFrame(columns=["tag", "status", "error"])
    t0 = time.time()
    res = []
    with mp.get_context("fork").Pool(workers) as pool:
        for i, r in enumerate(pool.imap_unordered(
                _worker, [(t, list(collections), do_random) for t in todo], chunksize=1), 1):
            res.append(r)
            if i % 25 == 0 or i == len(todo):
                el = time.time() - t0
                print(f"[{label}] {i}/{len(todo)}  {el/60:.1f} min  "
                      f"eta {(el/i)*(len(todo)-i)/60:.1f} min  "
                      f"fails={sum(1 for x in res if x[1]=='fail')}", flush=True)
    df = pd.DataFrame(res, columns=["tag", "status", "error"])
    for _, row in df[df.status == "fail"].iterrows():
        print(f"  FAILED {row.tag}: {row.error}", flush=True)
    return df


# --------------------------------------------------------------------------------------
def load_gsea(tag: str, collection: str) -> pd.DataFrame:
    return pd.read_parquet(gsea_path(tag, collection))


def pathway_pair_stats(dis: pd.DataFrame, rep: pd.DataFrame, fdr: float = 0.05) -> dict:
    """Pathway replication statistics, deliberately mirroring `metrics.pair_stats`.

    `sign_conc`     <-> gene_sign_conc        (sign agreement among discovery hits)
    `bothsig`       <-> gene_bothsig_fdr      (discovery hit also significant in replication)
    `bothsig_sign`  <-> gene_bothsig_sign     (also significant *and* same direction)
    """
    common = dis.index.intersection(rep.index)
    d, r = dis.loc[common], rep.loc[common]
    out = {"n_sets_common": len(common)}
    sd = d.index[d.fdr < fdr]
    out["n_dis_sig"] = len(sd)
    out["n_rep_sig"] = int((r.fdr < fdr).sum())
    if len(sd) >= 5:
        same = np.sign(d.loc[sd, "NES"]) == np.sign(r.loc[sd, "NES"])
        out["sign_conc"] = float(same.mean())
        out["bothsig"] = float((r.loc[sd, "fdr"] < fdr).mean())
        out["bothsig_sign"] = float(((r.loc[sd, "fdr"] < fdr) & same).mean())
    else:
        out["sign_conc"] = out["bothsig"] = out["bothsig_sign"] = np.nan
    if len(common) >= 10:
        out["all_sign_conc"] = float((np.sign(d.NES) == np.sign(r.NES)).mean())
        out["nes_pearson"] = float(np.corrcoef(d.NES.fillna(0), r.NES.fillna(0))[0, 1])
        from scipy import stats as _st
        out["nes_spearman"] = float(_st.spearmanr(d.NES.fillna(0), r.NES.fillna(0)).statistic)
    else:
        out["all_sign_conc"] = out["nes_pearson"] = out["nes_spearman"] = np.nan
    return out
