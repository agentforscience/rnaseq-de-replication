#!/usr/bin/env python3
"""Core differential-expression machinery for the replication meta-analysis.

One entry point, `run_de_task`, executes a single DESeq2 tumour-vs-normal contrast on a
(possibly subsampled, possibly label-shuffled) slice of one recount3 cohort and caches the
result to `results/de/<tag>.parquet`.  Everything else in the project consumes those files.

Design notes
------------
* Cohorts are loaded through `code/load_recount3.py::build_cohort`, which applies its
  expression filter *per cohort*.  Any statistic that compares two cohorts must therefore
  intersect their tested-gene sets first (see `metrics.py`); this module deliberately does
  not do that intersection so that each DE run stays a faithful standalone analysis.
* Loading is memoised per worker process, because the parallel driver hands contiguous
  blocks of same-cohort tasks to each worker.
* Every run is fully determined by its `DeTask`, and the task's `tag` is a hash-free,
  human-readable filename, so the cache is inspectable and reruns are idempotent.
"""
from __future__ import annotations

import os
import sys
import time
import json
import warnings
from dataclasses import dataclass, field, asdict
from functools import lru_cache
from typing import Optional, Sequence

import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")
# Keep BLAS single-threaded: the parallel driver owns the core budget, and nested
# threading makes 30 workers thrash.
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
           "NUMEXPR_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ.setdefault(_v, "1")

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "code"))

DE_DIR = os.path.join(ROOT, "results", "de")


@lru_cache(maxsize=3)
def load_cohort_cached(cohort: str):
    """Load and cache one cohort's (counts, coldata).  Memoised per worker process."""
    from load_recount3 import build_cohort
    counts, col = build_cohort(cohort)
    return counts, col


@dataclass(frozen=True)
class DeTask:
    """A single differential-expression run.

    Attributes
    ----------
    cohort : recount3 cohort name, e.g. "TCGA-LIHC" or "SRP068976".
    tag : output filename stem; also the cache key.
    n_per_group : if set, subsample to this many tumour and this many normal samples.
    seed : RNG seed for subsampling / label shuffling.
    runs : explicit list of run_ids to use (overrides subsampling). Used by split controls.
    shuffle : if True, permute the condition labels within the selected samples
              (negative control).
    """
    cohort: str
    tag: str
    n_per_group: Optional[int] = None
    seed: int = 0
    runs: Optional[tuple] = None
    shuffle: bool = False


def _select_samples(col: pd.DataFrame, task: DeTask) -> list:
    """Choose the sample (run) ids for a task, honouring explicit lists then subsampling."""
    if task.runs is not None:
        sel = [r for r in task.runs if r in col.index]
        if len(sel) != len(task.runs):
            raise ValueError(f"{task.tag}: {len(task.runs) - len(sel)} requested runs absent")
        return sel
    if task.n_per_group is None:
        return list(col.index)
    rng = np.random.default_rng(task.seed)
    sel = []
    for cond in ("tumor", "normal"):
        pool = col.index[col.condition == cond].to_numpy()
        if len(pool) < task.n_per_group:
            raise ValueError(f"{task.cohort}: only {len(pool)} {cond} samples, "
                             f"need {task.n_per_group}")
        sel += list(rng.choice(pool, task.n_per_group, replace=False))
    return sel


def run_de_task(task: DeTask, n_cpus: int = 1, overwrite: bool = False) -> str:
    """Execute one DE run and write `results/de/<tag>.parquet`.  Returns the path.

    The parquet holds DESeq2's results frame (baseMean, log2FoldChange, lfcSE, stat, pvalue,
    padj) indexed by gene symbol, restricted to genes with a non-missing padj (i.e. genes
    that survived independent filtering and are therefore genuinely *testable* in this run).
    """
    out = os.path.join(DE_DIR, f"{task.tag}.parquet")
    if os.path.exists(out) and not overwrite:
        return out

    from pydeseq2.dds import DeseqDataSet
    from pydeseq2.ds import DeseqStats

    t0 = time.time()
    counts, col = load_cohort_cached(task.cohort)
    sel = _select_samples(col, task)

    cs = counts.loc[sel]
    cl = col.loc[sel, ["condition"]].copy()

    if task.shuffle:
        # Negative control: permute labels while preserving the group sizes exactly.
        rng = np.random.default_rng(task.seed + 10_000)
        cl["condition"] = rng.permutation(cl["condition"].to_numpy())

    cl["condition"] = pd.Categorical(cl["condition"], categories=["normal", "tumor"])
    if cl["condition"].value_counts().min() < 2:
        raise ValueError(f"{task.tag}: a condition group has <2 samples")

    # Drop genes that are all-zero in this particular slice; DESeq2 cannot fit them and
    # they would otherwise inflate the multiple-testing burden for small subsamples.
    cs = cs.loc[:, (cs.sum(axis=0) > 0)]

    dds = DeseqDataSet(counts=cs, metadata=cl, design="~condition",
                       refit_cooks=True, quiet=True, n_cpus=n_cpus)
    dds.deseq2()
    st = DeseqStats(dds, contrast=["condition", "tumor", "normal"], quiet=True, n_cpus=n_cpus)
    st.summary()

    res = st.results_df.dropna(subset=["padj"]).copy()
    res.index.name = "gene"
    res.attrs = {}
    res.to_parquet(out)

    meta = dict(asdict(task), n_samples=len(sel), n_genes_tested=int(res.shape[0]),
                n_tumor=int((cl.condition == "tumor").sum()),
                n_normal=int((cl.condition == "normal").sum()),
                seconds=round(time.time() - t0, 1))
    meta["runs"] = list(task.runs) if task.runs else None
    with open(os.path.join(DE_DIR, f"{task.tag}.meta.json"), "w") as fh:
        json.dump(meta, fh)
    return out


def _worker(task_and_cpus):
    task, ncpu = task_and_cpus
    try:
        run_de_task(task, n_cpus=ncpu)
        return (task.tag, "ok", "")
    except Exception as exc:  # noqa: BLE001 - failures are recorded, never silently dropped
        return (task.tag, "fail", f"{type(exc).__name__}: {exc}")


def run_tasks(tasks: Sequence[DeTask], workers: int = 28, n_cpus: int = 1,
              label: str = "de") -> pd.DataFrame:
    """Run many DE tasks in parallel.  Tasks are sorted by cohort so that each worker
    receives contiguous same-cohort blocks and the per-process load cache actually hits."""
    import multiprocessing as mp

    os.makedirs(DE_DIR, exist_ok=True)
    todo = [t for t in tasks if not os.path.exists(os.path.join(DE_DIR, f"{t.tag}.parquet"))]
    print(f"[{label}] {len(tasks)} tasks, {len(tasks) - len(todo)} cached, {len(todo)} to run",
          flush=True)
    if not todo:
        return pd.DataFrame(columns=["tag", "status", "error"])

    todo = sorted(todo, key=lambda t: (t.cohort, t.tag))
    t0 = time.time()
    results = []
    with mp.get_context("fork").Pool(workers) as pool:
        for i, r in enumerate(pool.imap_unordered(_worker, [(t, n_cpus) for t in todo],
                                                  chunksize=3), 1):
            results.append(r)
            if i % 25 == 0 or i == len(todo):
                el = time.time() - t0
                nf = sum(1 for x in results if x[1] == "fail")
                print(f"[{label}] {i}/{len(todo)}  {el/60:.1f} min  "
                      f"eta {(el/i)*(len(todo)-i)/60:.1f} min  fails={nf}", flush=True)
    df = pd.DataFrame(results, columns=["tag", "status", "error"])
    for _, row in df[df.status == "fail"].iterrows():
        print(f"  FAILED {row.tag}: {row.error}", flush=True)
    return df


def load_de(tag: str) -> pd.DataFrame:
    """Read a cached DE result by tag."""
    return pd.read_parquet(os.path.join(DE_DIR, f"{tag}.parquet"))


def de_exists(tag: str) -> bool:
    return os.path.exists(os.path.join(DE_DIR, f"{tag}.parquet"))
