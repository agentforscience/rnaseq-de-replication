#!/usr/bin/env python3
"""Compute every pairwise replication statistic from the cached DE runs.

Four comparison arms, all written to `results/metrics/`:

  pairs_full.tsv    independent cohort pairs, native sample sizes   (E1, E2, E4)
  pairs_m10.tsv     independent cohort pairs, both at 10v10         (E3 - the primary arm)
  splits.tsv        within-cohort controls at 10v10                 (E7, E8, E9)

The matched-n arm is the primary one for every cross-arm claim: cohort sizes span 10v10 to
1135v114, so an unmatched comparison of "independent studies" against "random splits of one
study" would measure statistical power rather than replication.
"""
from __future__ import annotations

import os
import sys
import itertools

import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import design
import metrics
from de_core import load_de, de_exists

ROOT = design.ROOT
MET = os.path.join(ROOT, "results", "metrics")


def _job(args):
    """One (discovery tag, replication tag) comparison, returned as a flat record."""
    key, tag_a, tag_b = args
    try:
        a, b = load_de(tag_a), load_de(tag_b)
        rec = dict(key)
        rec.update(metrics.pair_stats(a, b))
        rec["status"] = "ok"
        return rec
    except Exception as exc:  # noqa: BLE001
        rec = dict(key)
        rec["status"] = f"fail: {type(exc).__name__}: {exc}"
        return rec


def _run(jobs, label, workers=28):
    import multiprocessing as mp
    import time
    print(f"[{label}] {len(jobs)} comparisons", flush=True)
    t0 = time.time()
    out = []
    with mp.get_context("fork").Pool(workers) as pool:
        for i, r in enumerate(pool.imap_unordered(_job, jobs, chunksize=4), 1):
            out.append(r)
            if i % 200 == 0 or i == len(jobs):
                el = time.time() - t0
                print(f"[{label}] {i}/{len(jobs)}  {el/60:.1f} min "
                      f"eta {(el/i)*(len(jobs)-i)/60:.1f} min", flush=True)
    return pd.DataFrame(out)


def cross_full(cohorts, pairs):
    jobs = []
    for _, p in pairs.iterrows():
        ta, tb = f"full__{p.discovery}", f"full__{p.replication}"
        if de_exists(ta) and de_exists(tb):
            jobs.append((dict(arm="cross_full", cancer_type=p.cancer_type,
                              discovery=p.discovery, replication=p.replication, seed=-1),
                         ta, tb))
    return _run(jobs, "cross_full")


def cross_m10(cohorts, pairs):
    """Independent cohort pairs with both studies subsampled to 10 tumour / 10 normal.

    Seed s of the discovery cohort is paired with seed s of the replication cohort, giving
    five independent power-matched replicates of every ordered pair.
    """
    jobs = []
    for _, p in pairs.iterrows():
        for s in range(design.N_SUBSAMPLE_SEEDS):
            ta, tb = f"m10__{p.discovery}__s{s}", f"m10__{p.replication}__s{s}"
            if de_exists(ta) and de_exists(tb):
                jobs.append((dict(arm="cross_m10", cancer_type=p.cancer_type,
                                  discovery=p.discovery, replication=p.replication, seed=s),
                             ta, tb))
    return _run(jobs, "cross_m10")


def splits():
    """Within-cohort controls.  Both directions of each split-pair are scored, because the
    A/B labelling is arbitrary and averaging both removes any asymmetry."""
    idx = pd.read_csv(os.path.join(MET, "split_index.tsv"), sep="\t")
    jobs = []
    for _, r in idx.iterrows():
        if not (de_exists(r.tag_a) and de_exists(r.tag_b)):
            continue
        ct = r.cohort.replace("TCGA-", "")
        for d, rep, dirn in ((r.tag_a, r.tag_b, "AB"), (r.tag_b, r.tag_a, "BA")):
            jobs.append((dict(arm=r.kind, cancer_type=ct,
                              discovery=r.cohort, replication=r.cohort,
                              cohort=r.cohort, rep=int(r.rep), direction=dirn, seed=int(r.rep)),
                         d, rep))
    return _run(jobs, "splits")


def main():
    os.makedirs(MET, exist_ok=True)
    cohorts = design.cohort_table()
    pairs = design.cohort_pairs(cohorts)

    for name, fn in (("pairs_full", lambda: cross_full(cohorts, pairs)),
                     ("pairs_m10", lambda: cross_m10(cohorts, pairs)),
                     ("splits", splits)):
        path = os.path.join(MET, f"{name}.tsv")
        if os.path.exists(path):
            print(f"[{name}] cached")
            continue
        df = fn()
        df.to_csv(path, sep="\t", index=False)
        nf = int((df.status != "ok").sum()) if "status" in df else 0
        print(f"wrote {path}  rows={len(df)}  failures={nf}")


if __name__ == "__main__":
    main()
