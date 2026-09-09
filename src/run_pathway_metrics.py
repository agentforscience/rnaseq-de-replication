#!/usr/bin/env python3
"""Pathway-level replication statistics, on exactly the same pairs as the gene-level arm.

Writes `results/metrics/pathways_<arm>.tsv` with one row per (pair, collection).  The
`random` collection is the size-matched random-gene-set null and is treated as just another
collection so that it goes through the identical code path -- no special-casing, which is
what makes it a fair null.
"""
from __future__ import annotations

import os
import sys

import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import design
import gsea

ROOT = design.ROOT
MET = os.path.join(ROOT, "results", "metrics")


def _job(args):
    key, ta, tb, coll = args
    rec = dict(key); rec["collection"] = coll
    try:
        a = gsea.load_gsea(ta, coll)
        b = gsea.load_gsea(tb, coll)
        rec.update(gsea.pathway_pair_stats(a, b))
        rec["status"] = "ok"
    except Exception as exc:  # noqa: BLE001
        rec["status"] = f"fail: {type(exc).__name__}: {exc}"
    return rec


def _run(jobs, label, workers=28):
    import multiprocessing as mp, time
    print(f"[{label}] {len(jobs)} comparisons", flush=True)
    t0 = time.time(); out = []
    with mp.get_context("fork").Pool(workers) as pool:
        for i, r in enumerate(pool.imap_unordered(_job, jobs, chunksize=8), 1):
            out.append(r)
            if i % 500 == 0 or i == len(jobs):
                print(f"[{label}] {i}/{len(jobs)}  {(time.time()-t0)/60:.1f} min", flush=True)
    return pd.DataFrame(out)


def build_jobs():
    cohorts = design.cohort_table()
    pairs = design.cohort_pairs(cohorts)
    jobs = []

    def add(key, ta, tb, colls):
        for c in colls:
            if os.path.exists(gsea.gsea_path(ta, c)) and os.path.exists(gsea.gsea_path(tb, c)):
                jobs.append((key, ta, tb, c))

    for _, p in pairs.iterrows():
        add(dict(arm="cross_full", cancer_type=p.cancer_type, discovery=p.discovery,
                 replication=p.replication, seed=-1),
            f"full__{p.discovery}", f"full__{p.replication}",
            ("hallmark", "random", "kegg", "reactome", "gobp"))
        for s in range(design.N_SUBSAMPLE_SEEDS):
            add(dict(arm="cross_m10", cancer_type=p.cancer_type, discovery=p.discovery,
                     replication=p.replication, seed=s),
                f"m10__{p.discovery}__s{s}", f"m10__{p.replication}__s{s}",
                ("hallmark", "random", "kegg", "reactome"))

    idx = pd.read_csv(os.path.join(MET, "split_index.tsv"), sep="\t")
    for _, r in idx.iterrows():
        ct = r.cohort.replace("TCGA-", "")
        for d, rep, dirn in ((r.tag_a, r.tag_b, "AB"), (r.tag_b, r.tag_a, "BA")):
            add(dict(arm=r.kind, cancer_type=ct, discovery=r.cohort, replication=r.cohort,
                     cohort=r.cohort, rep=int(r.rep), direction=dirn, seed=int(r.rep)),
                d, rep, ("hallmark", "random"))
    return jobs


def main():
    path = os.path.join(MET, "pathways.tsv")
    if os.path.exists(path):
        print("cached"); return
    jobs = build_jobs()
    df = _run(jobs, "pathways")
    df.to_csv(path, sep="\t", index=False)
    print(f"wrote {path}  rows={len(df)}  failures={int((df.status!='ok').sum())}")
    print(df.groupby(['arm', 'collection']).size())


if __name__ == "__main__":
    main()
