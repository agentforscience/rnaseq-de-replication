#!/usr/bin/env python3
"""Run pre-ranked GSEA over the cached DE results.

Coverage plan (chosen so the primary comparison is complete and the sensitivity analysis is
affordable):

  hallmark + random-set null   every DE run (cross-cohort, matched, and all controls)
  kegg + reactome              full-size and matched cross-cohort runs
  gobp                         full-size cross-cohort runs only

`hallmark` is the primary collection: 50 curated, deliberately non-redundant sets, which is
what most cancer papers actually report.  `random` is the size-matched null that decides
whether any pathway advantage is biological coherence or just aggregation.
"""
from __future__ import annotations

import os
import sys

import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import design
import gsea
from de_core import de_exists

ROOT = design.ROOT
MET = os.path.join(ROOT, "results", "metrics")


def tag_sets():
    cohorts = design.cohort_table()
    full = [f"full__{c}" for c in cohorts.cohort if de_exists(f"full__{c}")]
    m10 = [f"m10__{c}__s{s}" for c in cohorts.cohort
           for s in range(design.N_SUBSAMPLE_SEEDS) if de_exists(f"m10__{c}__s{s}")]
    idx = pd.read_csv(os.path.join(MET, "split_index.tsv"), sep="\t")
    split = sorted({t for col in ("tag_a", "tag_b") for t in idx[col] if de_exists(t)})
    return full, m10, split


def main(argv):
    full, m10, split = tag_sets()
    print(f"tags: full={len(full)} m10={len(m10)} split={len(split)}")
    stage = argv[1] if len(argv) > 1 else "all"
    logs = []
    if stage in ("all", "primary"):
        # Primary: hallmark + random-set null on everything.
        logs.append(gsea.run_gsea(full, ("hallmark",), do_random=True, label="h+rand:full"))
        logs.append(gsea.run_gsea(m10, ("hallmark",), do_random=True, label="h+rand:m10"))
        logs.append(gsea.run_gsea(split, ("hallmark",), do_random=True, label="h+rand:split"))
    if stage in ("all", "sensitivity"):
        logs.append(gsea.run_gsea(full + m10, ("kegg", "reactome"), label="kegg+react"))
        logs.append(gsea.run_gsea(full, ("gobp",), label="gobp:full"))
    log = pd.concat([x for x in logs if len(x)]) if logs else pd.DataFrame()
    if len(log):
        log.to_csv(os.path.join(ROOT, "logs", f"gsea_run_log_{stage}.tsv"), sep="\t", index=False)
        print(f"\nDone. ok={int((log.status=='ok').sum())} fail={int((log.status=='fail').sum())}")


if __name__ == "__main__":
    main(sys.argv)
