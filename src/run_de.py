#!/usr/bin/env python3
"""Execute every differential-expression run required by the study.

Four arms, all writing to `results/de/`:

  full      one DE per eligible cohort at its native sample size            (51 runs)
  matched   5 independent 10-vs-10 subsamples per cohort                    (255 runs)
  rsplit    15 disjoint random 10v10/10v10 split-pairs per large TCGA cohort (270 runs)
  tsplit    15 tissue-source-site-*disjoint* split-pairs per cohort          (270 runs)
  shuffle   label-permuted negative control on the rsplit sample sets        (144 runs)

Usage:  python src/run_de.py [arm ...]     (default: all arms, in the order above)
"""
from __future__ import annotations

import os
import sys
import json

import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import design
from de_core import DeTask, run_tasks, DE_DIR

ROOT = design.ROOT
N_SHUFFLE_REPS = 8


def full_tasks(cohorts):
    return [DeTask(cohort=c, tag=f"full__{c}") for c in cohorts.cohort]


def matched_tasks(cohorts):
    """10v10 subsamples.  Only cohorts with >=10 in both groups can participate."""
    out = []
    for _, r in cohorts.iterrows():
        if r.n_tumor < design.MATCHED_N or r.n_normal < design.MATCHED_N:
            continue
        for s in range(design.N_SUBSAMPLE_SEEDS):
            out.append(DeTask(cohort=r.cohort, tag=f"m10__{r.cohort}__s{s}",
                              n_per_group=design.MATCHED_N, seed=s))
    return out


def split_tasks(cohorts):
    """Within-cohort controls: random splits, TSS-disjoint splits, shuffled labels."""
    rs, ts, sh = [], [], []
    index = []
    for c in design.split_cohorts(cohorts):
        for rep in range(design.N_SPLIT_REPS):
            a, b = design.random_split_runs(c, rep)
            rs += [DeTask(cohort=c, tag=f"rsplit__{c}__r{rep}__A", runs=a, seed=rep),
                   DeTask(cohort=c, tag=f"rsplit__{c}__r{rep}__B", runs=b, seed=rep)]
            index.append(dict(kind="rsplit", cohort=c, rep=rep,
                              tag_a=f"rsplit__{c}__r{rep}__A", tag_b=f"rsplit__{c}__r{rep}__B"))
            if rep < N_SHUFFLE_REPS:
                sh += [DeTask(cohort=c, tag=f"shuf__{c}__r{rep}__A", runs=a, seed=rep, shuffle=True),
                       DeTask(cohort=c, tag=f"shuf__{c}__r{rep}__B", runs=b, seed=rep + 500, shuffle=True)]
                index.append(dict(kind="shuffle", cohort=c, rep=rep,
                                  tag_a=f"shuf__{c}__r{rep}__A", tag_b=f"shuf__{c}__r{rep}__B"))
            tt = design.tss_split_runs(c, rep)
            if tt is not None:
                ta, tb = tt
                ts += [DeTask(cohort=c, tag=f"tsplit__{c}__r{rep}__A", runs=ta, seed=rep),
                       DeTask(cohort=c, tag=f"tsplit__{c}__r{rep}__B", runs=tb, seed=rep)]
                index.append(dict(kind="tsplit", cohort=c, rep=rep,
                                  tag_a=f"tsplit__{c}__r{rep}__A", tag_b=f"tsplit__{c}__r{rep}__B"))
    pd.DataFrame(index).to_csv(os.path.join(ROOT, "results/metrics/split_index.tsv"),
                               sep="\t", index=False)
    return rs, ts, sh


def main(argv):
    arms = argv[1:] or ["full", "matched", "rsplit", "tsplit", "shuffle"]
    cohorts = design.cohort_table()
    cohorts.to_csv(os.path.join(ROOT, "results/metrics/cohorts_used.tsv"), sep="\t", index=False)
    design.cohort_pairs(cohorts).to_csv(
        os.path.join(ROOT, "results/metrics/pairs.tsv"), sep="\t", index=False)

    os.makedirs(DE_DIR, exist_ok=True)
    rs, ts, sh = split_tasks(cohorts)
    plan = {"full": full_tasks(cohorts), "matched": matched_tasks(cohorts),
            "rsplit": rs, "tsplit": ts, "shuffle": sh}

    logs = []
    for arm in arms:
        # Large TCGA cohorts benefit from >1 core; small subsamples do not, so the full-size
        # arm gets fewer, fatter workers and every other arm gets maximum breadth.
        workers, ncpu = (10, 3) if arm == "full" else (28, 1)
        logs.append(run_tasks(plan[arm], workers=workers, n_cpus=ncpu, label=arm))
    log = pd.concat(logs) if logs else pd.DataFrame()
    if len(log):
        log.to_csv(os.path.join(ROOT, "logs", "de_run_log.tsv"), sep="\t", index=False)
        print(f"\nDone. ok={int((log.status=='ok').sum())} fail={int((log.status=='fail').sum())}")


if __name__ == "__main__":
    main(sys.argv)
