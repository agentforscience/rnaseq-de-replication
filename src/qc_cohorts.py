#!/usr/bin/env python3
"""E0 - cohort label quality control.

The SRA tumour/normal labels in `datasets/manifest/samples.tsv` were derived by regular
expression from free-text sample descriptions, not by curation.  A cohort whose labels are
inverted or garbled would look exactly like a cohort that "failed to replicate", so the
replication numbers are only meaningful once each cohort is shown to produce a coherent
tumour-vs-normal contrast.

Two independent, pre-specified checks, both computed from the full-size DE run:

1. **Proliferation panel direction.**  Solid tumours over-express cell-cycle genes relative
   to adjacent normal tissue.  We require the median log2FC across a fixed proliferation
   panel to be positive.  This is direction-only and uses no external data.
2. **Agreement with the same-cancer-type TCGA reference.**  Spearman correlation of the Wald
   statistic against the TCGA cohort of the same cancer type over commonly tested genes.  A
   cohort with inverted labels gives a strongly *negative* correlation.  We require rho > 0.

A cohort fails QC if either check fails.  Failures are reported, excluded from the primary
analysis, and re-included in a sensitivity analysis.
"""
from __future__ import annotations

import os
import sys

import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import design
from de_core import load_de, de_exists

ROOT = design.ROOT

# Canonical proliferation / cell-cycle markers, up in essentially every solid tumour
# relative to matched normal tissue.  Fixed in advance; not tuned to the data.
PROLIFERATION = ["MKI67", "TOP2A", "CCNB1", "CDK1", "BUB1", "AURKA", "PLK1", "RRM2",
                 "TYMS", "PCNA", "CCNA2", "FOXM1", "UBE2C", "BIRC5", "TPX2"]

MIN_RHO = 0.0
MIN_PROLIF_LFC = 0.0


def main():
    cohorts = design.cohort_table()
    tags = {c: f"full__{c}" for c in cohorts.cohort}
    missing = [c for c, t in tags.items() if not de_exists(t)]
    if missing:
        raise SystemExit(f"missing full-size DE for: {missing}")

    de = {c: load_de(t) for c, t in tags.items()}
    ref = {ct: f"TCGA-{ct}" for ct in cohorts.cancer_type.unique()
           if f"TCGA-{ct}" in de}

    rows = []
    for _, r in cohorts.iterrows():
        d = de[r.cohort]
        panel = [g for g in PROLIFERATION if g in d.index]
        prolif = float(d.loc[panel, "log2FoldChange"].median()) if panel else np.nan
        rho, n_common = np.nan, 0
        refc = ref.get(r.cancer_type)
        if refc and refc != r.cohort:
            rd = de[refc]
            common = d.index.intersection(rd.index)
            n_common = len(common)
            if n_common > 500:
                rho = float(stats.spearmanr(d.loc[common, "stat"], rd.loc[common, "stat"]).statistic)
        n_sig = int(((d.padj < 0.05) & (d.log2FoldChange.abs() > 1)).sum())
        ok_prolif = bool(np.isfinite(prolif) and prolif > MIN_PROLIF_LFC)
        # If no TCGA reference exists for the type (LUNG), the correlation check is waived.
        ok_rho = bool(np.isnan(rho) or rho > MIN_RHO)
        rows.append(dict(cohort=r.cohort, cancer_type=r.cancer_type, source=r.source,
                         n_tumor=r.n_tumor, n_normal=r.n_normal,
                         n_genes_tested=len(d), n_sig_joint=n_sig,
                         prolif_lfc=prolif, ref_cohort=refc if refc != r.cohort else None,
                         ref_spearman=rho, n_common_ref=n_common,
                         qc_prolif=ok_prolif, qc_ref=ok_rho,
                         qc_pass=ok_prolif and ok_rho))
    qc = pd.DataFrame(rows).sort_values(["cancer_type", "cohort"])
    out = os.path.join(ROOT, "results/metrics/cohort_qc.tsv")
    qc.to_csv(out, sep="\t", index=False)

    print(qc.to_string(index=False, float_format=lambda x: f"{x:.3f}"))
    fails = qc[~qc.qc_pass]
    print(f"\nQC: {qc.qc_pass.sum()}/{len(qc)} cohorts pass")
    if len(fails):
        print("FAIL:")
        print(fails[["cohort", "cancer_type", "prolif_lfc", "ref_spearman"]].to_string(index=False))
    print(f"\nwrote {out}")


if __name__ == "__main__":
    main()
