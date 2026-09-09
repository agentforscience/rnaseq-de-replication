#!/usr/bin/env bash
# Full experimental pipeline, in order.  Every stage is cached: rerunning skips completed
# work, so the script is safe to restart.  Expected wall time on 32 cores: ~2.5 h from
# scratch (the DE stage dominates), or seconds if all intermediates are present.
set -euo pipefail
cd "$(dirname "$0")/.."
source .venv/bin/activate

echo "== 1/6  differential expression (all arms) =="
python src/run_de.py                      2>&1 | tee -a logs/run_de.log

echo "== 2/6  cohort label QC (E0) =="
python src/qc_cohorts.py                  2>&1 | tee    logs/qc.log

echo "== 3/6  gene-level pairwise metrics (E1-E4, E7-E9) =="
python src/run_metrics.py                 2>&1 | tee    logs/run_metrics.log

echo "== 4/6  pre-ranked GSEA + size-matched random-set null (E5, E6) =="
python src/run_gsea.py primary            2>&1 | tee    logs/run_gsea_primary.log
python src/run_gsea.py sensitivity        2>&1 | tee    logs/run_gsea_sensitivity.log

echo "== 5/6  pathway-level pairwise metrics =="
python src/run_pathway_metrics.py         2>&1 | tee    logs/run_pathway_metrics.log

echo "== 6/6  analysis, tables, figures =="
python src/analyze.py                     2>&1 | tee    logs/analyze.log

echo "done."
