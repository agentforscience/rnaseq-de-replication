# Research Plan

**Phase:** `experiment_runner` (Phases 0–6). This document begins with the Motivation &
Novelty assessment (Phase 0), then the full experimental plan (Phase 1). The **Direction
Budget** produced during resource gathering is retained verbatim at the end and is
unchanged — the three retained directions D1/D2/D3 are exactly what is implemented below.

---

## Motivation & Novelty Assessment

### Why This Research Matters

Differentially expressed (DE) gene lists are the single most common output of cancer
transcriptomics: a typical paper names a handful of genes from such a list and builds an
entire follow-up program — antibodies, knockdowns, prognostic panels, occasionally clinical
translation — on the assumption that those genes are a stable property of the disease rather
than of the cohort. The Reproducibility Project: Cancer Biology (Errington 2021) found a 46%
replication success rate with effect sizes a median 85% smaller, and Ioannidis (2009) could
not reproduce 10 of 18 published microarray analyses at all. If the underlying gene lists are
themselves unstable across cohorts, that is an upstream cause of downstream irreproducibility
that is cheap to measure and currently under-quantified for **RNA-seq** specifically (most of
the classic evidence — MAQC, Ein-Dor, Michiels, Shi — is microarray-era). Anyone who reads,
reviews, funds, or builds on a cancer DE paper benefits from a calibrated expectation of how
much of such a list should be believed.

### Gap in Existing Work

From `literature_review.md`:

1. **The canonical numbers are microarray-era and mostly technical.** Shi 2008 and MAQC
   measure *inter-site technical* replication of the same RNA. Schurch 2016 measures
   *within-experiment* replicate subsetting in yeast. Neither answers "two independent human
   cancer studies, different patients, different countries, different years".
2. **Pipeline heterogeneity is nearly always confounded with study heterogeneity.** Bigler
   2013 shows harmonising the pipeline moves concordance from 67% to >99% — so any
   cross-study number computed from author-supplied result tables is measuring pipelines as
   much as biology. recount3 removes that confound entirely by reprocessing every study with
   one pipeline (Monorail/STAR, GENCODE v26).
3. **The gene-vs-pathway comparison is almost always made unfairly.** The widely repeated
   claim that "pathway analysis is more robust" is typically supported by comparing a
   *thresholded* gene overlap against an *unthresholded* pathway direction agreement. Venet
   2011 showed random signatures are significant predictors of breast cancer outcome; the
   analogous null — do *random gene sets* also show >70% directional concordance? — is, to our
   knowledge, not reported anywhere in the cross-study replication literature.
4. **Nobody separates "did not replicate" from "was underpowered to replicate".** A
   replication rate is a single number that conflates genuine biological heterogeneity with
   the replication cohort simply having 10 samples.

### Our Novel Contribution

Three things that, together, are not available in the existing literature:

- **A modern, pipeline-controlled RNA-seq replication estimate.** 161 independent cohort
  pairs across 11 cancer types, all counts derived from raw reads by one pipeline, so the
  measured discordance is study/biology, not processing.
- **A like-for-like gene-vs-pathway comparison with a random-set null.** Every pathway
  statistic is reported next to its exact gene-level analogue (sign vs sign,
  both-significant vs both-significant) *and* next to size-matched random gene sets. This
  turns "pathways are more robust" from a rhetorical claim into a testable one that can fail.
- **A variance decomposition of non-replication.** Random 50/50 splits of a single cohort
  (pure sampling noise), tissue-source-site-**disjoint** splits of the same cohort (adds
  collection-centre/population/batch, pipeline and protocol held fixed), and genuinely
  independent cohorts — all at *identical* sample size (10 vs 10) — attribute the shortfall
  to power versus heterogeneity. Plus Storey's π₁, which measures how much signal is present
  in the replication cohort *below* the significance threshold.

### Experiment Justification

| Experiment | Why it is necessary |
|---|---|
| **E0 — Label QC** | SRA tumour/normal labels are regex-derived from free text (a stated risk). Before any replication number is trusted, each cohort must be shown to produce a coherent tumour-vs-normal contrast. Without this, a mislabelled cohort would masquerade as "non-replication". |
| **E1 — Full-size cross-cohort replication** | The headline number the hypothesis asks for: does the pooled gene-level replication rate fall below 40%? Must be computed on commonly-tested genes only. |
| **E2 — Selection-rule sweep** | Shi 2008 shows the *rule* used to pick genes moves overlap from 20–40% to ~90% on identical data. A single number without this axis is not interpretable; it also pre-empts the objection that the result is an artefact of the FDR+FC cutoff. |
| **E3 — Matched-n (10v10) replication** | Cohorts span 10v10 to 1135v114. Uncorrected, a cross-cohort comparison measures *power*, not replication. Fixing every cohort at 10v10 makes all pairs and all controls directly comparable. |
| **E4 — π₁ (Storey) in the replication cohort** | Separates "the effect is absent" from "the effect is present but sub-threshold". This is the difference between a crisis and a power problem, and it cannot be read off a replication rate. |
| **E5 — Pathway replication, matched comparators** | The second half of the hypothesis. Must use the same cohort pairs, the same DE results, and comparators matched to E1 (sign→sign, both-sig→both-sig). |
| **E6 — Size-matched random-gene-set null** | Determines whether any pathway advantage is *biological coherence* or merely *averaging 50–500 noisy genes*. This is the experiment that can falsify the "pathways are more robust" interpretation while leaving the descriptive number intact. |
| **E7 — Within-cohort random splits (ceiling)** | The pure-sampling-noise upper bound on replication at the same n. Whatever gap remains between this and cross-cohort is attributable to between-study differences. |
| **E8 — TSS-disjoint splits (intermediate control)** | Adds collection centre, population and batch while holding cancer type, protocol and pipeline fixed. Placed between E7 and E1, it partitions the gap into "site/batch" and "everything else". |
| **E9 — Shuffled-label negative control** | Establishes the floor. Confirms that the observed replication is not attainable by chance and that the metric behaves. |

---

## Research Question

For bulk RNA-seq studies of the same cancer type, processed through an identical pipeline:
what fraction of genes called differentially expressed (FDR < 0.05, |log2FC| > 1) in one
study are called differentially expressed at the same thresholds in an independent study;
and is directional concordance of GSEA normalized enrichment scores (NES) materially higher
than the matched gene-level comparator once a size-matched random-gene-set null is accounted
for?

## Background and Motivation

See *Motivation & Novelty Assessment* above and `literature_review.md`. Anchors we will be
measured against: 20–40% (Shi 2008, P-ranked inter-site), ~90% (Shi 2008, FC-ranked), 20–40%
(Schurch 2016, 3 vs 42 replicates), 46% (Errington 2021), 67%→>99% (Bigler 2013, pipeline
harmonisation).

## Hypothesis Decomposition

The stated hypothesis has two conjuncts and one interpretation. We pre-register them
separately, with directions of test fixed **before** looking at the full results.

| ID | Claim | Test | Pre-registered decision rule |
|---|---|---|---|
| **H1** | Gene-level replication rate < 40% | Pooled rate over independent cohort pairs, FDR<0.05 & \|log2FC\|>1, commonly-tested genes | One-sided Wilcoxon signed-rank of per-pair rates against 0.40; supported if pooled rate and its 95% CI lie below 0.40 |
| **H2** | Pathway NES directional concordance > 70% | Hallmark GSEA on the same pairs, sign of NES among discovery-significant sets | One-sided Wilcoxon against 0.70; supported if pooled concordance and CI lie above 0.70 |
| **H3** *(the interpretation)* | Pathway analysis is **more robust** than gene-level reporting | Matched comparators: pathway sign concordance vs **gene** sign concordance; pathway both-significant vs **gene** both-significant. Random-set null for the pathway numbers | Supported only if pathway > gene on **matched** comparators *and* pathway > size-matched random sets. Paired Wilcoxon across pairs, Bonferroni over the comparator families |

H1 and H2 can both be true while H3 is false. That is the outcome the design is built to be
able to detect, and it is the scientifically interesting one.

**Independent variables:** cohort pair identity; selection rule; sample size (full vs matched
10v10); split type (independent study / TSS-disjoint / random / shuffled); gene-set
collection; unit of analysis (gene vs pathway vs random set).
**Dependent variables:** replication rate, sign concordance, both-significant rate, Spearman
correlation of log2FC and of Wald statistic, π₁, NES sign concordance, NES correlation.

## Proposed Methodology

### Approach

Reanalyse from raw counts, never from published tables. All 73 cohorts come from recount3, so
alignment, annotation and quantification are constant by construction and the only thing that
varies across a "pair" is the study. Every cohort goes through the same DESeq2 (pydeseq2)
tumour-vs-normal contrast; every replication statistic is computed on the intersection of the
two cohorts' testable genes; every cross-cohort number has a within-cohort counterpart at the
same sample size.

### Experimental Steps

1. **Build cohorts** with `code/load_recount3.py::build_cohort` — coverage→read-count
   conversion, protein-coding genes expressed at ≥10 counts in ≥20% of samples, symbols
   collapsed by sum. *Rationale:* per-cohort filtering is correct for DE, and forces the
   gene-intersection step in step 5.
2. **DE at full size** for the 51 cohorts belonging to the 11 cancer types with ≥2 cohorts
   (`~condition`, Wald test, Cook's refitting, independent filtering). *Rationale:* DESeq2 is
   the field default and Soneson/Rapaport show tool choice is second-order next to cohort
   effects (a pruned direction, D12).
3. **Label QC (E0).** For every cohort, correlate its Wald-statistic vector against the
   same-cancer-type TCGA reference over common genes, and check a proliferation panel
   (MKI67, TOP2A, CCNB1, …) is up in "tumour". Cohorts failing QC are reported and excluded
   from primary analysis, with a sensitivity analysis including them.
4. **DE at matched n = 10 tumour vs 10 normal**, 5 independent random subsamples per cohort
   (seeds 0–4). *Rationale:* removes power as a confound and makes E1/E7/E8/E9 comparable.
5. **Pairwise gene metrics (E1–E4).** For every ordered pair within a cancer type: intersect
   testable genes; replication rate; sign concordance among discovery hits; Spearman of
   log2FC and of Wald statistic over all common genes; Storey π₁ on the replication cohort's
   p-values restricted to discovery hits.
6. **Selection-rule sweep (E2).** FDR-only; \|log2FC\|>1-only; joint (primary); top-100 /
   top-500 by p-value; top-100 / top-500 by \|log2FC\| among FDR<0.05.
7. **Pre-ranked GSEA (E5).** `gseapy.prerank` (fgsea algorithm) on the Wald statistic —
   hallmark (primary), KEGG, Reactome, GO:BP (sensitivity). Not GSEA-with-DESeq2-inside-the-
   permutation loop (Geistlinger 2021: costly, no change in result).
8. **Random-set null (E6).** 200 random gene sets drawn from the tested background with sizes
   matched to the hallmark size distribution, run through the identical prerank and identical
   concordance metrics.
9. **Within-cohort controls (E7–E9).** For each TCGA cohort with ≥25 normals: 20 disjoint
   random 10v10/10v10 split-pairs; 20 TSS-**disjoint** 10v10/10v10 split-pairs (no
   tissue-source-site appears on both sides); and a shuffled-condition-label control.
10. **Meta-analysis and inference.** Per-cancer-type summaries and a pooled estimate that
    weights cancer types equally, so LIHC/BRCA/PRAD (30 of 51 cohorts) cannot silently drive
    the headline number.

### Baselines

| Baseline | Role |
|---|---|
| Within-cohort **random** 10v10 splits | Sampling-noise **ceiling** — the best replication achievable when there is no study heterogeneity at all |
| Within-cohort **TSS-disjoint** 10v10 splits | Intermediate: adds centre/population/batch, holds protocol+pipeline fixed |
| **Shuffled-label** control | **Floor** — replication attainable by chance |
| **Size-matched random gene sets** | Null for every pathway concordance claim (Venet 2011) |
| **Gene-level sign / both-significant** | The matched comparator for the pathway claim (Geistlinger 2021) |
| Literature anchors | External calibration: 20–40%, ~90%, 46%, 67%→>99% |

### Evaluation Metrics

- **Replication rate** = |{discovery hits} ∩ {replication hits}| / |{discovery hits}|,
  over commonly-tested genes, both directions of every pair reported separately.
- **Sign concordance** of log2FC among discovery-significant genes (the matched analogue of
  pathway NES sign concordance).
- **Spearman ρ** of log2FC and of the Wald statistic over all commonly-tested genes —
  threshold-free, immune to cutoff artefacts.
- **Storey's π₁** on replication-cohort p-values restricted to discovery hits — the estimated
  fraction of discovery hits that are genuinely non-null in the replication cohort.
- **Pathway:** NES sign concordance, both-significant rate, Pearson/Spearman of NES.
- All reported with bootstrap 95% CIs over cohort pairs.

### Statistical Analysis Plan

- α = 0.05 throughout. Pair-level statistics are **not** independent (cohorts recur across
  pairs), so primary inference uses **cluster bootstrap resampling over cohorts**
  (10,000 resamples) rather than treating pairs as independent observations; Wilcoxon
  signed-rank tests are reported alongside as a secondary, anticonservative check and
  labelled as such.
- H1: one-sided test of pooled rate < 0.40. H2: one-sided test of pooled concordance > 0.70.
- H3: paired comparisons within cohort pair (pathway vs matched gene comparator; pathway vs
  random-set null), Wilcoxon signed-rank with Bonferroni correction across the comparator
  families, effect sizes as median paired differences with bootstrap CIs.
- Multiple testing inside DE is DESeq2's Benjamini–Hochberg; inside GSEA it is the
  permutation-based FDR q-value from fgsea.
- Reported effect sizes: median differences with CIs (rates are bounded and non-normal, so
  Cohen's d is not the right summary; rank-biserial correlation given for Wilcoxon tests).

## Expected Outcomes

- **Supports the hypothesis:** pooled joint-threshold replication < 40%; hallmark NES sign
  concordance > 70%; pathway comparators beat matched gene comparators *and* the random-set
  null.
- **Refutes / complicates:** (a) replication ≥ 40% once the pipeline is held constant —
  plausible given Bigler 2013; (b) gene-level *sign* concordance ≈ pathway sign concordance,
  making H3 an artefact of comparing unmatched quantities — flagged as likely by the
  preliminary smoke test (96% gene sign vs 95% pathway sign); (c) random gene sets showing
  high sign concordance, meaning the pathway advantage is aggregation, not biology.
- **The decomposition** will show whether the shortfall is power (random-split ceiling also
  low) or heterogeneity (random-split ceiling high, cross-cohort low).

## Timeline and Milestones

| Milestone | Estimate |
|---|---|
| Plan + scaffolding | 20 min |
| Full-size DE, 51 cohorts + label QC | 25 min |
| Matched-n DE (51 × 5 seeds) | 20 min |
| Within-cohort split DE (~900 runs) | 40 min |
| Pairwise gene metrics + selection sweep | 15 min |
| GSEA (all DE results × 4 collections) + random sets | 40 min |
| Analysis, statistics, figures | 40 min |
| Report + validation | 40 min |
| Buffer (~25%) | 60 min |

## Potential Challenges

| Risk | Mitigation |
|---|---|
| SRA labels wrong (regex-derived) | E0 label QC against TCGA reference + proliferation panel; failures excluded and reported |
| Cohort size confound | Matched 10v10 arm is the primary comparison for all cross-arm claims |
| Cancer-type imbalance (LIHC 13, BRCA 9, PRAD 8) | Type-equal-weighted pooling; per-type tables always shown |
| Pair non-independence | Cluster bootstrap over cohorts |
| "Normal" heterogeneous (adjacent vs healthy donor) | Recorded as a limitation; TSS-disjoint control isolates part of it |
| Tumour purity unmodelled (D11) | Stated limitation; promotion candidate if TSS split ≈ random split |
| GSEA runtime on large collections | prerank on Wald stat, 1000 permutations, cached; hallmark is primary |
| pydeseq2 non-convergence at 10v10 | Wrapped in try/except, failures logged and counted, never silently dropped |

## Success Criteria

The research succeeds if it delivers, **regardless of which way the hypothesis falls**:

1. A defensible pooled gene-level replication rate with a CI, per cancer type and pooled,
   at matched sample size, with the selection rule varied.
2. A pathway-level number computed on the same pairs with **matched** comparators and a
   random-set null, so that "pathways are more robust" is either supported or refuted with
   evidence rather than assumed.
3. A decomposition attributing non-replication to sampling versus heterogeneity via controls
   at identical n.
4. Every number reproducible from cached intermediates by rerunning `src/` end to end.

---
---

# Research Planning — Direction Budget

**Hypothesis under test.** Fewer than 40% of genes reported as differentially expressed
(FDR < 0.05, |log2FC| > 1) in one study replicate at the same thresholds in an independent
study of the same cancer type, but pathway-level enrichment (GSEA NES) shows > 70%
directional concordance — i.e. pathway analysis is more robust than gene-level reporting.

## Enumerated candidate directions and scoring

Scored 1–5 on: **Ev** = literature evidence that the direction is informative,
**Rel** = relevance to the hypothesis as stated, **IG** = expected information gain,
**Feas** = implementation feasibility with the resources gathered.

| # | Direction | Ev | Rel | IG | Feas | Total | Verdict |
|---|-----------|----|----|----|------|-------|---------|
| D1 | Cross-cohort gene-level replication rate on uniformly processed counts (recount3), TCGA + independent SRA cohorts, many cancer types | 5 | 5 | 5 | 5 | **20** | **KEEP** |
| D2 | Matched pathway-level replication (GSEA NES sign + significance) on the *same* cohort pairs, with a like-for-like gene-level comparator | 5 | 5 | 5 | 5 | **20** | **KEEP** |
| D3 | Decomposition of replication failure: sampling/power vs. study heterogeneity, via within-cohort resampled splits and TSS-stratified splits as a ceiling control | 5 | 4 | 5 | 5 | **19** | **KEEP** |
| D4 | Effect of the selection rule (FDR-only, FC-only, joint FDR+FC, top-N by rank) on replication rate | 5 | 3 | 4 | 5 | 17 | Fold into D1 as a sensitivity axis, not a separate direction |
| D5 | Pipeline heterogeneity: re-quantify raw FASTQs with different aligners/quantifiers to add pipeline variance | 4 | 3 | 4 | 1 | 12 | **PRUNE** — needs petabyte-scale FASTQ download + alignment; recount3 deliberately holds pipeline constant, which is the cleaner design |
| D6 | Microarray-vs-RNA-seq cross-platform replication | 4 | 2 | 3 | 2 | 11 | **PRUNE** — hypothesis is explicitly about *bulk RNA-seq*; adds a platform confound |
| D7 | Single-cell pseudobulk vs bulk replication | 2 | 1 | 3 | 2 | 8 | **PRUNE** — out of scope |
| D8 | Literature-mining of *reported* DEG lists from published cancer papers and checking them against fresh data | 3 | 4 | 4 | 2 | 13 | **PRUNE** — supplementary-table extraction is unreliable and non-uniform; reanalysis from counts is strictly better evidence |
| D9 | Prognostic/survival signature replication | 4 | 1 | 3 | 3 | 11 | **PRUNE** — a different claim (prediction, not differential expression) |
| D10 | Non-cancer disease replication as a contrast | 3 | 2 | 3 | 3 | 11 | **PRUNE** — hypothesis scopes to cancer |
| D11 | Cell-composition/tumour-purity confounding as a driver of discordance | 4 | 2 | 4 | 3 | 13 | **PRUNE** as a direction; keep as a covariate/caveat in D3 |
| D12 | DE-tool choice (DESeq2 vs edgeR vs limma-voom) as a source of discordance | 4 | 2 | 3 | 3 | 12 | **PRUNE** — Soneson & Delorenzi and Rapaport et al. already establish tool concordance is high relative to cohort effects |

## Retained directions (top 3)

### D1 — Gene-level replication across independent cohorts
For every ordered pair of independent cohorts of the same cancer type, run DESeq2
tumour-vs-normal on both, restrict to commonly tested genes, and measure the fraction of
discovery-cohort DE genes (FDR < 0.05, |log2FC| > 1) that are also DE at the same thresholds
in the replication cohort. Report per-pair rates and a meta-analytic pooled rate; test
against the pre-registered 40% threshold. Sensitivity axis (from D4): repeat under
FDR-only, FC-only, joint, and top-N selection rules, since Shi et al. (2008) show the
selection rule alone moves overlap from ~20–40% to ~90%.

### D2 — Pathway-level concordance on the same cohort pairs
Run pre-ranked GSEA (hallmark, KEGG, Reactome, GO-BP) on the same DE results and measure
directional (NES sign) concordance and both-significant replication. **Critical design
constraint** (from Geistlinger et al. 2021): pathway "concordance" is not comparable to gene
"replication" unless the comparators are matched. Report the pathway numbers alongside the
*gene-level sign concordance* and a *gene-set-sized random-set null*, so the 70% claim is
tested against a proper baseline rather than against a quantity that is trivially high.

### D3 — Decomposition: sampling noise vs. genuine study heterogeneity
Establish a replication ceiling by splitting single large cohorts (TCGA) into disjoint
halves — both at random (pure sampling noise) and stratified by tissue-source-site (adds
collection-centre, population, and batch differences while holding pipeline and protocol
fixed). The gap between random-split, TSS-split, and genuinely-independent-study replication
rates attributes the shortfall to statistical power versus between-study heterogeneity. This
is what turns a descriptive number into an explanation.

## Notes on pruning
Directions were pruned for one of three reasons: (a) infeasible at this compute/bandwidth
scale (D5), (b) outside the hypothesis as stated (D6, D7, D9, D10), or (c) already settled in
the literature or better handled as a covariate/sensitivity axis inside a retained direction
(D4, D11, D12). If experimentation shows the retained directions cannot separate the effects
— for instance if the TSS-stratified split turns out to be indistinguishable from a random
split — D11 (composition confounding) is the first candidate for promotion, and the ranking
change will be recorded in STATE.md.
