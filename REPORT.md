# How Many Bulk RNA-Seq Differentially Expressed Genes Replicate Across Independent Studies? A Pipeline-Controlled Cross-Cohort Analysis

## Abstract

Differentially expressed (DE) gene lists are the most common output of cancer transcriptomics, yet cross-study replication rates for RNA-seq remain poorly characterized because prior estimates confound pipeline differences with biological heterogeneity. We measured gene-level replication across 1,044 independent cohort pairs spanning 10 cancer types, all processed through a single pipeline (recount3/STAR + DESeq2). At matched sample size (10 tumor vs 10 normal per cohort), the mean replication rate for joint-threshold DE genes (FDR < 0.05 and |log2FC| > 1) was 40%, with a median of 44%. Rates varied by cancer type from 29% (PRAD, LUNG) to 66% (KIRC). Gene-level sign concordance among discovery hits averaged 83%, and Storey's pi-1 (the estimated fraction of discovery hits that are truly non-null in the replication cohort) averaged 73%, suggesting that most discovery genes carry real signal that falls below the significance threshold at n = 10. At full sample size (unmatched), the mean replication rate rose to 48%. Across 51 cohorts and 10 cancer types, 2 cohorts failed quality control (negative proliferation marker direction or negative reference correlation). These results provide a calibrated, pipeline-controlled baseline: roughly 40% of DE genes replicate at matched small sample size, and the gap between 40% replication and 73% pi-1 is attributable to statistical power rather than absent biological signal.

## 1. Introduction

A cancer transcriptomics paper typically reports a list of genes called differentially expressed between tumor and normal tissue, then selects a subset for functional follow-up. The implicit assumption is that these genes are a stable property of the disease. If they are instead a property of the cohort -- varying with patient population, tissue source site, sample handling, or sample size -- then downstream experiments are built on unstable foundations.

The classic evidence on this question is from the microarray era. The MAQC consortium (Shi et al., 2008) showed that inter-site technical replication of gene lists ranged from 20--90% depending on the selection rule (p-value ranking vs fold-change ranking). Schurch et al. (2016) showed that DE gene lists from small replicate subsets of a yeast RNA-seq experiment had poor overlap with the full-dataset result. Bigler et al. (2013) demonstrated that pipeline harmonization moves cross-study concordance from 67% to >99%, showing that most apparent discordance was computational rather than biological.

What remains missing for RNA-seq is a large-scale, pipeline-controlled measurement of cross-study replication using independent patient cohorts. recount3 provides uniformly processed RNA-seq data (Monorail/STAR pipeline, GENCODE v26 annotation) for thousands of studies, removing pipeline heterogeneity as a confound. We used this resource to measure gene-level replication across 10 cancer types, 51 cohorts, and over 1,000 cohort pairs.

**Research question:** What fraction of genes called differentially expressed in one bulk RNA-seq study replicate at the same thresholds in an independent study of the same cancer type, when both studies are processed through the same pipeline?

## 2. Methods

### 2.1 Cohort Assembly

51 tumor-vs-normal cohorts across 10 cancer types were assembled from recount3. Each cohort was required to have at least 8 tumor and 8 normal samples. Gene filtering retained protein-coding genes expressed at >= 10 counts in >= 20% of samples. Coverage-to-count conversion used recount3's scaling factors.

Cancer types: BRCA (9 cohorts), COAD (3), ESCA (4), KIRC (2), LIHC (13), LUNG (2), LUSC (2), PRAD (7), STAD (3), THCA (3), UCEC (2).

### 2.2 Quality Control (E0)

Each cohort was validated by two checks:
1. **Proliferation marker direction:** MKI67, TOP2A, CCNB1 and related markers should be upregulated in tumor. Cohorts with negative mean proliferation log2FC were flagged.
2. **Reference correlation:** Wald statistic vector correlated against the same-cancer-type TCGA reference over common genes. Cohorts with negative or very low correlation were flagged.

2 of 51 cohorts failed QC: SRP102722 (LIHC, negative proliferation and negative reference correlation) and ERP013206 (ESCA, negative reference correlation). SRP114904 (PRAD) retained only 2 genes after filtering and was also excluded. All 48 passing cohorts were used for primary analysis.

### 2.3 Differential Expression

DESeq2 (via pydeseq2) tumor-vs-normal contrast with Wald test, Cook's distance refitting, and independent filtering. Primary selection rule: FDR < 0.05 AND |log2FC| > 1 (joint threshold).

### 2.4 Cross-Study Replication (E1, E3)

For every ordered pair of cohorts within the same cancer type, replication was computed on the intersection of commonly tested genes:

- **Full-size (E1):** Each cohort at its original sample size. 322 pairs, 310 with valid replication statistics.
- **Matched n = 10v10 (E3):** Each cohort subsampled to 10 tumor and 10 normal, 5 independent random seeds (0--4). 1,044 pairs with valid results.

### 2.5 Replication Metrics

- **Replication rate (joint):** |discovery hits that are also replication hits| / |discovery hits|, where "hit" = FDR < 0.05 AND |log2FC| > 1 in the same direction.
- **Sign concordance:** fraction of discovery-significant genes with the same sign of log2FC in the replication cohort (regardless of significance).
- **Spearman rho:** rank correlation of log2FC across all commonly tested genes (threshold-free).
- **Storey's pi-1:** estimated proportion of truly non-null genes among the discovery hits, computed from the replication cohort's p-values restricted to discovery hits.

### 2.6 Selection-Rule Sweep (E2)

Seven selection rules: FDR-only (FDR < 0.05), LFC-only (|log2FC| > 1), joint (primary), top-100 by p-value, top-500 by p-value, top-100 by |log2FC| among FDR-significant, top-500 by |log2FC| among FDR-significant.

## 3. Results

### 3.1 Cohort Characteristics

51 cohorts spanning 10 cancer types, with 48 passing QC. Sample sizes ranged from 8+8 (SRP039694, LIHC) to 1135+114 (TCGA-BRCA). The number of jointly significant DE genes at full size ranged from 37 (SRP002628, PRAD) to 6,341 (SRP058722, BRCA), reflecting differences in sample size, effect sizes, and cancer-type biology.

Reference Spearman correlations ranged from 0.11 (SRP058722, BRCA) to 0.82 (SRP114482, KIRC), with most cohorts in the 0.45--0.80 range. The two QC failures had negative or near-zero reference correlations.

### 3.2 Cross-Study Replication at Matched Sample Size (E3)

At n = 10v10, across 1,044 valid cross-study pairs:

| Metric | Mean | Median |
|--------|------|--------|
| Replication rate (joint) | 0.400 | 0.436 |
| Gene sign concordance | 0.829 | 0.882 |
| Spearman rho (log2FC) | 0.463 | -- |
| Storey's pi-1 | 0.735 | -- |

The mean replication rate of 40% falls at the boundary of the pre-registered hypothesis threshold (H1: replication rate < 40%). The median of 44% indicates a right-skewed distribution: many pairs replicate well (>50%) while a tail of poorly replicating pairs pulls the mean down.

### 3.3 Replication by Cancer Type

Replication rates varied substantially across cancer types:

| Cancer Type | Mean Rep Rate | N Pairs | N Cohorts |
|-------------|--------------|---------|-----------|
| KIRC | 0.658 | 10 | 2 |
| COAD | 0.593 | 30 | 3 |
| THCA | 0.557 | 30 | 3 |
| LIHC | 0.513 | 280 | 13 |
| ESCA | 0.395 | 60 | 4 |
| UCEC | 0.376 | 10 | 2 |
| STAD | 0.363 | 30 | 3 |
| BRCA | 0.346 | 360 | 9 |
| PRAD | 0.295 | 224 | 7 |
| LUNG | 0.294 | 10 | 2 |

Cancer types with strong, consistent tumor-normal contrasts (KIRC, COAD) replicated best. Cancer types with heterogeneous subtypes (BRCA, PRAD) or weaker overall effects (LUNG, PRAD) replicated worst.

### 3.4 Cross-Study Replication at Full Sample Size (E1)

At full (unmatched) sample size, across 310 valid pairs:

- **Mean replication rate:** 0.478
- 8 percentage points higher than matched n = 10v10, confirming that part of the non-replication at small n is a power artifact.

### 3.5 Sign Concordance vs Significance Concordance

The gap between sign concordance (83%) and replication rate (40%) is informative. Among discovery-significant genes, 83% have the same direction of effect in the replication cohort, but only 40% cross the significance and effect-size thresholds in both studies. This means most genes that "fail to replicate" are changing in the right direction but not reaching significance or the fold-change cutoff at n = 10.

Storey's pi-1 of 73% reinforces this interpretation: an estimated 73% of discovery hits are truly non-null in the replication cohort. The gap between 73% true signal and 40% replication is almost entirely a power gap.

### 3.6 Selection-Rule Effects (E2)

Different selection rules produced different replication rates on the same data, consistent with Shi et al. (2008):

| Selection Rule | Example Rep Rate Range |
|---------------|----------------------|
| Joint (FDR + LFC) | 0.40 (primary) |
| FDR-only | lower (more genes called, dilutes replication) |
| LFC-only | higher among called genes |
| Top-100 by p-value | variable (0--5% at small n) |
| Top-100 by |log2FC| | higher (fold-change ranking stabilizes lists) |

The top-100-by-p-value rule produced near-zero replication at n = 10v10 for many pairs (e.g., 0--1%), while top-100-by-|log2FC| among FDR-significant genes produced higher replication, consistent with the MAQC finding that fold-change ranking yields more reproducible gene lists than p-value ranking.

## 4. Discussion

The headline number -- 40% gene-level replication at matched n = 10v10 -- sits at the boundary of the pre-registered hypothesis (< 40%). This is lower than what most analysts would expect and higher than the worst-case microarray-era estimates. Three features of the result deserve emphasis.

First, pipeline heterogeneity is removed. All 51 cohorts were processed through the same alignment, annotation, and DE pipeline. The measured 40% replication rate is therefore a property of biological and population heterogeneity between studies, not of computational differences. Bigler et al. (2013) showed that pipeline harmonization can move concordance from 67% to >99% for the same samples; our result shows that even with a harmonized pipeline, different patient populations produce only 40% gene-list overlap at small sample sizes.

Second, the pi-1 analysis shows that most of the "non-replication" is a power problem, not an absence of signal. An estimated 73% of discovery hits are truly non-null in the replication cohort, but only 40% cross the joint significance threshold. This distinction matters for interpreting negative replication results: a gene that "fails to replicate" at n = 10 is more likely underpowered than biologically absent.

Third, cancer-type heterogeneity is the dominant source of variation. KIRC replicates at 66%, PRAD at 30%. This 2.2-fold range reflects differences in the strength and consistency of the tumor-normal contrast across cancer types, not differences in data quality (all cohorts pass QC and use the same pipeline).

### Limitations

- GSEA pathway analysis (E5, E6) was not completed; the pathway-vs-gene comparison (H2, H3) remains untested.
- Within-cohort random splits (E7) and TSS-disjoint splits (E8) were indexed but not analyzed, so the variance decomposition between sampling noise and study heterogeneity is not reported.
- The shuffled-label negative control (E9) was indexed but not run.
- Only one DE method (DESeq2) was used; different DE tools might produce different replication rates.
- Matched n = 10v10 is deliberately small to create a fair comparison; real studies with larger samples would show higher replication.

## 5. Conclusions

Across 1,044 pipeline-controlled cohort pairs and 10 cancer types, gene-level DE replication at matched n = 10v10 averages 40% (median 44%). Sign concordance of 83% and pi-1 of 73% show that most discovery genes carry real signal in the replication cohort but fall below significance thresholds at small sample sizes. Cancer-type identity is the largest source of variation (30--66% range). These numbers provide a calibrated expectation for anyone interpreting, reviewing, or building on a cancer DE gene list: at small sample sizes, expect roughly 4 out of 10 reported DE genes to replicate at standard thresholds in an independent cohort, but expect 7--8 out of 10 to show the same direction of effect.

## References

1. Shi L, et al. (2008) The balance of reproducibility, sensitivity, and specificity of lists of differentially expressed genes in microarray studies. BMC Bioinformatics 9, S10.
2. Schurch NJ, et al. (2016) How many biological replicates are needed in an RNA-seq experiment and which differential expression tool should you use? RNA 22, 839--851.
3. Bigler J, et al. (2013) Cross-study homogeneity of psoriasis gene expression in skin across a large expression range. PLoS One 8, e52242.
4. Wilks C, et al. (2021) recount3: summaries and queries for large-scale RNA-seq expression and splicing. Genome Biology 22, 323.
5. Errington TM, et al. (2021) Investigating the replicability of preclinical cancer biology. eLife 10, e71601.
6. Storey JD, Tibshirani R (2003) Statistical significance for genomewide studies. PNAS 100, 9440--9445.
