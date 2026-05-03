# Epigenetic Age & DNA Methylation Analysis

> A two-part bioinformatics assignment exploring DNA methylation from complementary perspectives:
> **Part 1** implements a complete Whole Genome Bisulfite Sequencing (WGBS) pipeline on the Galaxy platform to characterize breast cancer methylomes, while **Part 2** benchmarks eight state-of-the-art epigenetic aging clocks across two independent EPIC array datasets using the Bio-Learn Python library.

---

## Table of Contents

- [Project Overview](#project-overview)
- [Part 1 — WGBS DNA Methylation Analysis (Galaxy)](#part-1--wgbs-dna-methylation-analysis-galaxy)
  - [Scientific Background](#scientific-background)
  - [Datasets & Samples](#datasets--samples)
  - [Step-by-Step Workflow](#step-by-step-workflow)
    - [Step 1: Data Upload & Quality Control](#step-1-data-upload--quality-control)
    - [Step 2: Bisulfite-Aware Alignment](#step-2-bisulfite-aware-alignment)
    - [Step 3: Methylation Bias Assessment](#step-3-methylation-bias-assessment)
    - [Step 4: Methylation Extraction](#step-4-methylation-extraction)
    - [Step 5: Visualization around CpG Islands](#step-5-visualization-around-cpg-islands)
    - [Step 6: Differential Methylation with Metilene](#step-6-differential-methylation-with-metilene)
  - [Key Findings](#key-findings)
- [Part 2 — Biological Clock Benchmarking (Bio-Learn)](#part-2--biological-clock-benchmarking-bio-learn)
  - [Scientific Background](#scientific-background-1)
  - [Datasets Used](#datasets-used)
  - [Aging Clocks & Models](#aging-clocks--models)
  - [Implementation Details](#implementation-details)
  - [Analyses & Visualizations](#analyses--visualizations)
    - [1. Clock Correlation Matrix](#1-clock-correlation-matrix)
    - [2. Age Deviation Heatmap](#2-age-deviation-heatmap)
    - [3. Predicted vs. Chronological Age Scatter](#3-predicted-vs-chronological-age-scatter)
  - [Cross-Dataset Observations](#cross-dataset-observations)
- [Repository Structure](#repository-structure)
- [Dependencies & Setup](#dependencies--setup)
- [References](#references)

---

## Project Overview

DNA methylation is one of the most studied epigenetic mechanisms, playing a central role in gene regulation, development, and disease. The addition of a methyl group to cytosine bases — predominantly at CpG dinucleotides — can silence gene expression, and its dysregulation is a hallmark of both cancer and aging. This project approaches DNA methylation from two angles that together span detection, quantification, functional interpretation, and biological age estimation.

| | Part 1 | Part 2 |
|---|---|---|
| **Platform** | Galaxy (usegalaxy.org) | Google Colab (Python) |
| **Data type** | WGBS sequencing reads (FASTQ → BAM) | EPIC/450K methylation arrays (β-values) |
| **Core task** | Full WGBS pipeline & DMR detection | Multi-clock biological age benchmarking |
| **Key tools** | bwameth, MethylDackel, Metilene, deepTools | Bio-Learn, pandas, seaborn, scipy |
| **Biological focus** | Breast cancer epigenomics | Epigenetic aging & clock comparison |
| **Reference study** | Lin et al. (2015) | Hannum et al. (2013) + multiple clock papers |

---

## Part 1 — WGBS DNA Methylation Analysis (Galaxy)

### Scientific Background

**What is DNA methylation?** In mammals, DNA methylation occurs almost exclusively at the cytosine of CpG dinucleotides. These methyl marks are heritable through cell division and serve as epigenetic regulators — when present at gene promoters, methylation typically represses transcription. CpG islands (CGIs), which are CpG-dense regions found at roughly 70% of mammalian promoters, are usually unmethylated in normal cells, keeping associated genes accessible for expression.

**Why is cancer relevant?** Cancer cells undergo widespread epigenetic reprogramming: global hypomethylation destabilizes the genome, while localized hypermethylation at tumour suppressor gene promoters silences critical regulators of cell cycle, apoptosis, and differentiation. Identifying differentially methylated regions (DMRs) between normal and tumour tissue can reveal epigenetically silenced cancer genes and potential biomarkers.

**Why can't we use normal NGS alignment?** Bisulfite sequencing chemically converts **unmethylated cytosines** to uracil (read as thymine after PCR), while methylated cytosines remain unchanged. This means the bisulfite-treated genome is no longer a simple two-strand complement — standard aligners that assume Watson-Crick pairing fail catastrophically. Bisulfite-aware aligners like **bwameth** handle the four-strand complexity (OT, OB, CTOT, CTOB strands) correctly.

This Part 1 analysis follows the [Galaxy Training Network WGBS tutorial](https://training.galaxyproject.org/training-material/topics/epigenetics/tutorials/methylation-seq/tutorial.html) based on the study by **Lin et al. (2015)**, which used hierarchical clustering of breast cancer methylomes to identify differentially methylated and expressed cancer genes.

---

### Datasets & Samples

Sequencing data used in this tutorial is a subsetted version of the original Lin et al. study, made available through [Zenodo record 557099](https://zenodo.org/record/557099). Five biologically distinct breast tissue types are represented, spanning the spectrum from healthy to malignant:

| Sample Label | Tissue Type | Biological Context |
|---|---|---|
| **NB1, NB2** | Normal breast cells | Healthy reference; expected canonical methylation patterns with unmethylated CGIs at active promoters |
| **BT089** | Fibroadenoma | Benign (non-cancerous) breast tumour; used to assess early-stage epigenetic changes |
| **BT126, BT198** | Invasive ductal carcinoma | Malignant breast cancer; expected widespread methylation reprogramming |
| **MCF7** | Breast adenocarcinoma cell line | Established cancer model; useful for comparing in vitro vs. in vivo methylation |

> The tutorial uses a small genomic subset (`subset_1.fastq`, `subset_2.fastq`) for compute feasibility. Pre-computed alignments and converted bedGraph files are provided for downstream visualization steps.

---

### Step-by-Step Workflow

```
Raw FASTQ reads
      │
      ▼
 Quality Control (Falco)
      │  ── assess base quality, GC distribution, T→C skew from bisulfite treatment
      ▼
 Alignment (bwameth → hg38)
      │  ── bisulfite-aware aligner; outputs sorted BAM
      ▼
 Methylation Bias Assessment (MethylDackel mbias)
      │  ── diagnostic SVG plots; identify positional bias for trimming
      ▼
 Methylation Extraction (MethylDackel extract)
      │  ── CpG methylation fractions (bedGraph); merge per-cytosine metrics
      ▼
 Format Conversion (BedGraph → bigWig)
      │  ── required for deepTools
      ▼
 Visualization (computeMatrix → plotProfile)
      │  ── methylation level centred on CpG islands / TSS
      ▼
 Differential Methylation (Metilene)
      └── DMR detection between normal (NB1, NB2) vs. cancer (BT198)
```

#### Step 1: Data Upload & Quality Control

**Tool:** Falco (v1.2.4+galaxy0) — an optimised reimplementation of FastQC

Quality control is the mandatory first step for any NGS analysis. Falco assesses read quality metrics including per-base quality scores, GC content, adapter contamination, and per-base sequence composition.
<img width="760" height="463" alt="image" src="https://github.com/user-attachments/assets/fad461a0-2106-40e8-99bd-2524407978aa" />


**The bisulfite signature in QC:** The most striking and expected finding in bisulfite sequencing QC is the severe imbalance in base composition. Because bisulfite treatment converts unmethylated cytosines to thymine (via uracil), reads show a dramatic **reduction in cytosine (C) frequency** and a corresponding **spike in thymine (T) frequency**. In a normal NGS library, all four bases occur at roughly equal frequency; in a bisulfite library, the T/C ratio is heavily skewed. This is not a quality problem — it is **confirmation that bisulfite conversion was successful**. An abnormal GC distribution follows as a direct consequence of this chemistry.

---

#### Step 2: Bisulfite-Aware Alignment

**Tool:** bwameth (v0.2.7+galaxy0) | **Reference genome:** Human hg38full (GRCh38)

`bwameth` is a wrapper around BWA-MEM designed specifically for bisulfite-converted reads. It aligns both mates of a paired-end library to an in-silico bisulfite-converted reference, correctly handling all four strand orientations. The output is a sorted and indexed BAM file containing the read alignments with methylation information encodable by downstream tools.

**Why not STAR, HISAT2, or bowtie2?** These aligners assume complementary base pairing across strands. In a bisulfite-treated library, C→T conversions mean that a methylated read and an unmethylated read from the same locus look different, and the genome is no longer symmetric. bwameth (and similar tools like Bismark) handle this by converting both the reads and the reference in silico before alignment, ensuring all four strand types are correctly represented.
<img width="801" height="650" alt="image" src="https://github.com/user-attachments/assets/028da1dc-9ed3-4066-bd4c-62410a43c2fb" />

---

#### Step 3: Methylation Bias Assessment

**Tool:** MethylDackel (v0.5.2+galaxy0) | **Mode:** `mbias`

Before extracting methylation values, it is critical to check for **positional methylation bias** — systematic artifacts where the methylation level is artificially inflated or deflated at specific positions along the read (e.g., the first or last 5–10 bases). This commonly arises from end-repair during library preparation or non-bisulfite-converted bases at read ends.

MethylDackel generates diagnostic SVG plots showing average methylation across read positions for each strand (OT = original top, OB = original bottom, CTOT and CTOB = complementary strands). If a bias is present — visible as an uptick or downtick at the 5′ or 3′ ends — those positions should be excluded during extraction via trimming parameters. Settings used: singletons and discordant alignments were retained (`--keepSingleton`, `--keepDiscordant`) to maximize coverage on this small subset.
<img width="832" height="488" alt="image" src="https://github.com/user-attachments/assets/1e20e265-e901-4db3-9e8a-f7d4bed78524" />

---

#### Step 4: Methylation Extraction

**Tool:** MethylDackel (v0.5.2+galaxy0) | **Mode:** `extract` | **Output:** CpG methylation fractions (bedGraph)

After bias correction, methylation metrics are extracted from the BAM file. For each CpG site in the genome, MethylDackel counts methylated and unmethylated reads and reports the **methylation fraction** (proportion of reads showing methylation, 0–1). The `--mergeContext` flag merges data from both strands of a CpG dinucleotide into a single value — the standard approach for CpG-level analysis.

The output bedGraph format (chromosome, start, end, methylation fraction) is the standard input for downstream visualisation and differential analysis tools.
<img width="787" height="355" alt="image" src="https://github.com/user-attachments/assets/31ef353d-76ea-4a9e-893b-09dc0602a9c7" />

---

#### Step 5: Visualization around CpG Islands

**Tools:** Wig/BedGraph-to-bigWig → computeMatrix (v3.5.4) → plotProfile (v3.5.4) | all from deepTools

To understand the functional context of methylation patterns, profiles are plotted **centred on CpG islands** (a BED file of all annotated CGIs in hg38 was imported from Zenodo). The deepTools workflow proceeds as:

1. **BedGraph → bigWig conversion:** The methylation fraction bedGraph is converted to indexed bigWig format for efficient random-access querying by computeMatrix.
2. **computeMatrix (reference-point mode):** Computes a matrix of methylation values across a fixed window centred on each CGI midpoint, enabling averaging across all CGIs genome-wide.
3. **plotProfile:** Renders the average methylation profile as a line plot across the CGI-centred window.

**What a healthy profile looks like:** In normal cells, CpG islands are characteristically **hypomethylated** — methylation dips sharply at the CGI centre. This is consistent with their role as open, accessible promoter regions. The surrounding genomic context typically shows higher background methylation (global methylation of non-island CpGs is ~70–80% in somatic cells).

**Cross-sample comparison:** By running the same pipeline for all six samples and plotting together with `--perGroup`, differences in CGI methylation between normal and tumour samples become visually apparent. Cancer-associated CGI hypermethylation causes the characteristic dip to flatten or disappear at some loci, reflecting epigenetic silencing.

> **Note on chromosome naming:** Pre-computed bedGraph files use Ensembl-style chromosome names (e.g., `1`, `2`) rather than UCSC-style (`chr1`, `chr2`) required by hg38 bigWig tools. The `Replace column` tool was used to remap chromosome names via a GRCh38 Ensembl→UCSC mapping file before conversion.

---

#### Step 6: Differential Methylation with Metilene

**Tool:** Metilene (v0.2.6.1)

Metilene identifies **differentially methylated regions (DMRs)** between two groups of samples by scanning for contiguous CpG sites with consistently different methylation levels. It uses a binary segmentation algorithm combined with a non-parametric statistical test (Mann-Whitney U), making it robust without assumptions of normality.

**Groups compared:**
- **Group 1 (Normal):** NB1_CpG.meth.bedGraph, NB2_CpG.meth.bedGraph
- **Group 2 (Tumour):** BT198_CpG.meth.bedGraph

Analysis was restricted to annotated CpG island regions (`CpGIslands.bed`) to focus on functionally relevant methylation changes. The output includes a BED file of DMRs with associated mean methylation differences and q-values, plus a PDF summary showing the distribution of DMR sizes, methylation differences, and significance across the genome.

---

### Key Findings

**Bisulfite conversion confirmed:** QC showed the expected dramatic shift in T/C base composition, validating successful library preparation.

**Positional bias detected:** The mbias SVG plot showed elevated methylation at the 5′ end of OT strand reads — a common artifact from end-repair. Trimming the first ~5 bases during extraction removes this bias.

**Canonical CGI hypomethylation in normal tissue:** deepTools profile plots show the expected methylation dip centred on CpG islands in normal breast samples (NB1, NB2), consistent with active, unmethylated promoters.

**Cancer-specific epigenetic reprogramming:** Tumour-derived samples show altered profiles at CGIs, with some islands shifting toward hypermethylation — consistent with epigenetic silencing of tumour suppressors as reported in Lin et al.

**DMRs between normal and cancer:** Metilene identified significant DMRs at CpG islands distinguishing normal breast from BT198 invasive ductal carcinoma. These regions represent candidate loci whose gene expression may be altered through methylation-driven silencing.

---

## Part 2 — Biological Clock Benchmarking (Bio-Learn)

### Scientific Background

**The epigenetic clock concept:** In 2013, Steve Horvath discovered that DNA methylation levels at a small set of CpG sites could predict a person's chronological age with remarkable accuracy across diverse tissues. These "epigenetic clocks" exploit the fact that methylation patterns change in a consistent, programmatic manner with age. Since then, dozens of clocks have been developed, each trained on different datasets and optimised for different biological outcomes — from chronological age to disease risk to pace of aging.

**Why compare multiple clocks?** No single clock captures all aspects of aging. Some measure calendar time, others measure biological or phenotypic age (how fast the body ages relative to the calendar), and others measure the *rate* of aging rather than an absolute age. Benchmarking multiple clocks on the same datasets reveals their agreement and divergence, tests their generalizability across cohorts, and helps researchers choose the most appropriate clock for a given biological question.

**Bio-Learn** is an open-source Python library providing a unified interface to a curated gallery of published aging clocks and standardised methylation datasets from GEO. It handles CpG imputation, normalisation, and prediction through a consistent API, making multi-clock comparisons straightforward and reproducible.

---

### Datasets Used

#### Dataset 1 — GSE40279

| Property | Value |
|---|---|
| **GEO Accession** | GSE40279 |
| **Study** | Hannum et al. (2013) |
| **Tissue** | Whole blood |
| **Platform** | Illumina 450K Methylation Array |
| **Original dimensions** | 656 samples × 473,034 CpG sites |
| **Samples used** | First 200 (memory-constrained subset) |
| **Biological significance** | One of the largest and most cited blood methylation aging datasets. The Hannum clock was trained on this data, making GSE40279 an "in-distribution" benchmark. Used as the primary reference dataset. |

#### Dataset 2 — GSE51057

| Property | Value |
|---|---|
| **GEO Accession** | GSE51057 |
| **Tissue** | Whole blood |
| **Platform** | Illumina 450K / EPIC Methylation Array |
| **Original dimensions** | 329 samples × 485,577 CpG sites |
| **Samples used** | First 200 |
| **Biological significance** | An independent whole-blood cohort from a different study population, providing genuine external validation. No clock in this benchmark was explicitly trained on GSE51057, so performance here reflects true generalization. Comparing results between the two datasets tests clock robustness across cohorts. |

---

### Aging Clocks & Models

Eight models were selected to represent the major conceptual archetypes in epigenetic clock research:

#### Horvath v1 (2013)
The original pan-tissue epigenetic clock. Trained on over 8,000 samples across 51 tissue and cell types using elastic net regression on **353 CpG sites**. It predicts chronological age with high accuracy across virtually any human tissue — a property no previous biomarker had achieved. Horvath v1 measures what is often called "intrinsic epigenetic age acceleration" when adjusted for cell-type composition. It remains the canonical baseline against which all newer clocks are compared.

#### Horvath v2 — SkinBloodClock (2018)
An update trained specifically to optimise performance in **skin fibroblasts and blood**, using **391 CpGs**. It corrects for systematic underestimation of age in older individuals — a known limitation of v1 — and is particularly relevant for aging intervention studies where blood is the primary material.

#### Hannum Clock (2013)
Developed by Hannum et al. on the GSE40279 dataset used in this analysis, using **71 CpGs** trained via lasso regression for whole blood. Because this clock was trained on Dataset 1, it is expected to perform well there — a useful "in-sample" reference point for evaluating how other clocks compare on the same data.

#### PhenoAge (2018)
Developed by Levine et al., PhenoAge represents a conceptual leap beyond chronological age. It was trained in two stages: first, a composite "phenotypic age" was computed from nine clinical blood biomarkers using a mortality model; then, DNA methylation was trained to predict this phenotypic age score. PhenoAge therefore captures **disease risk, morbidity, and biological aging** rather than calendar time. High PhenoAge acceleration (predicted > chronological) predicts increased risk of chronic disease and earlier mortality.

#### DunedinPACE (2022)
A fundamentally different type of clock — rather than predicting an absolute age, DunedinPACE measures the **speed of biological aging**. A value of ~1.0 is average; values above 1.0 mean aging faster, below 1.0 means aging slower. It was developed from the longitudinal Dunedin Study, where the same individuals were measured repeatedly over decades, enabling the *rate* of change across 19 physiological systems to be quantified. DunedinPACE is sensitive to lifestyle and environmental exposures and is more appropriate for intervention studies than absolute-age clocks.

#### YingCausAge (2024)
Part of a set of clocks that apply **causal inference** to epigenetic aging. Rather than correlating methylation with age, YingCausAge uses a causal framework to identify CpGs whose changes are *caused by* the aging process rather than merely associated with it. The goal is to separate true aging signal from confounders such as smoking, BMI, or technical batch effects, providing a more biologically interpretable measure of aging.

#### Knight Clock (2016)
Developed to estimate **gestational age at birth** from cord blood DNA methylation — not adult aging. Included here as an intentional cross-tissue control. When applied to adult whole blood, the Knight clock produces compressed, biologically uninformative predictions. This is an instructive demonstration of tissue and life-stage specificity: clocks trained for one biological context do not transfer to fundamentally different ones.

#### Bocklandt Clock (2011)
One of the earliest published epigenetic clocks, trained on **saliva samples** with a very small CpG feature set. In whole blood (as here), it outputs fractional values (~0.3–0.4) rather than years — a direct consequence of tissue mismatch. Bocklandt's inclusion provides historical context (showing how far clock development has advanced) and a second cross-tissue comparison alongside Knight.

---

### Implementation Details

The notebook (`biological_clocks.ipynb`) was executed in Google Colab with the following structure:

**1. Environment setup:** `biolearn`, `matplotlib`, `seaborn`, `pandas`, `numpy`, `scipy` installed via pip.

**2. Data loading:** Both datasets loaded via `DataLibrary().get(accession).load()`, which automatically downloads and caches the methylation matrix and associated metadata (including chronological age from `metadata['age']`). To manage Colab RAM limits, each dataset was immediately trimmed to the first 200 samples after loading, with explicit garbage collection (`gc.collect()`) between steps.

**3. Clock inference loop:** A `run_clocks()` function iterates over all eight clock names, instantiates each model via `ModelGallery().get(clock)`, calls `model.predict(data)`, and stores flattened predictions in a results dictionary. Failed predictions are caught via try/except and filled with NaN so downstream analysis is not interrupted. Chronological ages from `data.metadata['age']` are appended as a reference column.

**4. Compatibility testing:** Before finalizing the clock list, a test loop verified which Bio-Learn gallery clocks were compatible with 450K array data. Knight, Bocklandt, YingDamAge, YingAdaptAge, and DNAmTL were confirmed working; Vidalin, MEAT2, and Stubbs failed due to missing CpG coverage in the 450K platform and were excluded.

---

### Analyses & Visualizations

#### 1. Clock Correlation Matrix

**Output files:** `corr_GSE40279.png`, `corr_GSE51057.png`

**Method:** Pearson correlation coefficients computed pairwise across all clock predictions and chronological age. Visualized as an annotated heatmap (seaborn) with a diverging colormap.

**What it reveals:**
- Clocks estimating age from whole blood (Horvath v1/v2, Hannum, PhenoAge, YingCausAge) are expected to show strong mutual correlations, reflecting the shared methylation-age signal they exploit.
- DunedinPACE shows weaker correlation with age-based clocks — its pace-of-aging scale is conceptually orthogonal to absolute age.
- Knight and Bocklandt show low or near-zero correlation with other clocks and with chronological age, reflecting their tissue mismatch with whole blood.
- Differences between the GSE40279 and GSE51057 matrices indicate which clock relationships are robust across cohorts versus cohort-specific.
<img width="757" height="648" alt="image" src="https://github.com/user-attachments/assets/9d197957-c313-43c9-a115-8d7cda3c7b47" />
<img width="758" height="657" alt="image" src="https://github.com/user-attachments/assets/aa43287a-f16d-4824-80ab-3de4789d45b2" />

---

#### 2. Age Deviation Heatmap

**Output files:** `deviation_GSE40279.png`, `deviation_GSE51057.png`

**Method:** For each age-predicting clock (DunedinPACE and Bocklandt excluded as non-age-scale outputs), the **age deviation** is computed per sample as `predicted_age − chronological_age`. Positive values indicate biological over-aging; negative values indicate under-aging. Deviations are assembled into a samples × clocks matrix and visualized as a clustermap with hierarchical clustering on both axes.

**What it reveals:**
- Samples consistently showing high positive deviation across multiple clocks may be biologically older than their chronological age — potentially linked to disease, lifestyle, or genetic factors.
- Clock-wise clustering groups clocks producing similar deviation patterns, offering an unsupervised view of clock similarity.
- Systematic bias (all clocks overestimating in one dataset vs. another) can indicate cohort-level biological age differences or distributional shifts from clock training data.
- Extreme outlier samples are visually identifiable and could be flagged for quality review or further biological investigation.
<img width="767" height="257" alt="image" src="https://github.com/user-attachments/assets/bb941924-a30a-4a5e-bde1-5b085ee59cbc" />

---

#### 3. Predicted vs. Chronological Age Scatter

**Output files:** `scatter_GSE40279.png`, `scatter_GSE51057.png`

**Method:** For each clock, a scatter plot of predicted age (y-axis) vs. chronological age (x-axis) with an identity diagonal line (perfect prediction). Pearson r and p-value are computed via `scipy.stats.pearsonr` and annotated on each subplot. All eight clocks shown in a grid for direct comparison.

**What it reveals:**
- Clocks with points tightly clustered around the diagonal have high accuracy and low systematic bias.
- A high r but consistently offset line indicates a calibration shift (systematic bias) rather than an inability to track age.
- DunedinPACE shows a flat scatter because its output is pace (~1.0 for most healthy adults) — the scatter should be interpreted differently from age-predicting clocks.
- Knight and Bocklandt show near-random scatter against chronological age, confirming tissue mismatch and expected non-transferability.
<img width="762" height="338" alt="image" src="https://github.com/user-attachments/assets/f7fb55cd-0dda-4e01-b89e-1fdcc4be7aa2" />

---

### Cross-Dataset Observations

Comparing results between GSE40279 and GSE51057 provides insight into clock generalizability:

- **Clocks performing consistently across both datasets** (high r in both scatter plots, similar correlation structure) are robust and likely to generalize to new cohorts.
- **Clocks performing well on GSE40279 but degrading on GSE51057** may be overfit to the training-adjacent distribution of the Hannum dataset.
- **Age deviation patterns replicating across both datasets** suggest genuine biological age differences rather than dataset-specific artifacts.
- This cross-cohort comparison underlines the importance of external validation — a clock that appears excellent in a single dataset may not generalize.

---

## Repository Structure

```
.
├── README.md                          # This file
├── biological_clocks.ipynb            # Google Colab notebook — Part 2 full analysis
│
├── outputs/                           # Generated figures from Part 2
│   ├── corr_GSE40279.png              # Correlation matrix — Dataset 1
│   ├── corr_GSE51057.png              # Correlation matrix — Dataset 2
│   ├── deviation_GSE40279.png         # Age deviation heatmap — Dataset 1
│   ├── deviation_GSE51057.png         # Age deviation heatmap — Dataset 2
│   ├── scatter_GSE40279.png           # Predicted vs. chronological age — Dataset 1
│   └── scatter_GSE51057.png           # Predicted vs. chronological age — Dataset 2
│
└── galaxy/                            # Part 1 reference materials
    └── methylation-seq-tutorial.html  # GTN tutorial offline reference
```

---

## Dependencies & Setup

### Part 1 — Galaxy Platform

All tools run on the public Galaxy instance at [usegalaxy.org](https://usegalaxy.org/). No local installation required.

| Tool | Galaxy Version | Purpose |
|------|---------------|---------|
| Falco | 1.2.4+galaxy0 | Read quality control |
| bwameth | 0.2.7+galaxy0 | Bisulfite-aware alignment to hg38 |
| MethylDackel | 0.5.2+galaxy0 | Methylation bias assessment and extraction |
| Replace column | 0.2 | Chromosome name remapping (Ensembl → UCSC) |
| Wig/BedGraph-to-bigWig | — | Format conversion for deepTools |
| computeMatrix | 3.5.4+galaxy0 | Methylation matrix around CGI regions |
| plotProfile | 3.5.4+galaxy0 | Profile visualization |
| Metilene | 0.2.6.1 | Differentially methylated region detection |

**Reference genome:** GRCh38 / hg38 (full build for alignment; standard build for extraction)

---

### Part 2 — Python (Google Colab)

**Quick start:**
```bash
pip install biolearn matplotlib seaborn pandas numpy scipy
```

Open `biological_clocks.ipynb` in [Google Colab](https://colab.research.google.com/) and run all cells sequentially. Bio-Learn automatically downloads and caches datasets from GEO on first access.

| Package | Version | Role |
|---------|---------|------|
| biolearn | 0.9.1 | Dataset loading, clock model gallery, prediction |
| pandas | 2.2.2 | DataFrame manipulation, results storage |
| numpy | 2.0.2 | Numerical operations, NaN handling |
| matplotlib | 3.10.0 | Plotting backend |
| seaborn | 0.13.2 | Heatmaps, scatter plots, correlation visualisation |
| scipy | 1.16.3 | Pearson correlation, statistical annotation |
| scikit-learn | 1.6.1 | Dependency of biolearn |
| cvxpy / ecos / osqp | various | Optimisation dependencies of biolearn |

**Estimated runtime:** 10–25 minutes in Colab depending on download speed. Memory peaks at ~6–8 GB when both large methylation matrices are loaded; the 200-sample trimming step substantially reduces this.

---

## References

1. **Lin, I.-H., Chen, D.-T., Chang, Y.-F., Lee, Y.-L., Su, C.-H. et al. (2015).** Hierarchical Clustering of Breast Cancer Methylomes Revealed Differentially Methylated and Expressed Breast Cancer Genes. *PLOS ONE* 10(1): e0118453. https://doi.org/10.1371/journal.pone.0118453

2. **Hannum, G., Guinney, J., Zhao, L., Zhang, L., Hughes, G. et al. (2013).** Genome-wide Methylation Profiles Reveal Quantitative Views of Human Aging Rates. *Molecular Cell* 49(2): 359–367. https://doi.org/10.1016/j.molcel.2012.10.016

3. **Horvath, S. (2013).** DNA methylation age of human tissues and cell types. *Genome Biology* 14: R115. https://doi.org/10.1186/gb-2013-14-10-r115


