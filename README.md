# Epigenetics-Aging-biological-clock-Analysis
# Epigenetic Age & DNA Methylation Analysis

> A two-part bioinformatics assignment combining **WGBS methylation analysis** (Galaxy platform) with **multi-clock biological age benchmarking** (Python / Bio-Learn).

---

## Table of Contents

- [Project Overview](#project-overview)
- [Part 1 — WGBS DNA Methylation Analysis (Galaxy)](#part-1--wgbs-dna-methylation-analysis-galaxy)
  - [Background](#background)
  - [Datasets](#datasets)
  - [Workflow Summary](#workflow-summary)
  - [Key Results](#key-results)
- [Part 2 — Biological Clock Benchmarking (Bio-Learn)](#part-2--biological-clock-benchmarking-bio-learn)
  - [Datasets Used](#datasets-used)
  - [Aging Clocks & Models](#aging-clocks--models)
  - [Analyses Performed](#analyses-performed)
  - [Visualizations](#visualizations)
- [Repository Structure](#repository-structure)
- [Dependencies & Setup](#dependencies--setup)
- [References](#references)

---

## Project Overview

This assignment explores DNA methylation from two complementary angles:

| Part | Platform | Focus |
|------|----------|-------|
| **Part 1** | Galaxy (GTN Tutorial) | WGBS pipeline — QC → alignment → methylation extraction → DMR detection |
| **Part 2** | Google Colab (Bio-Learn) | Multi-clock biological age estimation across two independent methylation datasets |

---

## Part 1 — WGBS DNA Methylation Analysis (Galaxy)

### Background

DNA methylation — the addition of a methyl group to cytosine residues, primarily at CpG dinucleotides — is a key epigenetic mark. Because bisulfite treatment converts unmethylated cytosines to uracil (read as thymine), standard alignment tools fail; dedicated bisulfite-aware tools are required. This part follows the [Galaxy Training Network tutorial](https://training.galaxyproject.org/training-material/topics/epigenetics/tutorials/methylation-seq/tutorial.html) based on **Lin et al. (2015)**, a landmark study on hierarchical clustering of breast cancer methylomes.

### Datasets

Samples were drawn from five breast tissue types, representing a spectrum from normal to cancerous:

| Sample ID | Sample Type |
|-----------|-------------|
| NB1, NB2 | Normal breast cells |
| BT089 | Fibroadenoma (benign tumor) |
| BT126, BT198 | Invasive ductal carcinoma |
| MCF7 | Breast adenocarcinoma cell line |

Raw reads: `subset_1.fastq` / `subset_2.fastq` (subsetted for tutorial compute time). Full data available at [Zenodo record 557099](https://zenodo.org/record/557099).

### Workflow Summary

```
Raw FASTQ reads
      │
      ▼
 Quality Control (Falco / FastQC)
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
 Visualization (deepTools: computeMatrix → plotProfile)
      │  ── methylation level around CpG islands / TSS
      ▼
 Differential Methylation (Metilene)
      └── DMR detection between normal (NB1, NB2) vs. cancer (BT198)
```

### Key Results

**Quality Control:** The per-base sequence content showed a markedly reduced cytosine (C) frequency with a corresponding increase in thymine (T). This is the expected signature of bisulfite treatment — unmethylated cytosines are converted to uracil (sequenced as T), confirming successful library preparation.

**Methylation Bias:** The MethylDackel mbias plot revealed a positional bias at the 5′ end of reads on the OT (original top) strand, indicating that the first several bases should be excluded during extraction to avoid inflated methylation estimates.

**Methylation Profiles:** The `plotProfile` output showed a characteristic dip in CpG methylation centered on CpG islands, consistent with the known role of unmethylated CpG islands at active gene promoters. Normal breast samples (NB1, NB2) displayed relatively uniform profiles, while the cancer-derived samples exhibited altered patterns, reflecting epigenetic reprogramming.

**Differential Methylation (Metilene):** Comparison of normal breast vs. BT198 (invasive ductal carcinoma) identified differentially methylated regions (DMRs) concentrated at CpG islands. The Metilene PDF output illustrated the distribution of mean methylation differences and q-values across detected DMRs, highlighting genomic loci with significant epigenetic divergence between normal and tumor tissue.

---

## Part 2 — Biological Clock Benchmarking (Bio-Learn)

### Datasets Used

Two independent whole-blood EPIC array datasets were loaded via the [Bio-Learn](https://bio-learn.github.io/) `DataLibrary`:

#### GSE40279 — Large-Scale Blood Methylation Aging Study
- **Source:** Hannum et al. (2013) — one of the foundational aging methylome datasets
- **Platform:** Illumina 450K array
- **Original size:** 656 samples, 473,034 CpG sites
- **Used:** First 200 samples (memory-constrained subset)
- **Age range:** Adult blood donors, broad age span
- **Relevance:** The canonical reference dataset for training and validating methylation-based age clocks

#### GSE51057 — Blood Methylation Longitudinal / Cross-Sectional Study
- **Source:** Independent whole-blood methylation cohort
- **Platform:** Illumina 450K / EPIC array
- **Original size:** 329 samples, 485,577 CpG sites
- **Used:** First 200 samples
- **Relevance:** Serves as an external validation dataset; different cohort composition from GSE40279 enables cross-dataset benchmarking

### Aging Clocks & Models

Eight epigenetic clocks/models were selected to represent different biological concepts of aging:

| Clock | Type | Predicts | Key Feature |
|-------|------|----------|-------------|
| **Horvath v1** | Epigenetic clock | Chronological age | Pan-tissue; trained on 51 tissue/cell types; 353 CpGs |
| **Horvath v2** (SkinBloodClock) | Epigenetic clock | Chronological age | Optimized for blood and skin; 391 CpGs |
| **Hannum** | Epigenetic clock | Chronological age | Blood-specific; trained on GSE40279; 71 CpGs |
| **PhenoAge** | Biological age clock | Phenotypic/mortality age | Trained on clinical biomarkers + methylation; captures morbidity risk |
| **DunedinPACE** | Pace-of-aging | Rate of aging (0–2 scale) | Longitudinal; measures speed of biological aging, not absolute age |
| **YingCausAge** | Causal aging clock | Causal biological age | Uses causal inference framework to separate aging signal from confounders |
| **Knight** | Gestational age clock | Gestational age | Cord blood-derived; often shows flat/compressed predictions in adults |
| **Bocklandt** | Saliva-based clock | Chronological age | Trained on saliva; limited CpG set; outputs fractional values in adult blood |

> **Note on Knight & Bocklandt:** These clocks were trained on tissue types other than whole blood (cord blood and saliva, respectively). Their predictions in blood datasets are expected to be compressed or offset — this is an informative cross-tissue artifact rather than a failure.

### Analyses Performed

**1. Clock Inference Pipeline**

All 8 clocks were applied to each dataset using `ModelGallery().get(clock).predict(data)`. Results were compiled into a unified DataFrame alongside chronological age (`ChronologicalAge`) for each sample.

**2. Correlation Matrix** (`corr_GSE40279.png`, `corr_GSE51057.png`)

Pearson correlations were computed pairwise across all clock predictions and chronological age. The correlation heatmap (with annotated values) reveals which clocks agree with each other and which capture distinct biological signals.

**3. Age Deviation Heatmap** (`deviation_GSE40279.png`, `deviation_GSE51057.png`)

For each age-predicting clock, the deviation from chronological age (`predicted − chronological`) was computed per sample and visualized as a clustermap. DunedinPACE and Bocklandt were excluded from this plot as they operate on non-age scales. The heatmap reveals systematic over- or under-aging patterns and clusters samples by deviation profile.

**4. Predicted vs. Chronological Age Scatter** (`scatter_GSE40279.png`, `scatter_GSE51057.png`)

Scatter plots of predicted biological age against chronological age were generated for each clock, with a diagonal reference line (perfect accuracy). Pearson r and p-values are annotated per subplot, enabling direct assessment of clock accuracy and systematic bias across both datasets.

### Visualizations

| File | Description |
|------|-------------|
| `corr_GSE40279.png` | Clock-to-clock + chronological age correlation matrix — Dataset 1 |
| `corr_GSE51057.png` | Clock-to-clock + chronological age correlation matrix — Dataset 2 |
| `deviation_GSE40279.png` | Per-sample age deviation heatmap — Dataset 1 |
| `deviation_GSE51057.png` | Per-sample age deviation heatmap — Dataset 2 |
| `scatter_GSE40279.png` | Predicted vs. chronological age scatter (all clocks) — Dataset 1 |
| `scatter_GSE51057.png` | Predicted vs. chronological age scatter (all clocks) — Dataset 2 |

---

## Repository Structure

```
.
├── README.md
├── biological_clocks.ipynb        # Google Colab notebook — Part 2 analysis
│
├── outputs/                       # Generated figures
│   ├── corr_GSE40279.png
│   ├── corr_GSE51057.png
│   ├── deviation_GSE40279.png
│   ├── deviation_GSE51057.png
│   ├── scatter_GSE40279.png
│   └── scatter_GSE51057.png
│
└── galaxy_tutorial/               # Part 1 reference materials
    └── methylation-seq-tutorial.html
```

---

## Dependencies & Setup

### Part 1 (Galaxy)

All tools run on the [Galaxy public server](https://usegalaxy.org/). No local installation required.

| Tool | Version |
|------|---------|
| Falco | 1.2.4+galaxy0 |
| bwameth | 0.2.7+galaxy0 |
| MethylDackel | 0.5.2+galaxy0 |
| Wig/BedGraph-to-bigWig | — |
| computeMatrix (deepTools) | 3.5.4+galaxy0 |
| plotProfile (deepTools) | 3.5.4+galaxy0 |
| Metilene | 0.2.6.1 |

### Part 2 (Google Colab / Python)

```bash
pip install biolearn matplotlib seaborn pandas numpy scipy
```

| Package | Version used |
|---------|-------------|
| biolearn | 0.9.1 |
| pandas | 2.2.2 |
| numpy | 2.0.2 |
| matplotlib | 3.10.0 |
| seaborn | 0.13.2 |
| scipy | 1.16.3 |
| scikit-learn | 1.6.1 |

**Run the notebook:**

Open `biological_clocks.ipynb` in [Google Colab](https://colab.research.google.com/), run all cells top to bottom. Dataset download is handled automatically by Bio-Learn's `DataLibrary`. Estimated runtime: ~10–20 minutes depending on download speed.

---

## References

1. **Lin et al. (2015).** Hierarchical Clustering of Breast Cancer Methylomes Revealed Differentially Methylated and Expressed Breast Cancer Genes. *PLOS ONE* 10(1): e0118453. https://doi.org/10.1371/journal.pone.0118453

2. **Hannum et al. (2013).** Genome-wide Methylation Profiles Reveal Quantitative Views of Human Aging Rates. *Molecular Cell* 49(2): 359–367.

3. **Horvath S. (2013).** DNA methylation age of human tissues and cell types. *Genome Biology* 14: R115.

4. **Levine et al. (2018).** An epigenetic biomarker of aging for lifespan and healthspan. *Aging* 10(4): 573–591. [PhenoAge]

5. **Belsky et al. (2022).** DunedinPACE, a DNA methylation biomarker of the pace of aging. *eLife* 11: e73420.

6. **Ying et al. (2024).** Causally-inferred aging clocks. Bio-Learn model gallery. [YingCausAge]

7. **Knight et al. (2016).** An epigenetic clock for gestational age at birth. *Genome Biology* 17: 206.

8. **Bocklandt et al. (2011).** Epigenetic predictor of age. *PLOS ONE* 6(6): e14821.

9. **Wolff J., Ryan D., Moosmann V. (2017).** DNA Methylation data analysis. *Galaxy Training Materials*. https://training.galaxyproject.org/training-material/topics/epigenetics/tutorials/methylation-seq/tutorial.html

10. **Hiltemann S., Rasche H. et al. (2023).** Galaxy Training: A Powerful Framework for Teaching! *PLOS Computational Biology*. https://doi.org/10.1371/journal.pcbi.1010752
