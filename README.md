<div align="center">

# Integrated Bulk RNA-seq Meta-Analysis: Host Transcriptomic Response to ZIKV Infection
**A comprehensive, publication-quality reproducible workflow for elucidating Zika Virus (ZIKV) pathogenesis in A549 epithelial cells.**

[![Reproducible](https://img.shields.io/badge/Reproducible-Yes-success.svg)]()
[![R](https://img.shields.io/badge/Language-R_4.6-blue.svg)]()
[![DESeq2](https://img.shields.io/badge/Analysis-DESeq2-orange.svg)]()
[![License](https://img.shields.io/badge/License-MIT-green.svg)]()

</div>

---

## 📌 Visual Abstract

```mermaid
graph TD
    A1[GSE146423] --> B[Data Acquisition & QC]
    A2[GSE233049] --> B
    A3[GSE265922] --> B
    B --> C[Alignment & Gene Mapping]
    C --> D[Normalization & Filtering]
    D --> E[Differential Expression<br/>DESeq2 + apeglm]
    E --> F[Functional Enrichment<br/>GO / KEGG / GSEA]
    F --> G[Dataset Integration<br/>Overlapping Meta-Signature]
    G --> H[Biological Discovery<br/>ZIKV Antiviral Response]
```

---

## 📑 Table of Contents
1. [Project Highlights](#-project-highlights)
2. [Dataset Cards](#-dataset-cards)
3. [Scientific Results & Key Genes](#-scientific-results--key-genes)
4. [Pipeline Statistics Dashboard](#-pipeline-statistics-dashboard)
5. [Directory Structure](#-directory-structure)
6. [Computational Environment](#-computational-environment)
7. [Reproducibility Checklist](#-reproducibility-checklist)

---

## ✨ Project Highlights

- **Three independent transcriptomic datasets** integrated into a single pipeline.
- **Unified statistical preprocessing**, ensuring true cross-study comparability.
- **Advanced LFC Shrinkage (apeglm)** to aggressively reduce noise from low-count genes.
- **Comprehensive Functional Enrichment** (Gene Ontology, KEGG, and GSEA).
- **Publication-quality 600 DPI visualizations** generated uniformly across all data.
- **Fully reproducible R-scripted architecture**.

---

## 🗂 Dataset Cards

<details open>
<summary><b>Dataset 1: GSE146423</b></summary>
<blockquote>
<b>Model:</b> Human A549 Cells (ZIKV vs Mock) <br>
<b>Platform:</b> Illumina HiSeq 4000 <br>
<b>Format:</b> Pre-processed Entrez Count Matrix <br>
<b>Top 15 Upregulated Genes:</b> <i>PARP9, PARP14, DDX60, RIGI, IFIT2, IFIT3, OAS3, HELZ2, DDX60L, DTX3L, MX1, ISG15, TRANK1, IFIT1, HERC6</i>
</blockquote>
</details>

<details open>
<summary><b>Dataset 2: GSE233049</b></summary>
<blockquote>
<b>Model:</b> Human A549 Cells (ZIKV vs Mock) <br>
<b>Platform:</b> Illumina NovaSeq 6000 <br>
<b>Format:</b> Pre-processed Entrez Count Matrix <br>
<b>Top 15 Upregulated Genes:</b> <i>ISG15, IFI6, JUN, GBP1, IFI16, ATF3, IFIH1, SP110, TRANK1, NFKBIZ, PARP9, DTX3L, PARP14, PLSCR1, HERC5</i>
</blockquote>
</details>

<details open>
<summary><b>Dataset 3: GSE265922</b></summary>
<blockquote>
<b>Model:</b> Human A549 Cells (ZIKV vs Mock) <br>
<b>Platform:</b> Illumina NovaSeq 6000 <br>
<b>Format:</b> Raw STAR <code>ReadsPerGene.out.tab.gz</code> mapped to ENSEMBL <br>
<b>Top 15 Upregulated Genes:</b> <i>IFIT3, DHX58, PTGER4, PMAIP1, IFIT1, KLF4, KDM7A-DT, ATP4A, ACHE, TNFAIP3, IFIH1, IFIT2, TAPBPL, ATF3, PLEKHA4</i>
</blockquote>
</details>

---

## 🧬 Scientific Results & Key Genes

By standardizing the DESeq2 pipeline, we stripped away study-specific technical noise to identify a highly robust, mathematically conserved ZIKV host-response signature.

### 🥇 The "Universal 3/3" Conserved Core
Because individual transcriptomic studies suffer from distinct batch effects and varying viral MOIs, requiring a gene to be significantly perturbed in 3 out of 3 studies is an extremely restrictive filter. However, our rigorously corrected pipeline successfully resolved this noise to identify **66 universally upregulated genes** spanning the core antiviral architecture.

### 🥈 The "Robust 2/3" Meta-Signature (1,118 Genes)
Broadening the signature to genes perturbed in $\geq 2$ of the datasets completely reconstructed the **Type I Interferon Antiviral Response**, yielding **1,118 significantly upregulated** genes and **16 downregulated** genes, providing massive *in silico* biological validation of our pipeline's accuracy. 

Key gene families dominating this meta-signature:
- **Cytosolic RNA Sensors:** `IFIH1` (MDA5), `RIGI` (DDX58)
- **Viral Translation Inhibitors:** `IFIT1`, `IFIT2`, `IFIT3`, `IFIT5`, `IFITM1`
- **RNA Degradation Effectors:** `OAS2`, `OAS3`, `OASL`
- **Capsid Assembly Blockers:** `MX1`, `MX2`
- **Master Transcription Factors:** `IRF7`, `IRF9`

> [!IMPORTANT]
> The appearance of `IFIH1` (MDA5), `OAS3`, and `IFIT1` natively validates the analysis, as these are the primary cytosolic sensors and effectors for Flaviruses like Zika.

---

## 📊 Pipeline Statistics Dashboard

| Metric | Value |
| :--- | :--- |
| **Independent Datasets** | 3 |
| **Total Samples Analyzed** | 18 |
| **Genes Tested per Dataset** | ~16,700 - 20,600 |
| **Conserved Upregulated Genes ($\geq 2$ datasets)** | 1,118 |
| **Conserved Downregulated Genes ($\geq 2$ datasets)** | 16 |
| **Output Publication Figures** | 66 (22 per dataset) |

---

## 📁 Directory Structure

```text
Zika_wetlab/
├── meta-analysis/         # Integrated 2/3 overlap analysis & final gene lists
│   ├── results/
│   │   ├── Common_Upregulated_3of3.csv
│   │   ├── Common_Downregulated_3of3.csv
│   │   ├── Common_Upregulated_2of3.csv
│   │   └── Common_Downregulated_2of3.csv
│   └── meta_analysis.R    # Intersection aggregation script
│
├── GSE146423/             # Complete isolated pipeline for Dataset 1
│   ├── plots/             # 600 DPI PNGs (QC, Volcano, Heatmap, GSEA, Networks)
│   ├── results/           # Raw DESeq2 & Enrichment CSVs
│   └── GSE146423.R        # Execution script
│
├── GSE233049/             # Complete isolated pipeline for Dataset 2
│   ├── plots/
│   ├── results/
│   └── GSE233049.R
│
└── GSE265922/             # Complete isolated pipeline for Dataset 3
    ├── plots/
    ├── results/
    └── GSE265922.R
```

---

## 💻 Computational Environment

- **OS:** Windows 
- **R Version:** 4.6.0
- **Key Packages:** `DESeq2`, `apeglm`, `clusterProfiler`, `enrichplot`, `tidyverse`, `EnhancedVolcano`, `pheatmap`

---

## ✅ Reproducibility Checklist

- [x] Raw data matrix processing included
- [x] Strict statistical significance cutoffs applied uniformly
- [x] High-resolution 600 DPI plotting architecture standardized
- [x] Final pipeline executed successfully without manual intervention
- [x] Meta-analysis intersects thoroughly documented


## 🖼️ Figure Gallery

*A selection of the publication-quality visualizations automatically generated by this workflow.*

<details open>
<summary><b>Meta-Analysis Integrations</b></summary>
<br>
<img src="meta-analysis/plots/01_Venn_Up.png" width="45%" />
<img src="meta-analysis/plots/03_Common_Signature_Heatmap.png" width="45%" />
</details>

<details open>
<summary><b>Dataset Quality Control & Differential Expression (Example: GSE146423)</b></summary>
<br>
<img src="GSE146423/plots/QC/01_PCA.png" width="45%" />
<img src="GSE146423/plots/DEG/07_EnhancedVolcano.png" width="45%" />
</details>

<details open>
<summary><b>Functional Enrichment & Pathway Networks (Example: GSE146423)</b></summary>
<br>
<img src="GSE146423/plots/Enrichment/11_GO_BP_Up.png" width="45%" />
<img src="GSE146423/plots/GSEA/21_GSEA_CnetPlot.png" width="45%" />
</details>

