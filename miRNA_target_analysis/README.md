# miRNA Target Analysis — linking the 66-gene ZIKV core to miR-134 / 146a / 146b

This folder connects the bulk-RNA-seq **66-gene conserved ZIKV core** (from the
meta-analysis) to the three miRNAs tested in the wet-lab study
(`../miRNA_Antiviral_Main_Findings.md`): **miR-134-5p, miR-146a-5p, miR-146b-5p**.

## How it was made
Validated miRNA→target interactions were taken from **miRTarBase**
(`../00_Raw_miRNA_Data/hsa_MTI.csv`, human `hsa_MTI`) and intersected with the
66 conserved-up core genes (`../meta-analysis/results/Common_Upregulated_3of3.csv`).
Only human interactions for the exact mature miRNAs were used.

## Files
| File | Contents |
| :--- | :--- |
| **`miRNA_target_analysis.py`** | the analysis code (Python) — reads miRTarBase + the 66-core and writes all tables below. Run: `python miRNA_target_analysis.py` |
| **`miRNA_target_analysis.R`** | same analysis in R (uses `data.table::fread`). Run: `Rscript miRNA_target_analysis.R` |
| `66_core_genes_miRNA_annotated.csv` | all 66 core genes, per-dataset log2FC, and yes/no flags for miR-146a / 146b / 134 targeting + panel membership |
| `miRNA_target_interactions.csv` | the 21 validated miRNA→core-gene interactions, with evidence type (Functional / Weak) and PMIDs |
| `5_gene_validation_panel.csv` | the 5 chosen panel genes with log2FC in each dataset |
| `07_miRNA_Panel_5genes.png` | the panel figure (regenerate in R via `../meta-analysis/mirna_panel_figure.R`) |

**Raw input:** `../00_Raw_miRNA_Data/hsa_MTI.csv` — the miRTarBase human MTI master
file (miRNA→target validated interactions), used by both scripts above.

## Key numbers
- **miR-146a-5p** targets **15 / 66** core genes (incl. the master TF **IRF7**)
- **miR-146b-5p** targets **6 / 66**
- **miR-134-5p** targets **0 / 66** (it acts elsewhere — NOLC1, indirect)

## The 5-gene validation panel
**IRF7, IFIT3, RSAD2, CCL5, SAMD9L** — core genes that are validated miR-146a/146b
targets (IFIT3 & SAMD9L are dual 146a+146b targets). Present these as
**"genes of interest / miR-146-regulated core genes"**, NOT "top 5 by significance"
(all are strongly significant, but not the 5 lowest p-values).

## ⚠️ Interpretation caution (important)
miR-146a/146b are **antiviral** in the wet-lab study, yet they **repress** these
antiviral ISGs. So their antiviral mechanism is **ISG-independent** — do NOT write
"146a is antiviral because it boosts ISGs" (the direction is the opposite). In
natural infection the reciprocal pattern (miRNA down → ISG targets up = de-repression)
is what supports genuine regulation.

## Evidence caveat
Most miR-146a→ISG interactions trace to one high-throughput study (PMID 18057241);
CCL5 has independent evidence (PMID 24996260); the miR-146b entries have no PMID in
miRTarBase (weakest). The wet-lab qPCR/luciferase is what upgrades these to confirmed.
