# Transcriptomic Profiling of Vitiligo Skin (GSE311770)

## 🧬 Project Overview
This repository contains an end-to-end RNA-Seq differential expression analysis pipeline built in R. The objective is to uncover molecular mechanisms and biomarker signatures associated with vitiligo by comparing lesional skin biopsies against healthy controls.

- **Dataset:** [GSE311770](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE311770) (NCBI GEO)
- **Biological Insight:** Validated the hallmark autoimmune destruction of melanocytes via severe downregulation of pigmentary markers (`TYR`, `TRPM1`, `GPR143`).

---

## 🛠️ Pipeline Architecture
1. **Data Ingestion & QC:** Loaded raw RSEM count matrices and filtered low-count genes (`rowSums >= 10`).
2. **Metadata Engineering:** Programmatically parsed sample identifiers using regular expressions to assign clinical experimental groups (Lesional vs. Control).
3. **Differential Expression Analysis:** Implemented a Negative Binomial generalized linear model using `DESeq2`.
4. **Gene Annotation:** Mapped Ensembl Gene IDs to human Gene Symbols via `org.Hs.eg.db`.
5. **Visualization:** Generated publication-grade Volcano plots utilizing `EnhancedVolcano`.

---

## 📁 Repository Structure
```text
├── scripts/
    └── vitiligo_rnaseq_analysis.R    # Master R analysis script
├── results/
    ├── Vitiligo_Lesional_vs_Control_SigDEGs.csv  # Filtered significant DEGs
    └── Vitiligo_Volcano_Plot.pdf                 # Publication-quality volcano plot
