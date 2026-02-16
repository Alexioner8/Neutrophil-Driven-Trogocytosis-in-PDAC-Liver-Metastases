# Neutrophil-Driven Trogocytosis in PDAC Liver Metastases

Single-cell RNA-seq analysis of neutrophil functional states and intercellular interactions  
in pancreatic ductal adenocarcinoma (PDAC) liver metastases.

This repository contains scripts used for quality control, integration, clustering,
differential expression, pathway analysis, and visualization of single-cell RNA-seq data.

---

## Repository Structure

scripts/    - Analysis scripts  
QC/         - Quality control outputs  
de/         - Differential expression results  
results/    - Generated figures and tables (ignored by git)

---

## Main Script

scripts/process_Alejandro_final.R

---

## Computational Environment

- R version: 4.4.3  
- OS: Linux (SUSE Linux Enterprise Server 15 SP7)  
- BLAS/LAPACK: OpenBLAS  

---

## Core R Packages

Single-cell analysis:
- Seurat (v5.3.0)  
- SeuratObject  

Visualization:
- ggplot2  
- cowplot  
- patchwork  
- ComplexHeatmap  
- circlize  
- EnhancedVolcano  
- ggrepel  

Data manipulation:
- tidyverse  
- data.table  

Functional enrichment:
- clusterProfiler  
- enrichplot  

Parallelization:
- future  
- future.apply  

Export utilities:
- writexl  
- svglite  

---

