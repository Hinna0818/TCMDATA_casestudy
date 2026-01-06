# TCMDATA Case Study: Scorpion for Diabetic Nephropathy

This repository demonstrates the application of the **[TCMDATA](https://github.com/Hinna0818/TCMDATA)** R package to explore the molecular mechanism of **Scorpion (全蝎)** in treating **Diabetic Nephropathy (DN)**.

## Project Overview

We integrated **Network Pharmacology** with **Single-Cell RNA sequencing (scRNA-seq)** analysis to validate TCM targets in the cellular microenvironment.

**Analysis Workflow:**
1.  **Data Retrieval**: Used `TCMDATA` to fetch Scorpion ingredients and targets.
2.  **PPI Network**: Identified Hub Genes using topological metrics and clustering methods.
3.  **scRNA-seq Mapping**: Validated that core targets are predominantly expressed in **Macrophages, Neutrophils, and NK/T Cells**.
4.  **Cell-Cell Interaction (CCI)**: Visualized communication alterations using `CellChat` R package to support the findings in bulk RNA-seq.

## Key Biological Findings

Our analysis proposes that Scorpion alleviates DN by blocking **Immune-Stromal Crosstalk** through three key axes:

* **Anti-Adhesion (VCAM1 Axis)**: Targeting **VCAM1** to disrupt the interaction between immune cells and endothelial cells, reducing infiltration.
* **Anti-Inflammation (CXCL/TNF Axis)**: Inhibiting **CXCL1/8** and **TNF** to break the inflammatory recruitment loop driven by neutrophils and macrophages.
* **Anti-Fibrosis (Macrophage-Fibroblast Axis)**: Blocking pro-fibrotic signaling (e.g., **SPP1**) from macrophages to fibroblasts.
