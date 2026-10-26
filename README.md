# Targeting High-Risk Macrophage Phenotypes in the Breast Cancer Microenvironment: A Single-Cell Transcriptomic and Drug Repurposing Analysis
---
**Author(s):** Md Abu Sadat Soash, Raphael Kabir Niloy
*Department of Genetic Engineering and Biotechnology*<br>
*Shahjalal University of Science and Technology, Sylhet, Bangladesh*

---

## 📄 Abstract

The tumor microenvironment (TME) is a complex ecosystem where distinct cellular states dictate disease progression and therapeutic resistance. In this study, we performed a comprehensive single-cell transcriptomic dissection of the breast cancer TME using publicly available scRNA-seq data (12 paired samples across 6 patients; 98,410 high-quality cells). Unsupervised clustering identified 39 distinct cellular subpopulations. By applying the Single-Cell Phenotype-Associated Signatures (**scPAS**) risk-scoring algorithm mapped against TCGA-BRCA survival metrics, we discovered a highly immunosuppressive macrophage subpopulation that serves as a primary driver of high-risk clinical phenotypes. 

Pathway enrichment analysis revealed that these `scPAS+` macrophages are hyper-metabolic, relying heavily on **Oxidative Phosphorylation (OXPHOS)**, and are actively engaged in antigen presentation, phagocytosis, and chemotaxis to orchestrate an immunosuppressive ecosystem. To identify therapeutic vulnerabilities, we deployed dual *in silico* drug repurposing pipelines (**DrugReSC** and **ASGARD**). These analyses consistently identified targeted kinase inhibitors (e.g., Dabrafenib, Vemurafenib, Imatinib, Pazopanib) and proteasome inhibitors (e.g., Ixazomib) as top candidates capable of reversing the high-risk transcriptomic signature. Our findings provide a robust computational framework for identifying and targeting prognostic macrophage states to improve breast cancer immunotherapy outcomes.

---

## 📌 Key Findings & Study Highlights

* **Primary Prognostic Culprit:** Tumor-associated macrophages (TAMs) are the dominant drivers of poor clinical outcomes, with **52.38%** of annotated macrophages exhibiting the high-risk (`scPAS+`) survival phenotype, followed by unspecified myeloid cells (36.31%), IgG plasma cells (35.15%), inflammatory macrophages (28.74%), and pericytes (27.48%).
* **Scale of Analysis:** 124,373 raw single cells filtered down to **98,410 high-fidelity singlets** across 12 paired samples (6 primary breast tumor, 6 adjacent normal).
* **Survival Modeling (scPAS):** Model trained on TCGA-BRCA bulk RNA-seq and survival data identified **686 survival-associated features** (310 risk-associated, 376 protective). Overall, **8,151 cells (8.28%)** exhibited the high-risk phenotype.
* **Molecular Hallmarks:** Significant upregulation of **Oxidative Phosphorylation (OXPHOS)**, **Proteasome**, **Phagosome**, and **Lysosome biogenesis** in KEGG pathways, coupled with **antigen processing and presentation**, **phagocytosis**, and **chemotaxis** in Gene Ontology (GO) biological processes.
* **Dual Repurposing Concordance:** 
  - **DrugReSC** (AUCell Random Forest classification against LINCS L1000) identified BRAF/MEK inhibitors (**Dabrafenib**, **Vemurafenib**) and tyrosine kinase inhibitors (**Imatinib**).
  - **ASGARD** (cluster-specific Tumor vs. Normal perturbation matching) independently identified **Ixazomib** (Proteasome inhibitor, directly targeting the identified proteasome hyperactivation pathway) and **Pazopanib** (multi-kinase inhibitor).

---

## 📊 Visual Highlights

<p align="center">
  <table>
    <tr>
      <td align="center"><b>Phase 2: TME Integration & Clustering (Harmony)</b></td>
      <td align="center"><b>Phase 4: scPAS Risk Phenotype Distribution (UMAP)</b></td>
    </tr>
    <tr>
      <td><img src="result/2_integration/plots/Phase2_Integration_UMAP.png" alt="Integrated UMAP" width="450"/></td>
      <td><img src="result/4_scPAS/plots/Figure2A_Risk_UMAP.png" alt="Risk UMAP" width="450"/></td>
    </tr>
    <tr>
      <td align="center"><b>Phase 4: Risk Proportion by Cell Type</b></td>
      <td align="center"><b>Phase 7: KEGG Pathway Activation (scPAS+)</b></td>
    </tr>
    <tr>
      <td><img src="result/4_scPAS/plots/Figure2B_Risk_BarChart.png" alt="Risk Bar Chart" width="450"/></td>
      <td><img src="result/7_pathway/plots/Figure3A_KEGG_Dotplot.png" alt="KEGG Dotplot" width="450"/></td>
    </tr>
    <tr>
      <td colspan="2" align="center"><b>Phase 7: Gene Ontology (GO) Biological Processes (scPAS+)</b></td>
    </tr>
    <tr>
      <td colspan="2" align="center"><img src="result/7_pathway/plots/Figure3B_GO_Dotplot.png" alt="GO Dotplot" width="600"/></td>
    </tr>
  </table>
</p>

---

## 💊 Top In Silico Drug Repurposing Candidates

As detailed in **Table 1** of the manuscript, dual *in silico* profiling against the Broad Institute LINCS L1000 perturbational database yielded highly concordant therapeutic candidates:

| Pipeline | Rank | Drug Candidate | Primary Drug Class / Target | Predictive Score / Significance |
| :--- | :---: | :--- | :--- | :--- |
| **DrugReSC** | 1 | **Dabrafenib** | BRAF Inhibitor | Score: 1.000 |
| **DrugReSC** | 2 | **Vemurafenib** | BRAF Inhibitor | Score: 0.988 *(6 of top 10 signatures)* |
| **DrugReSC** | 5, 7 | **Imatinib** | Tyrosine Kinase Inhibitor | Top 10 Hit |
| **ASGARD** | 1 | **Mitoxantrone** | Topoisomerase II Inhibitor | Score: Inf ($FDR < 10^{-29}$) |
| **ASGARD** | 2 | **Vorinostat** | HDAC Inhibitor | Score: Inf ($FDR < 10^{-80}$) |
| **ASGARD** | 3 | **Mebendazole** | Microtubule Inhibitor | Score: 0.562 |
| **ASGARD** | 4 | **Ixazomib** | Proteasome Inhibitor | Score: 0.367 *(Validates Proteasome KEGG hit)* |
| **ASGARD** | 5 | **Pazopanib** | Tyrosine Kinase Inhibitor (VEGFR/PDGFR/KIT) | Score: 0.303 |

---

## 📋 Supplementary Quality Control Metrics

Summary of the dynamic 3-MAD filtering and `scDblFinder` doublet detection across all 12 paired samples (**Supplementary Table 1**):

| Sample ID | Tissue Status | Raw Cell Count | High-Quality Singlets | Doublets Removed | Max MT% Applied |
| :--- | :---: | :---: | :---: | :---: | :---: |
| **GSM8061824** | Tumor | 13,651 | 10,218 | 1,554 | 18.26% |
| **GSM8061825** | Tumor | 13,240 | 10,621 | 1,464 | 13.16% |
| **GSM8061826** | Tumor | 17,328 | 13,991 | 1,934 | 14.32% |
| **GSM8061827** | Tumor | 10,064 | 8,353 | 911 | 20.00% |
| **GSM8061828** | Tumor | 14,807 | 11,878 | 1,592 | 16.29% |
| **GSM8061829** | Tumor | 14,450 | 8,903 | 1,107 | 20.00% |
| **GSM8061830** | Normal | 5,706 | 4,846 | 449 | 8.24% |
| **GSM8061831** | Normal | 2,751 | 2,487 | 145 | 20.00% |
| **GSM8061832** | Normal | 4,917 | 4,360 | 325 | 16.10% |
| **GSM8061833** | Normal | 12,421 | 9,699 | 1,349 | 20.00% |
| **GSM8061834** | Normal | 9,512 | 8,171 | 810 | 19.41% |
| **GSM8061835** | Normal | 5,526 | 4,883 | 438 | 20.00% |
| **Total / Overall** | — | **124,373** | **98,410** | **12,078** | — |

---

## 📂 Repository Structure

```text
BRCA-scRNA-Drug-Repurposing/
├── manuscript.docx       # Full manuscript draft with figures and tables
├── data/                 # Raw 10X Genomics count matrices
│   ├── cancer/           # 6 Primary breast cancer samples (GSM8061824–GSM8061829)
│   └── normal/           # 6 Matched adjacent normal samples (GSM8061830–GSM8061835)
├── database/             # Reference databases and generation scripts
│   ├── DrugReSC/         # DrugReSC perturbagen signature references
│   ├── L1000/            # Broad Institute LINCS Level 5 databases (.gctx)
│   ├── TCGA_BRCA/        # TCGA Bulk RNA & clinical survival reference data
│   └── THBCA/            # The Human Breast Cancer Atlas (THBCA) reference
├── R/                    # Core Analytical Pipeline Scripts
│   ├── 1_qc.R            # Dynamic 3-MAD filtering & scDblFinder doublet removal
│   ├── 2_integration.R   # SCTransform v2 & Harmony integration across 12 samples
│   ├── 3_annotation.R    # SingleR cell-type annotation against THBCA
│   ├── 4_scPAS.R         # Survival risk scoring & prognostic cell stratification
│   ├── 4_graph.R         # Risk distribution plotting scripts
│   ├── 5_DrugResc.R      # AUCell Random Forest drug repurposing
│   ├── 6_Asgard.R        # Orthogonal cluster-specific mono-drug validation
│   ├── 7_pathway_enrichment.R # clusterProfiler GO/KEGG functional enrichment
│   └── 7_graph.R         # Functional pathway dotplot scripts
└── result/               # Checkpoints, tables, and publication-ready figures
    ├── 1_qc/             # Post-QC cell counts, doublet stats, metadata
    ├── 2_integration/    # Integrated Seurat object (Harmony) & UMAPs
    ├── 3_annotation/     # Annotated SingleR assignments & THBCA labels
    ├── 4_scPAS/          # scPAS risk scores, high-risk cell tables & plots
    ├── 5_drug_discovery/ # DrugReSC top candidate perturbagens
    ├── 6_Asgard/         # ASGARD drug scores and FDR tables
    └── 7_pathway/        # Enriched KEGG & GO pathway tables and plots
```

---

## 💻 Computational Environment & Requirements

The pipeline processes high-dimensional matrices (98,410 cells $\times$ 39 clusters) against high-throughput perturbation libraries (LINCS L1000). To avoid memory bottlenecks, execution was conducted on **Google Cloud Platform (GCP)**:

* **Machine Type:** `n1-standard-16`
* **vCPUs:** 16
* **Memory (RAM):** 60 GB
* **GPU Acceleration:** 1 $\times$ NVIDIA Tesla T4
* **R Version:** $\ge$ 4.3.0
* **Key Packages:** `Seurat v5`, `Harmony`, `scDblFinder`, `SingleR`, `scPAS`, `clusterProfiler`, `qs2`

---

## 🚀 Reproduction Guide

### 1. Build Reference Databases
Generate the required `.rds` and `.h5ad` reference matrices:
```bash
Rscript database/TCGA_BRCA/main.R
Rscript database/THBCA/main.R
Rscript database/DrugReSC/main.R
```

### 2. Execute Primary Pipeline Sequentially
```bash
# Step 1: Dynamic 3-MAD QC and doublet filtering
Rscript R/1_qc.R

# Step 2: SCTransform v2 normalization and Harmony batch integration
Rscript R/2_integration.R

# Step 3: SingleR automated cell-type annotation against THBCA
Rscript R/3_annotation.R

# Step 4: scPAS prognostic risk scoring & feature selection
Rscript R/4_scPAS.R
Rscript R/4_graph.R

# Step 5: DrugReSC perturbational drug repurposing
Rscript R/5_DrugResc.R

# Step 6: Orthogonal single-cell mono-drug repurposing with ASGARD
Rscript R/6_Asgard.R

# Step 7: Pathway and functional enrichment analysis (GO / KEGG)
Rscript R/7_pathway_enrichment.R
Rscript R/7_graph.R
```

---

## 📚 References

1. **Dataset Source:** Wen, Z., Liu, T., Xu, X., Acharya, N. et al. (2025). *Interleukin-16 enhances anti-tumor immune responses by establishing a Th1 cell-macrophage crosstalk through reprogramming glutamine metabolism in mice.* **Nature Communications**, 16(1), 2362. [GEO Accession: GSE254991](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE254991).
2. **scPAS Algorithm:** Zhao, X., et al. (2024). *scPAS: single-cell phenotype-associated signatures for identifying vulnerable subpopulations and drug repurposing.* **Briefings in Bioinformatics**, 25(2), bbae073.
3. **ASGARD Pipeline:** Zhang, X., et al. (2023). *ASGARD: A single-cell guided pipeline to aid repurposing of drugs.* **Nature Communications**, 14, 5327.
4. **Integration & Single-Cell Framework:** Hao, Y., et al. (2024). *Dictionary learning for integrative, multimodal and scalable single-cell analysis (Seurat v5).* **Nature Biotechnology**, 42, 293–304.
5. **Harmony Batch Correction:** Korsunsky, I., et al. (2019). *Fast, sensitive and accurate integration of single-cell data with Harmony.* **Nature Methods**, 16(12), 1289–1296.
6. **Enrichment Analysis:** Wu, T., et al. (2021). *clusterProfiler 4.0: A universal enrichment tool for interpreting omics data.* **The Innovation**, 2(3), 100141.
