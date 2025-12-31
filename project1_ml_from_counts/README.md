
## How to Run This Project

### 1. Environment setup
```bash
conda create -n cancer-rnaseq python=3.11 -y
conda activate cancer-rnaseq
pip install pandas numpy scikit-learn matplotlib


2. Differential expression & pathway analysis (R)

Rscript R/02_deseq2_brca.R
Rscript R/03_pathway_enrichment.R
Rscript R/04_marker_boxplots.R

3. Tumor vs normal ML classification

python python/01_build_ml_matrix.py
python python/02_train_logreg_auc.py
python python/03_export_logreg_coeffs.py

4. ER+ vs ER− classification

python python/04_define_er_status.py
python python/05_build_er_ml_matrix.py
python python/06_train_er_logreg.py
python python/07_export_er_logreg_coeffs.py

Outputs

All results (tables and figures) are written to the results/ directory.



## Key differentially expressed genes

DESeq2 identified 5,795 significant genes (padj < 0.05, |log2FC| > 1). Among the most strongly tumor-upregulated genes were:

CSAG1 (log2FC = 8.57), MAGEA6 (7.54), MAGEA3 (7.16) — Cancer-Testis Antigens (CTAs) that are frequently activated in tumors via epigenetic dysregulation and are well-known tumor-associated expression programs.

COL10A1 (7.10) — an extracellular matrix / stromal-associated gene often elevated in aggressive tumor microenvironments and remodeling contexts.

DSCAM-AS1 (6.89) — a long non-coding RNA reported as highly expressed in subsets of breast cancer and associated with hormone receptor biology.

In contrast, the most strongly normal-enriched genes included:

ADIPOQ (−5.32), LEP (−6.19), PLIN1 (−5.32), CIDEC (−5.31), GPD1 (−5.31) — a coherent adipogenesis/lipid storage program consistent with normal breast tissue’s adipose/stromal composition.

APOB (−5.42) — lipid transport metabolism, also consistent with a stronger metabolic program in non-tumor tissue.

Overall, downregulated genes form a coherent adipose/lipid metabolic signature, whereas tumor-upregulated genes reflect tumor-associated antigen expression and microenvironment remodeling.

## Pathway-level interpretation (Hallmark GSEA)

Hallmark GSEA on the ranked gene list showed a dominant tumor enrichment for proliferation and cell-cycle control:

HALLMARK_E2F_TARGETS (NES = 2.41, padj = 0.00147)

HALLMARK_G2M_CHECKPOINT (NES = 2.33, padj = 0.00147)

HALLMARK_MYC_TARGETS_V1/V2 (NES = 1.77 / 1.69, padj ≤ 0.00388)

HALLMARK_MITOTIC_SPINDLE (NES = 1.56, padj = 0.00261)

The leading-edge genes driving these enrichments include classic proliferation regulators such as MKI67, TOP2A, CDC20, CDK1, CCNB2, AURKA, AURKB, PLK1, confirming a strong tumor proliferative signature.

Breast-cancer–relevant hormone signaling also appeared strongly tumor-enriched:

HALLMARK_ESTROGEN_RESPONSE_LATE (NES = 1.65, padj = 0.00147)

HALLMARK_ESTROGEN_RESPONSE_EARLY (NES = 1.43, padj = 0.0146)
with leading-edge genes including TFF1, TFF3, CCND1, XBP1, consistent with estrogen-responsive transcriptional programs commonly observed in ER+ breast cancers.

Normal tissue was enriched for metabolic and tissue-identity programs:

HALLMARK_ADIPOGENESIS (NES = −2.48, padj = 0.00261)

HALLMARK_FATTY_ACID_METABOLISM (NES = −1.89, padj = 0.00261)

HALLMARK_XENOBIOTIC_METABOLISM (NES = −1.66, padj = 0.00261)


Together, these results indicate that BRCA tumor samples are characterized by strong proliferation/cell-cycle activity and estrogen signaling, while normal samples retain adipose/lipid metabolic features typical of breast tissue composition.

Some top-ranked genes (e.g., CGA/MUC2) likely reflect subtype heterogeneity and cellular composition effects within a large, multi-cohort TCGA-derived dataset; therefore, pathway-level enrichments and consistent program-level signatures (E2F/G2M/MYC and estrogen response) were used as the main biological anchors.

Classic proliferation markers (MKI67/TOP2A/CDC20) are elevated in tumors, while hormone receptor markers (ESR1/PGR) show subtype-dependent tumor expression, consistent with ER+ vs ER− heterogeneity.

### Summary
This analysis demonstrates a robust end-to-end RNA-seq workflow combining DESeq2 and pathway-level GSEA to recover canonical breast cancer biology, including proliferation-driven cell-cycle activation, estrogen-responsive transcriptional programs, and loss of normal adipose tissue identity. These results validate both the analytical pipeline and the biological interpretability of the findings.


### Machine Learning Classification (Tumor vs Normal)

Using the top 200 DE genes (ranked by adjusted p-value) as features, a Logistic Regression classifier achieved a 5-fold cross-validated ROC-AUC of 0.9999 using out-of-fold predicted probabilities. The resulting ROC curve is saved in results/ml_roc_curve.png.
This performance is consistent with the strong transcriptomic separation observed in DESeq2 and Hallmark GSEA (dominant E2F/G2M/MYC proliferation programs and breast-cancer–relevant estrogen response signatures).


Note: Tumor vs normal classification is an “easy” task in bulk RNA-seq because the transcriptomes differ strongly; high AUC is expected. More challenging extensions include subtype classification (e.g., ER+ vs ER−) or prediction tasks tied to clinical outcomes.

### Model interpretability (feature weights)

To interpret the classifier, we examined the fitted Logistic Regression coefficients (`results/ml_top25_logreg_genes.csv`).  
Genes with the strongest positive weights (tumor-associated) included **SEMA5B, COL10A1, MMP11, WISP1, and CCL11**, reflecting extracellular matrix remodeling, stromal interaction, and invasive tumor programs.

Conversely, genes with strong negative weights (normal-associated) included **HBB, ADIPOQ, LYVE1, CLEC3B, and C1QTNF9**, consistent with adipose, vascular, and normal tissue composition.

Canonical proliferation genes (**MKI67, TOP2A, CDC20, CDK1**) also contributed to classification but with more moderate coefficients, reflecting redundancy among tightly correlated cell-cycle features.

Overall, the model integrates both **cell-intrinsic proliferation signals** and **tissue/microenvironmental context**, yielding near-perfect tumor–normal discrimination (ROC-AUC = 0.9999).


## ER+ vs ER− Classification (Harder ML Task)

To demonstrate a more clinically relevant prediction task, we trained a classifier to distinguish ER-positive from ER-negative breast tumors. ER status was inferred from ESR1 expression by selecting extreme groups (bottom 40% = ER−, top 40% = ER+), excluding the middle 20% to reduce ambiguity (n = 493 ER−, n = 493 ER+). Direct receptor genes (ESR1/PGR) were excluded from features to avoid leakage.

Using a Logistic Regression model with 5-fold cross-validation and out-of-fold prediction, the classifier achieved a **ROC-AUC of 0.966**. The ROC curve is saved in `results/ml_er_roc_curve.png`.

This task is substantially more challenging than tumor–normal classification and indicates that the model captures broader ER-associated transcriptional programs beyond trivial marker genes.

### ER subtype model interpretability

To interpret the ER subtype classifier, we examined Logistic Regression coefficients (`results/ml_er_top25_logreg_genes.csv`).  
Genes with strong positive weights (ER+ associated) included **NUP210, CDC25C, CKAP2L, DTL, and FZD4**, reflecting nuclear organization, transcriptional regulation, and hormone-responsive cell-cycle control characteristic of ER-positive tumors.

In contrast, genes with strong negative weights (ER− associated) included **SIK2, AURKB, KIF23, ADAMTS5, and GSN**, capturing high mitotic activity, cytoskeletal remodeling, and extracellular matrix reorganization commonly observed in ER-negative (basal-like) breast cancers.

Importantly, direct hormone receptor genes (ESR1/PGR) were excluded from features, indicating that the classifier captures broader ER-associated transcriptional programs rather than trivial marker expression.

## Reproducibility

To recreate the analysis environment:

```bash
conda env create -f environment.yml
conda activate cancer-rnaseq-ml
```