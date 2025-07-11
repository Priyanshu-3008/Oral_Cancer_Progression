---
title: 'Differential Expression and Enrichment Analysis: Cancer vs Control'
author: "Priyanshu Sharma"
date: "`r Sys.Date()`"
output:
  html_document:
    toc: true
    toc_depth: 3
    theme: united
  pdf_document:
    toc: true
    toc_depth: '3'
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(
  echo = TRUE,
  message = FALSE,
  warning = FALSE,
  fig.width = 8,
  fig.height = 6
)
```

# Introduction

Oral cancer, a major subtype of head and neck cancers, remains one of the most aggressive and disfiguring malignancies worldwide. Despite advances in treatment, the 5-year survival rate remains low, partly due to delayed diagnosis and limited understanding of early-stage molecular changes. Oral squamous cell carcinoma (OSCC) often develops from premalignant lesions such as leukoplakia or hyperkeratosis through a multistep progression.

In this study, we analyze transcriptomic data from GEO dataset GSE227919, which includes samples from healthy controls, premalignant lesions, and OSCC cases. Our project specifically focuses on comparing cancer vs control groups to identify differentially expressed genes and enriched biological pathways. Through differential expression analysis and enrichment methods such as g:Profiler and GSEA, we aim to uncover molecular signatures that distinguish cancer from normal tissue, potentially aiding early detection and therapeutic targeting.

# 1. Load packages and data

```{r load-data}
library(tidyverse)
library(limma)
library(edgeR)
library(pheatmap)
library(GSVA)
library(gprofiler2)
library(clusterProfiler)
library(msigdbr)
library(enrichplot)
library(plotly)
library(DT)

# Read expression matrix and metadata
data_file <- "project_data.txt"
meta_file <- "project_metadata.csv"
expr <- read.table(data_file, header=TRUE, row.names=1) %>% as.matrix()
meta <- read.csv(meta_file)
```

# 2. Data preprocessing

```{r preprocess}
# Log2 transform (pseudocount)
log_expr <- log2(expr + 1)
# Subset Cancer vs Control
targets <- meta %>% select(sample, group) %>% filter(group %in% c("Control","Cancer"))
common <- intersect(colnames(log_expr), targets$sample)
log_expr <- log_expr[, common]
targets <- targets[match(common, targets$sample), ]
rownames(targets) <- targets$sample
annotation_col <- data.frame(Group=targets$group)
rownames(annotation_col) <- targets$sample
```

# 3. Differential Expression (limma-voom)

```{r de-analysis, message=FALSE, warning=FALSE}
# Design matrix
design <- model.matrix(~0 + factor(annotation_col$Group))
colnames(design) <- c("Cancer","Control")

# voom + linear model
v <- voom(expr[, common], design, plot=FALSE)
fit <- lmFit(v, design)
cont <- makeContrasts(Cancer_vs_Control = Cancer - Control, levels=design)
fit2 <- contrasts.fit(fit, cont) %>% eBayes()

# Extract results
deg <- topTable(fit2, coef="Cancer_vs_Control", adjust.method="BH", number=Inf)
deg <- tibble::rownames_to_column(deg, var="geneID")

# Classify as up/downregulated
deg <- deg %>%
  mutate(Direction = case_when(
    adj.P.Val < 0.05 & logFC > 1  ~ "Upregulated",
    adj.P.Val < 0.05 & logFC < -1 ~ "Downregulated",
    TRUE                          ~ "Not Significant"
  ))

# Top 20 DEGs (Up or Down by abs(logFC))
top_deg_table <- deg %>%
  filter(Direction %in% c("Upregulated", "Downregulated")) %>%
  arrange(desc(abs(logFC))) %>%
  slice(1:20)

# Show as interactive table
datatable(top_deg_table,
          options = list(pageLength = 10),
          rownames = FALSE,
          caption = "Top 20 Differentially Expressed Genes (Cancer vs Control)")
```

# 4. Heatmap of top 50 DEGs

```{r heatmap}
# Select top genes
top50 <- deg %>% arrange(adj.P.Val) %>% slice(1:50) %>% pull(geneID)
mat <- log_expr[top50, ]
mat_scaled <- t(scale(t(mat)))
# Plot
pheatmap(mat_scaled,
         annotation_col=annotation_col,
         cluster_rows=TRUE, cluster_cols=TRUE,
         scale="none",
         fontsize_row=6,
         show_rownames=TRUE,
         show_colnames=TRUE,
         color=colorRampPalette(c("navy","white","firebrick3"))(100))
```

# 5. Volcano Plot

```{r volcano}
# Add direction
deg$Direction <- ifelse(deg$adj.P.Val<0.05 & deg$logFC>1, "Up",
                        ifelse(deg$adj.P.Val<0.05 & deg$logFC< -1, "Down","NS"))

# Interactive
interactive_volcano <- ggplot(deg, aes(x=logFC, y=-log10(adj.P.Val),
               text=paste("Gene:", geneID,
                          "<br>logFC:", round(logFC,2),
                          "<br>adj.P.Val:", signif(adj.P.Val,3),
                          "<br>Direction:", Direction))) +
  geom_point(aes(color=Direction), alpha=0.7) +
  scale_color_manual(values=c(Up="red", Down="blue", NS="grey")) +
  geom_vline(xintercept=c(-1,1), linetype="dashed") +
  geom_hline(yintercept=-log10(0.05), linetype="dashed") +
  theme_minimal()

ggplotly(interactive_volcano, tooltip="text")
```

# 6. GO/KEGG/Reactome Enrichment (g:Profiler)

```{r gprofiler}
# Significant genes list
sig_genes <- deg %>% filter(adj.P.Val<0.05) %>% pull(geneID)
# Run g:Profiler
gostres <- gost(query=sig_genes, organism="hsapiens", correction_method="fdr")
# Interactive Manhattan
gostplot(gostres, interactive=TRUE, capped=TRUE)
```

# 7. GSEA (MSigDB C2)

```{r gsea}
hs_c2 <- msigdbr(species="Homo sapiens", category="C2") %>% select(gs_name, gene_symbol)
# Prepare ranked list
ranked <- deg %>% select(geneID, logFC) %>% deframe() %>% sort(decreasing=TRUE)
# Run GSEA
myGSEA <- GSEA(ranked, TERM2GENE=hs_c2, verbose=FALSE)
# Plot top pathways
gseaplot2(myGSEA, geneSetID=c(9,2,3,5), pvalue_table=FALSE)

# Bubble plot of top pathways
gsea_df <- as_tibble(myGSEA@result) %>%
  mutate(phenotype = if_else(NES > 0, "Cancer", "Control"),
         logp = -log10(p.adjust),
         label = stringr::str_wrap(ID, width = 50))

gsea_df_plot <- gsea_df %>%
  filter(!is.na(setSize), !is.na(logp), !is.na(NES)) %>%
  slice(1:20)

ggplot(gsea_df_plot, aes(x = phenotype, y = reorder(label, NES))) +
  geom_point(aes(size = setSize, color = NES, alpha = logp)) +
  scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  scale_alpha(range = c(0.6, 1)) +
  theme_bw() +
  theme(axis.text.y = element_text(size = 6)) +
  labs(title = "GSEA: Top Pathways", y = "Pathway", x = "Phenotype")
```

# 8. Discussion

In this transcriptome-based analysis of oral cancer using GEO dataset **GSE227919**, we focused on comparing gene expression profiles between oral squamous cell carcinoma (Cancer) and healthy control samples.

1. **Differential Gene Expression**
   - A total of thousands of genes were differentially expressed between cancer and control samples.
   - Notable upregulated genes included those involved in extracellular matrix remodeling and invasion, such as `MMP1` and `MMP10`.
   - Prominent downregulated genes included structural and immune-regulatory genes like `CGNL1`.
   - A volcano plot visualized the significant separation of up- and downregulated genes.

2. **Top DEGs**
   - Genes such as `PLAU`, `MMP1`, `IL24`, and `HOXC8` were among the most significantly altered, suggesting roles in tumor progression and epithelial plasticity.

3. **Functional Enrichment (GO/KEGG/Reactome)**
   - Gene Ontology (GO) analysis showed enrichment in processes like immune response, epithelial cell migration, and regulation of cell adhesion.
   - Reactome pathways such as *Scavenging by Class A Receptors* and *CD22-mediated BCR regulation* were significantly enriched.
   - KEGG pathways revealed altered drug metabolism, hinting at chemoresistance mechanisms in cancer cells.

4. **GSEA (Gene Set Enrichment Analysis)**
   - Gene sets associated with partial EMT, inflammatory signaling, and head and neck cancer signatures (e.g., *Rickman tumorigenesis*) were enriched in the cancer group.
   - A bubble plot highlighted the top enriched pathways, distinguishing cancer-promoting versus normal-protective gene sets.
