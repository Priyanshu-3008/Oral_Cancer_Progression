# Oral Cancer Progression Analysis

This project investigates transcriptomic changes across the progression of oral cancer using the publicly available **GEO dataset GSE227919**. It performs multi-stage **differential expression analysis** and **gene set enrichment analysis (GSEA)** to identify key molecular pathways driving the transition from normal tissue to dysplasia, HkNR (hyperkeratosis with no recurrence), and cancer.

## 🔬 Methods

- Differential Expression: `limma-voom`  
- Enrichment: `g:Profiler`, `GSEA` with MSigDB C2  
- Visualization: Volcano plots, heatmaps, interactive GSEA, NES trajectory plots  
- Tools: `R`, `clusterProfiler`, `gprofiler2`, `msigdbr`, `plotly`, `BiocParallel`

## 📊 Dataset

- **GEO Accession**: [GSE227919](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE227919)
- Groups: Control, Dysplasia, HkNR, Cancer

## 📈 Results

View the complete HTML report:  
👉 [Click to View RNAseq Analysis Report](https://yourusername.github.io/oral_cancer_progression)

*(Set up GitHub Pages to make this link work — see below.)*

## 📁 Files

| File | Description |
|------|-------------|
| `RNAseq_analysis.Rmd` | Main RMarkdown script |
| `RNAseq_analysis.html` | Knitted HTML report |
| `project_data.txt` | Gene expression matrix |
| `project_metadata.csv` | Sample metadata |

## 📄 License

This project is shared under the [MIT License](LICENSE).

---
