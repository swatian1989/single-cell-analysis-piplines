# Single-Cell Analysis by Dr Ali

A universal, tested, step-by-step guide to single-cell RNA-seq analysis, from a count matrix to annotated cell
types, differential expression, pathways and a ready-to-edit methods paragraph. It comes in two versions with
the same steps, settings and output names:

* **Python (Scanpy): [Single_Cell_Master_Guide.md](Single_Cell_Master_Guide.md)**
* **R (Seurat v5 + Bioconductor): [Single_Cell_Master_Guide_R.md](Single_Cell_Master_Guide_R.md)**

**R or Python?** Both give comparable results. Python scales better to very large datasets and has the
deep-learning tools (scVI, CellTypist); R has best-in-class doublet, ambient-RNA and bulk-statistics tools
and reads collaborators' Seurat objects directly. Section 0 of the R guide compares them step by step.

## What it covers

* **Any input format**: 10x folders and `.h5`, STARsolo, GEO matrix triplets and CSV/TSV/TXT tables (either
  orientation), `.h5ad` (finds raw counts, converts Ensembl IDs), `.loom`; `.zarr` in Python; Seurat and
  SingleCellExperiment `.rds` in R.
* **Any common data type**: scRNA-seq, snRNA-seq, CITE-seq, plate-based Smart-seq2, public atlases; pointers
  for scATAC/Multiome, spatial, TCR/BCR and multiplexed samples.
* **The full pipeline**: QC with per-sample adaptive thresholds → doublets → normalisation → highly variable
  genes → PCA → Harmony batch correction → UMAP and Leiden clustering → marker genes → cell-type annotation →
  composition → pseudobulk differential expression (DESeq2) → pathway and transcription-factor activity
  (decoupler) → saved results, settings, software versions and an auto-filled methods paragraph.
* **Optional modules**: automatic annotation (CellTypist / SingleR), sub-clustering, trajectories (PAGA /
  slingshot), cell-cell communication (LIANA), CITE-seq proteins, and pointers for scVI, RNA velocity,
  inferCNV, spatial, scATAC, TCR.

Every step has the code, a plain-English explanation, what to look at, and when to change the default.
You edit **one settings block**; the rest runs as written.

## Application: breast cancer (`breast_cancer_2027/`)

The Python pipeline applied, unchanged, to the 100,064-cell breast cancer atlas of Wu et al. 2021 (26 tumours), plus
a biomarker-axis map, single-cell copy-number inference and TCGA copy-number analysis of six metabolite–protein
"switch" axes. Report, code, figures and tables: [breast_cancer_2027/README.md](breast_cancer_2027/README.md).
Public data only; hypothesis-generating.

## Quick start

Python:

```bash
conda env create -f environment.yml
conda activate scanalysis
python -m ipykernel install --user --name scanalysis --display-name "Python (scanalysis)"
jupyter lab
```

R (or use your own R / RStudio: see Section 1 of the R guide):

```bash
conda env create -f environment_r.yml
conda activate scanalysis-r
Rscript -e 'options(repos = c(CRAN = "https://cloud.r-project.org")); BiocManager::install(c("SingleR", "scrapper", "slingshot"), update = FALSE, ask = FALSE)'
```

Then copy the guide's code blocks in order, starting with Step 0 (settings).

## Files

| File | What it is |
|---|---|
| `Single_Cell_Master_Guide.md` | The Python guide |
| `Single_Cell_Master_Guide_R.md` | The R guide |
| `environment.yml` / `environment_r.yml` | Exact software the guides were tested with |
| `tests/test_guide.py` | Re-runs both guides' own code blocks on public data |
| `breast_cancer_2027/` | The breast cancer application: report, scripts, figures and summary tables |
| `.gitignore` | Keeps data and results (`*.h5ad`, `*.rds`, `figures/`, `tables/`, ...) off GitHub |

## Testing

```bash
conda activate scanalysis
python tests/test_guide.py          # both guides (downloads ~25 MB of public data once; ~10 minutes)
python tests/test_guide.py python   # only the Python guide
python tests/test_guide.py r        # only the R guide (needs the scanalysis-r environment)
```

The script takes the code blocks straight from the guides, so what is tested is exactly what you read.
The tests and their results are described in Section 25 of the Python guide and Section 22 of the R guide.

## Author

Akbar Ali, Universiti Sains Malaysia. If you use these guides, please cite the tools they rely on
(Section 24 of the Python guide, Section 21 of the R guide).
