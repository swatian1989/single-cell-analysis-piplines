# Single-Cell RNA-seq Analysis: Universal Master Guide (R version)

**By Akbar Ali** · Universiti Sains Malaysia · Version 1.0 · 31 October 2026

The same pipeline as the Python guide ([Single_Cell_Master_Guide.md](Single_Cell_Master_Guide.md)), written in
R with **Seurat v5** and **Bioconductor**. Steps, settings and output files have the same names in both
guides, so you can move between them or check one against the other. This guide focuses on the R code; for
the longer background on *why* each step is done, it links to the matching section of the Python guide.

> **Tested.** Every core code block was run, top to bottom, on the same public datasets as the Python guide
> (see [Section 22](#22-how-this-guide-was-tested)). Blocks marked *untested here* need software that was not
> installed on the test machine.

---

## Contents

0. [R or Python?](#0-r-or-python)
1. [Install R and the packages (once)](#1-install-r-and-the-packages-once)
2. [Step 0: Settings](#2-step-0-settings)
3. [Step 1: Load the data (any format)](#3-step-1-load-the-data-any-format)
4. [Step 2: Quality control](#4-step-2-quality-control)
5. [Step 3: Doublets and ambient RNA](#5-step-3-doublets-and-ambient-rna)
6. [Step 4: Normalisation and cell cycle](#6-step-4-normalisation-and-cell-cycle)
7. [Step 5: Highly variable genes](#7-step-5-highly-variable-genes)
8. [Step 6: Scaling and PCA](#8-step-6-scaling-and-pca)
9. [Step 7: Batch correction with Harmony](#9-step-7-batch-correction-with-harmony)
10. [Step 8: Neighbours, UMAP and clustering](#10-step-8-neighbours-umap-and-clustering)
11. [Step 9: Marker genes](#11-step-9-marker-genes)
12. [Step 10: Cell-type annotation](#12-step-10-cell-type-annotation)
13. [Step 11: Cell-type composition](#13-step-11-cell-type-composition)
14. [Step 12: Differential expression (pseudobulk + DESeq2)](#14-step-12-differential-expression-pseudobulk--deseq2)
15. [Step 13: Pathways and transcription factors](#15-step-13-pathways-and-transcription-factors)
16. [Step 14: Save, export and write the methods](#16-step-14-save-export-and-write-the-methods)
17. [Advanced modules](#17-advanced-modules)
18. [Other data types in R](#18-other-data-types-in-r)
19. [Moving between R and Python](#19-moving-between-r-and-python)
20. [Troubleshooting (R)](#20-troubleshooting-r)
21. [References (R tools)](#21-references-r-tools)
22. [How this guide was tested](#22-how-this-guide-was-tested)

---

## 0. R or Python?

Both are excellent and give comparable results. The steps and the statistics are the same. Choose by **who
will run it and what they already know**, not by the data.

| | Python (Scanpy) | R (Seurat + Bioconductor) |
|---|---|---|
| Speed and memory on big data (> 200,000 cells) | Better | Slower, needs more RAM |
| Doublets, ambient RNA, bulk-style statistics | Good | Best-in-class (scDblFinder, SoupX, DESeq2, edgeR are native R) |
| Deep-learning tools (scVI, scANVI, totalVI, CellTypist) | Native | Through Python only |
| Reading Seurat `.rds` objects shared by others | Needs conversion | Native |
| Published tutorials for clinical / cancer single-cell papers | Many | Even more (Seurat is the most used toolkit) |
| Fits with this lab's other work (metabolomics, docking in Python) | Yes | Separate environment |

**Recommendation for this lab:** use **Python as the main pipeline** (you already work in Python, it scales
better, and the deep-learning tools are native), and use **R when a collaborator or reviewer works in Seurat,
or for the R-only tools** (scDblFinder, SoupX, CellChat, monocle3). Because both guides use the same steps and
settings, running a key result in both is also a cheap robustness check for a paper.

**What is the same and what differs between the two guides**

| Step | Python guide | This R guide | Same result? |
|---|---|---|---|
| Reader | `read_any()` | `read_any()` (also reads Seurat/SCE `.rds` directly) | Same cells and genes |
| QC | per-sample MAD rule | identical rule (raw MAD, same cut-offs) | Same cells removed |
| Doublets | Scrublet | **scDblFinder** (best in the Xi & Li 2021 benchmark) | Similar, not identical |
| Normalisation | counts to 10,000 + log1p | `NormalizeData(LogNormalize)`: the same formula | Same values |
| HVGs | Seurat-v3 method (vst) | Seurat `vst`: the same method | Very similar |
| PCA | no scaling | Seurat default: genes scaled first | Similar structure |
| Batch correction | Harmony | Harmony | Similar |
| Graph / clustering | kNN + Leiden | SNN graph + Leiden | Similar clusters, different numbers |
| Markers | Wilcoxon | Wilcoxon | Very similar |
| Annotation | marker-hit rule + CellTypist | marker-hit rule + **SingleR** | Similar |
| DE | pseudobulk + PyDESeq2 | pseudobulk + **DESeq2** (the original) | Very similar |
| Pathways | decoupler | decoupleR (same authors, same resources) | Very similar |

### How to read the code blocks

* **Core `r` blocks**: run them in order, top to bottom, in RStudio (or an R notebook / Quarto document).
* **Blocks that start with `# OPTIONAL`**: only if you need that analysis.
* `bash` blocks go in the Terminal.
* Figures go to `figures/`, tables to `tables/`, objects to `objects/` in your project folder.

The six golden rules, the glossary and the description of the upstream tools (FASTQ → count matrix) are the
same as in the Python guide: [Section 0](Single_Cell_Master_Guide.md#0-read-this-first) and
[Section 1](Single_Cell_Master_Guide.md#1-from-sequencer-to-count-matrix).

---

## 1. Install R and the packages (once)

**Option A: conda (what this guide was tested with; no administrator rights needed).** In the Terminal,
inside this repository folder:

```bash
conda env create -f environment_r.yml        # creates "scanalysis-r" (10-20 minutes)
conda activate scanalysis-r
# SingleR, scrapper and slingshot are not built by bioconda for Apple-Silicon Macs: compile them once (~10 min)
Rscript -e 'options(repos = c(CRAN = "https://cloud.r-project.org")); BiocManager::install(c("SingleR", "scrapper", "slingshot"), update = FALSE, ask = FALSE)'
R                                            # or: open RStudio from this environment
```

**Option B: your normal R / RStudio installation (R ≥ 4.4).** In R:

```r
# OPTIONAL (untested here): install from CRAN and Bioconductor instead of conda
install.packages(c("Seurat", "harmony", "hdf5r", "data.table", "ggplot2", "patchwork", "jsonlite", "R.utils", "leidenbase"))
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("scDblFinder", "scuttle", "DESeq2", "SingleR", "scrapper", "celldex", "decoupleR",
                       "OmnipathR", "zellkonverter", "slingshot", "SingleCellExperiment"))
```

On Windows and on a Mac with R from CRAN, Bioconductor installs ready-made binaries, so this takes minutes.
`scrapper` is only "suggested" by SingleR but SingleR 2.12 stops without it, so it is listed explicitly.

*What each package does:* **Seurat** (the core toolkit) · **harmony** (batch correction) · **hdf5r** (reads
10x `.h5` and `.loom` files) · **data.table** + **R.utils** (fast reading of big, gzipped CSV/TSV tables) ·
**leidenbase** (Leiden clustering) · **scDblFinder** (doublets) · **DESeq2** (pseudobulk DE) · **SingleR** +
**celldex** (automatic labels from reference data) · **decoupleR** + **OmnipathR** (pathway / TF activity) ·
**zellkonverter** (reads `.h5ad` files) · **slingshot** (trajectories) · **ggplot2 / patchwork** (plots).

---

## 2. Step 0: Settings

**This is the only block you must edit.** The names and meanings are exactly those of the Python guide
([Section 4](Single_Cell_Master_Guide.md#4-step-0-settings)); R uses `NULL` where Python uses `None`.

```r
# ===== STEP 0 · SETTINGS: the only block you need to edit =====
PROJECT_DIR <- file.path(Sys.getenv("HOME"), "Desktop", "my_sc_project_R")   # results are written here

# One entry per sample:  sample name = path to its data (absolute, or relative to PROJECT_DIR).
# Formats: 10x folder, 10x .h5, .h5ad, .loom, (prefix_)matrix.mtx(.gz), .csv/.tsv/.txt(.gz),
#          and R objects: Seurat or SingleCellExperiment .rds
SAMPLES <- c(
  P01_tumour = "data/P01_tumour/filtered_feature_bc_matrix",
  P01_normal = "data/P01_normal/filtered_feature_bc_matrix.h5"
)
SAMPLE_COL   <- NULL      # one file holding MANY samples? name its metadata column that says which sample
SAMPLE_SHEET <- NULL      # optional CSV, one row per sample, first column = sample name

SPECIES  <- "human"       # "human" or "mouse"
MATERIAL <- "cells"       # "cells" (scRNA-seq) or "nuclei" (snRNA-seq)

# --- Quality control (Step 2) ---
MIN_GENES  <- 200         # a cell must express at least this many genes (lower it for shallow data)
NMADS      <- 5           # outlier = more than 5 MADs from ITS OWN sample's median
MT_NMADS   <- 3           # stricter for mitochondrial %
MAX_PCT_MT <- if (MATERIAL == "cells") 20 else 5   # hard ceiling on mitochondrial %
MIN_CELLS_PER_GENE <- 3   # drop genes detected in fewer cells than this

# --- Main analysis (Steps 5-8) ---
BATCH_KEY   <- "sample"   # what to correct for; NULL = no correction
N_HVG       <- 2000
N_PCS       <- 30
N_NEIGHBORS <- 15
RESOLUTIONS <- c(0.2, 0.4, 0.6, 0.8, 1.0, 1.5)
RESOLUTION  <- 0.6

# --- Comparing conditions (Steps 11-12). Leave CONDITION_KEY <- NULL to skip ---
CONDITION_KEY <- NULL     # metadata column with the groups, e.g. "condition"
CONTRAST      <- NULL     # c(test group, reference group), e.g. c("tumour", "normal")
REPLICATE_KEY <- "sample" # the biological replicate (patient / animal / sample)
DE_COVARIATES <- character(0)   # e.g. c("patient") for a paired design
MIN_CELLS_PSEUDOBULK <- 10

SEED <- 0
```

---

## 3. Step 1: Load the data (any format)

### 3.1 Packages and helper functions

**What it does.** Loads the packages, creates the output folders, and defines `save_plot()` (saves a plot as a
300-dpi PNG and shows it) and `genes_present()` (finds genes whatever the capital letters, so `CD3E` also
finds mouse `Cd3e`). Note: in R, matrices are **genes x cells** (the opposite of Python).

```r
suppressPackageStartupMessages({
  library(SingleCellExperiment)     # load BEFORE Seurat: both have a function called Assays(); Seurat's must win
  library(Seurat)
  library(Matrix)
  library(ggplot2)
  library(patchwork)
})
options(future.globals.maxSize = 8 * 1024^3)    # allow Seurat to handle large objects
set.seed(SEED)

FIG <- file.path(PROJECT_DIR, "figures")
TAB <- file.path(PROJECT_DIR, "tables")
OBJ <- file.path(PROJECT_DIR, "objects")
for (d in c(FIG, TAB, OBJ)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

save_plot <- function(p, name, width = 7, height = 5) {
  # Save a ggplot / patchwork plot as figures/<name>.png and show it
  ggsave(file.path(FIG, paste0(name, ".png")), p, width = width, height = height, dpi = 300, bg = "white")
  if (interactive()) print(p)
}

genes_present <- function(names, genes) {
  # The genes from `genes` that exist in `names`, matched case-blind, in the data's own spelling
  lookup <- setNames(names, toupper(names))
  unname(lookup[intersect(toupper(genes), names(lookup))])
}

sessionInfo()$otherPkgs$Seurat$Version            # the Seurat version in use
```

### 3.2 The universal reader

**What it does.** `read_any(path)` recognises the format, reads it and returns three things: the RNA count
matrix (genes x cells), any extra modalities (CITE-seq antibodies → `ADT`, hashtags → `HTO`, CRISPR guides),
and any per-cell metadata stored in the file. It understands the same formats as the Python reader, plus
Seurat and SingleCellExperiment `.rds` objects (it takes their raw counts and metadata). `.zarr` is the one
format R cannot read well: convert it to `.h5ad` in Python first.

```r
# ---- Universal reader: (almost) any single-cell file -> list(rna = genes x cells, other = list(), meta) ----
GENE_HINTS <- c("ACTB", "GAPDH", "B2M", "MALAT1", "EEF1A1", "TMSB4X", "FTL", "FTH1", "RPLP0", "MT-CO1")
SYMBOL_COLUMNS <- c("gene_symbols", "gene_symbol", "gene_name", "gene_names", "feature_name", "symbol", "name")

sample_values <- function(m, n = 2000) {
  sub <- m[, seq_len(min(n, ncol(m))), drop = FALSE]
  if (inherits(sub, "sparseMatrix")) as(sub, "CsparseMatrix")@x else as.vector(sub)
}

looks_like_counts <- function(m) {
  # TRUE if the values are raw counts: whole numbers and never negative
  v <- sample_values(m)
  length(v) == 0 || (min(v) >= 0 && all(abs(v - round(v)) < 1e-6))
}

data_state <- function(m) {
  # "counts", "log" (already log-normalised) or "normalised" (e.g. TPM / CPM, not logged)
  if (looks_like_counts(m)) return("counts")
  v <- sample_values(m)
  if (min(v) < 0) stop("Negative values: these data are scaled/z-scored. Find the raw counts.")
  if (max(v) < 50) "log" else "normalised"
}

gene_like <- function(labels) {
  # How many labels look like gene symbols or Ensembl gene IDs
  up <- toupper(as.character(labels))
  sum(up %in% GENE_HINTS) + sum(grepl("^ENS[A-Z]*G[0-9]{6,}", up))
}

first_existing <- function(folder, names) {
  hits <- file.path(folder, names)
  hits <- hits[file.exists(hits)]
  if (length(hits)) hits[1] else NULL
}

as_sparse <- function(m) as(as(m, "CsparseMatrix"), "generalMatrix")

split_types <- function(m, types) {
  # Keep Gene Expression as RNA; other feature types become extra modalities (ATAC peaks are left out)
  if (!"Gene Expression" %in% types) return(list(rna = m, other = list()))
  other <- list()
  for (t in setdiff(unique(types), "Gene Expression")) {
    if (t == "Peaks") { message("  ", sum(types == t), " ATAC peaks left out (analyse ATAC separately)"); next }
    name <- switch(t, "Antibody Capture" = "ADT", "Multiplexing Capture" = "HTO",
                   "CRISPR Guide Capture" = "CRISPR", gsub("[^A-Za-z0-9]", "", t))
    other[[name]] <- m[types == t, , drop = FALSE]
    message("  ", sum(types == t), " '", t, "' features moved to assay ", name)
  }
  list(rna = m[types == "Gene Expression", , drop = FALSE], other = other)
}

use_symbols <- function(m, gene_table) {
  # If row names are Ensembl IDs and a symbol column exists, switch to symbols
  if (mean(grepl("^ENS[A-Z]*G[0-9]{6,}", rownames(m))) < 0.5) return(m)
  col <- intersect(SYMBOL_COLUMNS, colnames(gene_table))[1]
  if (is.na(col)) {
    message("  NOTE: genes are Ensembl IDs and no symbol column was found; mito genes and markers will not match.")
    return(m)
  }
  sym <- as.character(gene_table[[col]])
  ok <- !is.na(sym) & sym != "" & sym != "nan"
  rownames(m) <- make.unique(ifelse(ok, sym, rownames(m)))
  m
}

read_mtx_folder <- function(folder, prefix = "") {
  # matrix.mtx + features.tsv/genes.tsv + barcodes.tsv (10x, STARsolo, GEO; gzipped or not)
  mtx  <- first_existing(folder, paste0(prefix, c("matrix.mtx.gz", "matrix.mtx")))
  feat <- first_existing(folder, paste0(prefix, c("features.tsv.gz", "features.tsv", "genes.tsv.gz", "genes.tsv")))
  bcs  <- first_existing(folder, paste0(prefix, c("barcodes.tsv.gz", "barcodes.tsv")))
  if (is.null(mtx) || is.null(feat) || is.null(bcs))
    stop(folder, ": need ", prefix, "matrix.mtx, ", prefix, "features.tsv (or genes.tsv) and ", prefix, "barcodes.tsv")
  m <- as_sparse(Matrix::readMM(mtx))
  feat <- read.delim(feat, header = FALSE, colClasses = "character")
  bcs <- read.delim(bcs, header = FALSE, colClasses = "character")[[1]]
  if (all(dim(m) == c(length(bcs), nrow(feat)))) {
    m <- t(m)                                   # stored cells x genes: turn it round
  } else if (!all(dim(m) == c(nrow(feat), length(bcs)))) {
    stop("matrix is ", paste(dim(m), collapse = " x "), " but there are ", length(bcs), " barcodes and ",
         nrow(feat), " features")
  }
  rownames(m) <- make.unique(if (ncol(feat) > 1) feat[[2]] else feat[[1]])
  colnames(m) <- bcs
  types <- if (ncol(feat) > 2) feat[[3]] else rep("Gene Expression", nrow(feat))
  c(split_types(m, types), list(meta = NULL))
}

read_table_any <- function(path, genes_are_rows = NULL) {
  # A count table in CSV / TSV / TXT (optionally .gz), as often shared on GEO
  df <- data.table::fread(path, data.table = FALSE)          # separator detected automatically
  if (names(df)[1] %in% c("V1", "") || !is.numeric(df[[1]])) {
    rownames(df) <- make.unique(as.character(df[[1]]))     # first column = gene (or cell) names
    df[[1]] <- NULL
  }
  text <- !vapply(df, is.numeric, logical(1))
  extra <- df[, text, drop = FALSE]                         # e.g. a gene-description or cluster column
  df <- df[, !text, drop = FALSE]
  if (is.null(genes_are_rows)) genes_are_rows <- gene_like(rownames(df)) >= gene_like(colnames(df))
  m <- as_sparse(as.matrix(df))
  if (genes_are_rows) list(rna = m, other = list(), meta = NULL)
  else list(rna = t(m), other = list(), meta = if (ncol(extra)) extra else NULL)
}

read_h5ad_raw <- function(path) {
  # The .raw part of an .h5ad file as genes x cells (zellkonverter's R reader does not read .raw)
  f <- hdf5r::H5File$new(path, mode = "r")
  on.exit(f$close_all())
  if (!f$exists("raw") || !f$exists("raw/X")) return(NULL)
  X <- f[["raw/X"]]
  genes <- f[[paste0("raw/var/", hdf5r::h5attr(f[["raw/var"]], "_index"))]]$read()
  if (inherits(X, "H5D")) return(`rownames<-`(as_sparse(t(X$read())), genes))   # dense raw
  shape <- hdf5r::h5attr(X, "shape")                                          # n_cells, n_genes
  csr <- hdf5r::h5attr(X, "encoding-type") == "csr_matrix"
  m <- Matrix::sparseMatrix(i = X[["indices"]]$read() + 1, p = X[["indptr"]]$read(),
                            x = as.numeric(X[["data"]]$read()), dims = if (csr) rev(shape) else shape)
  m <- if (csr) m else t(m)              # a CSR cells x genes matrix is a CSC genes x cells matrix
  rownames(m) <- make.unique(genes)
  as_sparse(m)
}

read_h5ad_any <- function(path) {
  # An .h5ad file, read natively in R (no Python). Finds the raw counts even when X holds normalised values.
  sce <- suppressWarnings(zellkonverter::readH5AD(path, reader = "R", verbose = FALSE))   # harmless HDF5 warnings
  meta <- as.data.frame(colData(sce))
  genes <- as.data.frame(rowData(sce))
  choice <- assayNames(sce)[vapply(assayNames(sce), function(a) looks_like_counts(assay(sce, a)), logical(1))]
  choice <- choice[order(choice != "X", !grepl("count", choice, ignore.case = TRUE))]   # X first, then "counts"
  if (length(choice)) {
    if (choice[1] != "X") message("  X is not raw counts; using layer '", choice[1], "'")
    return(list(rna = use_symbols(as_sparse(assay(sce, choice[1])), genes), other = list(), meta = meta))
  }
  raw <- read_h5ad_raw(path)
  if (!is.null(raw) && looks_like_counts(raw)) {
    message("  X is not raw counts; using .raw")
    colnames(raw) <- colnames(sce)
    return(list(rna = raw, other = list(), meta = meta))
  }
  list(rna = use_symbols(as_sparse(assay(sce, "X")), genes), other = list(), meta = meta)   # processed values only
}

read_loom_any <- function(path) {
  # A .loom file (velocyto or anndata style)
  f <- hdf5r::H5File$new(path, mode = "r")
  on.exit(f$close_all())
  m <- f[["matrix"]]$read()
  ra <- names(f[["row_attrs"]]); ca <- names(f[["col_attrs"]])
  genes <- f[["row_attrs"]][[intersect(c("Gene", "var_names", "gene_names"), ra)[1]]]$read()
  cells <- f[["col_attrs"]][[intersect(c("CellID", "obs_names", "cell_names"), ca)[1]]]$read()
  if (nrow(m) == length(cells)) m <- t(m)       # make it genes x cells
  m <- as_sparse(m)
  rownames(m) <- make.unique(as.character(genes)); colnames(m) <- as.character(cells)
  list(rna = m, other = list(), meta = NULL)
}

read_rds_any <- function(path) {
  # A Seurat, SingleCellExperiment or plain matrix .rds file: take the raw counts and metadata
  x <- readRDS(path)
  if (inherits(x, "Seurat")) {
    if (inherits(x[[DefaultAssay(x)]], "Assay5")) x <- JoinLayers(x)
    return(list(rna = as_sparse(LayerData(x, layer = "counts")), other = list(), meta = x@meta.data))
  }
  if (inherits(x, "SingleCellExperiment")) {
    m <- if ("counts" %in% assayNames(x)) counts(x) else assay(x, 1)
    return(list(rna = as_sparse(m), other = list(), meta = as.data.frame(colData(x))))
  }
  if (inherits(x, c("Matrix", "matrix"))) return(list(rna = as_sparse(x), other = list(), meta = NULL))
  stop("Unknown object in ", basename(path), ": ", class(x)[1])
}

read_any <- function(path, genes_are_rows = NULL) {
  # Read one sample (or one multi-sample file) in any common single-cell format
  path <- path.expand(path)
  if (!grepl("^/", path)) path <- file.path(PROJECT_DIR, path)
  if (!file.exists(path)) stop("Not found: ", path)
  name <- tolower(basename(path))

  if (grepl("\\.zarr$", name)) {
    stop(basename(path), ": convert .zarr to .h5ad in Python first: anndata.read_zarr(p).write_h5ad('x.h5ad')")
  } else if (dir.exists(path)) {
    if (!is.null(first_existing(path, c("matrix.mtx.gz", "matrix.mtx")))) {
      x <- read_mtx_folder(path)
    } else {                                     # Cell Ranger / STARsolo output folders
      for (sub in c("filtered_feature_bc_matrix.h5", "outs/filtered_feature_bc_matrix.h5",
                    "sample_filtered_feature_bc_matrix.h5", "filtered_feature_bc_matrix",
                    "outs/filtered_feature_bc_matrix", "Solo.out/Gene/filtered", "Gene/filtered")) {
        if (file.exists(file.path(path, sub))) return(read_any(file.path(path, sub), genes_are_rows))
      }
      stop(path, ": no matrix.mtx or filtered_feature_bc_matrix inside")
    }
  } else if (grepl("\\.h5ad$", name)) {
    x <- read_h5ad_any(path)
  } else if (grepl("\\.h5$", name)) {
    h5 <- Read10X_h5(path)                       # Cell Ranger, Space Ranger, CellBender
    if (!is.list(h5)) h5 <- list("Gene Expression" = h5)
    x <- c(split_types(do.call(rbind, h5), rep(names(h5), vapply(h5, nrow, integer(1)))), list(meta = NULL))
  } else if (grepl("\\.loom$", name)) {
    x <- read_loom_any(path)
  } else if (grepl("matrix\\.mtx(\\.gz)?$", name)) {   # GEO style: GSM123_matrix.mtx.gz + GSM123_barcodes...
    x <- read_mtx_folder(dirname(path), prefix = sub("matrix\\.mtx(\\.gz)?$", "", basename(path)))
  } else if (grepl("\\.(csv|tsv|txt)(\\.gz)?$", name)) {
    x <- read_table_any(path, genes_are_rows)
  } else if (grepl("\\.rds$", name)) {
    x <- read_rds_any(path)
  } else {
    stop("Unknown format: ", basename(path))
  }
  rownames(x$rna) <- make.unique(rownames(x$rna))
  x
}
```

### 3.3 Load every sample and merge

**What it does.** Reads each sample, prefixes every cell name with its sample, keeps only genes measured in
**every** sample, and builds one Seurat object. Extra modalities present in every sample (e.g. CITE-seq
antibodies) become extra assays. Then the sample sheet is added, the values are checked (raw counts?), and
the settings are checked against the metadata columns.

*Why one object instead of Seurat's `merge()`?* `merge()` keeps the union of genes (missing genes become zeros,
which looks like "not expressed") and splits the data into one layer per sample. Joining on shared genes is
the honest choice and keeps the object simple.

```r
rna <- list(); metas <- list(); others <- list()
for (name in names(SAMPLES)) {
  message("Reading ", name, " ...")
  x <- read_any(SAMPLES[[name]])
  cells <- paste0(name, "_", colnames(x$rna))
  colnames(x$rna) <- cells
  meta <- if (is.null(x$meta)) data.frame(row.names = cells) else { rownames(x$meta) <- cells; x$meta }
  meta$sample <- if (is.null(SAMPLE_COL)) name else as.character(meta[[SAMPLE_COL]])
  for (k in names(x$other)) colnames(x$other[[k]]) <- cells
  rna[[name]] <- x$rna; metas[[name]] <- meta; others[[name]] <- x$other
  message("  ", ncol(x$rna), " cells x ", nrow(x$rna), " genes")
}

shared <- Reduce(intersect, lapply(rna, rownames))
all_genes <- Reduce(union, lapply(rna, rownames))
if (length(all_genes) > length(shared))
  message("NOTE: ", length(all_genes) - length(shared), " genes were not measured in every sample and were dropped.")
counts <- do.call(cbind, lapply(rna, function(m) m[shared, , drop = FALSE]))
meta <- as.data.frame(data.table::rbindlist(lapply(metas, function(d) cbind(.cell = rownames(d), d)), fill = TRUE))
rownames(meta) <- meta$.cell; meta$.cell <- NULL
rm(rna)

obj <- CreateSeuratObject(counts = counts, meta.data = meta[colnames(counts), , drop = FALSE])
rm(counts)
common_mod <- Reduce(intersect, lapply(others, names))
for (k in common_mod) {
  feats <- Reduce(intersect, lapply(others, function(o) rownames(o[[k]])))
  obj[[k]] <- CreateAssay5Object(counts = do.call(cbind, lapply(others, function(o) o[[k]][feats, , drop = FALSE])))
}
lost <- setdiff(Reduce(union, lapply(others, names)), common_mod)
if (length(lost)) message("NOTE: ", paste(lost, collapse = ", "), " not present in every sample, so dropped.")

if (!is.null(SAMPLE_SHEET)) {
  sheet_path <- if (grepl("^/", SAMPLE_SHEET)) SAMPLE_SHEET else file.path(PROJECT_DIR, SAMPLE_SHEET)
  sheet <- read.csv(sheet_path, colClasses = "character", check.names = FALSE)
  names(sheet)[1] <- "sample"
  missing <- setdiff(unique(obj$sample), sheet$sample)
  if (length(missing)) message("WARNING: samples missing from the sample sheet: ", paste(missing, collapse = ", "))
  extra <- sheet[match(obj$sample, sheet$sample), setdiff(names(sheet), "sample"), drop = FALSE]
  for (col in names(extra)) obj@meta.data[[col]] <- extra[[col]]
}

DATA_STATE <- data_state(LayerData(obj, layer = "counts"))
message(ncol(obj), " cells x ", nrow(obj), " genes from ", length(unique(obj$sample)), " sample(s); data = ", DATA_STATE)
message("Assays: ", paste(Assays(obj), collapse = ", "))
print(table(obj$sample))

for (key in c(BATCH_KEY, CONDITION_KEY, REPLICATE_KEY, DE_COVARIATES)) {   # catch typos in the settings now
  if (!key %in% colnames(obj@meta.data))
    stop("'", key, "' is not a metadata column. Columns are: ", paste(colnames(obj@meta.data), collapse = ", "))
}
```

**What to look at:** the same as in the Python guide
([Section 5.3](Single_Cell_Master_Guide.md#53-load-every-sample-and-merge)): cells per sample, number of genes,
and `data = counts`. Seurat may warn that feature names with underscores were changed to dashes: harmless.

---

## 4. Step 2: Quality control

**Why and what the metrics mean:** [Python guide, Section 6](Single_Cell_Master_Guide.md#6-step-2-quality-control)
(same metrics, same per-sample MAD rule, same tissue and shallow-data warnings). Seurat names the first two
metrics `nCount_RNA` (total counts) and `nFeature_RNA` (genes detected).

```r
counts <- LayerData(obj, layer = "counts")
up <- toupper(rownames(obj))
mt_genes   <- rownames(obj)[grepl("^MT[-:]", up)]                                   # mitochondrial
ribo_genes <- rownames(obj)[grepl("^(RP[SL][0-9]+[A-Z]?[0-9]*|RPLP[0-2]|RPSA)$", up)]  # ribosomal proteins
hb_genes   <- rownames(obj)[grepl("^HB[ABDEGMQZ]([0-9]|-|$)", up, perl = TRUE)]       # haemoglobins
message("genes found - mitochondrial: ", length(mt_genes), " | ribosomal: ", length(ribo_genes),
        " | haemoglobin: ", length(hb_genes))
if (!length(mt_genes)) message("WARNING: no mitochondrial genes found. Check gene names or the pattern above.")

total <- Matrix::colSums(counts)
pct_of <- function(genes) if (length(genes)) 100 * Matrix::colSums(counts[genes, , drop = FALSE]) / total else 0
top_n_pct <- function(m, n = 20) {
  # % of each cell's counts that come from its n most-expressed genes
  m <- as(m, "CsparseMatrix")
  vapply(seq_len(ncol(m)), function(j) {
    if (m@p[j + 1] == m@p[j]) return(NA_real_)
    v <- m@x[(m@p[j] + 1):m@p[j + 1]]
    100 * sum(sort(v, decreasing = TRUE)[seq_len(min(n, length(v)))]) / sum(v)
  }, numeric(1))
}
obj$pct_mt <- pct_of(mt_genes)
obj$pct_ribo <- pct_of(ribo_genes)
obj$pct_hb <- pct_of(hb_genes)
obj$pct_top20 <- top_n_pct(counts)
rm(counts)

p <- VlnPlot(obj, features = c("nFeature_RNA", "nCount_RNA", "pct_mt", "pct_hb"), group.by = "sample",
             pt.size = 0, ncol = 4) & theme(axis.title.x = element_blank())
save_plot(p, "02_qc_violin_before", width = 14, height = 4.5)
p <- ggplot(obj@meta.data, aes(nCount_RNA, nFeature_RNA, colour = pct_mt)) +
  geom_point(size = 0.4) + scale_x_log10() + scale_y_log10() + scale_colour_viridis_c() + theme_classic()
save_plot(p, "02_qc_counts_vs_genes")
```

```r
mad_outlier <- function(x, sample, nmads, upper_only = FALSE) {
  # Flag cells more than `nmads` MADs from the median of THEIR OWN sample (raw MAD, as in the Python guide)
  flag <- logical(length(x))
  for (s in unique(sample)) {
    i <- which(sample == s)
    med <- median(x[i]); d <- mad(x[i], constant = 1)
    if (d == 0) next                              # no spread (e.g. 0 % mito in nuclei): nothing to flag
    flag[i] <- x[i] > med + nmads * d | (!upper_only & x[i] < med - nmads * d)
  }
  flag
}

if (DATA_STATE == "counts") {
  obj$outlier <- mad_outlier(log1p(obj$nCount_RNA), obj$sample, NMADS) |
    mad_outlier(log1p(obj$nFeature_RNA), obj$sample, NMADS) |
    mad_outlier(obj$pct_top20, obj$sample, NMADS) |
    obj$nFeature_RNA < MIN_GENES
  obj$mt_outlier <- mad_outlier(obj$pct_mt, obj$sample, MT_NMADS, upper_only = TRUE) | obj$pct_mt > MAX_PCT_MT
} else {                                          # processed values: only the gene-number rule is meaningful
  obj$outlier <- obj$nFeature_RNA < MIN_GENES
  obj$mt_outlier <- FALSE
}

samples <- sort(unique(obj$sample))
count_by <- function(v, s = obj$sample) as.vector(tapply(v, factor(s, levels = samples), sum))
qc_table <- data.frame(sample = samples, cells_loaded = count_by(rep(1, ncol(obj))),
                       low_quality = count_by(obj$outlier), high_mito = count_by(obj$mt_outlier))
obj <- subset(obj, cells = colnames(obj)[!(obj$outlier | obj$mt_outlier)])
qc_table$after_qc <- as.vector(table(factor(obj$sample, levels = samples)))

keep <- Matrix::rowSums(LayerData(obj, layer = "counts") > 0) >= MIN_CELLS_PER_GENE
if (MATERIAL == "nuclei") keep <- keep & !rownames(obj) %in% mt_genes   # mito reads in nuclei = contamination
obj[["RNA"]] <- subset(obj[["RNA"]], features = rownames(obj)[keep])
print(qc_table)
message("Kept ", ncol(obj), " cells and ", nrow(obj), " genes.")
```

---

## 5. Step 3: Doublets and ambient RNA

### 5.1 Doublets with scDblFinder

**What it does.** scDblFinder (Germain et al. 2021) makes artificial doublets, then trains a classifier to
tell them from real cells; it was the most accurate method in an independent benchmark (Xi & Li 2021). It
runs **per sample**; samples with fewer than 100 cells after QC are skipped and named (report them). Background:
[Python guide, Section 7](Single_Cell_Master_Guide.md#7-step-3-doublets-and-ambient-rna).

```r
MIN_CELLS_DOUBLETS <- 100     # fewer cells than this: too few to model doublets, sample is skipped
skipped_doublets <- character(0)
obj$doublet_score <- NA_real_
obj$predicted_doublet <- FALSE
if (DATA_STATE == "counts") {
  for (s in unique(obj$sample)) {                       # doublets form within one sample
    cells <- colnames(obj)[obj$sample == s]
    if (length(cells) < MIN_CELLS_DOUBLETS) {
      message("  ", s, ": only ", length(cells), " cells, scDblFinder skipped for this sample (report it)")
      skipped_doublets <- c(skipped_doublets, s)
      next
    }
    set.seed(SEED)
    sce <- scDblFinder::scDblFinder(SingleCellExperiment(list(counts = LayerData(obj, layer = "counts")[, cells])))
    obj@meta.data[cells, "doublet_score"] <- sce$scDblFinder.score
    obj@meta.data[cells, "predicted_doublet"] <- sce$scDblFinder.class == "doublet"
  }
  print(round(100 * tapply(obj$predicted_doublet, obj$sample, mean), 1))
  qc_table$doublets <- count_by(obj$predicted_doublet)
  obj <- subset(obj, cells = colnames(obj)[!obj$predicted_doublet])
} else {
  message("Data are not raw counts: doublet detection skipped (scDblFinder needs counts).")
}
qc_table$final <- as.vector(table(factor(obj$sample, levels = samples)))
write.csv(qc_table, file.path(TAB, "02_qc_summary_per_sample.csv"), row.names = FALSE)
saveRDS(obj, file.path(OBJ, "02_after_qc.rds"))
print(qc_table)
```

### 5.2 Ambient RNA with SoupX

SoupX (Young & Behjati 2020) estimates the "soup" from empty droplets and subtracts it. It needs the Cell
Ranger `outs/` folder (raw **and** filtered matrices). Run it **before** this pipeline and save the corrected
matrix as a sample input. *Untested here (needs a raw matrix):*

```r
# OPTIONAL (untested here): ambient RNA removal with SoupX, one sample at a time
library(SoupX)
soup <- load10X("data/P01_tumour/outs")                 # reads raw + filtered matrices and Cell Ranger clusters
soup <- autoEstCont(soup)                               # estimates the contamination fraction
corrected <- adjustCounts(soup, roundToInt = TRUE)      # corrected counts, genes x cells
saveRDS(corrected, "data/P01_tumour_soupx.rds")         # then in Step 0:  P01_tumour = "data/P01_tumour_soupx.rds"
```

---

## 6. Step 4: Normalisation and cell cycle

**What it does.** `NormalizeData` scales every cell to 10,000 counts and takes `log(1 + x)`: the same formula
as the Python guide. The raw counts stay in the `counts` layer; the normalised values go to the `data` layer.
Then cells are scored for S and G2/M phase (same gene lists as the Python guide, old and new names).
Background: [Python guide, Section 8](Single_Cell_Master_Guide.md#8-step-4-normalisation-and-cell-cycle).

```r
if (DATA_STATE == "counts") {
  obj <- NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 1e4, verbose = FALSE)
} else if (DATA_STATE == "normalised") {            # e.g. TPM / CPM tables: only the log is missing
  LayerData(obj, layer = "data") <- log1p(LayerData(obj, layer = "counts"))
} else {                                            # already log-normalised
  LayerData(obj, layer = "data") <- LayerData(obj, layer = "counts")
}

S_GENES <- c("MCM5", "PCNA", "TYMS", "FEN1", "MCM2", "MCM4", "RRM1", "UNG", "GINS2", "MCM6", "CDCA7", "DTL",
             "PRIM1", "UHRF1", "MLF1IP", "CENPU", "HELLS", "RFC2", "RPA2", "NASP", "RAD51AP1", "GMNN", "WDR76",
             "SLBP", "CCNE2", "UBR7", "POLD3", "MSH2", "ATAD2", "RAD51", "RRM2", "CDC45", "CDC6", "EXO1",
             "TIPIN", "DSCC1", "BLM", "CASP8AP2", "USP1", "CLSPN", "POLA1", "CHAF1B", "BRIP1", "E2F8")
G2M_GENES <- c("HMGB2", "CDK1", "NUSAP1", "UBE2C", "BIRC5", "TPX2", "TOP2A", "NDC80", "CKS2", "NUF2", "CKS1B",
               "MKI67", "TMPO", "CENPF", "TACC3", "FAM64A", "PIMREG", "SMC4", "CCNB2", "CKAP2L", "CKAP2", "AURKB",
               "BUB1", "KIF11", "ANP32E", "TUBB4B", "GTSE1", "KIF20B", "HJURP", "CDCA3", "HN1", "JPT1", "CDC20",
               "TTK", "CDC25C", "KIF2C", "RANGAP1", "NCAPD2", "DLGAP5", "CDCA2", "CDCA8", "ECT2", "KIF23",
               "HMMR", "AURKA", "PSRC1", "ANLN", "LBR", "CKAP5", "CENPE", "CTCF", "NEK2", "G2E3", "GAS2L3",
               "CBX5", "CENPA")
s_genes <- genes_present(rownames(obj), S_GENES)
g2m_genes <- genes_present(rownames(obj), G2M_GENES)
if (length(s_genes) >= 10 && length(g2m_genes) >= 10) {
  obj <- CellCycleScoring(obj, s.features = s_genes, g2m.features = g2m_genes, set.ident = FALSE)
  print(table(obj$Phase))
} else {
  message("Too few cell-cycle genes found (gene names?): cell-cycle scoring skipped.")
}
```

---

## 7. Step 5: Highly variable genes

**What it does.** Seurat's `vst` method is the same Seurat-v3 method used in the Python guide: it fits a
smooth curve (LOESS) of variance against mean on the raw counts and picks the `N_HVG` genes that vary most
beyond it. If the curve fit fails (very small or sparse data), the code retries with a wider span and, as a
last resort, uses the dispersion method on the log-normalised values. The method used is stored for the
methods text. Background: [Python guide, Section 9](Single_Cell_Master_Guide.md#9-step-5-highly-variable-genes).

```r
find_hvgs <- function(obj) {
  # Mark N_HVG highly variable genes; the method used is stored in obj@misc$hvg_method
  if (DATA_STATE == "counts") {
    for (span in c(0.3, 0.5, 1.0)) {             # 0.3 = default; a wider span rescues the fit on sparse data
      res <- tryCatch(FindVariableFeatures(obj, selection.method = "vst", nfeatures = N_HVG,
                                           loess.span = span, verbose = FALSE),
                      error = function(e) { message("  vst fit failed with span ", span, ": ", conditionMessage(e)); NULL })
      if (!is.null(res)) {
        res@misc$hvg_method <- paste0("Seurat vst method on raw counts", if (span != 0.3) paste0(", LOESS span ", span) else "")
        return(res)
      }
    }
  }
  obj <- FindVariableFeatures(obj, selection.method = "dispersion", nfeatures = N_HVG, verbose = FALSE)
  obj@misc$hvg_method <- "dispersion method on log-normalised values"
  obj
}

obj <- find_hvgs(obj)
message(length(VariableFeatures(obj)), " highly variable genes; ", obj@misc$hvg_method)
p <- LabelPoints(VariableFeaturePlot(obj), points = head(VariableFeatures(obj), 10), repel = FALSE)
save_plot(p, "05_highly_variable_genes")
```

---

## 8. Step 6: Scaling and PCA

**What it does.** Following Seurat's standard workflow, the HVGs are first scaled (each gene to mean 0 and
variance 1, so every gene has equal weight), then PCA turns them into 50 components. This is the one
deliberate difference from the Python guide, which follows the current Scanpy tutorials and does not scale;
both are accepted practice. Choose `N_PCS` from the elbow plot. Background:
[Python guide, Section 10](Single_Cell_Master_Guide.md#10-step-6-pca).

```r
obj <- ScaleData(obj, features = VariableFeatures(obj), verbose = FALSE)
n_comps <- min(50, ncol(obj) - 1, length(VariableFeatures(obj)) - 1)
obj <- RunPCA(obj, npcs = n_comps, seed.use = SEED, verbose = FALSE)
save_plot(ElbowPlot(obj, ndims = n_comps), "06_pca_elbow")
save_plot(DimPlot(obj, reduction = "pca", group.by = "sample"), "06_pca_by_sample")
n_pcs <- min(N_PCS, n_comps)
```

---

## 9. Step 7: Batch correction with Harmony

**When to use it, and the danger when batch = condition:**
[Python guide, Section 11](Single_Cell_Master_Guide.md#11-step-7-batch-correction-integration). The R
`harmony` package works directly on the Seurat object and writes a new reduction called `harmony`; a UMAP of
the uncorrected data is made first so you can compare.

```r
n_batches <- if (is.null(BATCH_KEY)) 1 else length(unique(obj@meta.data[[BATCH_KEY]]))
REP <- "pca"
if (n_batches > 1) {
  # map BEFORE correction, for comparison. Its name must not contain "umap": DimPlot() would then draw it
  # instead of the corrected map, because it picks the first reduction whose name contains "umap".
  obj <- RunUMAP(obj, reduction = "pca", dims = 1:n_pcs, n.neighbors = N_NEIGHBORS, seed.use = SEED,
                 reduction.name = "before_correction", reduction.key = "BEFORE_", verbose = FALSE)
  obj <- harmony::RunHarmony(obj, group.by.vars = BATCH_KEY, reduction.use = "pca", dims.use = 1:n_pcs,
                             verbose = FALSE)
  REP <- "harmony"
}
message("Reduction used from now on: ", REP)
```

---

## 10. Step 8: Neighbours, UMAP and clustering

**What it does.** Builds the neighbour graph on the corrected components, draws the UMAP, and clusters with
**Leiden** at several resolutions (`algorithm = 4`; if Leiden is not available, it falls back to Louvain and
says so). Seurat clusters on a *shared*-nearest-neighbour (SNN) graph, Scanpy on a kNN graph, so cluster
numbers will differ between the guides even when the groups agree. How to read a UMAP honestly and how to
choose the resolution: [Python guide, Section 12](Single_Cell_Master_Guide.md#12-step-8-neighbours-umap-and-clustering).

```r
obj <- FindNeighbors(obj, reduction = REP, dims = 1:n_pcs, k.param = N_NEIGHBORS, verbose = FALSE)
obj <- RunUMAP(obj, reduction = REP, dims = 1:n_pcs, n.neighbors = N_NEIGHBORS, seed.use = SEED, verbose = FALSE)

show <- unique(c("sample", BATCH_KEY, CONDITION_KEY))
if ("before_correction" %in% Reductions(obj)) {
  p <- wrap_plots(lapply(show, function(g) DimPlot(obj, reduction = "before_correction", group.by = g) +
                           ggtitle(paste(g, "- before correction"))), ncol = length(show))
  save_plot(p, "08_umap_before_correction", width = 6 * length(show), height = 5)
}
p <- wrap_plots(lapply(show, function(g) DimPlot(obj, group.by = g)), ncol = length(show))
save_plot(p, "08_umap_after_correction", width = 6 * length(show), height = 5)

res_all <- sort(unique(c(RESOLUTIONS, RESOLUTION)))
clustered <- tryCatch(FindClusters(obj, resolution = res_all, algorithm = 4, random.seed = SEED, verbose = FALSE),
                      error = function(e) { message("Leiden not available (", conditionMessage(e), "); using Louvain"); NULL })
if (is.null(clustered)) {
  obj <- FindClusters(obj, resolution = res_all, algorithm = 1, random.seed = SEED, verbose = FALSE)
  obj@misc$cluster_method <- "Louvain"
} else {
  obj <- clustered
  obj@misc$cluster_method <- "Leiden"
}
rm(clustered)
res_col <- function(r) paste0(DefaultAssay(obj), "_snn_res.", r)
print(sapply(res_all, function(r) length(unique(obj@meta.data[[res_col(r)]]))))
p <- wrap_plots(lapply(RESOLUTIONS, function(r) DimPlot(obj, group.by = res_col(r), label = TRUE) +
                         NoLegend() + ggtitle(paste("resolution", r))), ncol = 3)
save_plot(p, "08_umap_resolutions", width = 15, height = 5 * ceiling(length(RESOLUTIONS) / 3))

obj$cluster <- factor(obj@meta.data[[res_col(RESOLUTION)]])
Idents(obj) <- "cluster"
save_plot(DimPlot(obj, label = TRUE) + NoLegend(), "08_umap_clusters")
```

```r
qc_cols <- intersect(c("nCount_RNA", "nFeature_RNA", "pct_mt", "doublet_score"), colnames(obj@meta.data))
qc_cols <- qc_cols[vapply(qc_cols, function(c) any(!is.na(obj@meta.data[[c]])), logical(1))]
p <- FeaturePlot(obj, features = qc_cols, ncol = 2)
if ("Phase" %in% colnames(obj@meta.data)) p <- p | DimPlot(obj, group.by = "Phase")
save_plot(p, "08_umap_qc_metrics", width = 14, height = 8)

mix <- round(100 * prop.table(table(obj$cluster, obj$sample), 1), 1)
write.csv(as.data.frame.matrix(mix), file.path(TAB, "08_cluster_by_sample_percent.csv"))
one_sample <- apply(mix, 1, max) > 80
message("Clusters where one sample makes up > 80% of cells: ",
        if (any(one_sample)) paste(rownames(mix)[one_sample], collapse = ", ") else "none")
saveRDS(obj, file.path(OBJ, "08_clustered.rds"))
```

---

## 11. Step 9: Marker genes

**What it does.** For every cluster, Wilcoxon test of each gene against all other cells (`FindAllMarkers`).
Columns: `avg_log2FC`, `pct.1` (fraction of cells expressing inside the cluster), `pct.2` (outside),
`p_val_adj`. How to read markers: [Python guide, Section 13](Single_Cell_Master_Guide.md#13-step-9-marker-genes).
*Speed:* for large datasets install the `presto` package (from GitHub, `immunogenomics/presto`); Seurat then
uses it automatically and the test becomes many times faster.

```r
markers <- FindAllMarkers(obj, only.pos = TRUE, test.use = "wilcox", min.pct = 0.1, logfc.threshold = 0.25,
                          verbose = FALSE)
write.csv(markers, file.path(TAB, "09_markers_all_clusters.csv"), row.names = FALSE)

good <- markers[markers$p_val_adj < 0.05 & markers$avg_log2FC > 1 & markers$pct.1 > 0.25, ]
top <- do.call(rbind, lapply(split(good, good$cluster), head, 10))
write.csv(top, file.path(TAB, "09_markers_top10_per_cluster.csv"), row.names = FALSE)
print(sapply(split(top$gene, top$cluster), function(g) paste(head(g, 6), collapse = ", ")))

top4 <- unique(unlist(lapply(split(good$gene, good$cluster), head, 4)))
save_plot(DotPlot(obj, features = top4) + RotatedAxis(), "09_markers_dotplot",
          width = max(8, 0.3 * length(top4)), height = 5)
```

---

## 12. Step 10: Cell-type annotation

### 12.1 Known markers: which clusters carry which marker sets

**What it does.** The same rule as the Python guide: for every cluster and cell type, the share of that cell
type's marker genes that are clear markers of the cluster (log2FC > 1, adjusted p < 0.05, expressed in > 25%
of the cluster's cells). The best-supported cell type is suggested; fewer than half its markers → "Unclear";
equal support for the runner-up → both names (e.g. "Monocyte / Macrophage").
Judging markers *relative to other clusters* stops genes that are high everywhere (e.g. ambient haemoglobin
in blood) from labelling every cluster; see [Python guide, Section 14](Single_Cell_Master_Guide.md#14-step-10-cell-type-annotation).

```r
MARKERS <- list(
  # immune
  "T cell"         = c("CD3D", "CD3E", "CD2", "TRAC"),
  "NK cell"        = c("NKG7", "GNLY", "KLRD1", "KLRF1"),
  "B cell"         = c("MS4A1", "CD79A", "CD19", "BANK1"),
  "Plasma cell"    = c("JCHAIN", "MZB1", "TNFRSF17", "DERL3"),
  "Monocyte"       = c("CD14", "LYZ", "S100A8", "S100A9", "FCN1"),
  "Macrophage"     = c("CD68", "C1QA", "C1QB", "C1QC"),
  "Dendritic cell" = c("FCER1A", "CD1C", "CLEC10A", "CLEC9A"),
  "pDC"            = c("LILRA4", "CLEC4C", "IL3RA"),
  "Mast cell"      = c("TPSAB1", "TPSB2", "CPA3", "MS4A2"),
  "Neutrophil"     = c("CSF3R", "FCGR3B", "CXCR2"),
  "Platelet"       = c("PPBP", "PF4", "GP9"),
  "Erythrocyte"    = c("HBB", "HBA1", "HBA2"),
  # structural
  "Epithelial"     = c("EPCAM", "KRT8", "KRT18", "KRT19", "CDH1"),
  "Endothelial"    = c("PECAM1", "VWF", "CDH5", "CLDN5"),
  "Fibroblast"     = c("COL1A1", "COL1A2", "DCN", "LUM", "PDGFRA"),
  "Pericyte/SMC"   = c("RGS5", "PDGFRB", "ACTA2", "MYH11", "TAGLN")
)
TISSUE_MARKERS <- list(
  # Example for human breast (check the markers for YOUR tissue in recent atlases before use):
  # "Luminal HR+"         = c("ESR1", "PGR", "ANKRD30A", "AR"),
  # "Luminal secretory"   = c("KIT", "ELF5", "SLPI", "LTF"),
  # "Basal/myoepithelial" = c("KRT5", "KRT14", "KRT17", "TP63"),
  # "Adipocyte"           = c("ADIPOQ", "PLIN1", "LEP")
)
MARKERS <- c(MARKERS, TISSUE_MARKERS)

present <- lapply(MARKERS, function(g) genes_present(rownames(obj), g))
present <- present[lengths(present) > 0]
save_plot(DotPlot(obj, features = present, group.by = "cluster") + RotatedAxis() +
            theme(strip.text.x = element_text(angle = 90, size = 7)),
          "10_known_markers_dotplot", width = 18, height = 6)

support <- do.call(rbind, lapply(levels(obj$cluster), function(cl) {
  m <- markers[markers$cluster == cl, ]
  rownames(m) <- m$gene
  do.call(rbind, lapply(names(present), function(ct) {
    g <- m[present[[ct]], ]                                   # rows of NA = gene is not a marker of this cluster
    hit <- !is.na(g$gene) & g$avg_log2FC > 1 & g$p_val_adj < 0.05 & g$pct.1 > 0.25
    data.frame(cluster = cl, cell_type = ct, hits = mean(hit),
               median_log2FC = suppressWarnings(median(pmin(pmax(g$avg_log2FC, -5), 5), na.rm = TRUE)))
  }))
}))
write.csv(support, file.path(TAB, "10_marker_support_cluster_x_celltype.csv"), row.names = FALSE)

lfc_rank <- ifelse(is.na(support$median_log2FC), -Inf, support$median_log2FC)
support <- support[order(support$cluster, -support$hits, -lfc_rank), ]
best <- support[!duplicated(support$cluster), ]
rest <- support[duplicated(support$cluster), ]
second <- rest[!duplicated(rest$cluster), ]
second_hits <- second$hits[match(best$cluster, second$cluster)]
second_lfc <- second$median_log2FC[match(best$cluster, second$cluster)]
tie <- !is.na(second_hits) & second_hits == best$hits & !is.na(second_lfc) & abs(best$median_log2FC - second_lfc) < 0.5
label <- ifelse(tie, paste(best$cell_type, "/", second$cell_type[match(best$cluster, second$cluster)]), best$cell_type)
suggestion <- data.frame(cluster = best$cluster,
                         suggested = ifelse(best$hits >= 0.5, label, "Unclear"),   # at least half must hit; ties show both
                         marker_hits = round(best$hits, 2), median_log2FC = round(best$median_log2FC, 2),
                         runner_up = second$cell_type[match(best$cluster, second$cluster)],
                         runner_up_hits = round(second$hits[match(best$cluster, second$cluster)], 2))
suggestion <- suggestion[match(levels(obj$cluster), suggestion$cluster), ]
suggestion$n_cells <- as.vector(table(obj$cluster)[suggestion$cluster])
write.csv(suggestion, file.path(TAB, "10_annotation_suggestion.csv"), row.names = FALSE)
print(suggestion, row.names = FALSE)
```

### 12.2 Your decision

Read the suggestion with the top markers and the dot plot; overwrite or refine names in `CLUSTER_TO_CELLTYPE`.

```r
CLUSTER_TO_CELLTYPE <- c(
  # "0" = "CD4 T cell",
  # "5" = "Doublets (remove)"
)
labels <- setNames(suggestion$suggested, suggestion$cluster)
if (length(CLUSTER_TO_CELLTYPE)) labels[names(CLUSTER_TO_CELLTYPE)] <- CLUSTER_TO_CELLTYPE
obj$cell_type <- factor(unname(labels[as.character(obj$cluster)]))

save_plot(DimPlot(obj, group.by = "cell_type", label = TRUE, repel = TRUE) + NoLegend(), "10_umap_cell_types")
save_plot(DotPlot(obj, features = present, group.by = "cell_type") + RotatedAxis() +
            theme(strip.text.x = element_text(angle = 90, size = 7)),
          "10_known_markers_by_cell_type", width = 18, height = 5)
print(table(obj$cell_type))
```

### 12.3 Automatic annotation with SingleR (optional, recommended as a second opinion)

**What it does.** SingleR (Aran et al. 2019) compares each cluster with labelled reference expression data
from **celldex** and gives the best-matching label: the R counterpart of CellTypist. References:
`MonacoImmuneData()` (blood immune cells), `HumanPrimaryCellAtlasData()` (many human cell types),
`MouseRNAseqData()` and `ImmGenData()` (mouse). The reference is downloaded once.

```r
# OPTIONAL: SingleR labels per cluster (downloads the celldex reference the first time)
ref <- celldex::MonacoImmuneData()
pred <- SingleR::SingleR(test = LayerData(obj, layer = "data"), ref = ref, labels = ref$label.main,
                         clusters = obj$cluster)
singler_label <- setNames(pred$pruned.labels, rownames(pred))
obj$singler <- factor(unname(singler_label[as.character(obj$cluster)]))
write.csv(data.frame(cluster = rownames(pred), ours = labels[rownames(pred)], singler = pred$pruned.labels),
          file.path(TAB, "10_cluster_vs_singler.csv"), row.names = FALSE)
save_plot(DimPlot(obj, group.by = "cell_type", label = TRUE, repel = TRUE) + NoLegend() |
            DimPlot(obj, group.by = "singler", label = TRUE, repel = TRUE) + NoLegend(),
          "10_umap_singler", width = 14, height = 6)
```

---

## 13. Step 11: Cell-type composition

Same table, plot and exploratory test as the Python guide (Mann-Whitney / Wilcoxon on per-sample %, BH
corrected), with the same warning about compositional data
([Section 15](Single_Cell_Master_Guide.md#15-step-11-cell-type-composition)). In R, **propeller** (package
`speckle`) is a good compositional test for the paper.

```r
tab <- table(obj@meta.data[[REPLICATE_KEY]], obj$cell_type)
props <- 100 * prop.table(tab, 1)
write.csv(as.data.frame.matrix(round(props, 2)), file.path(TAB, "11_cell_type_percent_per_sample.csv"))
long <- as.data.frame(props)
names(long) <- c("replicate", "cell_type", "percent")
p <- ggplot(long, aes(replicate, percent, fill = cell_type)) + geom_col(width = 0.85) + theme_classic() +
  labs(x = NULL, y = "% of cells", fill = NULL) + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
save_plot(p, "11_composition_per_sample", width = max(6, 0.45 * nrow(props) + 3))

if (!is.null(CONDITION_KEY)) {
  pairs <- unique(obj@meta.data[, c(REPLICATE_KEY, CONDITION_KEY)])
  if (any(duplicated(pairs[[REPLICATE_KEY]])))
    stop("Some ", REPLICATE_KEY, " values belong to more than one ", CONDITION_KEY, ": fix the sample sheet.")
  group_of <- setNames(as.character(pairs[[CONDITION_KEY]]), as.character(pairs[[REPLICATE_KEY]]))
  lv <- sort(unique(group_of))
  test_grp <- if (is.null(CONTRAST)) lv[2] else CONTRAST[1]
  ref_grp  <- if (is.null(CONTRAST)) lv[1] else CONTRAST[2]
  message("Comparison: ", test_grp, " (test) vs ", ref_grp, " (reference)")
  a_s <- names(group_of)[group_of == test_grp]
  b_s <- names(group_of)[group_of == ref_grp]
  comp <- do.call(rbind, lapply(colnames(props), function(ct) data.frame(
    cell_type = ct, mean_test = mean(props[a_s, ct]), mean_ref = mean(props[b_s, ct]),
    p_value = suppressWarnings(wilcox.test(props[a_s, ct], props[b_s, ct])$p.value))))
  comp$p_adj <- p.adjust(comp$p_value, method = "BH")
  names(comp)[2:3] <- paste0("mean_pct_", c(test_grp, ref_grp))
  write.csv(comp, file.path(TAB, "11_composition_test.csv"), row.names = FALSE)
  print(comp, digits = 3, row.names = FALSE)

  long$group <- group_of[as.character(long$replicate)]
  p <- ggplot(long, aes(cell_type, percent, fill = group)) + geom_boxplot(outlier.shape = NA) +
    geom_point(position = position_dodge(width = 0.75), size = 1) + theme_classic() +
    labs(x = NULL, y = "% of cells per sample") + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
  save_plot(p, "11_composition_by_group", width = max(6, 0.6 * ncol(props) + 2))
}
```

---

## 14. Step 12: Differential expression (pseudobulk + DESeq2)

**Why pseudobulk, and how to read the results:**
[Python guide, Section 16](Single_Cell_Master_Guide.md#16-step-12-differential-expression-pseudobulk). Here
the counts are summed per sample and cell type with one sparse matrix product, then analysed with the
original **DESeq2** (Love et al. 2014). Same filters, same design, same thresholds as the Python guide.
`log2FoldChange > 0` = higher in the test group.

```r
if (!is.null(CONDITION_KEY) && DATA_STATE == "counts") {
  suppressPackageStartupMessages(library(DESeq2))

  pseudobulk <- function(obj, group_key, sample_key, min_cells = 10) {
    # Add up the raw counts of all cells of one group (cell type) in one sample
    key <- paste(obj@meta.data[[sample_key]], obj@meta.data[[group_key]], sep = "|")
    n <- table(key)
    keep <- names(n)[n >= min_cells]
    use <- key %in% keep
    f <- factor(key[use], levels = keep)
    design <- Matrix::sparse.model.matrix(~ 0 + f)
    colnames(design) <- levels(f)
    sums <- LayerData(obj, layer = "counts")[, use] %*% design
    meta <- data.frame(row.names = keep, sample = sub("\\|.*$", "", keep), cell_type = sub("^[^|]*\\|", "", keep),
                       n_cells = as.vector(n[keep]))
    list(counts = as.matrix(sums), meta = meta)
  }

  pb <- pseudobulk(obj, "cell_type", REPLICATE_KEY, MIN_CELLS_PSEUDOBULK)
  info <- obj@meta.data[!duplicated(obj@meta.data[[REPLICATE_KEY]]), c(REPLICATE_KEY, CONDITION_KEY, DE_COVARIATES), drop = FALSE]
  extra <- info[match(pb$meta$sample, as.character(info[[REPLICATE_KEY]])), c(CONDITION_KEY, DE_COVARIATES), drop = FALSE]
  rownames(extra) <- NULL
  pb$meta <- cbind(pb$meta, extra)
  pb$meta <- pb$meta[as.character(pb$meta[[CONDITION_KEY]]) %in% c(test_grp, ref_grp), ]
  write.csv(pb$meta, file.path(TAB, "12_pseudobulk_samples.csv"))

  design_formula <- as.formula(paste("~", paste(c(DE_COVARIATES, CONDITION_KEY), collapse = " + ")))
  de_results <- list()
  for (ct in sort(unique(pb$meta$cell_type))) {
    m <- pb$meta[pb$meta$cell_type == ct, ]
    n <- table(factor(as.character(m[[CONDITION_KEY]]), levels = c(test_grp, ref_grp)))
    if (any(n < 2)) {
      message(ct, ": skipped (needs >= 2 samples per group, has ", paste(names(n), n, sep = "=", collapse = ", "), ")")
      next
    }
    cnt <- round(pb$counts[, rownames(m), drop = FALSE])
    cnt <- cnt[rowSums(cnt >= 10) >= min(n), , drop = FALSE]
    storage.mode(cnt) <- "integer"
    for (v in c(CONDITION_KEY, DE_COVARIATES)) m[[v]] <- factor(as.character(m[[v]]))
    m[[CONDITION_KEY]] <- relevel(m[[CONDITION_KEY]], ref = ref_grp)
    dds <- DESeqDataSetFromMatrix(cnt, m, design_formula)
    dds <- DESeq(dds, quiet = TRUE)
    res <- as.data.frame(results(dds, contrast = c(CONDITION_KEY, test_grp, ref_grp)))
    res <- res[order(res$padj), ]
    write.csv(res, file.path(TAB, paste0("12_DE_", gsub("[^A-Za-z0-9]+", "_", ct), ".csv")))
    de_results[[ct]] <- res
    sig <- !is.na(res$padj) & res$padj < 0.05 & abs(res$log2FoldChange) > 1
    message(ct, ": ", nrow(m), " samples, ", nrow(cnt), " genes tested, ",
            sum(sig & res$log2FoldChange > 0), " up / ", sum(sig & res$log2FoldChange < 0), " down")
  }

  if (length(de_results)) {
    de_all <- do.call(rbind, lapply(names(de_results), function(ct)
      cbind(gene = rownames(de_results[[ct]]), cell_type = ct, de_results[[ct]])))
    write.csv(de_all, file.path(TAB, "12_DE_all_cell_types.csv"), row.names = FALSE)

    de_all$significant <- !is.na(de_all$padj) & de_all$padj < 0.05 & abs(de_all$log2FoldChange) > 1
    lab <- do.call(rbind, lapply(split(de_all[de_all$significant, ], de_all$cell_type[de_all$significant]), head, 8))
    p <- ggplot(de_all, aes(log2FoldChange, -log10(pvalue), colour = significant)) +
      geom_point(size = 0.5) + scale_colour_manual(values = c(`FALSE` = "grey80", `TRUE` = "firebrick")) +
      geom_text(data = lab, aes(label = gene), size = 2, colour = "black", check_overlap = TRUE, vjust = -0.5) +
      facet_wrap(~ cell_type, scales = "free") + theme_classic() + theme(legend.position = "none") +
      labs(x = paste0("log2 FC (", test_grp, " / ", ref_grp, ")"), y = "-log10 p")
    save_plot(p, "12_volcano_per_cell_type", width = 12, height = 4 * ceiling(length(de_results) / 3))
  }
} else {
  message("Pseudobulk DE skipped (needs CONDITION_KEY and raw counts).")
}
```

---

## 15. Step 13: Pathways and transcription factors

**What it does.** decoupleR (Badia-i-Mompel et al. 2022; the R original of Python's decoupler, same authors
and resources) scores PROGENy pathways, CollecTRI transcription factors and MSigDB Hallmarks, (a) per cell and
then averaged per cell type, and (b) on the DE statistics of Step 12. Background and how to read it:
[Python guide, Section 17](Single_Cell_Master_Guide.md#17-step-13-pathways-and-transcription-factors). The
resources are downloaded from OmniPath the first time (internet needed). When the OmniPath web service is
down, decoupleR falls back to static copies; a resource that still cannot be downloaded is skipped with a
message instead of stopping the run. (On the test day, 31 Oct 2026, the CollecTRI fallback failed in
OmnipathR 4.0.0 while Python's decoupler downloaded CollecTRI fine: if that happens to you, run Step 13 in the
Python guide or try again later.)

```r
# OPTIONAL: pathway (PROGENy), transcription-factor (CollecTRI) and Hallmark activity with decoupleR
suppressPackageStartupMessages(library(decoupleR))
organism <- if (SPECIES == "human") "human" else "mouse"
fetch <- function(label, expr) tryCatch(expr, error = function(e) {
  message("  ", label, " could not be downloaded (", conditionMessage(e), "): skipped"); NULL })
progeny <- fetch("PROGENy", get_progeny(organism = organism, top = 500))
collectri <- fetch("CollecTRI", get_collectri(organism = organism, split_complexes = FALSE))
msig <- fetch("MSigDB", get_resource("MSigDB"))
hallmark <- if (is.null(msig)) NULL else
  unique(data.frame(source = msig$geneset[msig$collection == "hallmark"],
                    target = msig$genesymbol[msig$collection == "hallmark"], mor = 1))
nets <- Filter(Negate(is.null), list(PROGENy = progeny, Hallmark = hallmark, CollecTRI = collectri))
message("Resources available: ", paste(names(nets), collapse = ", "))

heatmap_plot <- function(long, value, title) {
  ggplot(long, aes(source, condition, fill = .data[[value]])) + geom_tile() +
    scale_fill_gradient2(low = "steelblue", mid = "white", high = "firebrick") + theme_minimal() +
    labs(x = NULL, y = NULL, title = title, fill = NULL) + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
}

# (a) per cell, then averaged per cell type (only the network genes are taken, to save memory)
if (!is.null(progeny)) {
  genes <- intersect(rownames(obj), progeny$target)
  acts <- run_ulm(mat = as.matrix(LayerData(obj, layer = "data")[genes, ]), net = progeny,
                  .source = "source", .target = "target", .mor = "weight", minsize = 5)
  acts$cell_type <- obj$cell_type[match(acts$condition, colnames(obj))]
  per_type <- aggregate(score ~ source + cell_type, data = acts, FUN = mean)
  names(per_type)[2] <- "condition"
  save_plot(heatmap_plot(per_type, "score", "PROGENy pathway activity (mean per cell type)"),
            "13_pathways_per_cell_type", width = 8, height = 5)
}

# (b) on the pseudobulk DE statistics: what changes between the conditions, per cell type
if (!is.null(CONDITION_KEY) && exists("de_all")) {
  stat <- tapply(de_all$stat, list(de_all$gene, de_all$cell_type), function(x) x[1])
  stat[is.na(stat)] <- 0
  for (res_name in names(nets)) {
    net <- nets[[res_name]]
    mor <- if (res_name == "PROGENy") "weight" else "mor"
    change <- run_ulm(mat = stat, net = net, .source = "source", .target = "target", .mor = mor, minsize = 5)
    write.csv(change, file.path(TAB, paste0("13_", res_name, "_change_", test_grp, "_vs_", ref_grp, ".csv")), row.names = FALSE)
    top_sets <- names(sort(tapply(abs(change$score), change$source, max), decreasing = TRUE))[1:min(25, length(unique(change$source)))]
    save_plot(heatmap_plot(change[change$source %in% top_sets, ], "score", paste(res_name, "change:", test_grp, "vs", ref_grp)),
              paste0("13_", res_name, "_change_heatmap"), width = 10, height = 4)
  }
}
```

---

## 16. Step 14: Save, export and write the methods

**What it does.** Saves the final Seurat object (`readRDS()` to open it again), CSV tables, every setting,
the package versions (`sessionInfo.txt`), and a draft methods paragraph with your numbers filled in.

```r
saveRDS(obj, file.path(OBJ, "final.rds"))
write.csv(obj@meta.data, file.path(TAB, "14_cell_metadata.csv"))
write.csv(Embeddings(obj, "umap"), file.path(TAB, "14_umap_coordinates.csv"))

settings <- list(SAMPLES = as.list(SAMPLES), SAMPLE_COL = SAMPLE_COL, SAMPLE_SHEET = SAMPLE_SHEET, SPECIES = SPECIES,
                 MATERIAL = MATERIAL, MIN_GENES = MIN_GENES, NMADS = NMADS, MT_NMADS = MT_NMADS, MAX_PCT_MT = MAX_PCT_MT,
                 MIN_CELLS_PER_GENE = MIN_CELLS_PER_GENE, BATCH_KEY = BATCH_KEY, N_HVG = N_HVG, N_PCS = N_PCS,
                 N_NEIGHBORS = N_NEIGHBORS, RESOLUTIONS = RESOLUTIONS, RESOLUTION = RESOLUTION,
                 CONDITION_KEY = CONDITION_KEY, CONTRAST = CONTRAST, REPLICATE_KEY = REPLICATE_KEY,
                 DE_COVARIATES = DE_COVARIATES, MIN_CELLS_PSEUDOBULK = MIN_CELLS_PSEUDOBULK, SEED = SEED,
                 DATA_STATE = DATA_STATE, reduction = REP)
jsonlite::write_json(settings, file.path(TAB, "run_settings.json"), auto_unbox = TRUE, pretty = TRUE, null = "null")
writeLines(capture.output(sessionInfo()), file.path(TAB, "sessionInfo.txt"))

pkgs <- c("Seurat", "SeuratObject", "harmony", "scDblFinder", "DESeq2", "SingleR", "celldex", "decoupleR",
          "zellkonverter", "Matrix", "ggplot2")
versions <- data.frame(package = c("R", pkgs), version = c(paste(R.version$major, R.version$minor, sep = "."),
  vapply(pkgs, function(p) if (requireNamespace(p, quietly = TRUE)) as.character(packageVersion(p)) else "not installed", "")))
write.csv(versions, file.path(TAB, "software_versions.csv"), row.names = FALSE)

counts_ok <- DATA_STATE == "counts"
methods <- c(
  sprintf("Single-cell data from %d samples were analysed in R %s.%s with Seurat %s (Hao et al. 2024).",
          length(unique(obj$sample)), R.version$major, R.version$minor, packageVersion("Seurat")),
  if (counts_ok) sprintf(paste("Cells more than %s median absolute deviations (MADs) from their sample's median in log total counts,",
    "log genes detected or the fraction of counts in the top 20 genes, cells with fewer than %s genes, and cells with a",
    "mitochondrial fraction more than %s MADs above the sample median or above %s%% were removed. Genes detected in fewer",
    "than %s cells were removed."), NMADS, MIN_GENES, MT_NMADS, MAX_PCT_MT, MIN_CELLS_PER_GENE)
  else sprintf("Only processed (non-count) values were available; cells with fewer than %s genes were removed.", MIN_GENES),
  if (counts_ok) paste0("Doublets were identified with scDblFinder (Germain et al. 2021) run separately for each sample and removed",
    if (length(skipped_doublets)) paste0(" (not screened, fewer than ", MIN_CELLS_DOUBLETS, " cells: ",
                                          paste(skipped_doublets, collapse = ", "), ").") else "."),
  sprintf("In total %s of %s cells passed quality control.", format(ncol(obj), big.mark = ","),
          format(sum(qc_table$cells_loaded), big.mark = ",")),
  if (counts_ok) "Counts were normalised to 10,000 per cell and log-transformed (LogNormalize).",
  sprintf("%s highly variable genes (%s) were scaled and used for principal component analysis.",
          format(length(VariableFeatures(obj)), big.mark = ","), obj@misc$hvg_method),
  if (REP == "harmony") sprintf("Batch effects between %s were corrected with Harmony (Korsunsky et al. 2019).", BATCH_KEY),
  sprintf(paste("A shared-nearest-neighbour graph (k = %s) on the first %s components was used for UMAP and %s clustering",
    "(resolution %s, %d clusters). Clusters were annotated from their Wilcoxon marker genes and canonical marker sets",
    "(share of each cell type's markers with log2 fold change > 1, adjusted p < 0.05 and detection in > 25%% of the",
    "cluster's cells), checked by eye on dot plots."), N_NEIGHBORS, n_pcs, obj@misc$cluster_method, RESOLUTION,
    length(levels(obj$cluster))),
  if (!is.null(CONDITION_KEY) && counts_ok) sprintf(paste("Differential expression between %s and %s was tested per cell type",
    "on pseudobulk counts (sum per %s and cell type, at least %s cells) with DESeq2 (Love et al. 2014; design %s);",
    "genes with adjusted p < 0.05 and |log2 fold change| > 1 were called differentially expressed."),
    test_grp, ref_grp, REPLICATE_KEY, MIN_CELLS_PSEUDOBULK, paste(deparse(design_formula), collapse = ""))
)
writeLines(paste(methods, collapse = " "), file.path(TAB, "methods_draft.txt"))
cat(readLines(file.path(TAB, "methods_draft.txt")), sep = "\n")
```

---

## 17. Advanced modules

All optional; run them after Step 10. When and why to use each:
[Python guide, Section 19](Single_Cell_Master_Guide.md#19-advanced-modules).

### 17.1 Sub-clustering one lineage

```r
# OPTIONAL: sub-cluster one cell type
SUBSET_TYPE <- "T cell"
SUB_RESOLUTION <- 0.5

n_sub <- sum(obj$cell_type == SUBSET_TYPE)
if (n_sub < 50) stop("Only ", n_sub, " cells labelled '", SUBSET_TYPE, "'. Labels are: ",
                     paste(levels(obj$cell_type), collapse = ", "))
sub <- subset(obj, cells = colnames(obj)[obj$cell_type == SUBSET_TYPE])
sub <- find_hvgs(sub)
message("Sub-cluster HVGs: ", sub@misc$hvg_method)
sub <- ScaleData(sub, features = VariableFeatures(sub), verbose = FALSE)
sub <- RunPCA(sub, npcs = min(50, ncol(sub) - 1), seed.use = SEED, verbose = FALSE)
sub_rep <- "pca"
if (!is.null(BATCH_KEY) && length(unique(sub@meta.data[[BATCH_KEY]])) > 1 && min(table(sub@meta.data[[BATCH_KEY]])) >= 30) {
  sub <- harmony::RunHarmony(sub, group.by.vars = BATCH_KEY, reduction.use = "pca", verbose = FALSE)
  sub_rep <- "harmony"
}
sub_dims <- 1:min(n_pcs, ncol(Embeddings(sub, "pca")))
sub <- FindNeighbors(sub, reduction = sub_rep, dims = sub_dims, verbose = FALSE)
sub <- RunUMAP(sub, reduction = sub_rep, dims = sub_dims, seed.use = SEED, verbose = FALSE)
sub <- FindClusters(sub, resolution = SUB_RESOLUTION, algorithm = if (obj@misc$cluster_method == "Leiden") 4 else 1,
                    random.seed = SEED, verbose = FALSE)
sub$subcluster <- Idents(sub)
sub_markers <- FindAllMarkers(sub, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.25, verbose = FALSE)
tag <- gsub("[^A-Za-z0-9]+", "_", SUBSET_TYPE)
write.csv(sub_markers, file.path(TAB, paste0("17_subcluster_markers_", tag, ".csv")), row.names = FALSE)
save_plot(DimPlot(sub, group.by = "subcluster", label = TRUE) | DimPlot(sub, group.by = "sample"),
          paste0("17_subclusters_", tag), width = 13, height = 5)

obj$cell_subtype <- as.character(obj$cell_type)
obj@meta.data[colnames(sub), "cell_subtype"] <- paste(SUBSET_TYPE, sub$subcluster)
obj$cell_subtype <- factor(obj$cell_subtype)
```

### 17.2 Trajectories and pseudotime with slingshot

slingshot (Street et al. 2018) fits smooth curves through the clusters in reduced space, starting from a
cluster **you** choose from biology.

```r
# OPTIONAL: slingshot lineages and pseudotime
ROOT_CLUSTER <- levels(obj$cluster)[1]     # set to the cluster where the process STARTS

emb <- Embeddings(obj, REP)[, 1:min(10, n_pcs)]
sds <- slingshot::slingshot(emb, clusterLabels = obj$cluster, start.clus = ROOT_CLUSTER)
pt <- slingshot::slingPseudotime(sds)
obj$pseudotime <- rowMeans(pt, na.rm = TRUE)              # average over lineages
message(ncol(pt), " lineage(s) found")
save_plot(FeaturePlot(obj, "pseudotime") | DimPlot(obj, group.by = "cluster", label = TRUE) + NoLegend(),
          "17_pseudotime", width = 13, height = 5)
```

### 17.3 CITE-seq proteins (ADT)

If the data had antibody features, they are already in the `ADT` assay. CLR normalisation across cells
(`margin = 2`) is Seurat's standard for proteins:

```r
# OPTIONAL: CITE-seq protein (ADT) normalisation and plots
obj <- NormalizeData(obj, assay = "ADT", normalization.method = "CLR", margin = 2, verbose = FALSE)
proteins <- head(rownames(obj[["ADT"]]), 6)
DefaultAssay(obj) <- "ADT"
save_plot(FeaturePlot(obj, features = proteins, ncol = 3), "17_cite_seq_proteins", width = 12, height = 7)
DefaultAssay(obj) <- "RNA"
```

For joint RNA + protein clustering use Seurat's weighted-nearest-neighbour method (`FindMultiModalNeighbors`,
Hao et al. 2021).

### 17.4 Cell-cell communication (LIANA or CellChat)

The R versions of LIANA and CellChat are installed from GitHub, not CRAN/Bioconductor. Either run the Python
LIANA block ([Python guide, Section 19.4](Single_Cell_Master_Guide.md#194-cell-cell-communication-liana)) on
the exported object, or in R. *Untested here (needs `remotes::install_github("saezlab/liana")`):*

```r
# OPTIONAL (untested here): LIANA consensus ranking in R
library(liana)
lr <- liana_wrap(obj, idents_col = "cell_type")
lr <- liana_aggregate(lr)
write.csv(lr, file.path(TAB, "17_ligand_receptor_liana.csv"), row.names = FALSE)
```

### 17.5 Tumour copy-number inference (inferCNV)

The original inferCNV (Bioconductor `infercnv`) needs JAGS installed on the computer and a gene-position
file. *Untested here:*

```r
# OPTIONAL (untested here): inferCNV with immune/endothelial cells as the normal reference
library(infercnv)
annot <- data.frame(row.names = colnames(obj), cell_type = obj$cell_type)
icnv <- CreateInfercnvObject(raw_counts_matrix = LayerData(obj, layer = "counts"), annotations_file = annot,
                             gene_order_file = "gene_positions_GRCh38.txt", ref_group_names = c("T cell", "Endothelial"))
icnv <- infercnv::run(icnv, cutoff = 0.1, out_dir = file.path(PROJECT_DIR, "infercnv"), cluster_by_groups = TRUE,
                      denoise = TRUE, HMM = FALSE)
```

---

## 18. Other data types in R

What changes for each data type is explained in the
[Python guide, Section 20](Single_Cell_Master_Guide.md#20-other-data-types). The main R tools:

| Data type | R tool |
|---|---|
| snRNA-seq | This pipeline with `MATERIAL <- "nuclei"`; SoupX or CellBender for ambient RNA |
| CITE-seq | Seurat `ADT` assay (Section 17.3), weighted-nearest-neighbour analysis; dsb package for normalisation |
| scATAC / Multiome | **Signac** (Stuart et al. 2021) or **ArchR** (Granja et al. 2021) |
| Spatial (Visium, Xenium, CosMx) | Seurat `Load10X_Spatial()`, `LoadXenium()`; deconvolution with **spacexr (RCTD)** |
| Trajectories | slingshot (Section 17.2), **monocle3** (Cao et al. 2019) |
| TCR / BCR | **scRepertoire** (Borcherding et al. 2020) |
| Cell-cell communication | **CellChat** (Jin et al. 2021), LIANA |
| Empty droplets (raw matrix) | **DropletUtils** `emptyDrops()` (Lun et al. 2019) |

---

## 19. Moving between R and Python

* **Python `.h5ad` → R:** `read_any("objects/final.h5ad")` (this guide's reader), or directly
  `zellkonverter::readH5AD("final.h5ad", reader = "R")` → SingleCellExperiment → `as.Seurat(sce, counts = "counts", data = "X")`.
* **R → Python:** the Python reader cannot open `.rds`. Export with `zellkonverter::writeH5AD()` (it uses a
  small private Python environment that it downloads the first time), or write the counts as a 10x-style
  folder, which every tool reads:

```r
# OPTIONAL (untested here): write a Seurat object's counts + metadata as files the Python guide can read
dir.create(file.path(PROJECT_DIR, "export_10x"), showWarnings = FALSE)
Matrix::writeMM(LayerData(obj, layer = "counts"), file.path(PROJECT_DIR, "export_10x", "matrix.mtx"))
writeLines(rownames(obj), file.path(PROJECT_DIR, "export_10x", "features.tsv"))
writeLines(colnames(obj), file.path(PROJECT_DIR, "export_10x", "barcodes.tsv"))
write.csv(obj@meta.data, file.path(PROJECT_DIR, "export_10x", "cell_metadata.csv"))
```

---

## 20. Troubleshooting (R)

The general problems (too many cells removed, one label everywhere, batch problems, memory) are the same as
in the [Python guide, Section 22](Single_Cell_Master_Guide.md#22-troubleshooting). R-specific ones:

| Problem | Likely cause and fix |
|---|---|
| `there is no package called ...` | Install it (Section 1), or activate the conda environment first: `conda activate scanalysis-r`. |
| `Read10X_h5` error about `hdf5r` | `install.packages("hdf5r")` (conda: it is in `environment_r.yml`). |
| `fread` cannot read a `.gz` file | `install.packages("R.utils")`. |
| "Feature names cannot have underscores ('_'), replacing with dashes" | Harmless Seurat warning. |
| `FindAllMarkers` is very slow | Install `presto` from GitHub (`immunogenomics/presto`); Seurat uses it automatically. |
| "Leiden not available ... using Louvain" | Install `leidenbase`, or accept Louvain (report which one you used; the methods text does). |
| `future.globals.maxSize` error | Raise the limit in Step 1: `options(future.globals.maxSize = 16 * 1024^3)`. |
| An `.rds` from a collaborator fails in `JoinLayers` | It is an old Seurat v3/v4 object: `obj <- UpdateSeuratObject(obj)` first. |
| `.h5ad` cannot be read | Very new `.h5ad` versions may not yet be supported by the native R reader: open it in Python and re-save, or use `readH5AD(..., reader = "python")`. |
| Memory error on a big dataset | R needs more RAM than Python for the same data; use the Python guide for > 200,000 cells, or a bigger machine. |

---

## 21. References (R tools)

All references of the Python guide ([Section 24](Single_Cell_Master_Guide.md#24-references)) apply. Additional,
R-specific ones:

* Satija R et al. (2015). Spatial reconstruction of single-cell gene expression data. *Nature Biotechnology* 33:495-502. doi:10.1038/nbt.3192
* Hao Y et al. (2021). Integrated analysis of multimodal single-cell data. *Cell* 184:3573-3587. doi:10.1016/j.cell.2021.04.048
* Hao Y et al. (2024). Dictionary learning for integrative, multimodal and scalable single-cell analysis. *Nature Biotechnology* 42:293-304. doi:10.1038/s41587-023-01767-y
* Amezquita RA et al. (2020). Orchestrating single-cell analysis with Bioconductor. *Nature Methods* 17:137-145. doi:10.1038/s41592-019-0654-x
* Aran D et al. (2019). Reference-based analysis of lung single-cell sequencing reveals a transitional profibrotic macrophage. *Nature Immunology* 20:163-172. doi:10.1038/s41590-018-0276-y
* Street K et al. (2018). Slingshot: cell lineage and pseudotime inference for single-cell transcriptomics. *BMC Genomics* 19:477. doi:10.1186/s12864-018-4772-0
* Lun ATL et al. (2019). EmptyDrops: distinguishing cells from empty droplets in droplet-based single-cell RNA sequencing data. *Genome Biology* 20:63. doi:10.1186/s13059-019-1662-y
* Stuart T et al. (2021). Single-cell chromatin state analysis with Signac. *Nature Methods* 18:1333-1341. doi:10.1038/s41592-021-01282-5
* Granja JM et al. (2021). ArchR is a scalable software package for integrative single-cell chromatin accessibility analysis. *Nature Genetics* 53:403-411. doi:10.1038/s41588-021-00790-6
* Cao J et al. (2019). The single-cell transcriptional landscape of mammalian organogenesis. *Nature* 566:496-502. doi:10.1038/s41586-019-0969-x
* Borcherding N, Bormann NL, Kraus G (2020). scRepertoire: an R-based toolkit for single-cell immune receptor analysis. *F1000Research* 9:47. doi:10.12688/f1000research.22139.2
* Jin S et al. (2021). Inference and analysis of cell-cell communication using CellChat. *Nature Communications* 12:1088. doi:10.1038/s41467-021-21246-9

Tools also cited in the Python guide: harmony (Korsunsky 2019), scDblFinder (Germain 2021), SoupX (Young &
Behjati 2020), DESeq2 (Love 2014), decoupleR (Badia-i-Mompel 2022), PROGENy, CollecTRI, MSigDB, scDblFinder
benchmark (Xi & Li 2021), pseudobulk (Squair 2021).

---

## 22. How this guide was tested

**How.** `python tests/test_guide.py r` takes every `r` block straight from this file, keeps the Step 0 block
followed by test settings, and runs the result with `Rscript` from the `scanalysis-r` environment. Blocks
marked *untested here* and `bash` blocks are not run. The data are the same public datasets as for the Python
guide ([Python guide, Section 25](Single_Cell_Master_Guide.md#25-how-this-guide-was-tested)).

**Where.** 31 October 2026, MacBook (Apple Silicon, macOS Sonoma), R 4.5.3, Bioconductor 3.22, the versions
listed in `environment_r.yml`. All five R tests: about 7 minutes.

| Test | What it covers | Result |
|---|---|---|
| `r_reader` | 13 cases: STARsolo folder, GEO prefixed `mtx.gz` + old `genes.tsv`, matrix stored cells x genes, Cell Ranger `outs/`, cells x genes TSV with a text column, R `write.table` TXT, `.h5ad` with Ensembl IDs + counts in a layer, `.h5ad` with counts in `.raw`, `.loom`, Seurat `.rds`, SingleCellExperiment `.rds`, refusal of `.zarr`, refusal of scaled data | all passed |
| `r_covid` | 7 blood samples in one `.h5ad`, a donor with 92 cells, Harmony, DESeq2, and the OPTIONAL blocks SingleR, decoupleR, sub-clustering, slingshot | passed |
| `r_formats` | 4 samples in 4 formats, sample sheet, DESeq2 with a covariate (`~run + condition`) | passed |
| `r_cite` | One CITE-seq sample: `ADT` assay with 17 proteins, CLR normalisation block | passed |
| `r_lognorm_nuclei` | Only log-normalised values, `MATERIAL <- "nuclei"`, no batch key | passed |

**What the testing changed.**

1. zellkonverter's native R reader **ignores the `.raw` part** of `.h5ad` files → a small HDF5 reader for `raw/X` was added.
2. Checking for `raw/X` in a file without any `raw` group raised an HDF5 error, so scaled data were "refused" for the wrong reason → fixed, and the test now demands the right message.
3. `SingleCellExperiment` and Seurat both define `Assays()`; loading SingleCellExperiment second broke the merge step → it is now loaded first.
4. SingleR 2.12 needs **scrapper**, which it only "suggests" → added to the install instructions.
5. **Seurat's `DimPlot()` drew the uncorrected map** in every plot after Step 7, because it picks the first reduction whose name contains "umap" → the before-correction map is now called `before_correction`, and the tests check that plots use the corrected map.
6. OmnipathR 4.0.0 could not download **CollecTRI** on the test day (its backup path failed) → each resource is now downloaded separately and skipped with a message if it fails; PROGENy and Hallmarks ran.
7. bioconda does not build SingleR, scrapper and slingshot for Apple Silicon → `environment_r.yml` plus a one-line source install.

**Do R and Python agree?** On the blood dataset, 3,169 cells passed QC in both guides. Where both gave a
label (not "Unclear"), they agreed on **1,289 of 1,365 cells (94%)**: platelets 201/201, monocytes, NK cells and
neutrophils closely. The differences are in the shallow lymphocytes: Python called a T-cell cluster that R left
"Unclear", and R found a B-cell cluster inside Python's "Unclear" group. SingleR (independent, reference-based)
called both of R's "Unclear" clusters CD4+ T cells, as CellTypist did in Python, which confirms these are T cells
too shallow for the cautious marker rule. SingleR called platelets "Progenitors" because its blood reference
contains no platelets: a reference can only give labels it contains. On the 10x PBMC data, R found monocytes,
T, B and NK cells, dendritic cells, platelets and a small JCHAIN/MZB1-high cluster labelled "Plasma cell" with
pDC as an equal runner-up (`runner_up_hits = 1.0`): a case for your judgement in Step 10.2.

**Not tested here.** Blocks marked *untested here* (Option B install, SoupX, LIANA, inferCNV, export to
10x files) and the Windows / CRAN-R installation route; datasets larger than ~5,000 cells.
