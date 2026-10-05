## ============================================================================
## Darmanis adult human brain scRNA-seq preprocessing
## Two-class benchmark: Neuron vs Astrocyte
## ============================================================================
##
## Dataset:
##   scRNAseq::DarmanisBrainData()
##
## Goal:
##   Construct a compact sparse-omics benchmark using two discrete adult
##   human brain cell classes: neurons and astrocytes.
##
## Important metadata:
##   age    = developmental stage / donor age
##   tissue = brain region (e.g., cortex, hippocampus)
##
## Adult definition used here:
##   all samples with age labels beginning with "postnatal"
##
## Main preprocessing:
##   1. Load DarmanisBrainData().
##   2. Keep postnatal cells only.
##   3. Keep neurons and astrocytes only.
##   4. Remove zero-library cells.
##   5. Remove genes detected in fewer than 5 retained cells.
##   6. Log-normalize.
##   7. Select top 200 HVGs.
##
## ============================================================================


## ---------------------------------------------------------------------------
## 0. Install/load required packages
## ---------------------------------------------------------------------------

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

bioc_pkgs <- c(
  "scRNAseq",
  "SingleCellExperiment",
  "scater",
  "scran"
)

for (pkg in bioc_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    BiocManager::install(
      pkg,
      ask = FALSE,
      update = FALSE
    )
  }
}

if (!requireNamespace("radviz3d", quietly = TRUE)) {
  install.packages("radviz3d")
}

library(scRNAseq)
library(SingleCellExperiment)
library(scater)
library(scran)


## ---------------------------------------------------------------------------
## 1. Load Darmanis human brain data
## ---------------------------------------------------------------------------

sce <- scRNAseq::DarmanisBrainData()

cat("Original dimensions (genes x cells):\n")
print(dim(sce))

cat("\nAvailable assays:\n")
print(assayNames(sce))

cat("\nAvailable metadata fields:\n")
print(colnames(colData(sce)))


## ---------------------------------------------------------------------------
## 2. Inspect age, tissue, and cell-type annotations
## ---------------------------------------------------------------------------

if (!"age" %in% colnames(colData(sce))) {
  stop("Metadata field 'age' was not found.")
}

if (!"cell.type" %in% colnames(colData(sce))) {
  stop("Metadata field 'cell.type' was not found.")
}

if (!"tissue" %in% colnames(colData(sce))) {
  stop("Metadata field 'tissue' was not found.")
}

age <- as.character(colData(sce)$age)
cell_type <- as.character(colData(sce)$cell.type)
tissue <- as.character(colData(sce)$tissue)

cat("\nAge counts:\n")
print(
  sort(
    table(age, useNA = "ifany"),
    decreasing = TRUE
  )
)

cat("\nTissue counts:\n")
print(
  sort(
    table(tissue, useNA = "ifany"),
    decreasing = TRUE
  )
)

cat("\nCell-type counts:\n")
print(
  sort(
    table(cell_type, useNA = "ifany"),
    decreasing = TRUE
  )
)


## ---------------------------------------------------------------------------
## 3. Restrict to adult/postnatal cells
## ---------------------------------------------------------------------------
##
## In this dataset:
##   tissue = cortex / hippocampus
## and should NOT be used to define adulthood.
##
## Age labels include:
##   postnatal 21 years
##   postnatal 22 years
##   ...
##   postnatal 63 years
##   prenatal 16-18 W
##
## We therefore define the adult benchmark as all cells with age beginning
## with "postnatal".

adult_keep <-
  !is.na(age) &
  grepl(
    "^postnatal",
    age,
    ignore.case = TRUE
  )

sce_adult <- sce[, adult_keep]

cat(
  "\nNumber of postnatal/adult cells:",
  ncol(sce_adult),
  "\n"
)

cat("\nAdult/postnatal age counts:\n")
print(
  sort(
    table(
      colData(sce_adult)$age,
      useNA = "ifany"
    ),
    decreasing = TRUE
  )
)

cat("\nAdult/postnatal cell-type counts:\n")
print(
  sort(
    table(
      colData(sce_adult)$cell.type,
      useNA = "ifany"
    ),
    decreasing = TRUE
  )
)


## ---------------------------------------------------------------------------
## 4. Identify neuron and astrocyte labels
## ---------------------------------------------------------------------------
##
## Exact spelling may vary slightly across package versions, so labels are
## matched case-insensitively.

adult_cell_type <-
  as.character(
    colData(sce_adult)$cell.type
  )

observed_adult_types <-
  unique(
    adult_cell_type[
      !is.na(adult_cell_type)
    ]
  )

cat("\nObserved adult cell-type labels:\n")
print(sort(observed_adult_types))


neuron_hits <-
  observed_adult_types[
    grepl(
      "^neuron",
      observed_adult_types,
      ignore.case = TRUE
    )
  ]

astro_hits <-
  observed_adult_types[
    grepl(
      "^astro",
      observed_adult_types,
      ignore.case = TRUE
    )
  ]


if (length(neuron_hits) != 1) {
  stop(
    "Expected exactly one neuron label, but found: ",
    paste(neuron_hits, collapse = ", ")
  )
}

if (length(astro_hits) != 1) {
  stop(
    "Expected exactly one astrocyte label, but found: ",
    paste(astro_hits, collapse = ", ")
  )
}


neuron_label <- neuron_hits
astro_label <- astro_hits

candidate_types <- c(
  neuron_label,
  astro_label
)

cat("\nSelected benchmark labels:\n")
print(candidate_types)


## ---------------------------------------------------------------------------
## 5. Keep adult neurons and astrocytes only
## ---------------------------------------------------------------------------

keep_type <-
  !is.na(adult_cell_type) &
  adult_cell_type %in% candidate_types

sce_sub <-
  sce_adult[, keep_type]

type <-
  droplevels(
    factor(
      colData(sce_sub)$cell.type,
      levels = candidate_types
    )
  )

cat("\nAdult neuron/astrocyte class counts:\n")
print(table(type))

cat(
  "\nTotal retained cells:",
  ncol(sce_sub),
  "\n"
)


## Verify both classes are present
if (nlevels(type) != 2) {
  stop(
    "Expected exactly 2 retained classes, but found ",
    nlevels(type),
    "."
  )
}


## Protect against unexpectedly small classes
min_class_n <- 20

class_counts <- table(type)

if (any(class_counts < min_class_n)) {
  stop(
    "At least one selected class contains fewer than ",
    min_class_n,
    " cells. Counts: ",
    paste(
      names(class_counts),
      as.integer(class_counts),
      sep = "=",
      collapse = ", "
    )
  )
}


## ---------------------------------------------------------------------------
## 6. Extract raw counts
## ---------------------------------------------------------------------------

if (!"counts" %in% assayNames(sce_sub)) {
  stop("Raw assay 'counts' was not found.")
}

counts_mat <-
  assay(
    sce_sub,
    "counts"
  )


## ---------------------------------------------------------------------------
## 7. Raw-count sparsity before filtering
## ---------------------------------------------------------------------------
##
## Sparsity = proportion of zero entries in the raw-count matrix.

sparsity_before <-
  mean(
    counts_mat == 0
  )

cat(
  "\nRaw-count sparsity before filtering:",
  round(sparsity_before, 4),
  "\n"
)


## ---------------------------------------------------------------------------
## 8. Basic cell filtering
## ---------------------------------------------------------------------------
##
## Remove zero-library cells as a minimal QC safeguard.

lib_sizes <-
  colSums(
    counts_mat
  )

keep_cells <-
  lib_sizes > 0

sce_sub <-
  sce_sub[, keep_cells]

type <-
  droplevels(
    type[keep_cells]
  )

cat(
  "\nCells after zero-library filtering:",
  ncol(sce_sub),
  "\n"
)

cat("\nClass counts after cell filtering:\n")
print(table(type))


## ---------------------------------------------------------------------------
## 9. Basic gene filtering
## ---------------------------------------------------------------------------
##
## Remove genes detected in fewer than 5 retained cells.
##
## This removes extremely rare genes while preserving the sparse nature of
## the scRNA-seq expression matrix.

counts_mat <-
  assay(
    sce_sub,
    "counts"
  )

keep_genes <-
  rowSums(
    counts_mat > 0
  ) >= 5

sce_sub <-
  sce_sub[
    keep_genes,
  ]

cat(
  "\nDimensions after gene filtering (genes x cells):\n"
)

print(
  dim(sce_sub)
)


## ---------------------------------------------------------------------------
## 10. Log-normalization
## ---------------------------------------------------------------------------
##
## logNormCounts() performs library-size normalization followed by log
## transformation and stores the result in the "logcounts" assay.

sce_sub <-
  scater::logNormCounts(
    sce_sub
  )


## ---------------------------------------------------------------------------
## 11. Select top 200 highly variable genes
## ---------------------------------------------------------------------------
##
## Primary analysis:
##   p = 200 HVGs
##
## This keeps the final clustering problem computationally manageable while
## remaining high-dimensional.

n_hvg <- 200

dec <-
  scran::modelGeneVar(
    sce_sub
  )

n_hvg_eff <-
  min(
    n_hvg,
    nrow(sce_sub)
  )

top_hvgs <-
  scran::getTopHVGs(
    dec,
    n = n_hvg_eff
  )

sce_hvg <-
  sce_sub[
    top_hvgs,
  ]

cat(
  "\nNumber of selected HVGs:",
  nrow(sce_hvg),
  "\n"
)

cat("\nFirst 20 selected HVGs:\n")

print(
  head(
    rownames(sce_hvg),
    20
  )
)

## Biological variance cutoff for top 200 HVGs
hvg_cutoff <- min(
  dec[top_hvgs, "bio"],
  na.rm = TRUE
)

## Percentile of that cutoff among all genes
hvg_cutoff_percentile <- mean(
  dec$bio <= hvg_cutoff,
  na.rm = TRUE
)

cat(
  "The top 200 HVGs correspond to genes above approximately the",
  round(100 * hvg_cutoff_percentile, 1),
  "th percentile of biological variance.\n"
)


## ---------------------------------------------------------------------------
## 12. Sparsity after HVG selection
## ---------------------------------------------------------------------------
##
## Sparsity is measured on the raw-count scale.

counts_hvg <-
  assay(
    sce_hvg,
    "counts"
  )

sparsity_after <-
  mean(
    counts_hvg == 0
  )

cat(
  "\nRaw-count sparsity before filtering:",
  round(sparsity_before, 4),
  "\n"
)

cat(
  "Raw-count sparsity after gene filtering/HVG selection:",
  round(sparsity_after, 4),
  "\n"
)


## ---------------------------------------------------------------------------
## 13. Construct cells x genes log-normalized matrix
## ---------------------------------------------------------------------------

X_log <-
  t(
    as.matrix(
      logcounts(sce_hvg)
    )
  )

storage.mode(X_log) <-
  "double"

cat(
  "\nLog-normalized matrix dimensions (cells x HVGs):\n"
)

print(
  dim(X_log)
)

type <-
  droplevels(
    type
  )

## Check consistency

stopifnot(
  nrow(X_log) == length(type)
)

stopifnot(
  nlevels(type) == 2
)

stopifnot(
  ncol(X_log) == nrow(sce_hvg)
)

## Create final dataset

darmanis_adult_2class_df <- data.frame(cell_type = type, X_log)
