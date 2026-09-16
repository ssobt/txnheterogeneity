# gene_retention_pct()

Calculate the percentage of genes that would be retained after filtering
based on a specified median expression cutoff. It is recommended to find
a cutoff that eliminates low expressing genes that are too noisy for
analysis (usually ~10-20% of genes are retained).

## Usage

``` r
gene_retention_pct(seurat_obj, cutoff)
```

## Arguments

- seurat_obj:

  A Seurat object containing the scRNA-seq data. The object should be
  pre-processed for cell quality control only.

- cutoff:

  A numeric value indicating the minimum median expression threshold for
  genes to be included in the analysis. Default is 0.1.

## Value

A numeric value indicating the percentage of genes that would be
retained after filtering based on the specified cutoff.

## Examples

``` r
gene_retention_pct(seurat_obj = CRISPRa_seurat, cutoff = 0.1)
```
