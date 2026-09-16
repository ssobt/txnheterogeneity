# sc_het()

Calculate transcriptional heterogeneity in scRNA-seq samples compared to
a control sample.

## Usage

``` r
sc_het(
  seurat_obj,
  cutoff = 0.1,
  seed = 42,
  sample_cells_per_guide_cutoff = NULL,
  meta_data_sample_column = NULL,
  sample_names = NULL,
  control_sample_name = NULL
)
```

## Arguments

- seurat_obj:

  A Seurat object containing the scRNA-seq data. The object should be
  pre-processed for cell quality control only.

- cutoff:

  A numeric value indicating the minimum median expression threshold for
  genes to be included in the analysis. Default is 0.1. It is
  recommended to adjust this cutoff using the gene_retention_pct()
  function to retain ~10-20% of genes.

- seed:

  An integer value to set the random seed for reproducibility. Default
  is 42.

- sample_cells_per_guide_cutoff:

  Number of cells to sample per guide. It is highly recommended to
  choose a number above 50 for a representative outlook of
  heterogeneity, even at the expense of losing some samples in the
  analysis. Default uses cell count of the sample with the minimum
  number of cells.

- meta_data_sample_column:

  A string indicating the column name in the Seurat object's metadata
  that contains the sample information for cells.

- sample_names:

  A character vector of sample names to include in the analysis that
  fall within samples in meta_data_sample_column.

- control_sample_name:

  A string indicating the name of the control sample in
  meta_data_sample_column. This sample will be used as the reference for
  heterogeneity comparisons.

  It is advised that you save the output of this function as an RDS file
  using saveRDS() immediately after running it for easy loading in
  future sessions and because memory failure can occur when running
  downstream plotting functions on large datasets.

## Value

A list with the following components:

- guide_subsetted_dataA list of matrices with the expression data of
  cells from each sample as well as samples constructed from randomizing
  cell identity to mitigate issues expected from noisy single cell
  expression data.

- master_df_listA list that contains CV values of each gene within a
  sample for all samples analyzed.

- mean_shifts_from_NTA list that contains gene expression mean of each
  gene for all samples analyzed.

- asymp_test_p_valsA dataframe containing the p-values indicating
  signiciant change in CV from the cvequality package's asymptotic test
  for each gene in sample vs control.

- order_of_guidesA reference for order of samples analyzed for use in
  plot_sc_het().

- significant_CV_gene_countCount of genes in each sample with
  significant change in CV vs control.

- control_sample_nameA reference for use in plot_sc_het().

- CV_pvals_adjAdjusted p-values for global CV changes of sample vs
  control using t.test().

- CV_order_of_guidesA reference for use in plot_sc_het().

- CV_ratios_dfRatios of CV values of genes in sample vs control for each
  sample.

- mean_pvals_adjAdjusted p-values for global mean changes of sample vs
  control using t.test().

- mean_order_of_guidesA reference for use in plot_sc_het().

- mean_ratios_dfRatios of mean values of genes in sample vs control for
  each sample.

- sample_cells_per_guide_cutoffThe number of cells sampled per guide for
  the analysis.

- cc_meta_dataMetadata of analyzed cells in Seurat object for use in
  plot_sc_het().

## Examples

``` r
sc_het_out = sc_het(seurat_obj = CRISPRa_seurat, cutoff = 0.1, seed = 42, sample_cells_per_guide_cutoff = 50, meta_data_sample_column = 'guide', sample_names = c('NT', 'RNF8', 'MIS18A'), control_sample_name = 'NT')
```
