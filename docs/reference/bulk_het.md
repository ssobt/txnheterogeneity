# bulk_het()

Calculate transcriptional heterogeneity in bulk RNA-seq data samples

## Usage

``` r
bulk_het(
  data,
  g1 = NULL,
  g2 = NULL,
  quant = NULL,
  stratifiers = NULL,
  cores = 1
)
```

## Arguments

- data:

  A matrix containing FPKM data with rows as genes and columns as
  samples.

- g1:

  A numeric vector with column names in data.

- g2:

  A numeric vector with column names in data.

- quant:

  If g1 and g2 are NULL, a number between 0.5 and 1 indicating the
  quantile to use for stratifying samples based on expression of
  stratifiers.

- stratifiers:

  A character vector of gene names to use for stratification.

- cores:

  Number of cores to use for parallel processing. Default is 1.

  It is advised that you save the output of this function as an RDS file
  using saveRDS() immediately after running it for easy loading in
  future sessions and because memory failure can occur when running
  downstream plotting functions on large datasets.

## Value

A list with the following components:

- reference_dfA data frame with gene names, CVs, and means for each
  group.

- genes_with_significant_CV_changeA data frame with counts of genes with
  significant CV changes between groups.

- z_scoreZ-score comparing observed significant CV changes to background
  distribution.

- z_score_pvalueP-value associated with the z-score.

- PCA_coordinates_g1PCA coordinates for group 1 samples.

- PCA_coordinates_g2PCA coordinates for group 2 samples.

- distance_data_g1Weighted distances to centroid for group 1 samples.

- distance_data_g2Weighted distances to centroid for group 2 samples.

OR

- reference_dfA data frame with gene names, CVs, and means for each
  group.

- z_scoreZ-score comparing observed significant CV changes to background
  distribution and assoiciated p-value.

- upper_quant_PCA_coordinatesPCA coordinates for upper quantile samples
  for each provided stratified gene.

- lower_quant_PCA_coordinatesPCA coordinates for lower quantile samples
  for each provided stratified gene.

- upper_quant_distance_dataWeighted distances to centroid for upper
  quantile samples for each provided stratified gene.

- lower_quant_distance_dataWeighted distances to centroid for lower
  quantile samples for each provided stratified gene.

## Examples

``` r
## Example 1: Compare heterogeneity between two known groups of samples
## input_mtx is a matrix with rows as genes and columns as samples
## g1 and g2 are character vectors with column names in input_mtx
# set seed for reproducibility of bkg distribution sampling
set.seed(123)
samples = sample(colnames(input_mtx), 200)
g1 = samples[1:100]
g2 = samples[101:200]
out = bulk_het(data = input_mtx, g1 = g1, g2 = g2, cores = 1)


## Example 2: Compare heterogeneity between samples stratified by expression of certain genes
## input_mtx is a matrix with rows as genes and columns as samples
# set seed for reproducibility of bkg distribution sampling
set.seed(123)
## will stratify samples based on top and bottom 25% expression of each stratifier gene
## top 25% (q75) will be compared to bottom 25% (q25)
out = bulk_het(data = input_mtx, quant = 0.75, stratifiers = rownames(input_mtx)[1:10], cores = 10) 
```
