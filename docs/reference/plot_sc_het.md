# plot_sc_het()

Plot output of sc_het()

## Usage

``` r
plot_sc_het(
  sc_het_output,
  plot_type = NULL,
  sample_names = NULL,
  seed = 20,
  y_label_position = NULL
)
```

## Arguments

- sc_het_output:

  A list containing the output of sc_het().

- plot_type:

  A string indicating the type of plot to generate. Options are 'CV
  violin', 'mean violin', 'cell cycle', or 'heatmap'. If 'heatmap',
  please provide sample_name that you want individual genes heatmapped
  for.

- sample_names:

  A string indicating the sample names to plot for CV violin, mean
  violin or heatmap (default is all samples). If plot_type = 'heatmap',
  a single sample name is required.

- seed:

  An integer to set the seed for reproducibility when generating
  heatmaps (default is 20).

- y_label_position:

  A numeric value indicating the y-axis position to place the adjusted
  p-values on the violin plots (default is 1.2).

## Value

A ggplot object.

## Examples

``` r
plot_sc_het(sc_het_output, plot_type = 'CV violin')
plot_sc_het(sc_het_output, plot_type = 'mean violin', sample_names = c('sgRNA_ABC', 'sgRNA_XYZ'))
plot_sc_het(sc_het_output, plot_type = 'cell cycle')
plot_sc_het(sc_het_output, plot_type = 'heatmap', sample_names = 'sgRNA_XYZ')
```
