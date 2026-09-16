# plot_bulk_het()

Plot output of bulk_het()

## Usage

``` r
plot_bulk_het(bulk_out, plot_type, components_to_plot, stratifier_gene)
```

## Arguments

- bulk_out:

  A list containing the output of bulk_het().

- plot_type:

  A string indicating the type of plot to generate. Options are 'PC
  scatter' or 'distance violin'.

- components_to_plot:

  A numeric vector indicating the PCA components to plot if plot_type =
  'PC scatter'.

- stratifier_gene:

  A string indicating the stratifier gene to plot. Required if using
  quantile-based grouping.

## Value

A ggplot object.

## Examples

``` r
set.seed(123)
out = bulk_het(data = input_mtx, quant = 0.75, stratifiers = rownames(input_mtx)[1:10], cores = 10) 
plot_bulk_het(bulk_out = out, components_to_plot = c(1,2), plot_type = 'PC scatter', stratifier_gene = 'UBE2Q2P3') + ggplot2::scale_color_brewer(palette = "Dark2")

plot_bulk_het(bulk_out = out, plot_type = 'distance violin', stratifier_gene = 'UBE2Q2P3')
```
