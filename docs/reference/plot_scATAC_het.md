# plot_scATAC_het()

Plot output of scATAC_het()

## Usage

``` r
plot_scATAC_het(scATAC_het_out, plot_type = NULL, sample_names = NULL)
```

## Arguments

- scATAC_het_out:

  A list containing the output of scATAC_het().

- plot_type:

  A string indicating the type of plot to generate. Options are 'CV
  ratios', 'mean ratios', 'CV between samples', or 'mean between
  samples'.

- sample_names:

  A string indicating the sample names to plot, do not include control
  sample. Default is all samples.

## Value

A ggplot object.

## Examples

``` r
plot_scATAC_het(scATAC_het_out, plot_type = 'CV ratios', sample_names = c('RNF8-Ci', 'MIS18A-Ci'))
plot_scATAC_het(scATAC_het_out, plot_type = 'CV between samples', sample_names = c('RNF8-Ci', 'MIS18A-Ci'))
```
