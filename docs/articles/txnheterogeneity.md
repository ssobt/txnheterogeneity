# Quantifying transcriptional heterogeneity

> Code is shown but not executed — every function here needs a real
> dataset (a Seurat object, an expression matrix, or 10x fragment
> files), none of which is small enough to ship inside the package.

## What this package measures

Most differential analysis compares **means**: is this gene higher in
condition A than B? This package compares **dispersion**: are the cells
or samples in condition A more *variable* than those in condition B?

That distinction matters in cancer, where intratumoral heterogeneity
predicts poor outcome, and where the interesting question is often not
which genes moved but whether the population became more or less
uniform.

Three assays, one idea, one interface:

| Assay | Metric | Entry point |
|----|----|----|
| Bulk RNA-seq | Gene-level coefficient of variation; also dispersion in PC space | [`bulk_het()`](https://ssobt.github.io/txnheterogeneity/reference/bulk_het.md) |
| scRNA-seq | Transcriptome-wide CV per condition | [`sc_het()`](https://ssobt.github.io/txnheterogeneity/reference/sc_het.md) |
| scATAC-seq | CV over called peaks, via [ArchR](https://www.archrproject.com) | [`scATAC_het()`](https://ssobt.github.io/txnheterogeneity/reference/scATAC_het.md) |

Each `*_het()` function has a matching `plot_*_het()`.

## Installation

The dependency set is heavy — ArchR, a BSgenome, and several
Bioconductor packages. A conda environment is strongly recommended; one
is shipped with the package:

``` bash
conda env create -f inst/txnheterogeneity.yml
conda activate txnheterogeneity
```

Then, from R:

\
`if`` ``(``!`[`require`](https://rdrr.io/r/base/library.html)`(`[`"devtools"`](https://devtools.r-lib.org/)`)``)`` `[`install.packages`](https://rdrr.io/r/utils/install.packages.html)`(``"devtools"``)`\
`devtools``::`[`install_github`](https://devtools.r-lib.org/reference/install-deprecated.html)`(``"ssobt/txnheterogeneity"``)`

If dependencies do not resolve automatically, install into a dedicated
library:

\
`devtools``::`[`install_github`](https://devtools.r-lib.org/reference/install-deprecated.html)`(``"ssobt/txnheterogeneity"``,`\
`                         upgrade ``=`` ``"never"``,`\
`                         lib ``=`` ``"/path/to/txnheterogeneity_lib/"``)`

## Bulk RNA-seq

[`bulk_het()`](https://ssobt.github.io/txnheterogeneity/reference/bulk_het.md)
works two ways. Either you already know your groups:

\
[`.libPaths`](https://rdrr.io/r/base/libPaths.html)`(``"/path/to/txnheterogeneity_lib/"``)`\
[`library`](https://rdrr.io/r/base/library.html)`(``txnheterogeneity``)`\
\
`# input_mtx: genes in rows, samples in columns`\
[`set.seed`](https://rdrr.io/r/base/Random.html)`(``123``)``                       ``# the background distribution is sampled`\
`samples`` ``<-`` `[`sample`](https://rdrr.io/r/base/sample.html)`(`[`colnames`](https://rdrr.io/r/base/colnames.html)`(``input_mtx``)``, ``200``)`\
`out`` ``<-`` `[`bulk_het`](https://ssobt.github.io/txnheterogeneity/reference/bulk_het.md)`(``data ``=`` ``input_mtx``,`\
`                g1 ``=`` ``samples``[``1``:``100``]``,`\
`                g2 ``=`` ``samples``[``101``:``200``]``,`\
`                cores ``=`` ``1``)`

…or you want samples stratified by how strongly they express candidate
genes — the mode used to screen TCGA for regulators of heterogeneity:

\
[`set.seed`](https://rdrr.io/r/base/Random.html)`(``123``)`\
`out`` ``<-`` `[`bulk_het`](https://ssobt.github.io/txnheterogeneity/reference/bulk_het.md)`(``data        ``=`` ``input_mtx``,`\
`                quant       ``=`` ``0.75``,          ``# top 25% vs bottom 25%`\
`                stratifiers ``=`` `[`rownames`](https://rdrr.io/r/base/colnames.html)`(``input_mtx``)``[``1``:``10``]``,`\
`                cores       ``=`` ``10``)`

Plot either the PC-space spread or the pairwise distance distribution:

\
[`plot_bulk_het`](https://ssobt.github.io/txnheterogeneity/reference/plot_bulk_het.md)`(``out``, components_to_plot ``=`` `[`c`](https://rdrr.io/r/base/c.html)`(``1``, ``2``)``,`\
`              plot_type ``=`` ``"PC scatter"``, stratifier_gene ``=`` ``"MIS18A"``)`` ``+`\
`  ``ggplot2``::`[`scale_color_viridis_d`](https://ggplot2.tidyverse.org/reference/scale_viridis.html)`(``)`\
\
[`plot_bulk_het`](https://ssobt.github.io/txnheterogeneity/reference/plot_bulk_het.md)`(``out``, plot_type ``=`` ``"distance violin"``, stratifier_gene ``=`` ``"MIS18A"``)`

![Output of
plot_bulk_het()](../reference/figures/bulk_scatter_violin.png)

## Single-cell RNA-seq

Requires a Seurat object with a metadata column identifying each cell’s
sample or guide.

\
`sc_out`` ``<-`` `[`sc_het`](https://ssobt.github.io/txnheterogeneity/reference/sc_het.md)`(`\
`  seurat_obj                   ``=`` ``CRISPRa_seurat``,`\
`  meta_data_sample_column      ``=`` ``"guide"``,`\
`  sample_names                 ``=`` `[`c`](https://rdrr.io/r/base/c.html)`(``"NT"``, ``"RNF8"``, ``"MIS18A"``)``,`\
`  control_sample_name          ``=`` ``"NT"``,`\
`  cutoff                       ``=`` ``0.1``,   ``# gene retention threshold`\
`  sample_cells_per_guide_cutoff ``=`` ``50``,   ``# equalize cells per condition`\
`  seed                         ``=`` ``42`\
`)`

Two arguments deserve attention, because they are what keep the result
from being an artifact:

- **`cutoff`** sets the minimum detection rate for a gene to be counted.
  CV is unstable for genes that are mostly zero, and without a floor the
  metric ends up reporting dropout rather than biology.
  [`gene_retention_pct()`](https://ssobt.github.io/txnheterogeneity/reference/gene_retention_pct.md)
  shows how many genes survive a given cutoff.
- **`sample_cells_per_guide_cutoff`** downsamples every condition to the
  same number of cells. Dispersion estimates scale with sample size, so
  unequal cell counts alone will produce an apparent heterogeneity
  difference.

\
[`plot_sc_het`](https://ssobt.github.io/txnheterogeneity/reference/plot_sc_het.md)`(``sc_out``, plot_type ``=`` ``"CV violin"``)`\
[`plot_sc_het`](https://ssobt.github.io/txnheterogeneity/reference/plot_sc_het.md)`(``sc_out``, plot_type ``=`` ``"mean violin"``, sample_names ``=`` `[`c`](https://rdrr.io/r/base/c.html)`(``"RNF8-Ci"``, ``"MIS18A-Ci"``)``)`\
[`plot_sc_het`](https://ssobt.github.io/txnheterogeneity/reference/plot_sc_het.md)`(``sc_out``, plot_type ``=`` ``"cell cycle"``)``     ``# is the shift just cycling?`\
[`plot_sc_het`](https://ssobt.github.io/txnheterogeneity/reference/plot_sc_het.md)`(``sc_out``, plot_type ``=`` ``"heatmap"``, sample_names ``=`` ``"RNF8-Ci"``)`

![Output of
plot_sc_het()](../reference/figures/sc_ratios.png)![Cell-cycle breakdown
from plot_sc_het()](../reference/figures/sc_cc_ht.png)

The `"cell cycle"` view is a control, not decoration: proliferating
cells are transcriptionally more variable, so a heterogeneity shift that
is really a shift in cycling fraction should be visible here before it
is interpreted as anything else.

## Single-cell ATAC-seq

Takes fragment files from Cell Ranger ATAC and calls peaks through
ArchR.

\
`# Combine replicates by mapping each sample name to a merged label`\
`rep_map`` ``<-`` `[`list`](https://rdrr.io/r/base/list.html)`(`\
`  `[`c`](https://rdrr.io/r/base/c.html)`(``"NTCi-1"``, ``"NTCi-2"``, ``"RNF8-Ci-1"``, ``"RNF8-Ci-2"``, ``"MIS18A-Ci-1"``, ``"MIS18A-Ci-2"``)``,`\
`  `[`c`](https://rdrr.io/r/base/c.html)`(``"NTCi"``,   ``"NTCi"``,   ``"RNF8-Ci"``,   ``"RNF8-Ci"``,   ``"MIS18A-Ci"``,   ``"MIS18A-Ci"``)`\
`)`\
\
`out`` ``<-`` `[`scATAC_het`](https://ssobt.github.io/txnheterogeneity/reference/scATAC_het.md)`(`\
`  scATAC_fragment_files ``=`` `[`paste0`](https://rdrr.io/r/base/paste.html)`(``folder``, ``"/"``, ``file_names``)``,`\
`  sample_names          ``=`` `[`names`](https://rdrr.io/r/base/names.html)`(``file_names``)``,`\
`  replicate_map         ``=`` ``rep_map``,`\
`  samples_to_analyze    ``=`` `[`c`](https://rdrr.io/r/base/c.html)`(``"NTCi"``, ``"RNF8-Ci"``, ``"MIS18A-Ci"``)``,`\
`  control_sample_name   ``=`` ``"NTCi"``,`\
`  savepath              ``=`` ``savep``,`\
`  threads               ``=`` ``35``,`\
`  genome                ``=`` ``"hg38"``,`\
`  TSS_bed_path          ``=`` ``"/path/to/GRCh38_transcriptsOnly.tss.bed"``,`\
`  macs2_path            ``=`` ``"/path/to/envs/archr/bin/macs2"`\
`)`\
\
[`plot_scATAC_het`](https://ssobt.github.io/txnheterogeneity/reference/plot_scATAC_het.md)`(``out``, plot_type ``=`` ``"CV between samples"``,`\
`                sample_names ``=`` `[`c`](https://rdrr.io/r/base/c.html)`(``"RNF8-Ci"``, ``"MIS18A-Ci"``)``)`

![Output of plot_scATAC_het()](../reference/figures/scatac.png)

Pass `macs2_path` explicitly. ArchR searches the system for MACS2 and
picks the wrong copy when several exist, which is common once you have
more than one conda environment.

## Interpreting a heterogeneity change

Dispersion metrics are easy to fool. Before concluding that a
perturbation changed heterogeneity, rule out the alternatives:

1.  **Depth and cell number** — equalized via
    `sample_cells_per_guide_cutoff`.
2.  **Detection rate** — controlled by `cutoff`; check with
    [`gene_retention_pct()`](https://ssobt.github.io/txnheterogeneity/reference/gene_retention_pct.md).
3.  **Cell cycle** — inspect the `"cell cycle"` plot.
4.  **Mean–variance coupling** — CV falls as mean expression rises, so
    compare the `"CV violin"` against the `"mean violin"`. A CV shift
    with no mean shift is the interesting case.

Significance throughout comes from sampled background distributions
rather than a parametric assumption, because gene-level CV is not well
behaved.

## Related work

This package underpins:

> Woo BJ\*, Sobti S\*, Suh JM, Yousefi H, Garcia K, Zhou S, Borah A,
> Goodarzi H. *Systematic identification of chromatin organizers as
> tuners of intratumoral heterogeneity.* bioRxiv 2026.
> [doi:10.64898/2026.04.18.719392](https://doi.org/10.64898/2026.04.18.719392)

For pathway shifts in Perturb-seq data, see
[VariPath](https://ssobt.github.io/VariPath/).
