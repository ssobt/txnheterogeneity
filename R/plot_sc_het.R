## plot_sc_het() ##
# packages necessary: ggplot2, ggpubr, ggforce, ComplexHeatmap, grid, circlize, viridis, scales, dplyr, tidyr

#' plot_sc_het()
#' @description Plot output of sc_het()
#'
#' @import ggplot2
#' @import ggpubr
#' @import dplyr
#' @import tidyr
#'
#' @param sc_het_output A list containing the output of sc_het().
#' @param plot_type A string indicating the type of plot to generate. Options are 'CV violin', 'mean violin', 'cell cycle', or 'heatmap'. If 'heatmap', please provide sample_name that you want individual genes heatmapped for.
#' @param sample_names A string indicating the sample names to plot for CV violin, mean violin or heatmap (default is all samples). If plot_type = 'heatmap', a single sample name is required.
#' @param seed An integer to set the seed for reproducibility when generating heatmaps (default is 20).
#' @param y_label_position A numeric value indicating the y-axis position to place the adjusted p-values on the violin plots (default is 1.2).
#' @returns A ggplot object.
#' @examples
#' plot_sc_het(sc_het_output, plot_type = 'CV violin')
#' plot_sc_het(sc_het_output, plot_type = 'mean violin', sample_names = c('sgRNA_ABC', 'sgRNA_XYZ'))
#' plot_sc_het(sc_het_output, plot_type = 'cell cycle')
#' plot_sc_het(sc_het_output, plot_type = 'heatmap', sample_names = 'sgRNA_XYZ')
#' @export


plot_sc_het <- function(sc_het_output, plot_type = NULL, sample_names = NULL, seed = 20, y_label_position = NULL){

  if (is.null(plot_type)){
    stop("Please provide a plot type. Options are 'CV violin', 'mean violin', 'cell cycle', or 'heatmap'.")
  }
  if (plot_type == 'CV violin'){
      ratios_df = sc_het_output$CV_ratios_df
      pvals_adj = sc_het_output$CV_pvals_adj
      order_of_guides = sc_het_output$CV_order_of_guides
      master_df_list = sc_het_output$master_df_list
      control_sample_name = sc_het_output$control_sample_name


      plotted_guides = names(master_df_list)[!(names(master_df_list) %in% c(control_sample_name, paste0('random_', control_sample_name)))]
      if(is.null(sample_names)){
          sample_names = plotted_guides
      }else{
          sample_names = sample_names
      }
      if (control_sample_name %in% sample_names){
          sample_names = sample_names[sample_names != control_sample_name]
      }else{
          sample_names = sample_names
      }
      if (is.null(y_label_position)){
          y_label_position = 1.2
      }
      ratios_df = ratios_df[!(is.infinite(ratios_df$CV_gdCV_ctrlratio) | is.nan(ratios_df$CV_gdCV_ctrlratio) | is.na(ratios_df$CV_gdCV_ctrlratio)),]
      plot = ratios_df %>% dplyr::filter(guide %in% sample_names) %>% ggplot(aes(guide, CV_gdCV_ctrlratio, fill = guide)) + ylab('log2(CV_guide / CV_control_guide)') +
        geom_violin(position=position_dodge()) + annotate("text", x = 1:length(sample_names), y = y_label_position, size = 5, angle='45', hjust = -0.2, label = pvals_adj[sample_names]) +
        ggforce::geom_sina(position=position_dodge(), size = 0.5) + geom_hline(yintercept = 0) + coord_cartesian(clip = 'off', ylim = c(-2.5,2.5)) + xlab('') +
        theme_classic() + theme(plot.margin = margin(1,1,1,1, "cm"), plot.title = element_text(hjust = 0.5), legend.position = 'none', axis.text.y = element_text(face = 'bold', size = 14), axis.text.x = element_text(face = 'bold', size = 12, angle = 45, hjust = 1))
      
      return(plot)

  }

  if (plot_type == 'mean violin'){
      ratios_df = sc_het_output$mean_ratios_df
      pvals_adj = sc_het_output$mean_pvals_adj
      order_of_guides = sc_het_output$mean_order_of_guides
      master_df_list = sc_het_output$master_df_list
      control_sample_name = sc_het_output$control_sample_name

      plotted_guides = names(master_df_list)[!(names(master_df_list) %in% c(control_sample_name, paste0('random_', control_sample_name)))]
      if(is.null(sample_names)){
          sample_names = plotted_guides
      }else{
          sample_names = sample_names
      }
      if (control_sample_name %in% sample_names){
          sample_names = sample_names[sample_names != control_sample_name]
      }else{
          sample_names = sample_names
      }
      if (is.null(y_label_position)){
          y_label_position = 1.2
      }
      ratios_df = ratios_df[!(is.infinite(ratios_df$mean_gdmean_ctrlratio) | is.nan(ratios_df$mean_gdmean_ctrlratio) | is.na(ratios_df$mean_gdmean_ctrlratio)),]
      plot = ratios_df %>% dplyr::filter(guide %in% sample_names) %>% ggplot(aes(guide, mean_gdmean_ctrlratio, fill = guide)) + ylab('log2(mean_guide / mean_control_guide)') +
        geom_violin(position=position_dodge()) + annotate("text", x = 1:length(sample_names), y = y_label_position, size = 5, angle='45', hjust = -0.2, label = pvals_adj[sample_names]) +
        ggforce::geom_sina(position=position_dodge(), size = 0.5) + geom_hline(yintercept = 0) + coord_cartesian(clip = 'off', ylim = c(-2.5,2.5)) + xlab('') +
        theme_classic() + theme(plot.margin = margin(1,1,1,1, "cm"), plot.title = element_text(hjust = 0.5), legend.position = 'none', axis.text.y = element_text(face = 'bold', size = 14), axis.text.x = element_text(face = 'bold', size = 12, angle = 45, hjust = 1))
      
      return(plot)
  }

  if (plot_type == 'cell cycle'){
    ## cell cycle phase distribution across samples ##
    control_sample_name = sc_het_output$control_sample_name
    guide_subsetted_data = sc_het_output$guide_subsetted_data
    if (control_sample_name %in% sample_names){
        sample_names = sample_names
    }else{
        sample_names = c(control_sample_name, sample_names)
    }
    selected_guides_data = guide_subsetted_data[sample_names]
    col_sizes = sapply(selected_guides_data, ncol)
    selected_guides_names = unlist(mapply(FUN = function(x, y) rep(x, y), x = sample_names, y = col_sizes, SIMPLIFY = FALSE))
    select_guide_ids = lapply(selected_guides_data, colnames)
    select_guide_ids = unlist(select_guide_ids)
    cc_df = data.frame(guide = selected_guides_names, cell_names = select_guide_ids)
    cc_df = cbind(cc_df, sc_het_output$cc_meta_data[select_guide_ids, c('S.Score', 'G2M.Score')])
    rownames(cc_df) = NULL
    cc_df = cc_df %>% tidyr::pivot_longer(cols = c('S.Score', 'G2M.Score'), names_to = 'cycle', values_to = 'score') 
    cc_df$guide = factor(cc_df$guide, levels = sample_names)
    cc_df$cycle = factor(cc_df$cycle, levels = c('S.Score', 'G2M.Score'))
    
    ## get all combinations between control and other samples for statistical comparison ##
    v1 <- c(control_sample_name)
    v2 <- setdiff(sample_names, control_sample_name)
    # Cartesian product using expand.grid
    combo_df <- expand.grid(v1 = v1, v2 = v2, stringsAsFactors = FALSE)
    # Convert rows to list of vectors
    combo_list <- split(combo_df, seq(nrow(combo_df)))
    combo_list <- lapply(combo_list, function(row) c(row$v1, row$v2))
    names(combo_list) = NULL

    plot = cc_df %>% ggplot(aes(x = guide, y = score, fill = guide)) + geom_violin(position=position_dodge()) + 
      ggforce::geom_sina(position=position_dodge(), size = 0.5) + scale_fill_viridis_d(alpha = 0.7, name = 'Percentile') + xlab('') +
      facet_wrap(~cycle) + stat_compare_means(comparisons = combo_list, method = 'wilcox.test') + theme_pubr() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
    
    return(plot)
  }


  if (plot_type == 'heatmap'){
    ## heatmap as of all genes in the sample condition vs control ##
    if (is.null(sample_names) | !(sample_names %in% as.character(sc_het_output$order_of_guides) | length(sample_names) != 1)){
      stop("Please provide a single valid sample name to plot heatmap for.")
    }
    master_df_list = sc_het_output$master_df_list
    sample_cells_per_guide_cutoff = sc_het_output$sample_cells_per_guide_cutoff
    control_sample_name = sc_het_output$control_sample_name

    lower_sample = master_df_list[[sample_names]] %>% filter(CV_gdCV_ctrlratio < 1) %>% select(CV_ctrl, CV_gd) %>% arrange(CV_gd)
    upper_sample = master_df_list[[sample_names]] %>% filter(CV_gdCV_ctrlratio > 1) %>% select(CV_ctrl, CV_gd) %>% arrange(CV_ctrl)
    sample_df = rbind(lower_sample, upper_sample) %>% as.matrix()
    colnames(sample_df) = c(paste0(control_sample_name, " (n = ", sample_cells_per_guide_cutoff, ")"), paste0(sample_names, " (n = ", sample_cells_per_guide_cutoff, ")"))

    set.seed(seed)
    coloring_map <- circlize::colorRamp2(c(0.26, 1.6), scales::viridis_pal(option = 'D')(5)[c(1,5)])

    perturb_perturb_corr_mtx_htmp <- function(matrix, title) {
      ComplexHeatmap::draw(ComplexHeatmap::Heatmap(matrix, name = title, show_column_names = TRUE, show_row_names = FALSE, show_column_dend = FALSE,
                  col = coloring_map, show_row_dend = FALSE, use_raster = TRUE, row_title = 'Population\ntranscript CV', 
                  row_title_side = 'left', column_title_side = 'top', column_names_side = 'top',
                  column_names_gp = grid::gpar(fontsize = 10), row_names_gp = grid::gpar(fontsize = 1), column_names_rot = 45,
                  cluster_rows = FALSE, cluster_columns = FALSE, heatmap_legend_param = list(title = "CV", at = c(0.26, 2))), padding = unit(c(2, 2, 2, 12), "pt"))
    }

    return(perturb_perturb_corr_mtx_htmp(sample_df))

  }


}

