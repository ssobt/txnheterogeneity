## plot_scATAC_het ##

# packages necessary: ggplot2, ggpubr, ggforce, grid, scales, dplyr, tidyr

#' plot_scATAC_het()
#' @description Plot output of scATAC_het()
#'
#' @import ggplot2
#' @import ggpubr
#' @import dplyr
#' @import tidyr
#'
#' @param scATAC_het_out A list containing the output of scATAC_het().
#' @param plot_type A string indicating the type of plot to generate. Options are 'CV ratios', 'mean ratios', 'CV between samples', or 'mean between samples'.
#' @param sample_names A string indicating the sample names to plot, do not include control sample. Default is all samples.
#' @returns A ggplot object.
#' @examples
#' plot_scATAC_het(scATAC_het_out, plot_type = 'CV ratios', sample_names = c('RNF8-Ci', 'MIS18A-Ci'))
#' plot_scATAC_het(scATAC_het_out, plot_type = 'CV between samples', sample_names = c('RNF8-Ci', 'MIS18A-Ci'))
#' @export


plot_scATAC_het <- function(scATAC_het_out, plot_type = NULL, sample_names = NULL) {

    if (is.null(plot_type)) {
        stop("Please provide a plot type: 'CV ratios', 'mean ratios', 'CV between samples', or 'mean between samples'")
    }
    

    ## violin plot for CV ratios ##
    
    if (plot_type == "CV ratios") {
        CV_ratios_df = scATAC_het_out$CV_ratios_df
        pvals_adj = scATAC_het_out$CV_pvals_adj
        control_sample_name = scATAC_het_out$control_sample_name
        CV_ratios_df = CV_ratios_df[!(is.infinite(CV_ratios_df$CV_gdCV_ctrlratio) | is.nan(CV_ratios_df$CV_gdCV_ctrlratio) | is.na(CV_ratios_df$CV_gdCV_ctrlratio)),]

        if (is.null(sample_names)) {
            sample_names = unique(CV_ratios_df$guide)
            sample_names = sample_names[!(sample_names %in% c(control_sample_name, paste0('random_', control_sample_name)))]
        }

        #plotted_guides = names(master_df_list)[!(names(master_df_list) %in% c(control_sample_name, paste0('random_', control_sample_name)))]
        plot = CV_ratios_df %>% dplyr::filter(guide %in% sample_names) %>% ggplot(aes(guide, CV_gdCV_ctrlratio, fill = guide)) + ylab(paste0('log2(Gene accessibility CV normalized to ', control_sample_name, ')')) + xlab('') +
        geom_violin(position=position_dodge()) + annotate("text", x = (1:length(sample_names)), y = 3.5, size = 5, angle='45', label = pvals_adj[sample_names]) +
        ggforce::geom_sina(position=position_dodge(), size = 0.5) + geom_hline(yintercept = 0) + coord_cartesian(clip = 'off') +
        theme_classic() + theme(plot.margin = margin(1,1,1,1, "cm"), plot.title = element_text(hjust = 0.5), legend.position = 'none', axis.text.y = element_text(face = 'bold', size = 14), axis.text.x = element_text(face = 'bold', size = 12, angle = 45, hjust = 1))

        return(plot)
    }

    
    ## violin plot for mean ratios ##

    if (plot_type == "mean ratios") {
        mean_ratios_df = scATAC_het_out$mean_ratios_df
        pvals_adj = scATAC_het_out$mean_pvals_adj
        control_sample_name = scATAC_het_out$control_sample_name
        mean_ratios_df = mean_ratios_df[!(is.infinite(mean_ratios_df$mean_gdmean_ctrlratio) | is.nan(mean_ratios_df$mean_gdmean_ctrlratio) | is.na(mean_ratios_df$mean_gdmean_ctrlratio)),]

        if (is.null(sample_names)) {
            sample_names = unique(mean_ratios_df$guide)
            sample_names = sample_names[!(sample_names %in% c(control_sample_name, paste0('random_', control_sample_name)))]
        }

        #plotted_guides = names(master_df_list)[!(names(master_df_list) %in% c(control_sample_name, paste0('random_', control_sample_name)))]
        plot = mean_ratios_df %>% dplyr::filter(guide %in% sample_names) %>% ggplot(aes(guide, mean_gdmean_ctrlratio, fill = guide)) + ylab(paste0('log2(Gene accessibility mean normalized to ', control_sample_name, ')')) + xlab('') +
        geom_violin(position=position_dodge()) + annotate("text", x = (1:length(sample_names)), y = 3.5, size = 5, angle='45', label = pvals_adj[sample_names]) +
        ggforce::geom_sina(position=position_dodge(), size = 0.5) + geom_hline(yintercept = 0) + coord_cartesian(clip = 'off') +
        theme_classic() + theme(plot.margin = margin(1,1,1,1, "cm"), plot.title = element_text(hjust = 0.5), legend.position = 'none', axis.text.y = element_text(face = 'bold', size = 14), axis.text.x = element_text(face = 'bold', size = 12, angle = 45, hjust = 1))

        return(plot)
    }



    ## violin plot for CV between samples ##

    if (plot_type == "CV between samples") {
        CVs_df = scATAC_het_out$CVs_df
        control_sample_name = scATAC_het_out$control_sample_name
        CVs_df = CVs_df[!(is.infinite(CVs_df$CV_gd) | is.nan(CVs_df$CV_gd) | is.na(CVs_df$CV_gd)),]

        if (is.null(sample_names)) {
            sample_names = unique(CVs_df$guide)
            sample_names = sample_names[!(sample_names %in% c(control_sample_name, paste0('random_', control_sample_name)))]
        }

        ## get all combinations between control and other samples for statistical comparison ##
        v1 <- c(control_sample_name)
        v2 <- setdiff(sample_names, control_sample_name)
        # Cartesian product using expand.grid
        combo_df <- expand.grid(v1 = v1, v2 = v2, stringsAsFactors = FALSE)
        # Convert rows to list of vectors
        combo_list <- split(combo_df, seq(nrow(combo_df)))
        combo_list <- lapply(combo_list, function(row) c(row$v1, row$v2))
        names(combo_list) = NULL
        
        comparisons = combo_list
        plot = ggplot(CVs_df, aes(guide, CV_gd)) + ylab('log2(Gene accessibility CV)') + xlab('') +
        geom_violin(fill = 'black', color = 'black', position=position_dodge()) + #annotate("text", x = 1:length(plotted_guides), y = 3.5, size = 5, angle='45', hjust = -0.2, label = pvals_adj) +
        ggforce::geom_sina(position=position_dodge(), size = 0.5) + coord_cartesian(clip = 'off') + ggpubr::stat_compare_means(comparisons = comparisons, method = 't.test', na.rm = T, step.increase = 0.1) + 
        theme_classic() + theme(plot.margin = margin(1,1,1,1, "cm"), plot.title = element_text(hjust = 0.5), legend.position = 'none', axis.text.y = element_text(face = 'bold', size = 14), axis.text.x = element_text(face = 'bold', size = 12, angle = 45, hjust = 1))

        return(plot)

    }


    ## violin plot for mean between samples ##

    if (plot_type == "mean between samples") {
        means_df = scATAC_het_out$means_df
        control_sample_name = scATAC_het_out$control_sample_name
        means_df = means_df[!(is.infinite(means_df$mean_gd) | is.nan(means_df$mean_gd) | is.na(means_df$mean_gd)),]

        if (is.null(sample_names)) {
            sample_names = unique(means_df$guide)
            sample_names = sample_names[!(sample_names %in% c(control_sample_name, paste0('random_', control_sample_name)))]
        }

        ## get all combinations between control and other samples for statistical comparison ##
        v1 <- c(control_sample_name)
        v2 <- setdiff(sample_names, control_sample_name)
        # Cartesian product using expand.grid
        combo_df <- expand.grid(v1 = v1, v2 = v2, stringsAsFactors = FALSE)
        # Convert rows to list of vectors
        combo_list <- split(combo_df, seq(nrow(combo_df)))
        combo_list <- lapply(combo_list, function(row) c(row$v1, row$v2))
        names(combo_list) = NULL
        
        comparisons = combo_list
        plot = ggplot(means_df, aes(guide, mean_gd)) + ylab('log2(Gene accessibility mean)') + xlab('') +
        geom_violin(fill = 'black', color = 'black', position=position_dodge()) + #annotate("text", x = 1:length(plotted_guides), y = 3.5, size = 5, angle='45', hjust = -0.2, label = pvals_adj) +
        ggforce::geom_sina(position=position_dodge(), size = 0.5) + coord_cartesian(clip = 'off') + ggpubr::stat_compare_means(comparisons = comparisons, method = 't.test', na.rm = T, step.increase = 0.1) + 
        theme_classic() + theme(plot.margin = margin(1,1,1,1, "cm"), plot.title = element_text(hjust = 0.5), legend.position = 'none', axis.text.y = element_text(face = 'bold', size = 14), axis.text.x = element_text(face = 'bold', size = 12, angle = 45, hjust = 1))

        return(plot)

    }
}