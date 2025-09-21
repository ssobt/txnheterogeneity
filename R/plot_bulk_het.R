## plot_bulk_het() ##

#' plot_bulk_het()
#' @description Plot output of bulk_het()
#'
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 geom_point
#' @importFrom ggplot2 geom_violin
#' @importFrom ggplot2 labs
#' @importFrom ggplot2 ggtitle
#' @importFrom ggplot2 theme
#'
#' @param bulk_out A list containing the output of bulk_het().
#' @param plot_type A string indicating the type of plot to generate. Options are 'PC scatter' or 'distance violin'.
#' @param components_to_plot A numeric vector indicating the PCA components to plot if plot_type = 'PC scatter'.
#' @param stratifier_gene A string indicating the stratifier gene to plot. Required if using quantile-based grouping.
#' @returns A ggplot object.
#' @examples
#' set.seed(123)
#' out = bulk_het(data = input_mtx, quant = 0.75, stratifiers = rownames(input_mtx)[1:10], cores = 10) 
#' plot_bulk_het(bulk_out = out, components_to_plot = c(1,2), plot_type = 'PC scatter', stratifier_gene = 'UBE2Q2P3') + ggplot2::scale_color_brewer(palette = "Dark2")
#'
#' plot_bulk_het(bulk_out = out, plot_type = 'distance violin', stratifier_gene = 'UBE2Q2P3')
#' @export

######## use devtools::document() to make the actual .Rd file in man folder for this function ########

plot_bulk_het <- function(bulk_out, plot_type, components_to_plot, stratifier_gene) {
    if (plot_type == 'PC scatter' & length(components_to_plot) == 2 & ('PCA_coordinates_g1' %in% names(bulk_out)) & ('PCA_coordinates_g2' %in% names(bulk_out))) {
        # Create a combined data frame for plotting
        PCA_coords_g1 = bulk_out$PCA_coordinates_g1
        PCA_coords_g2 = bulk_out$PCA_coordinates_g2
        combined_PCA_coords <- data.frame(
            PC_first = c(PCA_coords_g1[, components_to_plot[1]], PCA_coords_g2[, components_to_plot[1]]),
            PC_second = c(PCA_coords_g1[, components_to_plot[2]], PCA_coords_g2[, components_to_plot[2]]),
            Group = factor(c(rep("g1", nrow(PCA_coords_g1)), rep("g2", nrow(PCA_coords_g2))))
        )
    # Generate the plot
    plot = ggplot(combined_PCA_coords, aes(x = PC_first, y = PC_second, color = Group)) +
        geom_point(alpha = 1) + # facet_wrap(~Grouping_gene) +
        ggpubr::theme_pubr() + theme(legend.position = "right")
        labs(x = paste0("PC", components_to_plot[1]), y = paste0("PC", components_to_plot[2]))
    return(plot)
    }

    if (plot_type == 'PC scatter' & length(components_to_plot) == 2 & ('upper_quant_PCA_coordinates' %in% names(bulk_out)) & ('lower_quant_PCA_coordinates' %in% names(bulk_out)) & is.character(stratifier_gene)) {
        PCA_coords_upper = bulk_out$upper_quant_PCA_coordinates[[stratifier_gene]]
        PCA_coords_lower = bulk_out$lower_quant_PCA_coordinates[[stratifier_gene]]

        # Create a combined data frame for plotting
        combined_PCA_coords <- data.frame(
            PC_first = c(PCA_coords_upper[,components_to_plot[1]], PCA_coords_lower[,components_to_plot[1]]),
            PC_second = c(PCA_coords_upper[,components_to_plot[2]], PCA_coords_lower[,components_to_plot[2]]),
            Group = factor(c(rep(paste0('q', 100*out$quant), nrow(PCA_coords_upper)), rep(paste0('q', 100*(1-out$quant)), nrow(PCA_coords_lower))))
        )
        # Generate the plot
        plot = ggplot(combined_PCA_coords, aes(x = PC_first, y = PC_second, color = Group)) +
            geom_point(alpha = 1) + 
            ggpubr::theme_pubr() +
            labs(x = paste0("PC", components_to_plot[1]), y = paste0("PC", components_to_plot[2])) + ggtitle(stratifier_gene) + theme(plot.title = element_text(hjust = 0.5), legend.position = "right")
        return(plot)
    } 

    if (plot_type == 'distance violin' & length(components_to_plot) == 2 & ('distance_data_g1' %in% names(bulk_out)) & ('distance_data_g2' %in% names(bulk_out))) {
        # Create a combined data frame for plotting
        combined_distance_data <- data.frame(
            Distance = c(bulk_out$distance_data_g1, bulk_out$distance_data_g2),
            Group = factor(c(rep("g1", length(bulk_out$distance_data_g1)), rep("g2", length(bulk_out$distance_data_g2))))
        )
        # Generate the plot
        plot = ggplot(combined_distance_data, aes(x = Group, y = Distance, fill = Group)) +
            geom_violin(position=position_dodge(), alpha = 0.8) +
            ggforce::geom_sina(position=position_dodge(), size = 0.5) +
            ggpubr::theme_pubr() +
            labs(x = "Group", y = "Weighted distance") + ggtitle('Weighted distance from centroid in PC space') + theme(plot.title = element_text(hjust = 0.5), legend.position = "right")
        return(plot)
        ## check this violin plot code for accuracy
    } 

    if (plot_type == 'distance violin' & length(components_to_plot) == 2 & ('upper_quant_distance_data' %in% names(bulk_out)) & ('lower_quant_distance_data' %in% names(bulk_out)) & is.character(stratifier_gene)) {
        # Create a combined data frame for plotting
        combined_distance_data <- data.frame(
            Distance = c(bulk_out$upper_quant_distance_data[[stratifier_gene]], bulk_out$lower_quant_distance_data[[stratifier_gene]]),
            Group = factor(c(rep(paste0('q', 100*out$quant), length(bulk_out$upper_quant_distance_data[[stratifier_gene]])), rep(paste0('q', 100*(1-out$quant)), length(bulk_out$lower_quant_distance_data[[stratifier_gene]]))))
        )
        # Generate the plot
        plot = ggplot(combined_distance_data, aes(x = Group, y = Distance, fill = Group)) +
            geom_violin(position=position_dodge(), alpha = 0.8) +
            ggforce::geom_sina(position=position_dodge(), size = 0.5) +
            ggpubr::theme_pubr() +
            labs(x = "Group", y = "Weighted distance") + ggtitle('Weighted distance from centroid in PC space') + theme(plot.title = element_text(hjust = 0.5), legend.position = "right")
        return(plot)
        ## check this violin plot code for accuracy
    } 


}

