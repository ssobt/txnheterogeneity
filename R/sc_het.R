## Input bulk scRNA-seq data in a Seurat object ##
## User needs to perform their own cell QC, but provide pre-normalized data ##
## Used to compare heterogeneity between different guides or conditions. An unperturbed control is required. ##

# packages necessary: cvequality, qvalue, Seurat, pbapply, sparseMatrixStats, dplyr, tidyr

#' @noRd
mtx_subsetter = function(guide_name, meta_data_sample_column, matrix, meta_data){return(matrix[, as.character(meta_data[[meta_data_sample_column]]) == guide_name])}

#' @noRd
CV_calculator = function(x){return((sparseMatrixStats::rowSds(x))/(sparseMatrixStats::rowMeans2(x)))}

#' gene_retention_pct()
#' @description Calculate the percentage of genes that would be retained after filtering based on a specified median expression cutoff. It is recommended to find a cutoff that eliminates low expressing genes that are too noisy for analysis (usually ~10-20% of genes are retained).
#' @param seurat_obj A Seurat object containing the scRNA-seq data. The object should be pre-processed for cell quality control only.
#' @param cutoff A numeric value indicating the minimum median expression threshold for genes to be included in the analysis. Default is 0.1.
#' @returns A numeric value indicating the percentage of genes that would be retained after filtering based on the specified cutoff.
#' @examples
#' gene_retention_pct(seurat_obj = CRISPRa_seurat, cutoff = 0.1)
#'
#' @export
gene_retention_pct = function(seurat_obj, cutoff){
    filtered_raw_mtx <- seurat_obj@assays$RNA@data
    medians = sparseMatrixStats::rowMedians(filtered_raw_mtx)
    filtered_genes = as.numeric(medians) >= cutoff
    pct = 100*sum(filtered_genes)/length(filtered_genes)
    return(pct)
}

##### Heterogeneity functions #####

## sc_het() ##

#' sc_het()
#' @description Calculate transcriptional heterogeneity in scRNA-seq samples compared to a control sample.
#' @import dplyr
#' @import tidyr
#' @param seurat_obj A Seurat object containing the scRNA-seq data. The object should be pre-processed for cell quality control only.
#' @param sample_cells_per_guide_cutoff Number of cells to sample per guide. It is highly recommended to choose a number above 50 for a representative outlook of heterogeneity, even at the expense of losing some samples in the analysis. Default uses cell count of the sample with the minimum number of cells.
#' @param cutoff A numeric value indicating the minimum median expression threshold for genes to be included in the analysis. Default is 0.1. It is recommended to adjust this cutoff using the gene_retention_pct() function to retain ~10-20% of genes.
#' @param seed An integer value to set the random seed for reproducibility. Default is 42.
#' @param meta_data_sample_column A string indicating the column name in the Seurat object's metadata that contains the sample information for cells.
#' @param sample_names A character vector of sample names to include in the analysis that fall within samples in meta_data_sample_column.
#' @param control_sample_name A string indicating the name of the control sample in meta_data_sample_column. This sample will be used as the reference for heterogeneity comparisons.
#'
#' It is advised that you save the output of this function as an RDS file using saveRDS() immediately after running it for easy loading in future sessions and because memory failure can occur when running downstream plotting functions on large datasets.
#' @returns A list with the following components:
#' \itemize{
#'  \item{guide_subsetted_data}{A list of matrices with the expression data of cells from each sample as well as samples constructed from randomizing cell identity to mitigate issues expected from noisy single cell expression data.}
#'  \item{master_df_list}{A list that contains CV values of each gene within a sample for all samples analyzed.}
#'  \item{mean_shifts_from_NT}{A list that contains gene expression mean of each gene for all samples analyzed.}
#'  \item{asymp_test_p_vals}{A dataframe containing the p-values indicating signiciant change in CV from the cvequality package's asymptotic test for each gene in sample vs control.}
#'  \item{order_of_guides}{A reference for order of samples analyzed for use in plot_sc_het().}
#'  \item{significant_CV_gene_count}{Count of genes in each sample with significant change in CV vs control.}
#'  \item{control_sample_name}{A reference for use in plot_sc_het().}
#'  \item{CV_pvals_adj}{Adjusted p-values for global CV changes of sample vs control using t.test().}
#'  \item{CV_order_of_guides}{A reference for use in plot_sc_het().}
#'  \item{CV_ratios_df}{Ratios of CV values of genes in sample vs control for each sample.}
#'  \item{mean_pvals_adj}{Adjusted p-values for global mean changes of sample vs control using t.test().}
#'  \item{mean_order_of_guides}{A reference for use in plot_sc_het().}
#'  \item{mean_ratios_df}{Ratios of mean values of genes in sample vs control for each sample.}
#'  \item{sample_cells_per_guide_cutoff}{The number of cells sampled per guide for the analysis.}
#'  \item{cc_meta_data}{Metadata of analyzed cells in Seurat object for use in plot_sc_het().}
#' }
#'
#' @examples
#' sc_het_out = sc_het(seurat_obj = CRISPRa_seurat, cutoff = 0.1, seed = 42, sample_cells_per_guide_cutoff = 50, meta_data_sample_column = 'guide', sample_names = c('NT', 'RNF8', 'MIS18A'), control_sample_name = 'NT')
#' @export



sc_het <- function(seurat_obj, cutoff = 0.1, seed = 42, sample_cells_per_guide_cutoff = NULL, meta_data_sample_column = NULL, sample_names = NULL, control_sample_name = NULL){
    set.seed(seed)
    s.genes <- Seurat::cc.genes$s.genes
    g2m.genes <- Seurat::cc.genes$g2m.genes
    seurat_obj_with_cc_scores <- CellCycleScoring(seurat_obj, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
    seurat_obj <- NormalizeData(seurat_obj, normalization.method = 'RC')
    filtered_raw_mtx <- seurat_obj@assays$RNA@data
    medians = sparseMatrixStats::rowMedians(filtered_raw_mtx)
    median_df = data.frame(gene_medians = medians)

    ### removing genes with medians < cutoff -- low expressing genes have noisy expression and confound results ###
    filtered_meta_data <- seurat_obj@meta.data
    genes_to_keep = as.numeric(medians) >= cutoff
    message(paste0(100*sum(genes_to_keep)/length(genes_to_keep)), "% of genes kept after filtering out genes with median expression <", cutoff, ". Please adjust cutoff using gene_retention_pct() to find ideal cutoff for your dataset where low expressing genes are not retained. ~10-20% of genes kept is recommended.")
    filtered_raw_mtx <- seurat_obj@assays$RNA@data[genes_to_keep,]
    filtered_meta_data <- filtered_meta_data[colnames(filtered_raw_mtx),]

    used_genes = rownames(filtered_raw_mtx)
    if (is.null(sample_names)){sample_names = unique(as.character(filtered_meta_data[[meta_data_sample_column]]))}
    guide_subsetted_data = lapply(X = sample_names, FUN = mtx_subsetter, matrix = filtered_raw_mtx, meta_data = filtered_meta_data, meta_data_sample_column = meta_data_sample_column)
    names(guide_subsetted_data) = sample_names

    cell_counts = lapply(guide_subsetted_data, ncol)
    min_cells = min(unlist(cell_counts))
    sample_cells_per_guide_cutoff = max(sample_cells_per_guide_cutoff, min_cells)
    control_sample_cells = cell_counts[[control_sample_name]]
    if (sample_cells_per_guide_cutoff > control_sample_cells){
        stop(paste0("The control sample ", control_sample_name, " has only ", control_sample_cells, " cells. Please set sample_cells_per_guide_cutoff to be less than or equal to ", control_sample_cells, "."))
    }
    message(paste0("Sampling ", sample_cells_per_guide_cutoff, " cells per guide. If this is larger than specified, this is because all samples have at least ", min_cells, " cells."))
    ## remove samples with fewer than specified number of cells ##
    guide_subsetted_data = guide_subsetted_data[cell_counts >= sample_cells_per_guide_cutoff]
    if (length(guide_subsetted_data) != length(sample_names)){
        message(paste0(setdiff(sample_names, names(guide_subsetted_data)), "had cell counts below sample_cells_per_guide_cutoff and were filtered out."))
    }
    guide_subsetted_data = lapply(guide_subsetted_data, function(x) x[, sample(ncol(x), sample_cells_per_guide_cutoff)])


    ## add bkg distribution of guide subsetted mtxs to list
    ## this is just mtxs where the cells are randomly assigned to each guide
    ## the bkg distribution of mtxs is made to subtract out number of genes that are expected to have increased or decreased CV when cells assigned randomly to each guide
    ## Note: the number of cells assigned to each guide is kept the same

    randomized_filtered_raw_mtx = Reduce(cbind, guide_subsetted_data[names(guide_subsetted_data) != control_sample_name])
    randomized_cell_order = sample(colnames(randomized_filtered_raw_mtx), ncol(randomized_filtered_raw_mtx), replace = FALSE)
    guides = names(guide_subsetted_data)


    designation_vector = mapply(rep, guides[guides != control_sample_name], sample_cells_per_guide_cutoff, SIMPLIFY = TRUE)
    designation_vector = unlist(designation_vector)
    designation_vector = as.character(designation_vector)
    split_barcodes = split(randomized_cell_order, designation_vector)
    mtx_random_splitter = function(barcodes, mtx){return(mtx[,barcodes])}

    guide_random_subsetted_data = lapply(X = split_barcodes, FUN = mtx_random_splitter, mtx = randomized_filtered_raw_mtx)
    guide_random_subsetted_data = c(guide_subsetted_data[names(guide_subsetted_data) == control_sample_name], guide_random_subsetted_data)
    names(guide_random_subsetted_data) = paste('random', c(control_sample_name, guides[guides != control_sample_name]), sep = '_')
    guide_subsetted_data = c(guide_subsetted_data, guide_random_subsetted_data)


    ##### Calculate CVs #####
    CVs = lapply(X = guide_subsetted_data, FUN = CV_calculator)
    names(CVs) = names(guide_subsetted_data)


    ### creating first 4 columns of the table annotated above
    gene_means = lapply(guide_subsetted_data, sparseMatrixStats::rowMeans2)
    gene_sds = lapply(guide_subsetted_data, sparseMatrixStats::rowSds)
    names(gene_means) = names(guide_subsetted_data)
    names(gene_sds) = names(guide_subsetted_data)


    master_df_list = list()
    for (i in 1:length(guide_subsetted_data)){
        if (!startsWith(names(guide_subsetted_data)[i], 'random')){
            master_df_list[[i]] = data.frame(gene = rownames(filtered_raw_mtx), CV_ctrl = CVs[[control_sample_name]], 
                                            CV_gd = CVs[[i]], CV_gdCV_ctrlratio = CVs[[i]]/CVs[[control_sample_name]], mean_gdmean_ctrlratio = gene_means[[i]]/gene_means[[control_sample_name]])
            names(master_df_list)[i] <- names(guide_subsetted_data)[i]
        }
        if (startsWith(names(guide_subsetted_data)[i], 'random')){
            master_df_list[[i]] = data.frame(gene = rownames(filtered_raw_mtx), CV_ctrl = CVs[[control_sample_name]], 
                                            CV_gd = CVs[[i]], CV_gdCV_ctrlratio = CVs[[i]]/CVs[[control_sample_name]], mean_gdmean_ctrlratio = gene_means[[i]]/gene_means[[control_sample_name]])
            names(master_df_list)[i] <- names(guide_subsetted_data)[i]
        }
    }

    ### creating column 5 of the table annotated above
    for (i in 1:length(master_df_list)){
        master_df_list[[i]]$gene_status = 'NA'
        master_df_list[[i]]$gene_status[master_df_list[[i]]$CV_gdCV_ctrlratio == 1] = 'No Change'
        master_df_list[[i]]$gene_status[master_df_list[[i]]$CV_gdCV_ctrlratio > 1] = 'Increasing'
        master_df_list[[i]]$gene_status[master_df_list[[i]]$CV_gdCV_ctrlratio < 1] = 'Decreasing'               
    }


    ## performing CV equality aysmptotic test and adding its pval to master_df_list
    cells_per_guide = sapply(guide_subsetted_data, ncol)
    #cells_per_guide = rep(cells_per_guide, 2)
    names(cells_per_guide) = names(guide_subsetted_data)
    asymp_test_p_vals = as.data.frame(matrix(0, nrow = nrow(filtered_raw_mtx), ncol = length(master_df_list)))


    for (i in 1:length(master_df_list)){
        if (!startsWith(names(guide_subsetted_data)[i], 'random')){
            for (j in 1:nrow(filtered_raw_mtx)){
                test = cvequality::asymptotic_test2(k = 2, n = c(cells_per_guide[[control_sample_name]], cells_per_guide[[i]]), s = c(gene_sds[[control_sample_name]][j], gene_sds[[i]][j]), 
                                        x = c(gene_means[[control_sample_name]][j], gene_means[[i]][j]))
                asymp_test_p_vals[j,i] = test$p_value
            }
        }
        if (startsWith(names(guide_subsetted_data)[i], 'random')){
            for (j in 1:nrow(filtered_raw_mtx)){
                test = cvequality::asymptotic_test2(k = 2, n = c(cells_per_guide[[control_sample_name]], cells_per_guide[[i]]), s = c(gene_sds[[control_sample_name]][j], gene_sds[[i]][j]), 
                                        x = c(gene_means[[control_sample_name]][j], gene_means[[i]][j]))
                asymp_test_p_vals[j,i] = test$p_value
            }    
        }
        master_df_list[[i]]$p_val = asymp_test_p_vals[,i]
    }

    ##### Find genes that see a significant change in mean from control sample #####
    pbo = pbapply::pboptions(type="txt")
    cells_to_keep_from_fixed_cell_count <- colnames(Reduce(cbind, guide_subsetted_data[guides]))


    Idents(object = CRISPRa_seurat) <- CRISPRa_seurat@meta.data$guide
    CRISPRa_seurat = CRISPRa_seurat[genes_to_keep, cells_to_keep_from_fixed_cell_count]

    find_markers_wrapper = function(perturbed_gene){
                                        Seurat::FindMarkers(CRISPRa_seurat, ident.1 = control_sample_name, ident.2 = perturbed_gene, test.use = 't', verbose = FALSE)
    }

    mean_shifts_from_NT = pbapply::pblapply(guides[guides != control_sample_name], find_markers_wrapper)
    names(mean_shifts_from_NT) = guides[guides != control_sample_name]

    ## Find genes that see signinficant change in mean from control in bkg (randomized cell labels)

    CRISPRa_seurat_random = CRISPRa_seurat
    meta_temp = CRISPRa_seurat_random@meta.data
    meta_temp[randomized_cell_order, 'guide'] = designation_vector
    CRISPRa_seurat_random@meta.data = meta_temp

    find_markers_wrapper_random = function(perturbed_gene){
                                        Seurat::FindMarkers(CRISPRa_seurat_random, ident.1 = control_sample_name, ident.2 = perturbed_gene, test.use = 't', verbose = FALSE)
    }

    mean_shifts_from_NT_bkg = pbapply::pblapply(guides[guides != control_sample_name], find_markers_wrapper_random)
    names(mean_shifts_from_NT_bkg) = paste('random', sep = '_', guides[guides != control_sample_name])

    mean_shifts_from_NT = c(mean_shifts_from_NT, mean_shifts_from_NT_bkg)

    ## run qvalue instead of p.adjust? ##
    asymp_test_q_vals = apply(X = asymp_test_p_vals, MARGIN = 2, FUN = function(x) qvalue::qvalue(x)$qvalues)
    colnames(asymp_test_q_vals) = names(master_df_list)
    num_signif_genes = colSums(asymp_test_q_vals < 0.05, na.rm = TRUE)
    num_signif_genes_df = data.frame(guide = names(master_df_list), num_signif_genes = num_signif_genes)




    ## pvalues for CV and mean ratios for plotting ##

    # CV #
    plotted_guides = names(master_df_list)[!(names(master_df_list) %in% c(control_sample_name, paste0('random_', control_sample_name)))]
    plotted_guides_with_NT = names(master_df_list)
    ratios = lapply(names(master_df_list), function(x){log2(master_df_list[[x]]$CV_gdCV_ctrlratio)})
    names(ratios) = plotted_guides_with_NT
    ratios_df <- as.data.frame(do.call(cbind, ratios))
    colnames(ratios_df) = names(master_df_list)
    ratios_df = ratios_df[, !(colnames(ratios_df) %in% c(control_sample_name, paste0('random_', control_sample_name)))]
    ratios_df = tidyr::pivot_longer(ratios_df, cols = colnames(ratios_df), values_to = 'CV_gdCV_ctrlratio', names_to = 'guide')
    order_of_guides = ratios_df %>% group_by(guide) %>% summarize(means = mean(CV_gdCV_ctrlratio)) %>% arrange(means) %>% pull(guide) %>% as.character()
    ratios_df$guide = factor(ratios_df$guide, levels = order_of_guides)
    
    ratios = ratios[order_of_guides]
    ratios = lapply(ratios, function(x) x[!(is.infinite(x) | is.nan(x) | is.na(x))])
    tests = lapply(ratios[!endsWith(names(ratios), control_sample_name)], t.test, mu = 0) ## could cause problem (endsWith) in plotting later
    pvals = as.numeric(lapply(tests, function(x) return(x$p.val)))
    pvals_adj = signif(p.adjust(pvals), 4)
    pvals_adj[pvals_adj == 0] = '< e-300'
    names(pvals_adj) = order_of_guides
    CV_pvals_adj = pvals_adj
    CV_order_of_guides = order_of_guides
    CV_ratios_df = ratios_df

    # mean #
    plotted_guides = names(master_df_list)[!(names(master_df_list) %in% c(control_sample_name, paste0('random_', control_sample_name)))]
    plotted_guides_with_NT = names(master_df_list)
    ratios = lapply(names(master_df_list), function(x){log2(master_df_list[[x]]$mean_gdmean_ctrlratio)})
    names(ratios) = plotted_guides_with_NT
    ratios_df <- as.data.frame(do.call(cbind, ratios))
    colnames(ratios_df) = names(master_df_list)
    ratios_df = ratios_df[, !(colnames(ratios_df) %in% c(control_sample_name, paste0('random_', control_sample_name)))]
    ratios_df = tidyr::pivot_longer(ratios_df, cols = colnames(ratios_df), values_to = 'mean_gdmean_ctrlratio', names_to = 'guide')
    ratios_df$guide = factor(ratios_df$guide, levels = order_of_guides)

    ratios = ratios[order_of_guides]
    ratios = lapply(ratios, function(x) x[!(is.infinite(x) | is.nan(x) | is.na(x))])
    tests = lapply(ratios[!endsWith(names(ratios), control_sample_name)], t.test, mu = 0) ## could cause problem (endsWith) in plotting later
    pvals = as.numeric(lapply(tests, function(x) return(x$p.val)))
    pvals_adj = signif(p.adjust(pvals), 4)
    pvals_adj[pvals_adj == 0] = '< e-300'
    names(pvals_adj) = order_of_guides
    mean_pvals_adj = pvals_adj
    mean_order_of_guides = order_of_guides
    mean_ratios_df = ratios_df



    return_list = list(guide_subsetted_data = guide_subsetted_data, master_df_list = master_df_list, mean_shifts_from_NT = mean_shifts_from_NT, 
    asymp_test_p_vals = asymp_test_p_vals, order_of_guides = order_of_guides, significant_CV_gene_count = num_signif_genes_df, 
    control_sample_name = control_sample_name, CV_pvals_adj = CV_pvals_adj, CV_order_of_guides = CV_order_of_guides, CV_ratios_df = CV_ratios_df,
    mean_pvals_adj = mean_pvals_adj, mean_order_of_guides = mean_order_of_guides, mean_ratios_df = mean_ratios_df, sample_cells_per_guide_cutoff = sample_cells_per_guide_cutoff, cc_meta_data = seurat_obj_with_cc_scores@meta.data)
    
    return(return_list)

}