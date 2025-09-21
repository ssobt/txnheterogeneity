## Input bulk RNA-seq data
## Used to compare heterogeneity between patient groups or cell line groups (or conditions) in RNA seq data
# packages necessary: cvequality, qvalue, Seurat, ggplot2, parallel, ggforce, ggpubr, ggplot
##### 1. unit test bulk_het() PCA part for first half, 2 and 3. unit test latter half of bulk_het()
##### 4. for background distribution in bulk_het(), 
##### think about allowing user to specify columns in data that will serve as bkg distribution instead of randomly picking


##### Helper functions #####
## mtx_subsetter: Subsets the matrix into high and low expression groups based on the stratifying gene and quantile.

#' @noRd
mtx_subsetter = function(stratifying_gene, mtx, quant){
            high_cutoff_number = quantile(mtx[stratifying_gene,], quant)
            low_cutoff_number = quantile(mtx[stratifying_gene,], 1 - quant)
            group1_mtx = mtx[, mtx[stratifying_gene,] >= high_cutoff_number] ## higher expressing group
            group2_mtx = mtx[, mtx[stratifying_gene,] <= low_cutoff_number] ## lower expressing group
            return(list(high_exp_group = group1_mtx, low_exp_group = group2_mtx, high_cutoff_number = high_cutoff_number, low_cutoff_number = low_cutoff_number))
}

#' @noRd
CV_bulk_calculator = function(values){return(sd(values)/mean(values))}

#' @noRd
# two-tailed pvalue for z-score
z_pvalue <- function(z) {
    probability = pnorm(abs(z))
    if (probability < 0.5) {
        return(2 * probability)
    } else {
        return(2 * (1 - probability))
    }
}


##### Heterogeneity functions #####

## bulk_het() ##

#' bulk_het()
#' @description Calculate transcriptional heterogeneity in bulk RNA-seq data samples
#'
#' @param data A matrix containing FPKM data with rows as genes and columns as samples.
#' @param g1 A numeric vector with column names in data.
#' @param g2 A numeric vector with column names in data.
#' @param quant If g1 and g2 are NULL, a number between 0.5 and 1 indicating the quantile to use for stratifying samples based on expression of stratifiers.
#' @param stratifiers A character vector of gene names to use for stratification.
#' @param cores Number of cores to use for parallel processing. Default is 1.
#'
#' It is advised that you save the output of this function as an RDS file using saveRDS() immediately after running it for easy loading in future sessions and because memory failure can occur when running downstream plotting functions on large datasets.
#' @returns A list with the following components:
#' \itemize{
#'  \item{reference_df}{A data frame with gene names, CVs, and means for each group.}
#'  \item{genes_with_significant_CV_change}{A data frame with counts of genes with significant CV changes between groups.}
#'  \item{z_score}{Z-score comparing observed significant CV changes to background distribution.}
#'  \item{z_score_pvalue}{P-value associated with the z-score.}
#'  \item{PCA_coordinates_g1}{PCA coordinates for group 1 samples.}
#'  \item{PCA_coordinates_g2}{PCA coordinates for group 2 samples.}
#'  \item{distance_data_g1}{Weighted distances to centroid for group 1 samples.}
#'  \item{distance_data_g2}{Weighted distances to centroid for group 2 samples.}
#' }
#' OR
#' \itemize{
#'  \item{reference_df}{A data frame with gene names, CVs, and means for each group.}
#'  \item{z_score}{Z-score comparing observed significant CV changes to background distribution and assoiciated p-value.}
#'  \item{upper_quant_PCA_coordinates}{PCA coordinates for upper quantile samples for each provided stratified gene.}
#'  \item{lower_quant_PCA_coordinates}{PCA coordinates for lower quantile samples for each provided stratified gene.}
#'  \item{upper_quant_distance_data}{Weighted distances to centroid for upper quantile samples for each provided stratified gene.}
#'  \item{lower_quant_distance_data}{Weighted distances to centroid for lower quantile samples for each provided stratified gene.}
#' }
#'
#' @examples
#' ## Example 1: Compare heterogeneity between two known groups of samples
#' ## input_mtx is a matrix with rows as genes and columns as samples
#' ## g1 and g2 are character vectors with column names in input_mtx
#' # set seed for reproducibility of bkg distribution sampling
#' set.seed(123)
#' samples = sample(colnames(input_mtx), 200)
#' g1 = samples[1:100]
#' g2 = samples[101:200]
#' out = bulk_het(data = input_mtx, g1 = g1, g2 = g2, cores = 1)
#'
#'
#' ## Example 2: Compare heterogeneity between samples stratified by expression of certain genes
#' ## input_mtx is a matrix with rows as genes and columns as samples
#' # set seed for reproducibility of bkg distribution sampling
#' set.seed(123)
#' ## will stratify samples based on top and bottom 25% expression of each stratifier gene
#' ## top 25% (q75) will be compared to bottom 25% (q25)
#' out = bulk_het(data = input_mtx, quant = 0.75, stratifiers = rownames(input_mtx)[1:10], cores = 10) 
#' @export

bulk_het <- function(data, g1 = NULL, g2 = NULL, quant = NULL, stratifiers = NULL, cores = 1) {
    message('Use set.seed() to ensure reproducibility')
    ## data must be in matrix format with column names resembling sample/patient names
    if (is.matrix(data) & is.character(g1) & is.character(g2) & is.null(quant) & is.null(stratifiers)) {
        if (!all(g1 %in% colnames(data)) | !all(g2 %in% colnames(data))) {
            stop("g1 and g2 must be column names of the input matrix")
        }
        
        ##### CV analysis #####
        data = data[, c(g1, g2)]
        group1_mtx = data[, g1]
        group2_mtx = data[, g2]
        group1_exp_CV = apply(group1_mtx, 1, CV_bulk_calculator)
        group2_exp_CV = apply(group2_mtx, 1, CV_bulk_calculator)
        group1_exp_mean = apply(group1_mtx, 1, mean)
        group2_exp_mean = apply(group2_mtx, 1, mean)

        ## bkg sampling ##
        ## initially do 10, fix later with a more official answer based on a distribution or certain % of total genes in rows that could be stratified on or % of number of samples in total dataset given
        times_to_repeat_sampling = 10
        background_distribution = lapply(1:times_to_repeat_sampling, function(x) {
            g1_bkg = data[, sample(c(g1, g2), size = length(g1), replace = F)]
            g2_bkg = data[, sample(setdiff(c(g1, g2), colnames(g1_bkg)), size = length(g2), replace = F)]
            return(list(group1_bkg_mtx = g1_bkg, group2_bkg_mtx = g2_bkg))
        })       

        reference_df = data.frame(
            gene = rownames(data),
            group1_exp_CV = group1_exp_CV,
            group2_exp_CV = group2_exp_CV,
            group1_exp_mean = group1_exp_mean,
            group2_exp_mean = group2_exp_mean
        )

        ## perform CV equality test between sample/patient groups for each gene
        mtx = cbind(group1_mtx, group2_mtx)
        pvals = apply(mtx, 1, function(k) {cvequality::asymptotic_test(x = k, y = c(rep('g1', ncol(group1_mtx)), rep('g2', ncol(group2_mtx))))$p_value})
        pvals_bkg = parallel::mclapply(background_distribution, function(bkg) {
            mtx = cbind(bkg$group1_bkg_mtx, bkg$group2_bkg_mtx)
            apply(mtx, 1, function(k) {cvequality::asymptotic_test(x = k, y = c(rep('g1', ncol(bkg$group1_bkg_mtx)), rep('g2', ncol(bkg$group2_bkg_mtx))))$p_value})
        }, mc.cores = cores)

        pvals = ifelse(is.na(pvals), 1, pvals) ## convert NA pvals to 1 for qvalue calculation
        pvals_bkg = lapply(pvals_bkg, function(x) ifelse(is.na(x), 1, x)) ## convert NA pvals to 1 for qvalue calculation
        qvals = qvalue::qvalue(pvals, fdr.level = .05)
        qvals_bkg = lapply(pvals_bkg, function(p) qvalue::qvalue(p, fdr.level = .05))

        ## count number of genes with significantly different CVs between groups from CV equality test for each stratifier
        sig_gene_counts = sum(qvals$qvalues < 0.05 & !is.nan(qvals$qvalues))
        random_sig_counts = sapply(qvals_bkg, function(q) sum(q$qvalues < 0.05 & !is.nan(q$qvalues)))

        ## calculate z-score of stratifiers using randomly chosen stratifiers as background distribution
        z_score = (sig_gene_counts - mean(random_sig_counts))/sd(random_sig_counts)

        z_score_pval = z_pvalue(z_score)

        ##### PCA analysis #####
        
        ### This cell runs the PCA analysis on only the top 2000 DEGs ###
        ### edited out non-DEG genes here ###
        seurat_input_data = as.data.frame(cbind(genes = rownames(data), data))
        seurat_input_data = dplyr::distinct(seurat_input_data, genes, .keep_all = TRUE)
        seurat_input_data = as.matrix(seurat_input_data[, !colnames(seurat_input_data) %in% 'genes'])
        rna_seq_seurat <- Seurat::CreateSeuratObject(counts = seurat_input_data)
        rna_seq_seurat <- Seurat::NormalizeData(rna_seq_seurat)
        rna_seq_seurat <- Seurat::FindVariableFeatures(rna_seq_seurat)
        head(Seurat::VariableFeatures(rna_seq_seurat), 2000) -> var.genes
        var.genes.idx <- match(var.genes, row.names(rna_seq_seurat@assays$RNA))
        raw_counts_2000DEGs <- as.data.frame(rna_seq_seurat@assays$RNA@data[var.genes.idx,])

        ### note for PCA here we are using the SVD algorithm
        # The algorithm returns the same number of PCs as there are observations
        ## This is because SVD does not use the covariance matrix to compute PCs

        rownames(raw_counts_2000DEGs) <- row.names(rna_seq_seurat@assays$RNA)[var.genes.idx]
        pca1 = prcomp(t(raw_counts_2000DEGs), scale. = TRUE)
        PCA_coordinates = pca1$x
        PCA_coordinates_final = PCA_coordinates[, 1:50]
        PCA_coordinates_final = t(PCA_coordinates_final)
        eigs <- pca1$sdev^2
        eigen_values = eigs
        explained_variance = eigs/sum(eigs)
        explained_variance = explained_variance[1:50]

        # group PCA_coordinates_final based on either manually given g1 and g2 or quantile groupings
        PCA_coordinates_final = as.data.frame(PCA_coordinates_final)
        g1_PCA_coordinates = PCA_coordinates_final[, which(colnames(PCA_coordinates_final) %in% g1)]
        g2_PCA_coordinates = PCA_coordinates_final[, which(colnames(PCA_coordinates_final) %in% g2)]

        # transpose
        g1_PCA_coordinates = t(g1_PCA_coordinates)
        g2_PCA_coordinates = t(g2_PCA_coordinates)

        # find mean
        g1_PCA_coordinates_w_centroid = rbind(colMeans(g1_PCA_coordinates), g1_PCA_coordinates)
        g2_PCA_coordinates_w_centroid = rbind(colMeans(g2_PCA_coordinates), g2_PCA_coordinates)

        # find weighted distance to centroid
        g1_sample_weighted_distances_to_centroid = apply(g1_PCA_coordinates_w_centroid[-1,], 1, function(z, centroid, weights) {sqrt(sum(weights*(z - centroid)^2))}, g1_PCA_coordinates_w_centroid[1,], explained_variance)
        g2_sample_weighted_distances_to_centroid = apply(g2_PCA_coordinates_w_centroid[-1,], 1, function(z, centroid, weights) {sqrt(sum(weights*(z - centroid)^2))}, g2_PCA_coordinates_w_centroid[1,], explained_variance)

        # mean distance from centroid
        g1_mean_distance_to_centroid = mean(g1_sample_weighted_distances_to_centroid)
        g2_mean_distance_to_centroid = mean(g2_sample_weighted_distances_to_centroid)

        # mean distance ratio
        g1_g2_mean_distance_ratio = g1_mean_distance_to_centroid / g2_mean_distance_to_centroid

        message('If gene counts with significant changes in CV between groups is below randomized grouping, \nthere is no significant difference in heterogeneity between the groups')

        return(list(reference_df = reference_df, 
        genes_with_significant_CV_change = data.frame(genes_with_significant_CV_change = c(sig_gene_counts, mean(random_sig_counts)), type = c('given_sample_grouping', 'randomized_sample_grouping')), 
        z_score = z_score, z_score_pvalue = z_score_pval, PCA_coordinates_g1 = g1_PCA_coordinates, PCA_coordinates_g2 = g2_PCA_coordinates, distance_data_g1 = g1_sample_weighted_distances_to_centroid, distance_data_g2 = g2_sample_weighted_distances_to_centroid))

    }

    if (is.matrix(data) & is.null(g1) & is.null(g2) & is.numeric(quant) & is.character(stratifiers)){
        
        ##### CV analysis #####
        if (!(quant > 0.5 & quant < 1)){stop("quant must be a quantile between 0.5 and 1")}
        ## divide expression mtx into a list of mtxs subsetted by guide
        mtx = data
        subsetted_data = lapply(X = stratifiers, FUN = mtx_subsetter, 
                                        mtx = mtx, quant = quant)
        names(subsetted_data) = original_stratifiers = stratifiers

        ## remove any stratifiers that have 0 as lower quantile cutoff
        subsetted_data = subsetted_data[c(!sapply(subsetted_data, function(x) x$low_cutoff_number == 0))]
        stratifiers = modified_stratifiers = names(subsetted_data)

        ## add a random set of bkg distribtution on random genes to compare with (10% of total subset or minimum of 5)
        random_genes = sample(setdiff(rownames(mtx), stratifiers), max(ceiling(0.1*length(stratifiers)), 5))
        random_subsetted_data = lapply(X = random_genes, FUN = mtx_subsetter, mtx = mtx, quant = quant)
        names(random_subsetted_data) = paste0("random_", random_genes)
        ## remove any stratifiers that have 0 as lower quantile cutoff
        random_subsetted_data = random_subsetted_data[c(!sapply(random_subsetted_data, function(x) x$low_cutoff_number == 0))]

        ## combine real stratifiers and random stratifiers
        subsetted_data = c(subsetted_data, random_subsetted_data)

        if (all(original_stratifiers == modified_stratifiers)){
            message(paste0('All stratifiers were able to be used for analysis'))
        } else {
            message(paste0('The following stratifiers were removed from analysis because their lower quantile cutoff was 0: ', paste(setdiff(original_stratifiers, modified_stratifiers), collapse = ', ')))
        }

        reference_dfs = list()
        ## calculate CVs and means between sample/patient groups for reference
        for (i in 1:length(subsetted_data)){
            high_exp_CV = apply(subsetted_data[[i]][['high_exp_group']], 1, CV_bulk_calculator)
            low_exp_CV = apply(subsetted_data[[i]][['low_exp_group']], 1, CV_bulk_calculator)
            high_exp_mean = apply(subsetted_data[[i]][['high_exp_group']], 1, mean)
            low_exp_mean = apply(subsetted_data[[i]][['low_exp_group']], 1, mean)

            reference_dfs[[i]] = data.frame(
                gene = rownames(subsetted_data[[i]][['high_exp_group']]),
                high_exp_CV = high_exp_CV,
                low_exp_CV = low_exp_CV,
                high_exp_mean = high_exp_mean,
                low_exp_mean = low_exp_mean
            )
        }
        names(reference_dfs) = names(subsetted_data)

        ## perform CV equality test between sample/patient groups for each gene
        pvals = parallel::mclapply(subsetted_data, function(z){
                    apply(X = cbind(z[['high_exp_group']], z[['low_exp_group']]), MARGIN = 1, FUN = function(k){
                    cvequality::asymptotic_test(x = k, y = c(rep('high', length(k)/2), rep('low', length(k)/2)))$p_value
                })
        }, mc.cores = cores)
        names(pvals) = names(subsetted_data)
        
        ## adjust pvals to qvals
        pvals = lapply(pvals, function(x) ifelse(is.na(x), 1, x)) ## convert NA pvals to 1 for qvalue calculation
        qvals = lapply(X = pvals, FUN = qvalue::qvalue, fdr.level = .05)

        ## add pvals and qvals to reference_dfs
        for (i in 1:length(reference_dfs)){
            reference_dfs[[i]]$pval = pvals[[i]]
            reference_dfs[[i]]$qval = qvals[[i]]$qvalues
        }

        ## count number of genes with significantly different CVs between groups from CV equality test for each stratifier
        sig_gene_counts = lapply(reference_dfs, function(df) sum(df$qval < 0.05 & !is.nan(df$qval), na.rm = TRUE))
        names(sig_gene_counts) = names(reference_dfs)
        sig_gene_counts = data.frame(stratifier_genes = names(sig_gene_counts), genes_with_signigicant_CV_change = unlist(sig_gene_counts))

        ## calculate z-score of stratifiers using randomly chosen stratifiers as background distribution
        random_sig_counts = sig_gene_counts[startsWith(as.character(sig_gene_counts$stratifier_genes), 'random_'), ]
        stratifier_sig_counts = sig_gene_counts[!startsWith(as.character(sig_gene_counts$stratifier_genes), 'random_'), ]
        mean_random = mean(random_sig_counts$genes_with_signigicant_CV_change)
        sd_random = sd(random_sig_counts$genes_with_signigicant_CV_change)
        z_scores = (stratifier_sig_counts$genes_with_signigicant_CV_change - mean_random) / sd_random

        ##### PCA analysis #####
        
        ### This cell runs the PCA analysis on only the top 2000 DEGs ###
        ### edited out non-DEG genes here ###
        seurat_input_data = as.data.frame(cbind(genes = rownames(data), data))
        seurat_input_data = dplyr::distinct(seurat_input_data, genes, .keep_all = TRUE)
        seurat_input_data = as.matrix(seurat_input_data[, !colnames(seurat_input_data) %in% 'genes'])
        rna_seq_seurat <- Seurat::CreateSeuratObject(counts = seurat_input_data)
        rna_seq_seurat <- Seurat::NormalizeData(rna_seq_seurat)
        rna_seq_seurat <- Seurat::FindVariableFeatures(rna_seq_seurat)
        head(Seurat::VariableFeatures(rna_seq_seurat), 2000) -> var.genes
        var.genes.idx <- match(var.genes, row.names(rna_seq_seurat@assays$RNA))
        raw_counts_2000DEGs <- as.data.frame(rna_seq_seurat@assays$RNA@data[var.genes.idx,])

        ### note for PCA here we are using the SVD algorithm
        # The algorithm returns the same number of PCs as there are observations
        ## This is because SVD does not use the covariance matrix to compute PCs

        rownames(raw_counts_2000DEGs) <- row.names(rna_seq_seurat@assays$RNA)[var.genes.idx]
        pca1 = prcomp(t(raw_counts_2000DEGs), scale. = TRUE)
        PCA_coordinates = pca1$x
        PCA_coordinates_final = PCA_coordinates[, 1:50]
        PCA_coordinates_final = t(PCA_coordinates_final)
        eigs <- pca1$sdev^2
        eigen_values = eigs
        explained_variance = eigs/sum(eigs)
        explained_variance = explained_variance[1:50]

        # group PCA_coordinates_final based on either manually given g1 and g2 or quantile groupings
        PCA_coordinates_final = as.data.frame(PCA_coordinates_final)
        upper_quant_PCA_coordinates = lapply(subsetted_data, function(z) PCA_coordinates_final[, which(colnames(PCA_coordinates_final) %in% colnames(z$high_exp_group))])
        lower_quant_PCA_coordinates = lapply(subsetted_data, function(z) PCA_coordinates_final[, which(colnames(PCA_coordinates_final) %in% colnames(z$low_exp_group))])

        # transpose
        upper_quant_PCA_coordinates = lapply(upper_quant_PCA_coordinates, t)
        lower_quant_PCA_coordinates = lapply(lower_quant_PCA_coordinates, t)
        names(upper_quant_PCA_coordinates) = names(subsetted_data)
        names(lower_quant_PCA_coordinates) = names(subsetted_data)

        # find mean
        upper_quant_PCA_coordinates_w_centroid = lapply(upper_quant_PCA_coordinates, function(x) rbind(colMeans(x), x))
        lower_quant_PCA_coordinates_w_centroid = lapply(lower_quant_PCA_coordinates, function(x) rbind(colMeans(x), x))

        # find weighted distance to centroid
        upper_quant_sample_weighted_distances_to_centroid = lapply(upper_quant_PCA_coordinates_w_centroid, function(x) apply(x[-1,], 1, function(z, centroid, weights) {sqrt(sum(weights*(z - centroid)^2))}, x[1,], explained_variance))
        lower_quant_sample_weighted_distances_to_centroid = lapply(lower_quant_PCA_coordinates_w_centroid, function(x) apply(x[-1,], 1, function(z, centroid, weights) {sqrt(sum(weights*(z - centroid)^2))}, x[1,], explained_variance))

        # mean distance from centroid
        upper_quant_mean_distance_to_centroid = lapply(upper_quant_sample_weighted_distances_to_centroid, mean)
        lower_quant_mean_distance_to_centroid = lapply(lower_quant_sample_weighted_distances_to_centroid, mean)

        # mean distance ratio
        upper_lower_mean_distance_ratio = mapply(function(upper, lower) {upper / lower}, upper_quant_mean_distance_to_centroid, lower_quant_mean_distance_to_centroid)

        

        return(list(
            reference_dfs = reference_dfs,
            sig_gene_counts,
            z_scores = data.frame(stratifier_genes = stratifier_sig_counts$stratifier_genes, z_score = z_scores, pvalue = sapply(z_scores, z_pvalue)),
            upper_quant_PCA_coordinates = upper_quant_PCA_coordinates,
            lower_quant_PCA_coordinates = lower_quant_PCA_coordinates,
            upper_quant_distance_data = upper_quant_sample_weighted_distances_to_centroid,
            lower_quant_distance_data = lower_quant_sample_weighted_distances_to_centroid,
            quant = quant
        ))


    }

}




