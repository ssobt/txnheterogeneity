

# packages necessary: ArchR, BSgenome.Hsapiens.UCSC.hg38, plyr, Matrix, cvequality, sparseMatrixStats, dplyr, tidyr


#' @noRd
mtx_subsetter_scATAC = function(guide_nm, matrix, meta_data){return(matrix[, meta_data$guide == guide_nm])}


##### Heterogeneity functions #####

## scATAC_het() ##

#' scATAC_het()
#' @description Calculate transcriptional heterogeneity in scATAC-seq samples compared to a control sample.
#' @import ArchR
#' @import Matrix
#' @import BSgenome.Hsapiens.UCSC.hg38
#' @param scATAC_fragment_files A character vector with file paths to files containing quantified scATAC-seq fragment files (files must end in sample_name.fragments.tsv.gz) using 10x Genomics Cell Ranger ATAC pipeline or other software.
#' @param control_sample_name The name of the control sample.
#' @param sample_names A character vector with names of the samples corresponding to the fragment files provided.
#' @param savepath A file path to a directory where intermediate files and results will be saved.
#' @param threads Number of threads to use for ArchR processing. Default is 10.
#' @param genome The genome to use for ArchR processing. Default is "hg38".
#' @param TSS_bed_path A file path to a bed file containing TSS regions for the genome being used. This is used to create a TSS peak matrix for heterogeneity calculations. Please refer to documentation for ArchR::addPeakSet() for more details.
#' @param macs2_path A file path to the MACS2 executable in your current environment. This is used by ArchR to call peaks.
#' @param fixed_cell_count An integer specifying a fixed number of cells to use per sample for heterogeneity calculations. If NULL, the minimum number of cells across samples will be used. Default is NULL.
#' @param replicate_map A list of length 2 containing vectors for renaming samples if replicates are present. The first vector should contain the original sample names, and the second vector should contain the new names. If NULL, no renaming will be done. Default is NULL.
#' @param genes_to_use A character vector of gene names to use for heterogeneity calculations. If NULL, all genes will be used. Default is NULL.
#' @param samples_to_analyze A character vector of sample names to include in the heterogeneity calculations. If NULL, all samples will be used. Default is NULL.
#' @param seed An integer specifying the seed for random number generation. Default is 22.
#' @examples
#' out = scATAC_het(scATAC_fragment_files = paste0(folder, '/', file_names), 
#' sample_names = names(file_names), 
#' replicate_map = list(c('NTCi-1', 'NTCi-2', 'NTCa-1', 'NTCa-2', 'RNF8-Ci-1', 'RNF8-Ci-2', 'RNF8-Ca-1', 'RNF8-Ca-2', 'MIS18A-Ci-1', 'MIS18A-Ci-2', 'MIS18A-Ca-1', 'MIS18A-Ca-2'), c('NTCi', 'NTCi', 'NTCa', 'NTCa', 'RNF8-Ci', 'RNF8-Ci', 'RNF8-Ca', 'RNF8-Ca', 'MIS18A-Ci', 'MIS18A-Ci', 'MIS18A-Ca', 'MIS18A-Ca')), ## if replicates present, this replicate map must be done for ALL samples inputted into the function
#' samples_to_analyze = c('NTCi', 'RNF8-Ci', 'MIS18A-Ci'), control_sample_name = 'NTCi',
#' savepath = savep, 
#' threads = 35, 
#' genome = "hg38", 
#' TSS_bed_path = '/home/ssobti/projects/heterogeneity_brian/uploaded_data/scATACseq/genome/GRCh38_transcriptsOnly.tss.bed',
#' macs2_path = '/home/ssobti/miniconda3/envs/archr/bin/macs2') ## please make sure to provide this if using, ArchR gets confused if multiple copies exist
#' @return 
#' A list containing:
#' \item{filtered_gene_scores_mtx}{A filtered gene score matrix used for heterogeneity calculations.}
#' \item{filtered_gene_score_meta_data}{The metadata associated with the filtered gene score matrix.}
#' \item{control_sample_name}{The name of the control sample.}
#' \item{CV_ratios_df}{A dataframe containing log2(CV_guide/CV_control) values for each gene across all guides.}
#' \item{CV_order_of_guides}{The order of guides for CV ratio plotting.}
#' \item{CV_pvals_adj}{The adjusted p-values from t-tests of log2(CV_guide/CV_control) values for each guide.}
#' \item{mean_ratios_df}{A dataframe containing log2(mean_guide/mean_control) values for each gene across all guides.}
#' \item{mean_order_of_guides}{The order of guides for mean ratio plotting.}
#' \item{mean_pvals_adj}{The adjusted p-values from t-tests of log2(mean_guide/mean_control) values for each guide.}
#' \item{CVs_df}{A dataframe containing log2(CV_guide) values for each gene across all guides.}
#' \item{means_df}{A dataframe containing log2(mean_guide) values for each gene across all guides.}
#' @details
#' This function takes in fragment files from scATAC-seq data and calculates transcriptional heterogeneity for each gene in each sample 
#' compared to a control sample. The function uses the ArchR package to process the scATAC-seq data and obtain gene score matrices. 
#' It then calculates the coefficient of variation (CV) for each gene in each sample and compares it to the CV in the control sample. 
#' The function also calculates the mean expression for each gene in each sample and compares it to the mean in the control sample. 
#' The results are returned in a list containing dataframes for CV ratios, mean ratios, and their associated p-values.
#' @export

# input is fragment files from scATAC-seq data
scATAC_het <- function(scATAC_fragment_files = NULL, sample_names = NULL, control_sample_name = NULL, savepath = NULL, threads = 10, genome = "hg38", TSS_bed_path = NULL, macs2_path = NULL, fixed_cell_count = NULL, replicate_map = NULL, genes_to_use = NULL, samples_to_analyze = NULL, seed = 22) {

    ###### ArchR matrix extraction ###### ------------------------------------------------------------------
    message("BSgenome.Hsapiens.UCSC.hg38 and ArchR are required to be attached for scATAC_het(). Please load both using library() before running this function.")
    if (!("BSgenome.Hsapiens.UCSC.hg38" %in% .packages()) & ("ArchR" %in% .packages())) {
        stop("BSgenome.Hsapiens.UCSC.hg38 or ArchR is not loaded. Please load both using library() before running this function.")
    }
    
    set.seed(seed)
    setwd(savepath)
    addArchRThreads(threads = threads)
    addArchRGenome(genome)
    addArchRLocking(locking = TRUE)
    if (is.null(scATAC_fragment_files)) {
        stop("Please provide a vector of fragment files.")
    }
    if (is.null(sample_names)) {
        stop("Please provide a vector of sample names.")
    }
    if (length(scATAC_fragment_files) != length(sample_names)) {
        stop("Please provide a vector of sample names of the same length as the vector of fragment files.")
    }
    if (is.null(control_sample_name)) {
        stop("Please provide a control sample name.")
    }
    if (is.null(replicate_map) & !(control_sample_name %in% c(sample_names))) {
        stop("Control sample name not found in sample names.")
    }
    if (!is.null(replicate_map) & !(control_sample_name %in% c(replicate_map[[2]]))) {
        stop("Control sample name not found in renamed replicate map.")
    }
    if (is.null(savepath)) {                
        stop("Please provide a savepath.")
    }
    if (is.null(TSS_bed_path)) {
        stop("Please provide a TSS bed file path.")
    }
    if (is.null(macs2_path)) {
        stop('Please provide the path to MACS2 in your current environment. If you are using a conda environment for ArchR, you can activate it and run `which macs2` to find the path.')
    }
    if (is.null(fixed_cell_count)){
        message('No fixed_cell_count provided. The minimum number of cells across samples will be used for heterogeneity calculations.')
    }

    names(scATAC_fragment_files) = sample_names
    ArrowFiles <- createArrowFiles(
    inputFiles = scATAC_fragment_files,
    sampleNames = names(scATAC_fragment_files),
    minTSS = 4, #Dont set this too high because you can always increase later
    minFrags = 1000, 
    addTileMat = TRUE,
    addGeneScoreMat = TRUE
    )

    doubScores <- addDoubletScores(
    input = ArrowFiles,
    k = 10, #Refers to how many cells near a "pseudo-doublet" to count.
    knnMethod = "UMAP", #Refers to the embedding to use for nearest neighbor search.
    LSIMethod = 1
    )

    proj <- ArchRProject(
    ArrowFiles = ArrowFiles, 
    outputDirectory = savepath,
    copyArrows = TRUE #This is recommened so that you maintain an unaltered copy for later usage.
    )

    proj <- filterDoublets(ArchRProj = proj)
    proj <- saveArchRProject(ArchRProj = proj)

    ### Obtaining matrices ###
    proj <- addIterativeLSI(ArchRProj = proj, useMatrix = "TileMatrix", name = "IterativeLSI")
    proj <- addClusters(input = proj, reducedDims = "IterativeLSI")
    proj <- addGroupCoverages(proj, force = F, groupBy = 'Sample', threads = 1)

    proj <- addReproduciblePeakSet(proj, peakMethod = "macs2", groupBy = 'Sample', pathToMacs2 = macs2_path)
    #head(as.data.frame(proj@cellColData))
    txnhet_archr_macs2_project <- addPeakMatrix(proj, threads = 1)
    txnhet_archr_macs2_project <- saveArchRProject(ArchRProj = txnhet_archr_macs2_project)

    ## Global peak matrix ##
    # currently the MACS2 peaks for ALL regions loaded in the peak matrix
    txnhet_archr_macs2_peaks = getMatrixFromProject(txnhet_archr_macs2_project, useMatrix = 'PeakMatrix')
    txnhet_archr_macs2_peak_mtx = assays(txnhet_archr_macs2_peaks)$PeakMatrix
    txnhet_archr_macs2_peak_meta = as.data.frame(colData(txnhet_archr_macs2_peaks))
    saveRDS(txnhet_archr_macs2_peak_mtx, file = paste0(savepath,'/txnhet_archr_macs2_peak_mtx.rds'))
    saveRDS(txnhet_archr_macs2_peak_meta, file = paste0(savepath,'/txnhet_archr_macs2_peak_meta.rds'))

    TSS_granges_obj = genomation::readBed(TSS_bed_path, track.line = FALSE, remove.unusual = FALSE, zero.based = TRUE)
    seqlevelsStyle(TSS_granges_obj) <- "UCSC"

    ## TSS peak matrix ##
    txnhet_archr_macs2_TSS_project <- addPeakSet(ArchRProj = txnhet_archr_macs2_project, peakSet = TSS_granges_obj, force = TRUE)
    txnhet_archr_macs2_TSS_project <- addPeakMatrix(txnhet_archr_macs2_TSS_project)
    txnhet_archr_macs2_TSS_peaks = getMatrixFromProject(txnhet_archr_macs2_TSS_project, useMatrix = 'PeakMatrix')
    txnhet_archr_macs2_peak_TSS_mtx = assays(txnhet_archr_macs2_TSS_peaks)$PeakMatrix
    txnhet_archr_macs2_peak_TSS_meta = as.data.frame(colData(txnhet_archr_macs2_TSS_peaks))
    rownames(txnhet_archr_macs2_peak_TSS_mtx) = TSS_granges_obj$name
    saveRDS(txnhet_archr_macs2_peak_TSS_mtx, file = paste0(savepath,'/txnhet_archr_macs2_peak_TSS_mtx.rds'))
    saveRDS(txnhet_archr_macs2_peak_TSS_meta, file = paste0(savepath,'/txnhet_archr_macs2_peak_TSS_meta.rds'))

    ## Gene score matrix ##
    txnhet_archr_macs2_gene_scores = getMatrixFromProject(txnhet_archr_macs2_project, useMatrix = 'GeneScoreMatrix')
    txnhet_archr_macs2_gene_scores_mtx = assays(txnhet_archr_macs2_gene_scores)$GeneScoreMatrix
    txnhet_archr_macs2_gene_scores_meta = as.data.frame(colData(txnhet_archr_macs2_gene_scores))
    rownames(txnhet_archr_macs2_gene_scores_mtx) = rowData(txnhet_archr_macs2_gene_scores)$name
    saveRDS(txnhet_archr_macs2_gene_scores_mtx, file = paste0(savepath,'/txnhet_archr_macs2_gene_scores_mtx.rds'))
    saveRDS(txnhet_archr_macs2_gene_scores_meta, file = paste0(savepath,'/txnhet_archr_macs2_gene_scores_meta.rds'))
    


    ###### CV (heterogeneity) analysis ###### ------------------------------------------------------------------

    message('Starting CV calculations...(2/3)')
    if (is.null(genes_to_use)){
        genes_to_use = rownames(txnhet_archr_macs2_gene_scores_mtx)
    }

    if (is.null(samples_to_analyze)){
        samples_to_analyze = unique(txnhet_archr_macs2_gene_scores_meta$Sample)
    }

    filtered_meta_data = txnhet_archr_macs2_gene_scores_meta
    ## combine replicates if user specifies ##
    if (!is.null(replicate_map)){

        filtered_meta_data$Sample = plyr::mapvalues(filtered_meta_data$Sample, from = replicate_map[[1]], to = replicate_map[[2]])
    }
    filtered_meta_data = filtered_meta_data %>% dplyr::filter(Sample %in% samples_to_analyze)
    filtered_raw_mtx = txnhet_archr_macs2_gene_scores_mtx[rownames(txnhet_archr_macs2_gene_scores_mtx) %in% genes_to_use, rownames(filtered_meta_data)]
    colnames(filtered_meta_data)[colnames(filtered_meta_data) == 'Sample'] = 'guide'




    ## cells per guide count
    guides = unique(filtered_meta_data$guide) 

    cells_per_guide = vector()
    for (i in 1:length(guides)){
        cells_per_guide[i] = length(which(filtered_meta_data$guide == guides[i]))
    }

    cell_gd_count = data.frame(guide_name = guides, cell_count = cells_per_guide)

    if(is.null(fixed_cell_count)){
        fixed_cell_count = min(cell_gd_count$cell_count)
    }
    message(paste0('Using a fixed cell count of ', fixed_cell_count, ' per guide for heterogeneity calculations.'))



    guide_subsetted_data = lapply(X = guides, FUN = mtx_subsetter_scATAC, matrix = filtered_raw_mtx, meta_data = filtered_meta_data)

    names(guide_subsetted_data) = guides

    ## control number of cells per guide to be equivalent

    guide_subsetted_data = guide_subsetted_data[cells_per_guide >= fixed_cell_count]
    cells_to_discard = list()

    set.seed(seed)
    for (i in 1:length(guide_subsetted_data)){
        idx_to_keep = sample(1:ncol(guide_subsetted_data[[i]]), fixed_cell_count, replace = FALSE)
        idx_to_discard = setdiff(1:ncol(guide_subsetted_data[[i]]), idx_to_keep)
        cells_to_discard[[i]] = colnames(guide_subsetted_data[[i]])[idx_to_discard]
        guide_subsetted_data[[i]] = guide_subsetted_data[[i]][,idx_to_keep]
    }

    cells_to_discard = unlist(cells_to_discard)

    ## add bkg distribution of guide subsetted mtxs to list
    ## this is just mtxs where the cells are randomly assigned to each guide
    ## the bkg distribution of mtxs is made to subtract out number of genes that are expected to have increased or decreased CV when cells assigned randomly to each guide
    ## Note: the number of cells assigned to each guide is kept the same
    set.seed(seed)

    randomized_filtered_raw_mtx = Reduce(cbind, guide_subsetted_data[!startsWith(names(guide_subsetted_data), control_sample_name)])
    randomized_cell_order = sample(colnames(randomized_filtered_raw_mtx), ncol(randomized_filtered_raw_mtx), replace = FALSE)


    designation_vector = mapply(rep, guides[!startsWith(guides, control_sample_name)], fixed_cell_count, SIMPLIFY = TRUE)
    designation_vector = unlist(designation_vector)
    designation_vector = as.character(designation_vector)
    split_barcodes = split(randomized_cell_order, designation_vector)
    mtx_random_splitter = function(barcodes, mtx){return(mtx[,barcodes])}

    guide_random_subsetted_data = lapply(X = split_barcodes, FUN = mtx_random_splitter, mtx = randomized_filtered_raw_mtx)
    guide_random_subsetted_data = c(guide_subsetted_data[startsWith(names(guide_subsetted_data), control_sample_name)], guide_random_subsetted_data) ##
    names(guide_random_subsetted_data) = paste('random', c(control_sample_name, guides[!startsWith(guides, control_sample_name)]), sep = '_') ##
    guide_subsetted_data = c(guide_subsetted_data, guide_random_subsetted_data)

    ## CV calculator
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
                                            CV_gd = CVs[[i]], CV_gdCV_ctrlratio = CVs[[i]]/CVs[[control_sample_name]], mean_gd = gene_means[[i]], mean_gdmean_ctrlratio = gene_means[[i]]/gene_means[[control_sample_name]])
            names(master_df_list)[i] <- names(guide_subsetted_data)[i]
        }
        if (startsWith(names(guide_subsetted_data)[i], 'random')){
            master_df_list[[i]] = data.frame(gene = rownames(filtered_raw_mtx), CV_ctrl = CVs[[control_sample_name]], 
                                            CV_gd = CVs[[i]], CV_gdCV_ctrlratio = CVs[[i]]/CVs[[control_sample_name]], mean_gd = gene_means[[i]], mean_gdmean_ctrlratio = gene_means[[i]]/gene_means[[control_sample_name]])
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
    #cells_per_guide = rep(cells_per_guide, 2)
    #names(cells_per_guide) = names(guide_subsetted_data)
    #asymp_test_p_vals = as.data.frame(matrix(0, nrow = nrow(filtered_raw_mtx), ncol = length(master_df_list)))


    #for (i in 1:length(master_df_list)){
    #    if (!startsWith(names(guide_subsetted_data)[i], 'random')){
    #        for (j in 1:nrow(filtered_raw_mtx)){
    #            test = cvequality::asymptotic_test2(k = 2, n = c(cells_per_guide[control_sample_name], cells_per_guide[i]), s = c(gene_sds[[control_sample_name]][j], gene_sds[[i]][j]), 
    #                                    x = c(gene_means[[control_sample_name]][j], gene_means[[i]][j]))
    #            asymp_test_p_vals[j,i] = test$p_value
    #        }
    #    }
    #    if (startsWith(names(guide_subsetted_data)[i], 'random')){
    #        for (j in 1:nrow(filtered_raw_mtx)){
    #            test = cvequality::asymptotic_test2(k = 2, n = c(cells_per_guide[control_sample_name], cells_per_guide[i]), s = c(gene_sds[[control_sample_name]][j], gene_sds[[i]][j]), 
    #                                    x = c(gene_means[[control_sample_name]][j], gene_means[[i]][j]))
    #            asymp_test_p_vals[j,i] = test$p_value
    #        }    
    #    }
    #    master_df_list[[i]]$p_val = asymp_test_p_vals[,i]
    #}

    saveRDS(master_df_list, paste0(savepath, '/CVs.rds'))

    ### pvalues for CV and mean ratios for plotting ###

    # CV ratios #
    plotted_guides = names(master_df_list)[!(names(master_df_list) %in% c(control_sample_name, paste0('random_', control_sample_name)))]
    plotted_guides_with_NT = names(master_df_list)
    ratios = lapply(names(master_df_list), function(x){log2(master_df_list[[x]]$CV_gdCV_ctrlratio)})
    names(ratios) = plotted_guides_with_NT
    ratios_df <- as.data.frame(do.call(cbind, ratios))
    colnames(ratios_df) = names(master_df_list)
    ratios_df = ratios_df[, !(colnames(ratios_df) %in% c(control_sample_name, paste0('random_', control_sample_name)))]
    ratios_df = tidyr::pivot_longer(ratios_df, cols = colnames(ratios_df), values_to = 'CV_gdCV_ctrlratio', names_to = 'guide')
    order_of_guides = ratios_df %>% group_by(guide) %>% summarize(means = mean(CV_gdCV_ctrlratio, na.rm = T)) %>% arrange(means) %>% pull(guide) %>% as.character()
    ratios_df$guide = factor(ratios_df$guide, levels = order_of_guides)

    ratios = ratios[order_of_guides]
    ratios = lapply(ratios, function(x) x[!(is.infinite(x) | is.nan(x) | is.na(x))])
    tests = lapply(ratios, t.test, mu = 1)
    pvals = as.numeric(lapply(tests, function(x) return(x$p.val)))
    pvals_adj = signif(p.adjust(pvals), 4)
    pvals_adj[pvals_adj == 0] = '< e-300'
    names(pvals_adj) = order_of_guides
    CV_pvals_adj = pvals_adj
    CV_order_of_guides = order_of_guides
    CV_ratios_df = ratios_df

    # Direct CVs comparisons #
    plotted_guides = names(master_df_list)[!(names(master_df_list) %in% names(master_df_list)[startsWith(names(master_df_list), 'random_')])]
    plotted_guides_with_NT = names(master_df_list)
    CVs_list = lapply(names(master_df_list), function(x){log2(master_df_list[[x]]$CV_gd)})
    names(CVs_list) = plotted_guides_with_NT
    CVs_df <- as.data.frame(do.call(cbind, CVs_list))
    colnames(CVs_df) = names(master_df_list)
    CVs_df = CVs_df[, !(colnames(CVs_df) %in% names(master_df_list)[startsWith(names(master_df_list), 'random_')])]
    CVs_df = tidyr::pivot_longer(CVs_df, cols = colnames(CVs_df), values_to = 'CV_gd', names_to = 'guide')
    order_of_guides = CVs_df %>% group_by(guide) %>% summarize(means = mean(CV_gd, na.rm = T)) %>% arrange(means) %>% pull(guide) %>% as.character()
    order_of_guides = c(order_of_guides[order_of_guides == control_sample_name] ,order_of_guides[!order_of_guides == control_sample_name])
    CVs_df$guide = factor(CVs_df$guide, levels = order_of_guides)

    ###### Mean analysis ###### ------------------------------------------------------------------

    message('Starting mean calculations...(3/3)')

    guide_subsetted_data_pseudobulk = lapply(X = guide_subsetted_data, FUN = function(x) sparseMatrixStats::rowMeans2(x))
    guide_subsetted_data_means = guide_subsetted_data_pseudobulk
    guide_subsetted_data_means = lapply(X = guide_subsetted_data_means, FUN = function(x) data.frame(gene = rownames(guide_subsetted_data[[1]]), mean = x))
    saveRDS(guide_subsetted_data_means, paste0(savepath, '/means.rds'))


    # mean ratios #
    plotted_guides = names(master_df_list)[!(names(master_df_list) %in% c(control_sample_name, paste0('random_', control_sample_name)))]
    plotted_guides_with_NT = names(master_df_list)
    ratios = lapply(names(master_df_list), function(x){log2(master_df_list[[x]]$mean_gdmean_ctrlratio)})
    names(ratios) = plotted_guides_with_NT
    ratios_df <- as.data.frame(do.call(cbind, ratios))
    colnames(ratios_df) = names(master_df_list)
    ratios_df = ratios_df[, !(colnames(ratios_df) %in% c(control_sample_name, paste0('random_', control_sample_name)))]
    ratios_df = tidyr::pivot_longer(ratios_df, cols = colnames(ratios_df), values_to = 'mean_gdmean_ctrlratio', names_to = 'guide')
    order_of_guides = ratios_df %>% group_by(guide) %>% summarize(means = mean(mean_gdmean_ctrlratio, na.rm = T)) %>% arrange(means) %>% pull(guide) %>% as.character()
    ratios_df$guide = factor(ratios_df$guide, levels = order_of_guides)

    ratios = ratios[order_of_guides]
    ratios = lapply(ratios, function(x) x[!(is.infinite(x) | is.nan(x) | is.na(x))])
    tests = lapply(ratios, t.test, mu = 1)
    pvals = as.numeric(lapply(tests, function(x) return(x$p.val)))
    pvals_adj = signif(p.adjust(pvals), 4)
    pvals_adj[pvals_adj == 0] = '< e-300'
    names(pvals_adj) = order_of_guides
    mean_pvals_adj = pvals_adj
    mean_order_of_guides = order_of_guides
    mean_ratios_df = ratios_df

    # Direct means comparisons #
    plotted_guides = names(master_df_list)[!(names(master_df_list) %in% names(master_df_list)[startsWith(names(master_df_list), 'random_')])]
    plotted_guides_with_NT = names(master_df_list)
    means_list = lapply(names(master_df_list), function(x){log2(master_df_list[[x]]$mean_gd)})
    names(means_list) = plotted_guides_with_NT
    means_df <- as.data.frame(do.call(cbind, means_list))
    colnames(means_df) = names(master_df_list)
    means_df = means_df[, !(colnames(means_df) %in% names(master_df_list)[startsWith(names(master_df_list), 'random_')])]
    means_df = tidyr::pivot_longer(means_df, cols = colnames(means_df), values_to = 'mean_gd', names_to = 'guide')
    order_of_guides = means_df %>% group_by(guide) %>% summarize(means = mean(mean_gd, na.rm = T)) %>% arrange(means) %>% pull(guide) %>% as.character()
    order_of_guides = c(order_of_guides[order_of_guides == control_sample_name] ,order_of_guides[!order_of_guides == control_sample_name])
    means_df$guide = factor(means_df$guide, levels = order_of_guides)

    return(list(filtered_gene_scores_mtx = filtered_raw_mtx, filtered_gene_score_meta_data = filtered_meta_data, control_sample_name = control_sample_name,
        CV_ratios_df = CV_ratios_df, CV_order_of_guides = CV_order_of_guides, CV_pvals_adj = CV_pvals_adj,
        mean_ratios_df = mean_ratios_df, mean_order_of_guides = mean_order_of_guides, mean_pvals_adj = mean_pvals_adj,
        CVs_df = CVs_df, means_df = means_df))

}