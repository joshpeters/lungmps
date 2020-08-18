# ---
# Description: Functions used in this project
# *** NOT ALL FUNCTIONS ARE USED. SOME FUNCTIONS REPRESENT LEGACY OR RELATED IDEAS
# Author: Josh Peters
# ---

AddBatchFreq <- function(object, batch = "batch") {
  sample_df <- as.data.frame(table(object[[batch]]))
  colnames(sample_df) <- c("batch", "freq")
  merged <- merge(object[[]], sample_df, by = "batch")
  rownames(merged) <- merged$cell_barcode
  merged <- merged %>% select(freq)
  object <- AddMetaData(object, merged, col.name = "batch_frequency")
  return(object)
}

Scrub <- function(object, batch.var = "batch") {
  scrublet <- reticulate::import("scrublet")
  object$scrublet_score <- "NA"
  object$scrublet_label <- "NA"
  sample_df <- as.data.frame(table(object[[batch.var]]))
  colnames(sample_df) <- c("batch", "freq")
  sample_df$batch <- as.character(sample_df$batch)
  sample_df$freq <- as.numeric(sample_df$freq)

  for (i in seq_along(1:length(sample_df$batch))) {
    freq <- as.numeric(sample_df$freq[i])
    cells <- (object[[batch.var, drop = TRUE]] == sample_df$batch[i])
    if (freq < 100) {
      message(glue(">> Only {freq} cells, skipping doublet prediction"))
      object[["scrublet_score"]][cells, ] <- NA
      object[["scrublet_label"]][cells, ] <- NA
    } else {
      matrix <- as.matrix(GetAssayData(object, slot = "counts")[, cells])
      scrublet_object <- scrublet$Scrublet(t(matrix), expected_doublet_rate = 4.6e-06*freq)
      message(glue(">> Scrublet object created for iteration {i}/{length(sample_df$batch)}"))
      scores <- scrublet_object$scrub_doublets(min_counts = 3, min_cells = 3,
        min_gene_variability_pctl = 85, n_prin_comps = as.integer(30), verbose = TRUE)
      message(glue(">> Identified {sum(as.vector(scores[[2]]))}/{length(scores[[2]])} cells as doublets"))
      object[["scrublet_score"]][cells, ] <- scores[[1]]
      object[["scrublet_label"]][cells, ] <- scores[[2]]
    }
  }
  return(object)
}

#' Convert gene names based on Ensembl v86 and HGNC databases
#' @param genes List of genes to convert
#' @param verbose Verbosity toggle
#' @param hgnc HGNC database to use
#' @return Gene ids and conversions data frame
ConvertGeneNames <- function(genes, verbose = TRUE) {

  hgnc <- readRDS("data/hgnc.rds")

  # map IDs using the most recent ENSEMBL version
  num_genes <- length(genes)
  assertthat::assert_that(length(unique(genes)) == num_genes)
  ids <- data.frame(orig_name = genes)
  ensembl_mapped <- ensembldb::mapIds(EnsDb.Hsapiens.v86::EnsDb.Hsapiens.v86, keys = as.character(genes), keytype = "SYMBOL", column = "GENEID")
  ensembl_mapped <- ConvertNamedVecToDF(ensembl_mapped)
  ids <- merge(ids, ensembl_mapped, by.x = "orig_name", by.y = "name", all.x = TRUE, sort = FALSE)
  colnames(ids) <- c("orig_name", "id")
  if (verbose) usethis::ui_info("{sum(duplicated(ids$id))} duplicated Ensembl mappings")
  if (verbose) usethis::ui_oops("{sum(is.na(ids$id))} ({round(sum(is.na(ids$id))/num_genes*100, 2)}%) unidentified genes after Ensembl mapping")
  assertthat::assert_that(all.equal(as.character(ids$orig_name), as.character(genes)))
  ids <- ids %>% mutate_all(as.character)

  # check HGNC alias database
  unids <- ids[is.na(ids$id), ]
  unids$lower <- tolower(unids$orig_name)
  hgnc$lower <- tolower(hgnc$symbol)
  unids <- merge(unids, hgnc, by = "lower", all.x = TRUE, sort = FALSE)
  colnames(unids) <- c("lower", "orig_name", "ensembl_id", "hgnc_symbol", "hgnc_id")
  unids <- unids[!duplicated(unids$orig_name), ]
  unids <- unids %>% mutate_all(as.character)
  unids <- unids[, c("orig_name", "hgnc_symbol", "hgnc_id")]
  ids <- merge(ids, unids, by = "orig_name", all.x = TRUE, sort = FALSE)
  assertthat::assert_that(length(unique(ids$orig_name)) == num_genes)
  ids$combined_id <- as.character(ids$id)
  ids$combined_id[is.na(ids$combined_id)] <- as.character(ids$hgnc_id[is.na(ids$combined_id)])
  if (verbose) usethis::ui_oops("{sum(is.na(ids$combined_id))} ({round(sum(is.na(ids$combined_id))/num_genes*100, 2)}%) unidentified genes remain after alias lookup")
  colnames(ids)
  colnames(ids) <- c("orig_name", "ensembl_id", "hgnc_symbol", "hgnc_id", "id")

  # reidentify gene symbols
  symbols <- ensembldb::mapIds(EnsDb.Hsapiens.v86::EnsDb.Hsapiens.v86, keys = as.character(ids$id[!is.na(ids$id)]), keytype = "GENEID", column = "SYMBOL")
  symbols <- ConvertNamedVecToDF(symbols)
  ids <- merge(ids, symbols, by.x = "id", by.y = "name", all.x = TRUE, sort = FALSE)
  colnames(ids)[6] <- "ensembl_name"
  if (verbose) usethis::ui_done("{round(sum(ids$orig_name %in% ids$ensembl_name)/nrow(ids)*100, 2)}% agreement")
  ids <- ids[!duplicated(ids), ]
  ids <- ids %>% mutate_all(as.character)
  ids$name <- ids$ensembl_name
  ids$name[is.na(ids$name)] <- as.character(ids$orig_name[is.na(ids$name)])
  ids$name[ids$name %in% ids$name[duplicated(ids$name)] & ids$orig_name != ids$name] <- ids$orig_name[ids$name %in% ids$name[duplicated(ids$name)] & ids$orig_name != ids$name]

  assertthat::assert_that(nrow(ids) == num_genes)
  assertthat::assert_that(any(duplicated(ids$name)) == FALSE)

  rownames(ids) <- ids$orig_name
  return(ids)
}

#' Add percentage expression families
#' @param object Seurat object to add metadata to
#' @return Seurat object with metadata columns added
AddExprMeta <- function(object) {

  object[["percent_mit"]] <- PercentageFeatureSet(object, pattern = "^MT-|^MTRNR|^MTERF|^MTFR")

  # https://www.genenames.org/data/genegroup/#!/group/1054
  ribo <- intersect(rownames(object), c("MRPL1", "MRPL2", "MRPL3", "MRPL4", "MRPL10", "MRPL11", "MRPL12",
    "MRPL13", "MRPL14", "MRPL15", "MRPL16", "MRPL17", "MRPL18", "MRPL19",
    "MRPL20", "MRPL21", "MRPL22", "MRPL23", "MRPL24", "MRPL27", "MRPL28",
    "MRPL30", "MRPL32", "MRPL33", "MRPL34", "MRPL35", "MRPL36", "MRPL37",
    "MRPL38", "MRPL39", "MRPL40", "MRPL41", "MRPL42", "MRPL43", "MRPL44",
    "MRPL45", "MRPL46", "MRPL47", "MRPL48", "MRPL49", "MRPL50", "MRPL51",
    "MRPL52", "MRPL53", "MRPL54", "MRPL55", "MRPL57", "MRPL58", "RPLP0",
    "RPLP1", "RPLP2", "RPL3", "RPL3L", "RPL4", "RPL5", "RPL6", "RPL7",
    "RPL7A", "RPL7L1", "RPL8", "RPL9", "RPL10", "RPL10A", "RPL10L",
    "RPL11", "RPL12", "RPL13A", "RPL14", "RPL15", "RPL17", "RPL18A",
    "RPL19", "RPL21", "RPL22", "RPL23", "RPL23A", "RPL24", "RPL26",
    "RPL26L1", "RPL27", "RPL27A", "RPL28", "RPL29", "RPL30", "RPL31",
    "RPL32", "RPL34", "RPL35", "RPL35A", "RPL36", "RPL36A", "RPL36AL",
    "RPL37", "RPL37A", "RPL38", "RPL39", "RPL39L", "UBA52", "RPL41",
    "MRPL1", "MRPL2", "MRPL3", "MRPL4", "MRPL9", "MRPL10", "MRPL11",
    "MRPL12", "MRPL13", "MRPL14", "MRPL15", "MRPL16", "MRPL17", "MRPL18",
    "MRPL19", "MRPL20", "MRPL21", "MRPL22", "MRPL23", "MRPL24", "MRPL27",
    "MRPL28", "MRPL30", "MRPL32", "MRPL33", "MRPL34", "MRPL35", "MRPL36",
    "MRPL37", "MRPL38", "MRPL39", "MRPL40", "MRPL41", "MRPL42", "MRPL43",
    "MRPL44", "MRPL45", "MRPL46", "MRPL47", "MRPL48", "MRPL49", "MRPL50",
    "MRPL51", "MRPL52", "MRPL53", "MRPL54", "MRPL55", "MRPL57", "MRPS2",
    "MRPS5", "MRPS6", "MRPS7", "MRPS9", "MRPS10", "MRPS11", "MRPS12",
    "MRPS14", "MRPS15", "MRPS16", "MRPS17", "MRPS18A", "MRPS18B",
    "MRPS18C", "MRPS21", "MRPS22", "MRPS23", "MRPS24", "MRPS25",
    "MRPS26", "MRPS27", "MRPS28", "DAP3", "MRPS30", "MRPS31", "MRPS33",
    "MRPS34", "MRPS35", "MRPS36", "MRPS2", "MRPS10", "MRPS11", "MRPS12",
    "MRPS14", "MRPS15", "MRPS16", "MRPS17", "MRPS18A", "MRPS18B",
    "MRPS18C", "MRPS21", "MRPS22", "MRPS23", "MRPS24", "MRPS25",
    "MRPS26", "MRPS27", "MRPS28", "MRPS30", "MRPS31", "MRPS33", "MRPS34",
    "MRPS35", "MRPS36", "RPS2", "RPS3", "RPS3A", "RPS4X", "RPS4Y1",
    "RPS4Y2", "RPS5", "RPS6", "RPS7", "RPS8", "RPS9", "RPS10", "RPS11",
    "RPS12", "RPS13", "RPS14", "RPS15", "RPS15A", "RPS16", "RPS17",
    "RPS18", "RPS19", "RPS20", "RPS21", "RPS23", "RPS24", "RPS25",
    "RPS26", "RPS27", "RPS27A", "RPS27L", "RPS28", "RPS29", "FAU"))
  object[["percent_rib"]] <- PercentageFeatureSet(object, features = ribo)

  # https://www.genenames.org/data/genegroup/#!/group/582
  hsp <- intersect(rownames(object), c("BBS10", "BBS12", "TCP1", "CCT2", "CCT3", "CCT4", "CCT5", "CCT6A",
    "CCT6B", "CCT7", "CCT8", "CLPB", "HSPD1", "HSPE1", "MKKS", "DNAJA1",
    "DNAJA2", "DNAJA3", "DNAJA4", "DNAJB1", "DNAJB2", "DNAJB3", "DNAJB4",
    "DNAJB5", "DNAJB6", "DNAJB7", "DNAJB8", "DNAJB9", "DNAJB11",
    "DNAJB12", "DNAJB13", "DNAJB14", "DNAJC1", "DNAJC2", "DNAJC3",
    "DNAJC4", "DNAJC5", "DNAJC5B", "DNAJC5G", "DNAJC6", "DNAJC7",
    "DNAJC8", "DNAJC9", "DNAJC10", "DNAJC11", "DNAJC12", "DNAJC13",
    "DNAJC14", "DNAJC15", "DNAJC16", "DNAJC17", "DNAJC18", "DNAJC19",
    "HSCB", "DNAJC21", "DNAJC22", "SEC63", "DNAJC24", "DNAJC25",
    "GAK", "DNAJC27", "DNAJC28", "SACS", "DNAJC30", "HSPA1A", "HSPA1B",
    "HSPA1L", "HSPA2", "HSPA4", "HSPA4L", "HSPA5", "HSPA6", "HSPA7",
    "HSPA8", "HSPA9", "HSPA12A", "HSPA12B", "HSPA13", "HSPA14", "HSPH1",
    "HYOU1", "HSP90AA1", "HSP90AA3P", "HSP90AB1", "HSP90B1", "TRAP1",
    "HSPB1", "HSPB2", "HSPB3", "CRYAA", "CRYAB", "HSPB6", "HSPB7",
    "HSPB8", "HSPB9", "ODF1", "HSPB11"))

  object[["percent_hsp"]] <- PercentageFeatureSet(object, features = hsp)

  hb <- intersect(rownames(object),
    c("HBA1", "HBA2", "HBB", "HBD", "HBE1", "HBG1", "HBG2", "HBM", "HBQ1", "HBZ"))
  object[["percent_hgb"]] <- PercentageFeatureSet(object, features = hb)

  return(object)
}

#' Filter cells from Seurat object
#' @param object Seurat object
#' @param variables Variables used to drop cells based on MADs
#' @param batch Batch column
#' @param nmads MADS to use for batch-specific filtering
#' @param percent_mito_cutoff
#' @param percent_ribo_cutoff
#' @param percent_hb_cutoff
#' @param percent_hsp_cutoff
#' @param nCount_RNA_cutoff
#' @param nFeature_RNA_cutoff
#' @param remove Toggle to subset cells from object and return filtered object
#' @param qc_column Column name to store quality control results in
#' @return Seurat object with qc_column
FilterCells <- function(
  object,
  variables = c("nCount_RNA", "nFeature_RNA"),
  batch = "batch",
  nmads = 3,
  percent_mito_cutoff = 20,
  percent_ribo_cutoff = 50,
  percent_hb_cutoff = 5,
  percent_hsp_cutoff = 5,
  nCount_RNA_cutoff = 200,
  nFeature_RNA_cutoff = 100,
  remove = FALSE,
  qc_column = "passing_qc"
) {

  drop <- list()
  tryCatch(drop[[length(drop)+1]] <- object[["percent_mit"]][[1]] > percent_mito_cutoff,
    error = function(e) {message("No % mito threshold needed")})

  tryCatch(drop[[length(drop)+1]] <- object[["percent_rib"]][[1]] > percent_ribo_cutoff,
    error = function(e) {message("No % hemoglobin threshold needed")})

  tryCatch(drop[[length(drop)+1]] <- object[["percent_hgb"]][[1]] > percent_hb_cutoff,
    error = function(e) {message("No % hemoglobin threshold needed")})

  tryCatch(drop[[length(drop)+1]] <- object[["percent_hsp"]][[1]] > percent_hsp_cutoff,
    error = function(e) {message("No % hemoglobin threshold needed")})

  tryCatch(drop[[length(drop)+1]] <- object[["nCount_RNA"]][[1]] < nCount_RNA_cutoff,
    error = function(e) {message("No nCount threshold needed")})

  tryCatch(drop[[length(drop)+1]] <- object[["nFeature_RNA"]][[1]] < nFeature_RNA_cutoff,
    error = function(e) {message("No nCount threshold needed")})

  threshold_drop <- Reduce("|", drop)
  threshold_drop[is.na(threshold_drop)] <- FALSE
  ui_info("{sum(threshold_drop)} cells dropped by thresholds")
  thresholded_object <- object[, !(threshold_drop)]

  drop <- list()
  for (i in seq_along(variables)) {
    drop[[i]] <- scater::isOutlier(thresholded_object[[variables[i]]][[1]], nmads = nmads,
      log = TRUE, type = "both", batch = thresholded_object[[batch]][[1]])
  }
  all_drop <- c(drop, list(threshold_drop))
  all_drop <- Reduce("|", all_drop)
  all_drop[is.na(all_drop)] <- FALSE
  ui_info("{sum(all_drop)} cells dropped by all filters")
  ui_info("{sum(Reduce('|', drop))} cells dropped by MADs")

  if (remove) {
    object <- object[, !(all_drop)]
  }
  if (!remove) {
    object[[qc_column]] <- "TRUE"
    object[[qc_column]][all_drop, ] <- "FALSE"
  }
  return(object)
}

#' Process Seurat object utilizing a standard procedure
BaseProcess <- function(object, nfeatures = 3000) {
  object <- Seurat::NormalizeData(object)
  object <- Seurat::FindVariableFeatures(object, selection.method = "vst", nfeatures = nfeatures)
  if (ncol(object) > 2E4) {
    object <- Seurat::ScaleData(object, features = VariableFeatures(object), block.size = 250)
  } else {
    object <- Seurat::ScaleData(object, features = VariableFeatures(object), block.size = 1000)
  }
  object <- Seurat::RunPCA(object, npcs = 50, features = VariableFeatures(object),
    verbose = FALSE, seed.use = 1, weight.by.var = TRUE)
  usethis::ui_todo("PCA calculated... determining PCs")
  object <- CalculatePCHeuristics(object)
  usethis::ui_info("TP estimation: {unique(object@misc$pc_heuristics$tp_pc)}
    Max distance: {unique(object@misc$pc_heuristics$max_dist_pcs)}
    Derivative: {unique(object@misc$pc_heuristics$deriv_pcs)}
    Range: {unique(object@misc$pc_heuristics$percent_cutoff_pcs)}")
  return(object)
}

#' Calculate heuristic metrics to determine optimal number of principal components for downstream analyses
#'
#' @param object Seurat object
#' @param percent.stdev.range Percent of PCA standard deviation range to use for
#' @param set.pcs Number of PCs to use for maximum distance and derivative calculations
#' @param derivs.change.threshold Threshold used for second derivative changes
#'
#' @return Returns `object` with optimal PC information within the `object@misc` slot
#' @export
#'
CalculatePCHeuristics <- function(
  object,
  reduction = "pca",
  store.name = "pc_heuristics",
  percent.stdev.range = 0.05,
  set.pcs = NA,
  derivs.change.threshold = 0.02,
  force_tp = FALSE
) {

  deriv <- function(x, y) diff(y) / diff(x)
  middle_pts <- function(x) x[-1] - diff(x) / 2
  if (reduction == "harmony") {
    object@reductions$harmony@stdev <- as.numeric(apply(object@reductions$harmony@cell.embeddings, 2, sd))
  }

  object_stdev <- Seurat::Stdev(object, reduction = reduction)

  stdev_range <- range(object_stdev)[2] - range(object_stdev)[1]
  cutoff <- min(object_stdev) + stdev_range * percent.stdev.range
  if (is.na(set.pcs)) {
    pcs_to_use <- max(which(object_stdev > cutoff))
    usethis::ui_info("Using percent of range cutoff ({usethis::ui_value(pcs_to_use)}) for distance and derivative determinations")
  } else {
    pcs_to_use <- set.pcs
    usethis::ui_info("Using set number of PCs ({usethis::ui_value(pcs_to_use)}) for distance and derivative determinations")
  }

  pcs_stdev <- object_stdev[1:pcs_to_use]
  pcs <- seq(1:length(pcs_stdev))

  slope <- (pcs_stdev[1] - pcs_stdev[length(pcs)])/(pcs[1]-pcs[length(pcs)])
  intercept <- pcs_stdev[1]-slope*pcs[1]
  b <- c(pcs[1], pcs_stdev[1])
  c <- c(pcs[length(pcs_stdev)], pcs_stdev[length(pcs_stdev)])
  dist <- vector()
  for (i in seq_along(1:length(pcs_stdev))) {
    a <- c(pcs[i], pcs_stdev[i])
    v1 <- b - c
    v2 <- a - b
    m <- cbind(v1,v2)
    dist[i] <- abs(det(m))/sqrt(sum(v1*v1))
  }
  derivs <- data.frame(pcs = 1:pcs_to_use, std = object_stdev[1:pcs_to_use],
    d1 = c(0, deriv(1:pcs_to_use, object_stdev[1:pcs_to_use])),
    d2 = c(0, 0, deriv(middle_pts(1:pcs_to_use), deriv(1:pcs_to_use, object_stdev[1:pcs_to_use]))))

  derivs$max_dist_pcs <- which.max(dist)
  derivs$percent_cutoff_pcs <- max(which(object_stdev > cutoff))
  derivs$deriv_pcs <- min(which(derivs$d2[which.max(dist):pcs_to_use] < derivs.change.threshold &
      derivs$d2[which.max(dist):pcs_to_use] > -derivs.change.threshold)) + which.max(dist) - 1
  derivs$slope <- slope
  derivs$intercept <- intercept

  if (ncol(object) < 50000 | force_tp == TRUE) {
    usethis::ui_todo("Calculating TP")
    tp <- intrinsicDimension::maxLikGlobalDimEst(object@reductions$pca@cell.embeddings, k = 10)
    tp <- ceiling(tp$dim.est)
  } else {
    usethis::ui_oops("Object size too large, tp = 101")
    tp <- 101
  }
  derivs$tp_pc <- tp
  object@misc[[store.name]] <- derivs
  usethis::ui_done("Heuristic metrics stored in {ui_field(VarToString(object@misc))} slot named {ui_value(store.name)}.")
  return(object)
}

HarmonyCorrection <- function (
  object,
  batch
)
{
  # determine dims to use
  # harmony_dims <- tryCatch(
  #   {ifelse(plyr::round_any(nrow(object@misc$pc_heuristics), 5, f = ceiling) < 30,
  #     30, plyr::round_any(nrow(object@misc$pc_heuristics), 5, f = ceiling))},
  #   error = function(cond) {return(FALSE)})
  # if (harmony_dims == FALSE) harmony_dims <- 30
  harmony_dims <- 30
  usethis::ui_info("Using {harmony_dims} dimensions for integration")
  object$harmony_dims_used <- harmony_dims
  dims <- unique(object$harmony_dims_used)

  # run harmony
  corrected <- harmony::RunHarmony(object, group.by.vars = batch, dims.use = 1:dims,
    max.iter.harmony = 100, verbose = TRUE, assay.use = "RNA", plot_convergence = FALSE,
    thetha = 2, sigma = 0.1, lambda = 1, max.iter.cluster = 30)

  # return corrected object
  return(corrected)
}

BaseCluster <- function(
  object,
  dims = 20,
  graph.name = "RNA_snn",
  res.start = 0.001,
  res.end = 1,
  num.res = 30,
  num.clusters.lower = 5,
  num.clusters.upper = 30,
  mod.similarity = 0.99,
  #clustering.algorithm = 4,
  verbose = TRUE
) {

  # find optimal Leiden clustering based on maximal modularity
  g <- igraph::graph_from_adjacency_matrix(adjmatrix = object@graphs[[graph.name]], mode = "directed", weighted = TRUE)
  resolution_parameter <-  signif(exp(seq(log(res.start), log(res.end), length.out = num.res)), 3)
  table_results <- LeidenClustering(g, resolution_parameter)
  first_round_table_results <- table_results
  lower_res <- min(table_results$resolution_parameter[table_results$cluster_count >= num.clusters.lower])
  upper_res <- max(table_results$resolution_parameter[table_results$cluster_count <= num.clusters.upper])
  resolution_parameter <- signif(seq(lower_res, upper_res, length.out = num.res), 3)
  ui_done("First round of clustering complete.
    Lower resolution: {resolution_parameter[1]},
    Upper resolution: {resolution_parameter[length(resolution_parameter)]}")
  table_results <- LeidenClustering(g, resolution_parameter)
  max_modularity <- max(table_results$modularity)
  modularity_cutoff <- max_modularity*mod.similarity
  max_resolution <- table_results$resolution_parameter[table_results$modularity == max_modularity]
  ui_info("{sum(table_results$modularity >= modularity_cutoff)} alternative resolutions")
  if (any(table_results$modularity >= modularity_cutoff)) {
    ui_info("Choosing non-optimal modularity...")
    final_resolution_parameter <- max(table_results$resolution_parameter[table_results$modularity >= modularity_cutoff])
  } else if (!any(table_results$modularity >= modularity_cutoff)) {
    final_resolution_parameter <- table_results$resolution_parameter[table_results$modularity == max_modularity]
  }
  ui_done("Second round of clustering complete.
    Final resolution: {final_resolution_parameter} vs. {max_resolution}")
  final_result <- leidenbase::leiden_find_partition(g,
    partition_type = "CPMVertexPartition",
    initial_membership = NULL,
    edge_weights = NULL,
    node_sizes = NULL,
    seed = as.numeric(Sys.time()),
    resolution_parameter = final_resolution_parameter,
    num_iter = 30,
    verbose = verbose)
  out_result <- list(membership = final_result[['membership']],
    modularity = final_result[['modularity']])
  names(out_result$membership) <- colnames(object@graphs[[graph.name]])
  leiden_results <- list(optim_res = out_result, second_results = table_results, first_results = first_round_table_results)
  ids <- GroupSingletons(leiden_results$optim_res$membership, object@graphs[[graph.name]], min.size = 10, group.singletons = TRUE, verbose = TRUE)
  object <- AddMetaData(object, ids, col.name = "base_leiden")
  object@misc$leiden_results <- leiden_results
  table(object@active.ident)
  return(object)
}

GroupSingletons <- function(ids, SNN, min.size = 9, group.singletons = TRUE, verbose = TRUE) {
  # identify singletons
  singletons <- c()
  singletons <- names(x = which(x = table(ids) <= min.size))
  singletons <- intersect(x = unique(x = ids), singletons)
  if (!group.singletons) {
    ids[which(ids %in% singletons)] <- "singleton"
    return(ids)
  }
  # calculate connectivity of singletons to other clusters, add singleton
  # to cluster it is most connected to
  cluster_names <- as.character(x = unique(x = ids))
  cluster_names <- setdiff(x = cluster_names, y = singletons)
  connectivity <- vector(mode = "numeric", length = length(x = cluster_names))
  names(x = connectivity) <- cluster_names
  new.ids <- ids
  for (i in singletons) {
    i.cells <- names(which(ids == i))
    for (j in cluster_names) {
      j.cells <- names(which(ids == j))
      subSNN <- SNN[i.cells, j.cells]
      set.seed(1) # to match previous behavior, random seed being set in WhichCells
      if (is.object(x = subSNN)) {
        connectivity[j] <- sum(subSNN) / (nrow(x = subSNN) * ncol(x = subSNN))
      } else {
        connectivity[j] <- mean(x = subSNN)
      }
    }
    m <- max(connectivity, na.rm = T)
    mi <- which(x = connectivity == m, arr.ind = TRUE)
    closest_cluster <- sample(x = names(x = connectivity[mi]), 1)
    ids[i.cells] <- closest_cluster
  }
  if (length(x = singletons) > 0 && verbose) {
    message(paste(
      length(x = singletons),
      "singletons identified.",
      length(x = unique(x = ids)),
      "final clusters."
    ))
  }
  return(ids)
}

LLClassifier <- function(object, index) {
  scmatrix <- GetAssayData(object, assay = "RNA", slot = "counts")

  index_norm <- t(t(index) / colSums(index))
  index_norm = index_norm[rownames(index_norm) %in% rownames(scmatrix), ]
  index_norm = as.data.frame(index_norm)
  scmatrix = scmatrix[match(rownames(index_norm), rownames(scmatrix)),]

  loglikelihood = apply(scmatrix, 2, function(x) colSums(x * log(index_norm)))

  post_probs = apply(loglikelihood , 2, function(x) {exp(x - max(x) - log(sum(x / max(x))))})
  top_classified = apply(post_probs, 2, function(x) {
    rownames(post_probs)[order(x, decreasing = T)[1]]
  })
  return(list(ll = loglikelihood, pp = post_probs, call = top_classified))
}

WriteNMF <- function(
  object,
  rootpath,
  assay = "RNA",
  slot = "counts",
  zip = TRUE
) {

  # if (!is_empty(remove_genes)) {
  #   keep_genes <- setdiff(rownames(data), remove_genes)
  #   genes <- genes[genes$gene %in% keep_genes, , drop = FALSE]
  #   data <- data[keep_genes, ]
  # }

  matrix <- GetAssayData(object, assay = assay, slot = slot)
  fraction_neg <- sum(matrix < 0)/length(matrix)*100
  message(glue(">> {fraction_neg}% values are negative"))
  matrix[matrix < 0] <- 0

  DropletUtils::write10xCounts(path = rootpath,
    x = matrix,
    barcodes = colnames(matrix),
    gene.id = rownames(matrix),
    gene.symbol = rownames(matrix),
    gene.type = "Gene Expression",
    overwrite = TRUE,
    type = "sparse",
    version = "3")

  meta_path <- glue("{rootpath}/meta.tsv")
  write_tsv(object@meta.data, meta_path)

  var_genes <- as.data.frame(VariableFeatures(object, assay = assay))
  colnames(var_genes) <- "gene"
  var_features_path <- glue("{rootpath}/var_genes.tsv")
  write_tsv(var_genes, var_features_path, col_names = FALSE)

  if (zip) {
    invisible(lapply(list(meta_path), R.utils::gzip))
  }

}

#' Rename cells and stash identity
#' @param object Seurat object
#' @param markers Markers
#' @param slot Slot of object to choose
#' @return stm_logFC dataframe with second-to-max logFCs
#' @importFrom Seurat Idents

CalculateSTMRatio <- function(
  object,
  features,
  assay = "RNA",
  slot = "data",
  pseudocount.use = 1
) {

  idents_all <- sort(x = unique(x = Idents(object = object)))
  mean.fxn <- if (slot != "scale.data") {
    switch(
      EXPR = slot,
      "data" = function(x) {
        return(log(x = mean(x = expm1(x = x)) + pseudocount.use))
      },
      "counts" = function(x) {
        return(log(x = mean(x = x) + pseudocount.use))
      }
    )
  } else {
    mean
  }

  MaxN <- function(x, n) {
    order(x, decreasing = TRUE)[n]
  }

  GetSTM <- function(x, idx) {
    x <- as.numeric(x[1:idx])
    max1 <- MaxN(x[1:idx], n = 1)
    max2 <- MaxN(x[1:idx], n = 2)
    diff <- x[max1]/x[max2]
    return(diff)
  }

  mean_expr <- list()
  for (i in 1:length(x = idents_all)) {
    mean_expr[[i]] <- apply(
      X = GetAssayData(object, assay = assay, slot = slot)[features, WhichCells(object, idents = idents_all[i])],
      MARGIN = 1,
      FUN = mean.fxn
    )
  }

  mean_exprs <- as.data.frame(matrix(unlist(mean_expr),
    nrow = length(unlist(mean_expr[1]))))
  colnames(mean_exprs) <- idents_all
  mean_exprs$genes <- features
  mean_exprs$stm <- apply(mean_exprs, 1, GetSTM, idx = length(idents_all))
  return(mean_exprs)
}

PlotGeneFraction <- function(object, gene, cluster, dist.only = FALSE) {
  norm_gene <- object@assays$RNA@data[gene, ]
  norm_gene <- data.frame(gene = norm_gene, cluster = object[[cluster]])
  gene_frac <- norm_gene %>%
    group_by(!!sym(cluster)) %>%
    summarize(n = n(), n_above = sum(gene > 0), frac_above = n_above/n)

  if (dist.only) {
    dist <- ggplot(gene_frac, aes(x = reorder(!!sym(cluster), frac_above), y = frac_above)) +
      scale_y_continuous(breaks = seq(0.05, 1, 0.05)) +
      geom_point() + theme_classic() + labs(x = "Cluster", y = "Fraction cells expressing")
    return(dist)
  }

  if (!dist.only) {
    umap <- DimPlot(object, group.by = cluster, label = TRUE) + NoLegend()
    vln <- VlnPlot(object, features = gene, group.by = cluster, slot = "data") + NoLegend()
    feature <- FeaturePlot(object, features = gene, order = TRUE, min.cutoff = "q10")
    dist <- ggplot(gene_frac, aes(x = reorder(!!sym(cluster), frac_above), y = frac_above)) +
      scale_y_continuous(breaks = seq(0.05, 1, 0.05)) +
      geom_point() + theme_classic() + labs(x = "Cluster", y = "Fraction cells expressing")

    top_row <- cowplot::plot_grid(umap, feature, labels = c("A", "B"), label_size = 12, nrow = 1)
    bottom_row <- cowplot::plot_grid(vln, dist, labels = c("C", "D"), rel_widths = c(2, 1), nrow = 1)
    return(cowplot::plot_grid(top_row, bottom_row, ncol = 1))
  }

}

PlotFeatureFraction <- function(object, feature, cluster, dist.only = FALSE) {
  norm_feature <- object@meta.data[, feature]
  norm_feature <- data.frame(feature = norm_feature, cluster = object[[cluster]])
  feature_frac <- norm_feature %>%
    group_by(!!sym(cluster)) %>%
    summarize(n = n(), n_above = sum(feature > 0), frac_above = n_above/n)

  if (dist.only) {
    dist <- ggplot(feature_frac, aes(x = reorder(!!sym(cluster), frac_above), y = frac_above)) +
      scale_y_continuous(breaks = seq(0.05, 1, 0.05)) +
      geom_point() + theme_classic() + labs(x = "Cluster", y = "Fraction cells expressing")
    return(dist)
  }

  if (!dist.only) {
    umap <- DimPlot(object, group.by = cluster, label = TRUE, reduction = "harmony_umap") + NoLegend()
    vln <- VlnPlot(object, features = feature, group.by = cluster, slot = "data", sort = TRUE, pt.size = 0) + NoLegend()
    feature_plot <- FeaturePlot(object, features = feature, order = TRUE, min.cutoff = "q10", reduction = "harmony_umap")
    dist <- ggplot(feature_frac, aes(x = reorder(!!sym(cluster), frac_above), y = frac_above)) +
      scale_y_continuous(breaks = seq(0.05, 1, 0.05)) +
      geom_point() + theme_classic() + labs(x = "Cluster", y = "Fraction cells expressing")

    top_row <- cowplot::plot_grid(umap, feature_plot, labels = c("A", "B"), label_size = 12, nrow = 1)
    bottom_row <- cowplot::plot_grid(vln, dist, labels = c("C", "D"), rel_widths = c(2, 1), nrow = 1)
    return(cowplot::plot_grid(top_row, bottom_row, ncol = 1))
  }

}

LabelGeneFraction <- function(object, gene, cluster, threshold) {
  norm_gene <- object@assays$RNA@data[gene, ]
  norm_gene <- data.frame(gene = norm_gene, cluster = object[[cluster]])
  gene_frac <- norm_gene %>%
    group_by(!!sym(cluster)) %>%
    summarize(n = n(), n_above = sum(gene > 0), frac_above = n_above/n)
  gene_pos_clusters <- droplevels(as.factor(gene_frac[gene_frac$frac_above > threshold, ][[cluster]]))
  object[[gene]] <- FALSE
  gene_pos_clusters <- droplevels(as.factor(gene_frac[gene_frac$frac_above > threshold, ][[cluster]]))
  gene_pos_clusters
  object[[gene]] <- FALSE
  ident <- object@active.ident
  object <- SetIdent(object, value = cluster)
  gene_pos_cells <- WhichCells(object, idents = as.character(gene_pos_clusters))
  object[[gene]][gene_pos_cells, ]<- TRUE
  object@active.ident <- ident
  return(object)
}

DetermineUMAPModules <- function(
  cds,
  genes,
  res.start = 0.01,
  res.end = 1,
  num.res = 30,
  iterations,
  blacklist,
  sim.cutoff
) {

  if (!missing(blacklist)) {
    genes <- setdiff(genes, blacklist)
  }

  assertthat::assert_that(all(genes %in% rownames(cds)))

  # generate modules
  modules <- list()
  res_seq <- exp(seq(log(res.start), log(res.end), length.out = num.res))
  for (i in seq_along(1:iterations)) {
    usethis::ui_todo("Generating modules for iteration {i}")
    gene_module_df <- monocle3::find_gene_modules(cds[genes, ],
      resolution = res_seq,
      verbose = FALSE)
    modules[[i]] <- gene_module_df
  }

  # compare module similiarity
  modules_split <- lapply(modules, function(x) {
    x <- split(x$id, x$module)
    names(x) <- paste0("M", names(x))
    return(x)
  })
  names(modules_split) <- paste0("R", seq(iterations))

  comparisons <- pbapply::pblapply(names(modules_split), function(x) lapply(names(modules_split), function(y) {
    gom.obj <- GeneOverlap::newGOM(modules_split[[x]], modules_split[[y]], genome.size = length(genes))
    odds <- GeneOverlap::getMatrix(gom.obj, "odds.ratio")
    jacc <- GeneOverlap::getMatrix(gom.obj, "Jaccard")
    pvals <- GeneOverlap::getMatrix(gom.obj, "pval")
    total <- ConvertGOMatrix(pvals, jacc, odds, x, y)
    return(total)
  }))

  binded_comps <- do.call(rbind, do.call(rbind, comparisons))
  binded_comps <- binded_comps %>% tidyr::unite("set1_name", set1, dataset1, remove = FALSE)
  binded_comps <- binded_comps %>% tidyr::unite("set2_name", set2, dataset2, remove = FALSE)
  filtered_comps <- binded_comps %>% filter(dataset1 != dataset2 & jacc >= sim.cutoff & adj_pvals < 0.05)

  dedupped <- filtered_comps[!duplicated(t(apply(filtered_comps[c("set1_name", "set2_name")], 1, sort))), ]
  comp_matrix <- filtered_comps %>% dplyr::select(set1_name, set2_name, jacc)
  comp_matrix <- comp_matrix %>% tidyr::spread(set2_name, jacc)
  rownames(comp_matrix) <- comp_matrix$set1_name
  comp_matrix[, 1] <- NULL
  comp_matrix <- as.matrix(comp_matrix)
  mat <- comp_matrix
  mat[is.na(mat)] <- 0

  d <- dist(x = mat, method = "euclidean")
  hc <- hclust(d, method = "ward.D2")
  sub_grp <- data.frame(cluster = cutree(hc, k = round(mean(sapply(modules_split, length)))))
  sub_grp$module <- rownames(sub_grp)
  groups <- split(sub_grp$module, sub_grp$cluster)

  split_groups <- lapply(groups, stringr::str_split, pattern = "_")
  module_genes <- list()
  for (i in seq_along(1:length(split_groups))) {
    usethis::ui_todo("Aggregating module {i}
      {length(split_groups[[i]])} iterations found")
    genes <- list()
    for (j in seq_along(1:length(split_groups[[i]]))) {
      x <- split_groups[[i]]
      replicate <- x[[j]][2]
      module <- x[[j]][1]
      genes[[j]] <- modules_split[[replicate]][module][[1]]
    }
    module_genes[[i]] <- genes
  }

  return(module_genes)
}

ConvertGOMatrix <- function(pvals, jacc, odds, comp1, comp2) {
  pvals <- as.data.frame(pvals)
  pvals$set1 <- rownames(pvals)
  pvals <- gather(pvals, key = "set2", value = "pvals", -set1)
  pvals$adj_pvals <- p.adjust(pvals$pvals, method = "fdr")
  pvals$log_adj_pvals <- -log10(pvals$adj_pvals)

  jacc <- as.data.frame(jacc)
  jacc$set1 <- rownames(jacc)
  jacc <- gather(jacc, key = "set2", value = "jacc", -set1)

  odds <- as.data.frame(odds)
  odds$set1 <- rownames(odds)
  odds <- gather(odds, key = "set2", value = "odds", -set1)

  total <- pvals
  total$jacc <- jacc$jacc
  total$odds <- odds$odds
  total$comp1 <- comp1
  total$comp2 <- comp2
  return(total)
}



RemoveSinglets <- function(
  object,
  cell.limit = 10,
  cluster
){
  non_singlets <- names(table(object[[cluster]]))[table(object[[cluster]]) >= cell.limit]
  message(glue("Removing {sum(table(object[[cluster]])[table(object[[cluster]]) < cell.limit])} cells as small clusters"))
  object <- SetIdent(object, value = cluster)
  filtered <- object[, WhichCells(object, idents = non_singlets)]
  filtered@active.ident <- droplevels(filtered@active.ident)
  return(filtered)

}

ReplaceClusterAssignment <- function(cluster_to_replace, cluster_replacement, object) {
  object@active.ident[object@active.ident == cluster_to_replace] <- cluster_replacement
  return(object)
}

CellTypeEnrichment <- function (
  markers,
  cluster
) {
  dbs <- c("ARCHS4_Tissues",
    "ARCHS4_Cell-lines")
  enr <- enrichR::enrichr(markers %>% filter(padj < 1E-3 & logFC > log(1.25)) %>% filter(group == cluster) %>% pull(feature), dbs)
  print(enr$ARCHS4_Tissues[1:5, 1])
  print(enr$`ARCHS4_Cell-lines`[1:5, 1])
}

GOEnrichment <- function (
  markers,
  cluster
) {
  if (!any(grepl("^padj", colnames(markers)))) {
    colnames(markers)[grep("p_val_adj", colnames(markers))] <- "padj"
  }

  if (!any(grepl("^logFC", colnames(markers)))) {
    colnames(markers)[grep("logFC", colnames(markers))] <- "logFC"
  }

  if (!any(grepl("^group", colnames(markers)))) {
    colnames(markers)[grep("cluster", colnames(markers))] <- "group"
  }

  if (!any(grepl("^feature", colnames(markers)))) {
    colnames(markers)[grep("gene", colnames(markers))] <- "feature"
  }

  library(enrichR)
  dbs <- c("GO_Biological_Process_2018",
    "GO_Cellular_Component_2018", "GO_Molecular_Function_2018")
  enr <- enrichR::enrichr(markers %>% filter(padj < 1E-3 & logFC > log(1.25)) %>% filter(group == cluster) %>% pull(feature), dbs)
  return(enr)
}


SetCellTypes <- function (
  object,
  labels,
  meta.slot
) {
  names(x = labels) <- levels(x = object)
  object <- RenameIdents(object, labels)
  object <- BuildClusterTree(object, reorder = TRUE)
  object[[meta.slot]] <- Idents(object)
  return(object)

}

AddCellProps <- function (
  object,
  batch,
  type,
  total,
  meta.slot
) {

  type_counts <- object@meta.data %>%
    dplyr::group_by(!!sym(batch), !!sym(type)) %>%
    dplyr::summarize(freq = n())
  total_counts <- object@meta.data %>%
    dplyr::group_by(!!sym(batch), !!sym(total)) %>%
    dplyr::summarize(freq = n())

  props <- merge(type_counts, total_counts, by = batch)
  props <- props %>% dplyr::filter(!!sym(total) == TRUE & !!sym(type) == TRUE) %>%
    mutate(type_props = freq.x/freq.y) %>% dplyr::select(!!sym(batch), type_props)

  meta <- object@meta.data
  tryCatch(meta <- tibble::rownames_to_column(meta, "cell_barcode"), error = function(e) {message("cell_barcode already present")})
  merged <- merge(object@meta.data, props, by = batch)
  rownames(merged) <- merged$cell_barcode
  merged <- merged %>% dplyr::select(type_props)
  object <- AddMetaData(object, merged, col.name = meta.slot)
  return(object)

}

ConvertToCDS <- function(
  object,
  assay = "RNA",
  method_umap = "mnp_harmony_umap",
  grouping = "cell_subset_mnp"
) {

  if (assay == "RNA") {
    matrix <- object@assays[[assay]]@counts[rownames(object@reductions$pca@feature.loadings), ]
    pd <- object@meta.data
    fd <- as.data.frame(rownames(matrix))
    colnames(fd) <- c("gene_short_name")
    rownames(fd) <- fd$gene_short_name
    cds <- monocle3::new_cell_data_set(matrix,
      cell_metadata = pd, gene_metadata = fd)
    SingleCellExperiment::reducedDims(cds)[["UMAP"]] <- object@reductions[[method_umap]]@cell.embeddings
    SingleCellExperiment::reducedDims(cds)[["PCA"]] <- object@reductions$harmony@cell.embeddings
    cds_dims <- ncol(object@reductions$harmony)
    cds@preprocess_aux$gene_loadings <- object@reductions$harmony@feature.loadings[, 1:cds_dims]
    cds@preprocess_aux$prop_var_expl <- object@reductions$pca@stdev[1:cds_dims]^2 /
      sum(object@reductions$pca@stdev[1:cds_dims]^2)
  }

  if (assay == "integrated") {
    matrix <- object@assays[[assay]]@data[rownames(object@reductions$pca@feature.loadings), ]
    pd <- object@meta.data
    fd <- as.data.frame(rownames(matrix))
    colnames(fd) <- c("gene_short_name")
    rownames(fd) <- fd$gene_short_name
    cds <- new_cell_data_set(matrix,
      cell_metadata = pd, gene_metadata = fd)
    reducedDims(cds)[["UMAP"]] <- object@reductions[[method_umap]]@cell.embeddings
    reducedDims(cds)[["PCA"]] <- object@reductions$pca@cell.embeddings
    cds_dims <- ncol(object@reductions$pca)
    cds@preprocess_aux$gene_loadings <- object@reductions$pca@feature.loadings[, 1:cds_dims]
    cds@preprocess_aux$prop_var_expl <- object@reductions$pca@stdev[1:cds_dims]^2 /
      sum(object@reductions$pca@stdev[1:cds_dims]^2)
  }

  return(cds)
}

CollapseModules <- function (
  modules,
  keep.percent = 50
) {
  collapsed <- lapply(modules, function(x) Reduce(intersect, x))
  tallied <- lapply(modules, function(x) {
    top_genes <- as.data.frame(table(unlist(x)))
    top_genes_filtered <- as.character(top_genes[top_genes$Freq >= keep.percent/100*max(top_genes$Freq), ] %>% pull(Var1))
  })
  len <- sapply(tallied, length)
  len <- max(len) - len
  df <- as.data.frame(mapply(function(x, y) c(x, rep(NA, y)), tallied, len))
  colnames(df) <- paste0("M", seq(1:ncol(df)))
  return(df)
}



#' Plot proportions between 2 variables using stacked bar plot
#'
#' @export
PlotBarProportions <- function(object, xvar, fillvar) {
  res_props <- object[[]] %>%
    dplyr::group_by(!!sym(xvar), !!sym(fillvar)) %>%
    dplyr::summarise(n = n()) %>%
    dplyr::mutate(prop = n / sum(n))
  plot <- ggplot() +
    geom_bar(aes_string(x = xvar, y = "prop", fill = fillvar), width = 0.8,
      data = res_props, stat = "identity", color = "black") +
    scale_x_discrete(limits = rev(levels(res_props[, xvar, drop = TRUE]))) +
    scale_y_continuous(expand = c(0, 0))
  return(plot)
}

LoadModulesUMAP <- function(path, num.random.modules) {
  modules <- read_csv(path)
  modules_list <- split(modules$id, modules$module)
  num_modules <- length(modules_list)
  names(modules_list) <- paste0("M", names(modules_list))
  all_genes <- unname(unlist(modules_list))

  for (i in seq_along(1:num.random.modules)) {
    message(glue("Generating random gene list {i}..."))
    random_genes <- sample(all_genes, 50)
    modules_list[[num_modules + i]] <- random_genes
  }
  names(modules_list)[c((num_modules+1):(num_modules+num.random.modules))] <- paste0("C", as.character(c(1:num.random.modules)))
  return(modules_list)
}

ConvertConsensusModules <- function(df, num.random.modules = 3) {
  df <- unlist(lapply(split(df, seq_along(df[, 1])), as.list), recursive = FALSE)
  list <- lapply(df, function(x) x[!is.na(x)])
  names(list) <- paste0("M", seq(length(list)))
  all_genes <- unname(unlist(list))
  num_modules <- length(list)
  for (i in seq_along(1:num.random.modules)) {
    message(glue("Generating random gene list {i}..."))
    random_genes <- sample(all_genes, round(mean(sapply(list, length))))
    list[[num_modules + i]] <- random_genes
  }
  names(list)[c((num_modules+1):(num_modules+num.random.modules))] <- paste0("C", as.character(c(1:num.random.modules)))
  return(list)
}



GetNamedFactorLabels <- function(object, column) {
  cell_labels <- as.factor(object@meta.data[, column, drop = TRUE])
  names(cell_labels) <- rownames(object@meta.data)
  return(cell_labels)
}

GetNamedCharLabels <- function(object, column) {
  cell_labels <- as.character(object@meta.data[, column, drop = TRUE])
  names(cell_labels) <- rownames(object@meta.data)
  return(cell_labels)
}

Reclassify <- function(
  object,
  reduction,
  label_col,
  seed = 1,
  percent = 1,
  ensembles = 3,
  prob = FALSE,
  balance = FALSE,
  sampling.iter = 3
) {

  reduced.matrix = t(Embeddings(object, reduction = reduction))
  cell.labels = GetNamedCharLabels(object, column = label_col)

  models <- list()
  num_labels <- length(unique(cell.labels))

  for (l in 1:ensembles) {

    # initialize each ensemble
    message(glue::glue("Generating models for ensemble... {l}"))
    set.seed(seed + l)

    # set initial training and labels
    X <- reduced.matrix
    Y <- cell.labels
    model <- c()
    prob_mat <- c()

    for (i in 1:sampling.iter) {
      message(glue::glue("Sampling iteration... {i}"))

      # train SVM model
      tmp <- t(X)
      rownames(tmp) <- NULL
      model <- e1071::svm(tmp, factor(Y), probability = TRUE, kernel = "linear")
      prob_mat <- attr(predict(model, t(reduced.matrix), decision.values = FALSE,
        probability = TRUE), "probabilities")

      # reset training and label data
      X <- c()
      Y <- c()

      # for the number of labels
      for (j in 1:ncol(prob_mat)) {

        vote_class <- prob_mat[cell.labels == colnames(prob_mat)[j], ]
        idx <- c()

        if (balance == FALSE) {

          idx <- sample(1:nrow(vote_class),
            size = nrow(vote_class) * percent, replace = TRUE, prob = vote_class[, j])

        } else {

          sample_size <- round(median(table(cell.labels)))

          if (nrow(vote_class) > sample_size) {

            idx <- sample(1:nrow(vote_class),
              size = sample_size * percent, replace = TRUE, prob = vote_class[, j])

          } else {

            idx <- sample(1:nrow(vote_class),
              size = nrow(vote_class) * percent, replace = TRUE, prob = vote_class[, j])
          }
        }

        X <- cbind(X, reduced.matrix[, rownames(vote_class)[idx]])
        Y <- c(Y, cell.labels[rownames(vote_class)[idx]])
      }
    }

    models[[l]] <- model

  }

  message(glue::glue("Finalizing predictions..."))
  predict_mat <- matrix(0, nrow = ncol(reduced.matrix), ncol = length(table(cell.labels)))

  final <- c()
  for (l in 1:ensembles) {
    tmp <- attr(predict(models[[l]], newdata = t(reduced.matrix),
      probability = TRUE), "prob")[, names(table(cell.labels))]
    predict_mat <- predict_mat + tmp
  }

  if (prob == TRUE) {
    final <- apply(predict_mat, 1, max)
    names(final) <- names(table(cell.labels))[apply(predict_mat, 1, which.max)]
  }

  else {
    final <- names(table(cell.labels))[apply(predict_mat, 1, which.max)]
  }

  return(list(final = final, models = models, prob = predict_mat))
}

StdIntegrateHarmony <- function (object, batch, features, vars.to.regress = NULL, npcs = 50,
  ...)
{
  object <- object %>% NormalizeData() %>% FindVariableFeatures(selection.method = "vst",
    nfeatures = 3000) %>% ScaleData(vars.to.regress = vars.to.regress) %>%
    RunPCA(npcs = npcs, verbose = FALSE)
  integrated <- RunHarmony(object, "batch", ...)
  return(integrated)
}

GenerateCellXGeneFile <- function(object, object_path) {
  library(sceasy)
  library(reticulate)
  use_virtualenv("~/env")
  loompy <- reticulate::import('loompy')

  message(glue::glue("Writing object to \n{object_path}"))
  sceasy::convertFormat(object, from = "seurat", to = "anndata",
    outFile = object_path)
}

DeterminekBETBatchEffects <- function(
  object,
  pca_dims,
  iteration.size = 1000,
  iterations = 3,
  batch_reduction = "pca"
) {

  subset_size <- iteration.size/length(Cells(object))
  observed_mean <- c()
  silo <- c()

  for (j in seq(iterations)) {
    message(glue("Iteration {j}"))

    cells <- list()
    for (i in unique(object$batch)) {
      #message(glue("Sampling batch {i}"))
      batch_size <- length(object$batch[object$batch == i])
      subset_id <- sample.int(n = batch_size, size = floor(subset_size * batch_size))
      cells[[i]] <- names(object$batch[object$batch == i])[subset_id]
    }
    cells <- unlist(cells)
    message(glue("{length(cells)} cells sampled"))
    data <- Embeddings(object, reduction = batch_reduction)
    data <- data[cells, ]
    batch <- object$batch[cells]
    batched <- kBET::kBET(data, batch, plot = FALSE, verbose = FALSE,
      do.pca = FALSE, k0 = round(0.2*nrow(data)))

    observed_mean <- c(observed_mean, batched$summary$kBET.observed[1])
    silo <- c(silo, summary(cluster::silhouette(as.numeric(as.factor(batch)),
      dist(data[, seq_len(pca_dims)])))$avg.width)
  }
  return(list(observed_mean, silo))
}

ClusterModularity <- function(graph, clusters, get.weights=FALSE, get.values=NULL, as.ratio=FALSE) {
  by.clust <- split(seq_along(clusters), clusters)
  uclust <- names(by.clust)
  nclust <- length(uclust)

  mod.mat <- matrix(0, nclust, nclust)
  dimnames(mod.mat) <- list(uclust, uclust)
  clust.total <- numeric(nclust)
  names(clust.total) <- uclust

  # Calculating the observed weight within/between clusters.
  for (x in seq_along(by.clust)) {
    current <- by.clust[[x]]
    for (y in seq_len(x)) {
      other <- by.clust[[y]]
      grmat <- graph[current,other,drop=FALSE]
      grsum <- sum(grmat)

      if (x==y) {
        old.diag <- sum(diag(grmat))
        diag(grmat) <- 0

        self.sum <- sum(grmat)
        if (!igraph::is.directed(graph)) {
          # Only count undirected edges within a cluster once.
          self.sum <- self.sum/2
        } else {
          # Need to count directed edges between different nodes twice.
          grsum <- grsum + self.sum
        }

        # Self-edges contribute twice to total node weight,
        # according to igraph::modularity.
        grsum <- grsum + old.diag

        mod.mat[x,y] <- self.sum + old.diag
      } else {
        if (is.directed(graph)) {
          grsum <- grsum + sum(graph[other,current])
        }
        mod.mat[y,x] <- grsum
      }

      # If x==y, this is equivalent to adding edge weights for each node twice.
      # THIS IS DELIBERATE; the total per-node weight is that for all edges
      # involving each node, so edges between nodes in the same cluster are
      # allowed to be double-counted here when aggregating node weights per cluster.
      clust.total[x] <- clust.total[x] + grsum
      if (x!=y) {
        clust.total[y] <- clust.total[y] + grsum
      }
    }
  }

  # Calcuating the expected weight if they were randomly distributed.
  # Note some effort is involved in 'folding' it into an upper triangular.
  total.weight <- sum(mod.mat)
  clust.prop <- clust.total/sum(clust.total)
  expected.mat <- tcrossprod(clust.prop) * total.weight

  expected.mat[lower.tri(expected.mat)] <- 0
  old.diag <- diag(expected.mat)
  expected.mat <- expected.mat * 2
  diag(expected.mat) <- old.diag

  get.weights <- .switch_arg_names(get.values, get.weights)
  if (get.weights) {
    list(observed=mod.mat, expected=expected.mat)
  } else if (as.ratio) {
    mod.mat/expected.mat
  } else {
    1/total.weight * (mod.mat - expected.mat)
  }
}

DotPlot_ol <- function (object, assay = NULL, features, cols = c("lightgrey",
  "blue"), col.min = -2.5, col.max = 2.5, dot.min = 0, dot.scale = 6,
  group.by = NULL, split.by = NULL, scale.by = "radius", scale.min = NA,
  scale.max = NA)
{
  assay <- assay %||% DefaultAssay(object = object)
  DefaultAssay(object = object) <- assay
  scale.func <- switch(EXPR = scale.by, size = scale_size,
    radius = scale_radius, stop("'scale.by' must be either 'size' or 'radius'"))
  data.features <- FetchData(object = object, vars = features)
  data.features$id <- if (is.null(x = group.by)) {
    Idents(object = object)
  }
  else {
    object[[group.by, drop = TRUE]]
  }
  if (!is.factor(x = data.features$id)) {
    data.features$id <- factor(x = data.features$id)
  }
  id.levels <- levels(x = data.features$id)
  data.features$id <- as.vector(x = data.features$id)
  if (!is.null(x = split.by)) {
    splits <- object[[split.by, drop = TRUE]]
    if (length(x = unique(x = splits)) > length(x = cols)) {
      stop("Not enought colors for the number of groups")
    }
    cols <- cols[1:length(x = unique(x = splits))]
    names(x = cols) <- unique(x = splits)
    data.features$id <- paste(data.features$id, splits,
      sep = "_")
    unique.splits <- unique(x = splits)
    id.levels <- paste0(rep(x = id.levels, each = length(x = unique.splits)),
      "_", rep(x = unique(x = splits), times = length(x = id.levels)))
  }
  data.plot <- lapply(X = unique(x = data.features$id), FUN = function(ident) {
    data.use <- data.features[data.features$id == ident,
      1:(ncol(x = data.features) - 1), drop = FALSE]
    avg.exp <- apply(X = data.use, MARGIN = 2, FUN = function(x) {
      return(mean(x = expm1(x = x)))
    })
    pct.exp <- apply(X = data.use, MARGIN = 2, FUN = PercentAbove,
      threshold = 0)
    return(list(avg.exp = avg.exp, pct.exp = pct.exp))
  })
  names(x = data.plot) <- unique(x = data.features$id)
  data.plot <- lapply(X = names(x = data.plot), FUN = function(x) {
    data.use <- as.data.frame(x = data.plot[[x]])
    data.use$features.plot <- rownames(x = data.use)
    data.use$id <- x
    return(data.use)
  })
  data.plot <- do.call(what = "rbind", args = data.plot)
  if (!is.null(x = id.levels)) {
    data.plot$id <- factor(x = data.plot$id, levels = id.levels)
  }
  avg.exp.scaled <- sapply(X = unique(x = data.plot$features.plot),
    FUN = function(x) {
      data.use <- data.plot[data.plot$features.plot ==
          x, "avg.exp"]
      data.use <- scale(x = data.use)
      data.use <- MinMax(data = data.use, min = col.min,
        max = col.max)
      return(data.use)
    })
  avg.exp.scaled <- as.vector(x = t(x = avg.exp.scaled))
  if (!is.null(x = split.by)) {
    avg.exp.scaled <- as.numeric(x = cut(x = avg.exp.scaled,
      breaks = 20))
  }
  data.plot$avg.exp.scaled <- avg.exp.scaled
  data.plot$features.plot <- factor(x = data.plot$features.plot,
    levels = rev(x = features))
  data.plot$pct.exp[data.plot$pct.exp < dot.min] <- NA
  data.plot$pct.exp <- data.plot$pct.exp * 100
  if (!is.null(x = split.by)) {
    splits.use <- vapply(X = strsplit(x = as.character(x = data.plot$id),
      split = "_"), FUN = "[[", FUN.VALUE = character(length = 1L),
      2)
    data.plot$colors <- mapply(FUN = function(color, value) {
      return(colorRampPalette(colors = c("grey", color))(20)[value])
    }, color = cols[splits.use], value = avg.exp.scaled)
  }
  color.by <- ifelse(test = is.null(x = split.by), yes = "avg.exp.scaled",
    no = "colors")
  if (!is.na(x = scale.min)) {
    data.plot[data.plot$pct.exp < scale.min, "pct.exp"] <- scale.min
  }
  if (!is.na(x = scale.max)) {
    data.plot[data.plot$pct.exp > scale.max, "pct.exp"] <- scale.max
  }
  plot <- ggplot(data = data.plot, mapping = aes_string(x = "features.plot",
    y = "id")) + geom_point(mapping = aes_string(size = "pct.exp",
      fill = color.by), shape = 21, alpha = 0.9) + scale.func(range = c(0, dot.scale),
        limits = c(scale.min, scale.max)) + theme(axis.title.x = element_blank(),
          axis.title.y = element_blank()) + guides(size = guide_legend(title = "Percent Expressed")) +
    labs(x = "Features", y = ifelse(test = is.null(x = split.by),
      yes = "Identity", no = "Split Identity")) + cowplot::theme_cowplot()
  if (!is.null(x = split.by)) {
    plot <- plot + scale_fill_identity()
  }
  else if (length(x = cols) == 1) {
    plot <- plot + scale_fill_distiller(palette = cols)
  }
  else {
    plot <- plot + scale_fill_gradient(low = cols[1], high = cols[2])
  }
  if (is.null(x = split.by)) {
    plot <- plot + guides(fill = guide_colorbar(title = "Average Expression"))
  }
  return(plot)
}

LabelCyclingCells <- function(object, grouping_var = "cell_subset_mnp") {

  object@meta.data[, grep("s_score", colnames(object@meta.data))] <- NULL
  object@meta.data[, grep("g2m_score", colnames(object@meta.data))] <- NULL
  object@meta.data[, grep("phase", colnames(object@meta.data))] <- NULL

  # calculate cell cycle scores
  object <- Seurat::CellCycleScoring(object,
    s.features = cc.genes.updated.2019$s.genes,
    g2m.features = cc.genes.updated.2019$g2m.genes,
    set.ident = FALSE, search = FALSE)
  head(object[[]])
  colnames(object@meta.data)[grep("S.Score", colnames(object@meta.data))] <- "s_score"
  colnames(object@meta.data)[grep("G2M.Score", colnames(object@meta.data))] <- "g2m_score"
  colnames(object@meta.data)[grep("Phase", colnames(object@meta.data))] <- "phase"
  object$cc_diff <- object$s_score - object$g2m_score

  # summarize scores based on clusters
  summary <- object@meta.data %>%
    group_by(!!sym(grouping_var)) %>%
    summarize(median = median(cc_diff), mean = mean(cc_diff), sd = sd(cc_diff))

  # calculate outlier cell cutoffs
  cc_diff_norm_mean <- mean(summary$mean[  !(summary$mean > (median(summary$mean) + 2*mad(summary$mean))) | ! (summary$mean <(median(summary$mean) - 2*mad(summary$mean)))])
  cc_diff_norm_sd <- mean(summary$sd[  !(summary$mean > (median(summary$mean) + 2*mad(summary$mean))) | ! (summary$mean <(median(summary$mean) - 2*mad(summary$mean)))])
  upper <- cc_diff_norm_mean + 2*cc_diff_norm_sd
  lower <- cc_diff_norm_mean - 2*cc_diff_norm_sd

  # assign labels
  object$cc_diff_label <- ifelse(object$cc_diff > upper | object$cc_diff < lower, "Cycling", "Non-cycling")

  return(object)
}

ScoreLineages <- function(object, grouping_var = "base_clustering", genes) {

  # hla_genes <- c("HLA-DQB1", "HLA-DRB1", "HLA-DQA1", "HLA-DRB5", "HLA-DOA",
  #   "HLA-DRA", "HLA-C", "HLA-E", "HLA-G", "HLA-F", "HLA-H", "HLA-A",
  #   "HLA-F-AS1", "HLA-DPB1", "HLA-DRB6", "HLA-DPA1", "HLA-DQB2",
  #   "HLA-B", "HLA-DPA3", "HLA-DQA2", "HLA-DMB", "HLA-L", "HLA-V",
  #   "HLA-DRB9", "HLA-DMA", "HLA-J", "HLA-DPB2", "HLA-DOB", "HLA-K",
  #   "HLA-DPA2", "HLA-DQB1-AS1", "HLA-S", "HLA-W")
  # mac_genes <- c("C1QA", "C1QB", "C1QC", "CST3", "HLA-DRA", "LYZ", "TYROBP", "CD68", "CD33",
  #   "ITGAX", "CD14", "ITGAM", "MERTK", "CD40", "CD63", "CD80", "CD86",
  #   "SIGLEC1", "CD163", "TLR2", "TLR4", "MSR1", "MRC1", "MRC1", "CD209",
  #   "CD32", hla_genes)
  # mon_genes <- c("S100A8", "S100A9", "VCAN", "TIMP1", "CD14", "FCGR3A", "THBS1")
  # dc_genes <- c(hla_genes, "CLEC9A", "CLEC7A", "ITGAM", "ANPEP", "CD33", "CD80",
  #   "CD83", "CD86", "CCR7", "CD1C", "CD1B")
  # b_genes <- c("CD19", "MS4A1", "CD22", "CD70", "CD79A", "CD79B", "IGKC", "IGHA1",
  #   "IGL", "IGHM", "IGHE", "IGHG1")
  # t_genes <- c("CD2", "CD3D", "CD3E", "CD3G", "CD4", "CD5", "CD7", "FOXP3",
  #   "IL2RA", "IL7R", "CD8A", "CD8B", "CD27", "CD28", "SELL", "CD69",
  #   "CTLA4", "PD1", "CTLA4", "BTLA", "CD40LG", "TNFRSF4", "TNFRSF9")
  # nk_genes <- c("GZMB", "GZMA", "GNLY", "KLRK1", "NCR1", "NCR2", "NCR3", "KLRG1",
  #   "IL2RB", "KLRD1", "KLRB1", "CD244", "NCAM1", "CD226", "NKG7")
  # mas_genes <- c("CPA3", "TPSB2", "TPSAB1", "KIT")
  # neu_genes <- c("CEACAM8", "FUT8", "SIGLEC8", "CCR3", "FCGR3B", "CSF3R", "IL5RA",
  #   "CCR3", "FUT4", "CD66", "BST1", "FCAR", "CD93", "MME")
  # fib_genes <- c(rownames(object)[grep("^COL[[:digit:]]", rownames(object))],
  #   "FBLN2", "MGP", "CFD", "ACTA2", "ADAMTS1", "C1R", "C3", "CALD1", "CCDC80",
  #   "DCN", "FBLN1", "FILIP1L", "IGF1", "IGFBP7",
  #   "MYL9", "RGS5", "SPARCL1", "TAGLN", "TPM1")
  # epi_genes <- c("AKAP9", "AQP3", "CAPS", "CCDC146", "DHCR24", "DNAAF1", "HSP90AA1",
  #   "LRRIQ1", "MUC1", "NAPSA", "RSPH1", "SFTPA1", "SFTPA2", "SFTPB",
  #   "SFTPC", "SLC34A2", "TMEM190", "TPPP3", "TSPAN1")
  # endo_genes <- c("AQP1", "CLDN5", "GNG11", "PALMD", "PCDH17", "SPARCL1", "TCF4",
  #   "TM4SF1")
  # lineage_genes <- list(mac_genes, mon_genes, dc_genes,
  #   b_genes, t_genes, nk_genes,
  #   mas_genes, neu_genes,
  #   fib_genes, epi_genes, endo_genes)
  lineage_genes <- lapply(genes, function(x) intersect(x, rownames(object)))
  object@meta.data[, grep("_score", colnames(object@meta.data))] <- NULL
  object <- AddModuleScore(object, features = lineage_genes, nbin = 20, ctrl = 30,
    seed = 1, name = "lineage_scores", search = FALSE)
  colnames(object@meta.data)[grep("lineage_scores", colnames(object@meta.data))] <- names(genes)
  return(object)

  # # summarize scores based on clusters
  # summary <- object@meta.data %>%
  #   group_by(!!sym(grouping_var)) %>%
  #   summarize(median = median(mnp_diff_score), mean = mean(mnp_diff_score), sd = sd(mnp_diff_score))
  #
  # # calculate outlier cell cutoffs
  # median(summary$mean)
  # mad(summary$mean)
  # median(summary$mean) + 1.5*mad(summary$mean)
  # mnp_diff_norm_mean <- mean(summary$mean[!(summary$mean > (median(summary$mean) + 2*mad(summary$mean)))])
  # mnp_diff_norm_sd <- mean(summary$sd[!(summary$mean > (median(summary$mean) + 2*mad(summary$mean)))])
  # upper <- mnp_diff_norm_mean + 1*mnp_diff_norm_sd
  #
  # object$mnp_diff_label <- ifelse(object$mnp_diff_score > upper, "MNP", "Non-MNP")
  #
  # return(object)
}

AssignLabels <- function(object, scores) {
  for (i in seq_along(scores)) {
    score <- scores[i]
    label <- paste0(gsub( "_.*$", "", score), "_label")
    mm <- mixtools::normalmixEM(object[[score, drop = TRUE]])
    #plot(mm, which = 2)
    threshold <- mm$mu[1] + 2*mm$sigma[1]
    passing_clusters <- object@meta.data %>%
      select(base_clustering, !!sym(score)) %>%
      group_by(base_clustering) %>%
      summarise(median = median(!!sym(score))) %>%
      filter(median > threshold) %>%
      pull(base_clustering)
    passing_clusters <- as.character(passing_clusters)
    object[[label, drop = TRUE]] <- FALSE
    object[[label, drop = TRUE]][object$base_clustering %in% passing_clusters] <- TRUE
  }
  return(object)
}

DimPlotPlus <- function(object, ...) {
  dimplot_plus <- Seurat::DimPlot(object, ...) +
    colorspace::scale_color_discrete_qualitative("Dark 3", name = "") +
    labs(x = "UMAP1", y = "UMAP2") +
    theme_minimal(base_size = 18) +
    theme(
      axis.title = element_text(face = "bold", margin = margin(6, 6, 6, 6, "pt")),
      panel.grid = element_blank(),
      axis.text = element_blank(),
      plot.margin = unit(rep(12, 4), "pt"),
      plot.subtitle = element_text(face = "bold")
    )
  return(dimplot_plus)
}

RGB2Hex <- function(r,g,b) rgb(r, g, b, maxColorValue = 255)

gammaNormMix <- function(data, thresh = 1e-07, maxiter = 10000,
  removeZeroes = TRUE, plot = TRUE, hist = TRUE, hist_col = "light cyan",
  verbose = FALSE, forceExponential = FALSE, calculateAreaDifference = FALSE,
  minDataPoints = 5, onlyAddCurves = FALSE, addContextData = FALSE,
  contextData = NULL) {


  # fitting a 2 component normal and gamma mixture model

  # add other data to fit the model with as well, but only
  # return the classification for those we're interested in

  if (addContextData) {
    nOriginal = length(data)
    data <- c(data, contextData)
  }

  # assume all values exactly zero already belong to the gamma
  # comp and remove them from the EM algorithm

  if (removeZeroes) {
    nonZeroInd = which(data > 0)
    x = data[nonZeroInd]
  } else {
    x = data
  }

  if (length(x) < minDataPoints) {
    if (verbose)
      cat("Not enough data points to fit mixture model!")
    return(NA)
  }

  # initiate
  n = length(x)
  z = stats::rbinom(n, 1, 0.5)
  if(sum(z) == 0){z[1] = 1} ## Break out of a sequence of zeroes error
  z_iter = z
  mu = -100
  mu_iter = 10
  sig2 = -100
  sig2_iter = 0
  alpha = -100
  alpha_iter = 1
  beta = -100
  beta_iter = 1
  rho = -100
  rho_iter = 0.5
  niter = 0

  while (any(c(abs(mu - mu_iter) > thresh, abs(sig2 - sig2_iter) >
      thresh, abs(alpha - alpha_iter) > thresh, abs(beta -
          beta_iter) > thresh, abs(rho - rho_iter) > thresh)) &
      (niter < maxiter)) {

    # save old parameters
    mu = mu_iter
    sig2 = sig2_iter
    alpha = alpha_iter
    beta = beta_iter
    rho = rho_iter
    if (forceExponential)
      alpha_iter = 1

    niter = niter + 1

    # M step
    mu_iter = sum(z_iter * x)/sum(z_iter)
    sig2_iter = sum(z_iter * (x - mu_iter) * (x - mu_iter))/sum(z_iter)
    if (sig2_iter <= 0 | is.na(sig2_iter))
      sig2_iter = 1e-11
    beta_iter = alpha_iter * sum(1 - z_iter)/sum((1 - z_iter) *
        x)
    if (beta_iter <= 0 | is.na(beta_iter))
      beta_iter = 3
    if (!forceExponential) {
      alpha_iter = distr::igamma(sum((log(beta_iter) +
          log(x)) * (1 - z_iter))/sum(1 - z_iter))
    }
    if (alpha_iter > 150 | is.na(alpha_iter))
      alpha_iter = 150
    rho_iter = sum(z_iter)/n


    # E step
    eta_iter = -0.5 * log(2 * pi * sig2_iter) - ((x - mu_iter) *
        (x - mu_iter))/(2 * sig2_iter) - alpha_iter * log(beta_iter) +
      log(gamma(alpha_iter)) - (alpha_iter - 1) * log(x) +
      beta_iter * x + log(rho_iter/(1 - rho_iter))
    z_iter = 1/(1 + exp(-eta_iter))

    if (verbose)
      cat(niter, mu_iter, sqrt(sig2_iter), alpha_iter,
        beta_iter, rho_iter, "\n")
  }


  ll <- sum(log(rho_iter * stats::dnorm(x, mu_iter, sqrt(sig2_iter)) +
      (1 - rho_iter) * stats::dgamma(x, shape = alpha_iter,
        rate = beta_iter)))


  xg <- seq(0, max(x) + 1, length.out = 300)
  c1g <- rho_iter * stats::dnorm(xg, mu_iter, sqrt(sig2_iter))

  c2g <- (1 - rho_iter) * stats::dgamma(xg, shape = alpha_iter,
    rate = beta_iter)
  fg <- rho_iter * stats::dnorm(xg, mu_iter, sqrt(sig2_iter)) +
    (1 - rho_iter) * stats::dgamma(xg, shape = alpha_iter,
      rate = beta_iter)

  if (plot) {
    if (hist) {
      hist(x, probability = TRUE, col = hist_col, breaks = 50,
        main = NA, xlab = NA, ylab = "Density (zeroes removed)",
        ylim = c(0, 0.6), xlim = c(0, 20))
    }
    if (!onlyAddCurves) {
      graphics::lines(stats::density(x, from = 0), lty = 2,
        lwd = 2, col = scales::alpha("darkgrey", 0.6))
    }
    graphics::lines(xg, c1g, col = scales::alpha("red", 0.6), lwd = 2)  #Normal Lines
    graphics::lines(xg, c2g, col = scales::alpha("blue", 0.6), lwd = 2)  #Gamma lines
    graphics::lines(xg, fg, col = scales::alpha("black", 0.6), lwd = 2)  #Mixture model line

    if (onlyAddCurves)
      return(list(xg = xg, c1g = c1g, c2g = c2g, fg = fg))
  }
  if (calculateAreaDifference) {
    f1 <- stats::approxfun(xg, (stats::approxfun(stats::density(x,
      from = 0)))(xg) - fg)
    # piecewise linear function
    f2 <- function(x) abs(f1(x))
    # take the positive value
    AreaDifference = stats::integrate(f2, min(x[x != 0]),
      max(x))$value
  } else {
    AreaDifference = NULL
  }

  if (removeZeroes) {
    z = rep(0, length(data))
    z[nonZeroInd] <- z_iter
  } else {
    z = z_iter
  }

  # force prob expression values above the max to stay the same
  # value
  maxdata = data[which.max(z)]
  z[which(data > maxdata)] <- max(z)



  if (addContextData) {
    z <- z[seq_len(nOriginal)]
  }
  if (plot) {
    if (addContextData) {
      graphics::points(data[seq_len(nOriginal)], z * 0,
        pch = "|", cex = 1, col = scales::alpha(grDevices::rgb(z,
          0, 1 - z), 0.4))
    } else {
      graphics::points(data, z * 0, pch = "|", cex = 1,
        col = scales::alpha(grDevices::rgb(z, 0, 1 - z), 0.4))
    }
  }
  model_bic <- bic(ll, n, 5)
  model_aic <- aic(ll, 5)
  model_icl_bic <- icl_bic(ll, z, n, 5)
  return(list(probExpressed = z, propExpressed = n * rho_iter/length(data),
    numExpressed = length(which(z > 0.5)), mu = mu_iter,
    sd = sqrt(sig2_iter), alpha = alpha_iter, beta = beta_iter,
    rho = rho_iter, niter = niter, loglik = ll, BIC = model_bic,
    AIC = model_aic, ICL_BIC = model_icl_bic, AreaDifference = AreaDifference))
}

bic <- function(loglik, n, p) {
  return(-2 * loglik + p * log(n))
}


aic <- function(loglik, p) {
  return(-2 * loglik + 2 * p)
}  #Tend to fit more component


icl_bic <- function(loglik, postprob, n, p) {
  postprob <- postprob[postprob > 0]
  EN = -sum(postprob * log(postprob))
  return(-2 * loglik + 2 * EN + p * log(n))
}

GenerateAllModules <- function(
  cds,
  genes,
  res.start = 0.01,
  res.end = 1,
  num.res = 10,
  iterations,
  ...
) {
  iterations = as.list(seq(1, iterations, 1))
  res_seq <- exp(seq(log(res.start), log(res.end), length.out = num.res))
  modules <- future_lapply(iterations, FUN = function(x) {
    modules_iter <- monocle3::find_gene_modules(cds[genes, ],
      resolution = res_seq, umap.n_neighbors = 30, umap.fast_sgd = TRUE,
      louvain_iter = 1, cores = 4, random_seed = as.numeric(Sys.time()), ...)
    return(modules_iter)
  })
  return(modules)
}

ComputePartitions <- function(g, optim_res, qval_thresh = 0.05, verbose = FALSE)
{
  cell_membership <- as.factor(igraph::membership(optim_res))
  membership_matrix <- Matrix::sparse.model.matrix(~cell_membership +
      0)
  num_links <- Matrix::t(membership_matrix) %*% igraph::as_adjacency_matrix(g) %*%
    membership_matrix
  diag(num_links) <- 0
  louvain_modules <- levels(cell_membership)
  edges_per_module <- Matrix::rowSums(num_links)
  total_edges <- sum(num_links)
  theta <- (as.matrix(edges_per_module)/total_edges) %*% Matrix::t(edges_per_module/total_edges)
  var_null_num_links <- theta * (1 - theta)/total_edges
  num_links_ij <- num_links/total_edges - theta
  cluster_mat <- pnorm_over_mat(as.matrix(num_links_ij), var_null_num_links)
  num_links <- num_links_ij/total_edges
  cluster_mat <- matrix(stats::p.adjust(cluster_mat), nrow = length(louvain_modules),
    ncol = length(louvain_modules))
  sig_links <- as.matrix(num_links)
  sig_links[cluster_mat > qval_thresh] = 0
  diag(sig_links) <- 0
  cluster_g <- igraph::graph_from_adjacency_matrix(sig_links,
    weighted = T, mode = "undirected")
  list(cluster_g = cluster_g, num_links = num_links, cluster_mat = cluster_mat)
}

CalculateModuleRelativeFC <- function(object,
  score.features,
  ctrl.features,
  module.name,
  group.by ) {
  ctrl_features <- intersect(rownames(object), ctrl.features)
  if (length(ctrl_features) != length(ctrl.features)) {
    usethis::ui_warn("Removed {length(ctrl.features)-length(ctrl_features)} genes, not in object")
  }
  score_features <- intersect(rownames(object), score.features)
  if (length(score_features) != length(score.features)) {
    usethis::ui_warn("Removed {length(score.features)-length(score_features)} genes, not in object")
  }
  assay_data <- GetAssayData(object = object, slot = "data", assay = "RNA")
  ctrl_scores <- Matrix::colMeans(x = assay_data[ctrl_features, ])
  feature_scores <- Matrix::colMeans(x = assay_data[score_features, ])
  return_scores <- (feature_scores-ctrl_scores)/(ctrl_scores)
  return_scores <- data.frame(barcode = names(return_scores),
    score = return_scores,
    group = object@meta.data[, group.by, drop = TRUE],
    module = module.name,
    row.names = NULL, stringsAsFactors = FALSE)
  return(return_scores)
}

CalculateModuleExpression <- function(object,
  score.features,
  #ctrl.features,
  module.name,
  group.by,
  slot = "data") {

  score_features <- intersect(rownames(object), score.features)
  if (length(score_features) != length(score.features)) {
    usethis::ui_warn("Removed {length(score.features)-length(score_features)} genes, not in object")
  }
  assay_data <- GetAssayData(object = object, slot = slot, assay = "RNA")
  feature_scores <- Matrix::colMeans(x = assay_data[score_features, , drop = FALSE])
  return_scores <- feature_scores
  return_scores <- data.frame(barcode = names(return_scores),
    score = return_scores,
    group = object@meta.data[, group.by, drop = TRUE],
    module = module.name,
    row.names = NULL, stringsAsFactors = FALSE)
  return(return_scores)
}


GenerateNNs <- function(
  object_list,
  object_names,
  genes,
  projected = TRUE,
  distance = "cosine", k = 30, pruning = 1/15
){
  nn <- pblapply(object_names, function(x) {
    loadings <- Loadings(object_list[[x]], reduction = "harmony", projected = projected)
    loadings <- loadings[genes[genes %in% rownames(loadings)], ]
    if (distance == "cosine") {
      sim <- loadings / sqrt(rowSums(loadings * loadings))
      sim <- sim %*% t(sim)
      dist_matrix <- as.dist(1 - sim)
    }
    if (distance == "euclidean") {
      dist_matrix <- dist(loadings, method = "euclidean")
    }
    nn <- Seurat::FindNeighbors(dist_matrix, k.param = 30, prune.SNN = 1/15,
      nn.method = "rann", verbose = TRUE)
    return(nn)
  })
  return(nn)
}

LeidenClustering <- function(
  g,
  resolution_parameters) {
  table_results <- data.frame(
    resolution_parameter = double(),
    quality              = double(),
    modularity           = double(),
    significance         = double(),
    number_clusters      = integer())
  best_modularity <- -1
  best_result <- NULL
  best_resolution_parameter <- "No resolution"
  for(i in 1:length(resolution_parameters)) {
    cur_resolution_parameter <- resolution_parameters[i]
    cluster_result <- leidenbase::leiden_find_partition(g,
      partition_type = "CPMVertexPartition",
      initial_membership = NULL,
      edge_weights = NULL,
      node_sizes = NULL,
      seed = as.numeric(Sys.time()),
      resolution_parameter = cur_resolution_parameter,
      num_iter = 1,
      verbose = FALSE)
    table_results <- rbind(table_results, data.frame(
      resolution_parameter = cur_resolution_parameter,
      quality              = cluster_result[['quality']],
      modularity           = cluster_result[['modularity']],
      significance         = cluster_result[['significance']],
      cluster_count        = max(cluster_result[['membership']])))
  }
  return(table_results)
}

LeidenSNN <- function(
  snn,
  mode = "undirected",
  iterations = 10,
  res.start = 1E-5,
  res.end = 1E0,
  num.res = 20,
  num.modules.upper = 50,
  num.modules.lower = 5,
  verbose = TRUE,
  modularity.similarity = 0.99
) {

  use.genes <- rownames(snn)
  g <- igraph::graph_from_adjacency_matrix(snn, mode = mode, weighted = TRUE)
  resolution_parameter <-  signif(exp(seq(log(res.start), log(res.end), length.out = num.res)), 3)
  table_results <- LeidenClustering(g, resolution_parameter)
  first_round_table_results <- table_results
  lower_res <- min(table_results$resolution_parameter[table_results$cluster_count >= num.modules.lower])
  upper_res <- max(table_results$resolution_parameter[table_results$cluster_count <= num.modules.upper])
  resolution_parameter <- signif(seq(lower_res, upper_res, length.out = num.res), 3)
  ui_done("First round of clustering complete.
    Lower resolution: {resolution_parameter[1]},
    Upper resolution: {resolution_parameter[length(resolution_parameter)]}")

  iterations = as.list(seq(1, iterations, 1))
  leiden_results_total <- future.apply::future_lapply(iterations, FUN = function(x) {

    table_results <- LeidenClustering(g, resolution_parameter)
    max_modularity <- max(table_results$modularity)
    modularity_cutoff <- max_modularity*modularity.similarity
    max_resolution <- table_results$resolution_parameter[table_results$modularity == max_modularity]

    ui_info("{sum(table_results$modularity >= modularity_cutoff)} alternative resolutions")
    if (any(table_results$modularity >= modularity_cutoff)) {
      ui_info("Choosing non-optimal modularity...")
      final_resolution_parameter <- max(table_results$resolution_parameter[table_results$modularity >= modularity_cutoff])
    } else if (!any(table_results$modularity >= modularity_cutoff)) {
      final_resolution_parameter <- table_results$resolution_parameter[table_results$modularity == max_modularity]
    }
    ui_done("Second round of clustering complete.
    Final resolution: {final_resolution_parameter} vs. {max_resolution}")

    final_result <- leidenbase::leiden_find_partition(g,
      partition_type = "CPMVertexPartition",
      initial_membership = NULL,
      edge_weights = NULL,
      node_sizes = NULL,
      seed = as.numeric(Sys.time()),
      resolution_parameter = final_resolution_parameter,
      num_iter = 10,
      verbose = verbose)
    out_result <- list(membership = final_result[['membership']],
      modularity = final_result[['modularity']])
    names(out_result$membership) = use.genes
    leiden_results <- list(optim_res = out_result, second_results = table_results, first_results = first_round_table_results)
    return(leiden_results)
  })
  return(leiden_results_total)

}



FindGeneModules <- function(
  loadings,
  use.genes,
  k = 20,
  weight = FALSE,
  verbose = TRUE,
  res.start = 1E-5,
  res.end = 1E0,
  num.res = 30,
  num.modules.lower = 5,
  num.modules.upper = 50,
  iterations = 10,
  modularity.similarity = 0.95
) {

  # Check genes
  use.genes <- intersect(rownames(loadings), use.genes)
  ui_info("Utilizing {length(use.genes)} genes")

  loadings <- loadings[use.genes, ]
  tmp <- RANN::nn2(loadings, loadings, k, searchtype = "standard")
  neighbor_matrix <- tmp[[1]][, -1]
  dist_matrix <- tmp[[2]][, -1]
  links <- jaccard_coeff(neighbor_matrix, weight)
  links <- links[links[, 1] > 0, ]
  relations <- as.data.frame(links)
  colnames(relations) <- c("from", "to", "weight")
  relations$from <- use.genes[relations$from]
  relations$to <- use.genes[relations$to]
  g <- igraph::graph.data.frame(relations, directed = FALSE)
  resolution_parameter <-  signif(exp(seq(log(res.start), log(res.end), length.out = num.res)), 3)
  table_results <- LeidenClustering(g, resolution_parameter)
  first_round_table_results <- table_results
  lower_res <- min(table_results$resolution_parameter[table_results$cluster_count >= num.modules.lower])
  upper_res <- max(table_results$resolution_parameter[table_results$cluster_count <= num.modules.upper])
  resolution_parameter <- signif(seq(lower_res, upper_res, length.out = num.res), 3)
  ui_done("First round of clustering complete.
    Lower resolution: {resolution_parameter[1]},
    Upper resolution: {resolution_parameter[length(resolution_parameter)]}")

  iterations = as.list(seq(1, iterations, 1))
  leiden_results_total <- future.apply::future_lapply(iterations, FUN = function(x) {

    table_results <- LeidenClustering(g, resolution_parameter)
    max_modularity <- max(table_results$modularity)
    modularity_cutoff <- max_modularity*modularity.similarity
    max_resolution <- table_results$resolution_parameter[table_results$modularity == max_modularity]

    ui_info("{sum(table_results$modularity >= modularity_cutoff)} alternative resolutions")
    if (any(table_results$modularity >= modularity_cutoff)) {
      ui_info("Choosing non-optimal modularity...")
      final_resolution_parameter <- max(table_results$resolution_parameter[table_results$modularity >= modularity_cutoff])
    } else if (!any(table_results$modularity >= modularity_cutoff)) {
      final_resolution_parameter <- table_results$resolution_parameter[table_results$modularity == max_modularity]
    }
    ui_done("Second round of clustering complete.
    Final resolution: {final_resolution_parameter} vs. {max_resolution}")

    final_result <- leidenbase::leiden_find_partition(g,
      partition_type = "CPMVertexPartition",
      initial_membership = NULL,
      edge_weights = NULL,
      node_sizes = NULL,
      seed = as.numeric(Sys.time()),
      resolution_parameter = final_resolution_parameter,
      num_iter = 10,
      verbose = verbose)
    out_result <- list(membership = final_result[['membership']],
      modularity = final_result[['modularity']])
    names(out_result$membership) = use.genes
    leiden_results <- list(optim_res = out_result, second_results = table_results, first_results = first_round_table_results)
    return(leiden_results)
  })
  return(leiden_results_total)
}

SpreadMultipleColumns <- function(df, key, value) {
  # quote key
  keyq <- rlang::enquo(key)
  # break value vector into quotes
  valueq <- rlang::enquo(value)
  s <- rlang::quos(!!valueq)
  df %>% gather(variable, value, !!!s) %>%
    unite(temp, !!keyq, variable) %>%
    spread(temp, value)
}

#' Modifies plot axes text
#'
#' Increases font sizes and sets color to black for plot axes
#' @param size Size of text desired
#' @param text.color Color of axis text
#' @param title.color Color of title text
#' @param ... Other arguments to pass to theme()
#' @return ggplot2 theme
#' @usage IncreaseAxesSizes(size = 18, text.color = "black", title.color = "black", ...)
#' @export

IncreaseAxesSizes <- function(size = 18, text.color = "black", title.color = "black", ...) {
  AxesTheme <- theme(
    plot.title = element_text(size = size*1.25, color = title.color),
    axis.text.x = element_text(size = size, color = text.color),
    axis.text.y = element_text(size = size, color = text.color),
    axis.title.x = element_text(size = size*1.25, color = title.color, face = "bold", margin = margin(10, 0, 0, 0)),
    axis.title.y = element_text(size = size*1.25, color = title.color, face = "bold", margin = margin(0, 10, 0, 0)),
    # Validate the theme
    validate = TRUE,
    # Extra paramaters
    ...)
  return(AxesTheme)
}



#' Computes corrected expression ratios to shared vector of SEGs
CompareModuleExpr <- function(object, features, segs) {

  # set correct assay
  assay.old <- DefaultAssay(object = object)
  assay <- assay %||% assay.old
  DefaultAssay(object = object) <- assay
  assay.data <- GetAssayData(object = object)

  if (is.null(x = features)) {
    stop("Missing input feature list")
  }

  if (is.null(x = segs)) {
    stop("Missing input SEGs list")
  }

  features <- lapply(
    X = features,
    FUN = function(x) {
      missing.features <- setdiff(x = x, y = rownames(x = object))
      if (length(x = missing.features) > 0) {
        warning(
          "The following features are not present in the object: ",
          paste(missing.features, collapse = ", "),
          call. = FALSE,
          immediate. = TRUE
        )
      }
      return(intersect(x = x, y = rownames(x = object)))
    }
  )

  segs <- lapply(
    X = segs,
    FUN = function(x) {
      missing.segs <- setdiff(x = x, y = rownames(x = object))
      if (length(x = missing.segs) > 0) {
        warning(
          "The following features are not present in the object: ",
          paste(missing.segs, collapse = ", "),
          call. = FALSE,
          immediate. = TRUE
        )
      }
      return(intersect(x = x, y = rownames(x = object)))
    }
  )

  # calculate feature scores
  features.scores <- matrix(
    data = numeric(length = 1L),
    nrow = cluster.length,
    ncol = ncol(x = object)
  )
  for (i in 1:cluster.length) {
    features.use <- features[[i]]
    data.use <- assay.data[features.use, , drop = FALSE]
    features.scores[i, ] <- Matrix::colMeans(x = data.use)
  }

}

#' Computes features' correlation with technical covariates
ComputeTechCor <- function(object,
  features,
  grouping.var,
  tech.vars,
  method = "spearman"
) {

  groups <- FetchData(object = object, vars = grouping.var)
  tech <- FetchData(object, vars = tech.vars[1])

  cors <- pbmapply(FUN = function(x) {
    cells <- Cells(object)[which(x = groups == y)]
    feature <- FetchData(object, vars = x)
    cors <- cor.test(x = feature, y = data, alternative = "two.sided", method = method)

  }, x = features, y = levels(groups))


}

IterativeCluster <- function(
  object,
  graph_name,
  tree = "pca",
  #starting.res.addition = 0.5,
  dims = 20,
  res.start = 0.001,
  res.end = 1,
  num.res = 30,
  num.modules.lower = 5,
  num.modules.upper = 30,
  modularity.similarity = 0.99,
  #clustering.algorithm = 4,
  verbose = TRUE
) {

  # define a starting resolution by finding the Louvain resolution 3 steps higher than # Louvain clusters == # walktrap clusters
  # res_df <- object@meta.data[, grep("harmony_snn", colnames(object[[]]))]
  # num_walktrap <- length(unique(object[["base_walktrap", drop = TRUE]]))
  # num_louvain <- apply(res_df, 2, function(x) length(unique(x)))
  # num_louvain[which.min(abs(num_louvain - num_walktrap))]
  # start_res_name <- names(num_louvain)[which.min(abs(num_louvain - num_walktrap))]
  # start_res <- as.numeric(gsub("harmony_snn_res.", "", start_res_name)) + starting.res.addition
  # start_res_name <- paste0("harmony_snn_res.", as.character(start_res))
  # if (verbose) ui_done("Finding clusters at resolution {start_res}")

  # find optimal Leiden clustering based on maximal modularity
  g <- igraph::graph_from_adjacency_matrix(adjmatrix = object@graphs[[graph_name]], mode = "directed", weighted = TRUE)
  resolution_parameter <-  signif(exp(seq(log(res.start), log(res.end), length.out = num.res)), 3)
  table_results <- LeidenClustering(g, resolution_parameter)
  first_round_table_results <- table_results
  lower_res <- min(table_results$resolution_parameter[table_results$cluster_count >= num.modules.lower])
  upper_res <- max(table_results$resolution_parameter[table_results$cluster_count <= num.modules.upper])
  resolution_parameter <- signif(seq(lower_res, upper_res, length.out = num.res), 3)
  ui_done("First round of clustering complete.
    Lower resolution: {resolution_parameter[1]},
    Upper resolution: {resolution_parameter[length(resolution_parameter)]}")
  table_results <- LeidenClustering(g, resolution_parameter)
  max_modularity <- max(table_results$modularity)
  modularity_cutoff <- max_modularity*modularity.similarity
  max_resolution <- table_results$resolution_parameter[table_results$modularity == max_modularity]
  ui_info("{sum(table_results$modularity >= modularity_cutoff)} alternative resolutions")
  if (any(table_results$modularity >= modularity_cutoff)) {
    ui_info("Choosing non-optimal modularity...")
    final_resolution_parameter <- max(table_results$resolution_parameter[table_results$modularity >= modularity_cutoff])
  } else if (!any(table_results$modularity >= modularity_cutoff)) {
    final_resolution_parameter <- table_results$resolution_parameter[table_results$modularity == max_modularity]
  }
  ui_done("Second round of clustering complete.
    Final resolution: {final_resolution_parameter} vs. {max_resolution}")
  final_result <- leidenbase::leiden_find_partition(g,
    partition_type = "CPMVertexPartition",
    initial_membership = NULL,
    edge_weights = NULL,
    node_sizes = NULL,
    seed = as.numeric(Sys.time()),
    resolution_parameter = final_resolution_parameter,
    num_iter = 30,
    verbose = verbose)
  out_result <- list(membership = final_result[['membership']],
    modularity = final_result[['modularity']])
  names(out_result$membership) <- colnames(object@graphs[[graph_name]])
  leiden_results <- list(optim_res = out_result, second_results = table_results, first_results = first_round_table_results)

  ids <- GroupSingletons(leiden_results$optim_res$membership, object@graphs[[graph_name]], min.size = 5, group.singletons = TRUE, verbose = TRUE)

  object <- AddMetaData(object, ids, col.name = "optimal_leiden")
  object@misc$leiden_results <- leiden_results
  table(object@active.ident)

  # cluster at that resolution and build tree
  # object <- FindClusters(object, resolution = start_res, algorithm = clustering.algorithm,
  #   random.seed = 1, graph.name = glue("{method}_snn"))
  # object <- SetIdent(object, value = start_res_name)
  object <- SetIdent(object, value = "optimal_leiden")
  if (tree == "pca") object <- BuildClusterTree(object, verbose = FALSE, dims = 1:dims, reorder = TRUE, reorder.numeric = TRUE)
  if (tree == "expr") object <- BuildClusterTree(object, verbose = FALSE, reorder = TRUE, reorder.numeric = TRUE)

  # plot starting resolution
  a <- PlotAnnotatedClusters(object)
  if (verbose) a

  count <- 2
  repeat_count <- 0
  plots <- list(a)
  stop_signal = FALSE
  while (stop_signal == FALSE) {

    # identify marker genes for that resolution
    passing_markers <- IdentifyQualityMarkers(object)
    marker_counts <- table(passing_markers$group)[table(passing_markers$group) >= 3]
    if (verbose) print(marker_counts)
    if (exists("failed_clusters")) { previous_clusters <- failed_clusters }
    failed_clusters <- as.character(unique(Idents(object))[!(unique(Idents(object)) %in% names(marker_counts))])


    # merge failed clusters
    if (length(failed_clusters) > 0) {
      if (verbose) { ui_oops("Cluster {failed_clusters} failed to pass marker QC") }
      merged <- MergeUndefinedClusters(object, marker_counts)
      stop_signal <- merged$stop
      object <- merged$object
      if (tree == "pca") { object <- BuildClusterTree(object, verbose = FALSE, dims = 1:dims, reorder = TRUE, reorder.numeric = TRUE) }
      if (tree == "expr") { object <- BuildClusterTree(object, verbose = FALSE, reorder = TRUE, reorder.numeric = TRUE) }
      #if (verbose) PlotClusterTree(object)
      plots[[count]] <- PlotAnnotatedClusters(object)
      count <- count + 1
      if (exists("previous_clusters")) {
        if (identical(previous_clusters, failed_clusters)) { repeat_count <- repeat_count + 1 }
      }

    } else {
      ui_done("No failing clusters, skipping merge step")
      stop_signal = TRUE
    }

    if (repeat_count > 2) { stop_signal = TRUE }
    if (count > 11) { stop_signal = TRUE }
  }
  return(list(object = object, plots = plots))
}

IdentifyQualityMarkers <- function(
  object,
  tlogFC = log(1.25),
  tauc = 0.5,
  tpadj = 0.001,
  pseudocount = 1,
  tstm = 1.05
) {
  current_markers <- presto::wilcoxauc(object, assay = "data", seurat_assay = "RNA")
  current_markers <- current_markers %>% filter(logFC > tlogFC & auc > tauc & padj < tpadj)
  stm_markers <- CalculateSTMRatio(object, features = current_markers$feature, pseudocount.use = pseudocount, slot = "data")
  stm_markers <- stm_markers %>%
    gather(key = "cluster", value = "mean", -c(genes, stm)) %>%
    filter(stm >= tstm) %>%
    group_by(genes) %>%
    top_n(1, mean) %>%
    mutate(feature_group = paste0(genes, "_", cluster))
  passing_markers <- current_markers %>%
    mutate(feature_group = paste0(feature, "_", group)) %>%
    filter(feature_group %in% stm_markers$feature_group)
  return(passing_markers)
}

MergeUndefinedClusters <- function(
  object,
  marker_counts,
  marker_threshold = 20
) {
  clusters <- unique(object@active.ident)
  failed_clusters <- as.character(unique(Idents(object))[!(unique(Idents(object)) %in% names(marker_counts))])
  strong_clusters <- as.character(names(marker_counts)[marker_counts >= marker_threshold])
  missing_clusters <- as.character(sort(as.numeric(clusters[clusters %in% failed_clusters])))
  cluster_tree <- data.frame(node = object@tools$BuildClusterTree$edge[, 1], cluster = object@tools$BuildClusterTree$edge[, 2])
  nodes <- cluster_tree[cluster_tree$cluster %in% missing_clusters, "node"]
  merge_clusters <- cluster_tree[cluster_tree$node %in% nodes, ]

  if (!all(merge_clusters$cluster %in% clusters) | any(merge_clusters$cluster %in% strong_clusters)) {
    ui_oops("Removing singleton merging with node or singleton merging with strong cluster")
    merge_clusters <- merge_clusters[merge_clusters$cluster %in% clusters, ]
    merge_clusters <- merge_clusters[!(merge_clusters$cluster %in% strong_clusters), ]
  }
  if (nrow(merge_clusters) > 1) {
    rename_idents <- setNames(as.character(merge_clusters$node), merge_clusters$cluster)
    assertthat::assert_that(all(names(rename_idents) %in% levels(object@active.ident)))
    object <- RenameIdents(object, rename_idents)
    object@active.ident <- droplevels(object@active.ident)
    stop_signal = FALSE
  } else {
    usethis::ui_oops("No clusters to merge")
    stop_signal = TRUE
  }
  return(list(object = object, stop = stop_signal))
}

ExtractClusterCorr <- function(object,
  ident,
  features = Seurat::VariableFeatures(object),
  method = "spearman"
) {
  expr <- GetAssayData(object, slot = "data")
  stopifnot(ident %in% Seurat::Idents(object))
  cluster_expr <- t(expr[features, Seurat::WhichCells(object, idents = ident)])
  cluster_corr <- cor(as.matrix(cluster_expr), method = method)
  return(cluster_corr)
}

ExtractVariableFeatures <- function(name, object_ls, recalc = FALSE) {
  if (recalc) {
    object_ls[[name]] <- FindVariableFeatures(object_ls[[name]], selection.method = "vst", nfeatures = 4000, loess.span = 0.3, verbose = TRUE)
  }
  var_df <- object_ls[[name]]@assays$RNA@meta.features
  var_df$set <- name
  var_df <- rownames_to_column(var_df, var = "feature")
  return(var_df)
}


CalcDiff <- function(i, obj){
  first <- i
  #second <- combn(cluster_names, 2)[2, i]
  second <- "NonMP"
  comparison = paste0(first, "__", "NonMP")
  print(comparison)
  out <- FindMarkers(obj, ident.1 = first, ident.2 = second, logfc.threshold = log(1.5), verbose = TRUE,
    test.use = "LR", latent.vars = "batch", min.pct = 0.2)
  out$comparison <- comparison
  return(out)
}

CalcAUC <- function(i, obj) {
  comparison = paste0(i, "_", "NonMP")
  out <- presto::wilcoxauc(obj, assay = "data", seurat_assay = "RNA",
    verbose = FALSE, groups_use = c(i, "NonMP"))
  out$comparison <- comparison
  out <- out %>% filter(auc > 0.5 & logFC > log(1.1) & group != "NonMP")
  ui_done("Found markers for cluster {i}")
  Clean()
  return(out)
}

PlotFactorUMAP <- function(df, fill, label = TRUE, pal = "Dark 2") {
  if (label == TRUE) {
    centroids <- factor_umap %>% group_by(!!sym(fill)) %>% summarize(mean_d1 = mean(UMAP1), mean_d2 = mean(UMAP2))
    factor_plot <- ggplot(factor_umap, aes(x = UMAP1, y = UMAP2, fill = !!sym(fill))) +
      geom_point(shape = 21, color = "black", alpha = 0.8, size = 3, na.value = "gray90") +
      ggrepel::geom_label_repel(data = centroids, aes(x = mean_d1, y = mean_d2, label = !!sym(fill)),
        fill = "white", fontface = "bold", color = "black") +
      labs(x = "UMAP1", y = "UMAP2") +
      colorspace::scale_fill_discrete_qualitative(palette = "Dark 2", name = "Partitions") +
      labs(x = "UMAP 1", y = "UMAP 2", subtitle = glue("By variable: {snakecase::to_snake_case(fill)}")) +
      theme_void(base_size = 18) +
      theme(
        plot.title = element_text(hjust = 0, face = "plain", margin = margin(0, 0, 12, 0, unit = "pt")),
        plot.subtitle = element_text(hjust = 0, face = "plain", margin = margin(0, 0, 12, 0, unit = "pt")),
        plot.caption = element_text(hjust = 1, vjust = 1, face = "plain", margin = margin(0, 0, 0, 0, unit = "pt")),
        axis.title.x = element_text(hjust = 0, face = "plain", margin = margin(12, 0, 0, 0, unit = "pt")),
        axis.title.y = element_text(hjust = 0, angle = 90, face = "plain", margin = margin(0, 12, 0, 0, unit = "pt")),
        panel.background = element_rect(fill = "transparent", color = "black", size = 1), # bg of the panel
        plot.background = element_rect(fill = "transparent", color = "transparent"), # bg of the plot
        panel.grid.major = element_blank(), # get rid of major grid
        panel.grid.minor = element_blank(), # get rid of minor grid
        legend.background = element_rect(fill = "transparent", color = "transparent"), # get rid of legend bg
        legend.box.background = element_rect(fill = "transparent", color = "transparent"), # get rid of legend panel bg
        legend.position = "none"
      )
  } else {
    factor_plot <- ggplot(factor_umap, aes(x = UMAP1, y = UMAP2, fill = !!sym(fill))) +
      geom_point(shape = 21, color = "black", alpha = 0.8, size = 3) +
      labs(x = "UMAP1", y = "UMAP2") +
      colorspace::scale_fill_discrete_qualitative(palette = "Dark 2", name = "Partitions") +
      labs(x = "UMAP 1", y = "UMAP 2", subtitle = glue("By variable: {snakecase::to_snake_case(fill)}")) +
      theme_void(base_size = 18) +
      theme(
        plot.title = element_text(hjust = 0, face = "plain", margin = margin(0, 0, 12, 0, unit = "pt")),
        plot.subtitle = element_text(hjust = 0, face = "plain", margin = margin(0, 0, 12, 0, unit = "pt")),
        plot.caption = element_text(hjust = 1, vjust = 1, face = "plain", margin = margin(0, 0, 0, 0, unit = "pt")),
        axis.title.x = element_text(hjust = 0, face = "plain", margin = margin(12, 0, 0, 0, unit = "pt")),
        axis.title.y = element_text(hjust = 0, angle = 90, face = "plain", margin = margin(0, 12, 0, 0, unit = "pt")),
        panel.background = element_rect(fill = "transparent", color = "black", size = 1), # bg of the panel
        plot.background = element_rect(fill = "transparent", color = "transparent"), # bg of the plot
        panel.grid.major = element_blank(), # get rid of major grid
        panel.grid.minor = element_blank(), # get rid of minor grid
        legend.background = element_rect(fill = "transparent", color = "transparent"), # get rid of legend bg
        legend.box.background = element_rect(fill = "transparent", color = "transparent"), # get rid of legend panel bg
        legend.position = "none"
      )
  }
  return(factor_plot)
}

CompileFactorReps <- function(
  factors,
  kfactors = 20,
  num_reps = 30,
  rho = 0.3,
  quant_thresh = 0.90
) {

  factors <- Reduce(rbind, factors)
  ui_info("{nrow(factors)} factors by {ncol(factors)} genes")
  norm_w <- apply(factors, 1, function(x) x/L2Norm(x))
  norm_w <- t(norm_w)

  knn <- rho*num_reps
  nn <- RANN::nn2(data = norm_w, query = norm_w, k = knn)
  avg_dist_to_nns <- apply(nn$nn.dists, 1, function(x) sum(x)/knn)
  ui_info("Removing factors with average distance > {signif(quantile(avg_dist_to_nns, quant_thresh), 2)}")
  fnorm_w <- norm_w[avg_dist_to_nns < quantile(avg_dist_to_nns, quant_thresh), ]
  dist_w <- dist(fnorm_w)
  kfnorm_w <- kmeans(dist_w, centers = kfactors, nstart = 30, iter.max = 30)
  ui_todo("Compiling median factor loadings...")
  median_w <- map(unique(kfnorm_w$cluster), ~ {
    factors <- as.numeric(names(kfnorm_w$cluster)[kfnorm_w$cluster == .x])
    subset_w <- fnorm_w[factors, ]
    median_factors <- apply(subset_w, 2, median)
    return(median_factors)
  })
  median_w <- Reduce(rbind, median_w)
  median_w <- apply(median_w, 2, function(x) x/sum(x))
  rownames(median_w) <- paste0("F_", rep(1:nrow(median_w)))
  return(median_w)
}


DefineConsensusMarkers <- function(
  markers,
  object_names,
  k,
  logFC_threshold = log(1.25),
  auc_threshold = 0.6,
  padj_threshold = 1E-3,
  multi_marker_threshold = 0.6,
  num_genes = 50,
  min_degs = 5,
  min_jacc_remove = 0.8
) {

  # pull markers
  pb <- progress_estimated(length(object_names))
  def_markers <- map(.x = object_names, .f = ~ {
    fil_markers <- markers[markers$set == .x, ] %>%
      dplyr::filter(logFC >= logFC_threshold & auc >= auc_threshold & padj <= padj_threshold)
    num_group <- length(unique(fil_markers$group))
    bad_markers <- fil_markers %>% dplyr::group_by(feature) %>% dplyr::summarize(frac = n()/num_group) %>%
      filter(frac >= multi_marker_threshold) %>% dplyr::pull(feature)
    fil_markers <- fil_markers %>% dplyr::filter(!feature %in% bad_markers)
    fil_markers <- fil_markers %>% dplyr::group_by(group) %>% dplyr::top_n(num_genes, auc)
    pb$tick()$print()
    return(fil_markers)
  }, markers = markers)

  # split each dataset into a list of genes
  markers_ls <- map(.x = def_markers, .f = ~ {
    marker_ls <- split(.x$feature, .x$group)
    return(marker_ls)
  })
  num_clusters <- unlist(lapply(markers_ls, length))
  summary(num_clusters)
  markers_ls <- unlist(markers_ls, recursive = FALSE)

  # remove clusters with less than 5 markers
  markers_ls <- markers_ls[sapply(markers_ls, length) >=  min_degs]

  # tabulate and calculate jaccard distance
  markers_totab <- markers_ls
  markers_tab <- qdapTools::mtabulate(markers_totab)
  jacc_dist <- (ade4::dist.binary(markers_tab, method = 1, diag = FALSE, upper = FALSE))^2
  jacc_mtx <- as.matrix(jacc_dist)
  diag(jacc_mtx) <- 1
  org_jacc_mtx <- jacc_mtx

  # filter poor performing matrices, lists with no matches above 0.1
  thresholds <- (rowSums(jacc_mtx < min_jacc_remove))
  remove <- names(thresholds[thresholds <= 0])
  remove <- !(names(markers_totab) %in% remove)
  markers_totab <- markers_totab[remove]
  markers_tab <- qdapTools::mtabulate(markers_totab)
  jacc_dist <- (ade4::dist.binary(markers_tab, method = 1, diag = FALSE, upper = FALSE))^2
  jacc_mtx <- as.matrix(jacc_dist)
  diag(jacc_mtx) <- 1

  # perform clustering
  hc <- hclust(jacc_dist, method = "ward.D2")

  # define groups and generate cluster annotation dataframe
  clustering <- cutree(hc, k = k)
  anno <- data.frame(clustering)
  anno$id <- rownames(anno)
  anno$id <- factor(anno$id, levels = hc$labels[hc$order])

  # add dataset information to anno
  dataset_names <- stringr::str_match(pattern = "^(.*)[[:punct:]](\\d{1,2})$", string = anno$id)
  anno$dataset <- dataset_names[, 2, drop = TRUE]
  anno$dataset <- as.character(anno$dataset)
  anno$orig_cluster <- dataset_names[, 3, drop = TRUE]
  colnames(anno) <- c("Agglo_Cluster", "ID", "Dataset", "Orig_Cluster")
  #clusters <- split(anno[, c("Dataset", "ID", "Orig_Cluster", "Agglo_Cluster")], anno$Agglo_Cluster)

  #return
  return(list(orig_matrix = org_jacc_mtx, markers = def_markers, matrix = jacc_mtx, hclust = hc, anno = anno))
}

PlotHeatmap <- function(matrix, row_cluster, column_cluster, anno,
  filename = glue("plots/markers_simmtx_{k}_{format(Sys.time(), '%y%m%d%I%M')}.png")) {

  k <- length(unique(anno$Agglo_Cluster))
  # set matrix colors
  # blues <- colorspace::sequential_hcl("Blues 2", n = length(c(0.1, seq(0.6, 1, 0.05))))
  # blues[length(blues)] <- "#FFFFFF"
  # colors <- circlize::colorRamp2(c(0.1, seq(0.6, 1, 0.05)), blues)

  # set matrix colors
  # blues <- colorspace::sequential_hcl("Blues 2", n = length(seq(0, 1, 0.1)))
  # blues[length(blues)] <- "#FFFFFF"
  # colors <- circlize::colorRamp2(seq(0, 1, 0.1), blues)
  colfunc <- colorRampPalette(c("#0072B2", "white"))
  colors <- circlize::colorRamp2(c(0, seq(0, 1, 0.1)), colfunc(length(c(0, seq(0, 1, 0.1)))))

  # set dataset colors
  dataset_cols <- colorspace::sequential_hcl("Grays", n = length(unique(anno$Dataset)))
  names(dataset_cols) <- unique(anno$Dataset)

  # set dataset annotation
  dataset <- rowAnnotation(Dataset = anno$Dataset,
    col = list(Dataset = dataset_cols),
    gp = gpar(col = "black"),
    show_legend = FALSE, annotation_name_rot = 0, annotation_name_side = "top")

  #png(filename, width = 12, height = 10, units = "in", res = 400)
  pdf(filename, width = 12, height = 10)
  ht <- Heatmap(matrix,
    border = TRUE, col = colors,
    cluster_rows = row_cluster, cluster_columns = column_cluster, show_column_dend = FALSE,
    show_row_dend = TRUE,
    show_row_names = FALSE, show_column_names = FALSE,
    left_annotation = dataset,
    bottom_annotation = HeatmapAnnotation(
      Cluster = anno_block(
        gp = gpar(fill = colorspace::qualitative_hcl("Dark 3", n = k)),
        labels = as.character(unique(anno$Agglo_Cluster[order(anno$ID)])),
        labels_gp = gpar(fontsize = 8, col = "black", fontface = "bold")
      ),
      #Hierarchical = as.character(anno$C1),
      #KMeans = as.character(anno$C2),
      #col = list(KMeans = anno_colors$C2, Hierarchical = anno_colors$C1),
      #show_legend = c(FALSE, FALSE, FALSE)
      show_legend = c(FALSE)
    ),
    row_split = k, column_split = k, row_title = NULL,
    row_gap = unit(0, "mm"), column_gap = unit(0, "mm"),
    heatmap_legend_param = list(
      title = "Jaccard\nDistance",
      at = c(0, 0.2, 0.4, 0.6, 0.8, 1.0),
      border = "black",
      legend_height = unit(4, "cm"),
      legend_width = unit(1, "cm")
      #direction = "horizontal",
      #title_position = "lefttop"
    ),
    height = 4,
    use_raster = TRUE, raster_quality = 1)
  draw(ht, heatmap_legend_side = "right")
  dev.off()
}

# determine k
PlotAggloHeuristics <- function(
  matrix
) {
  wss <- factoextra::fviz_nbclust(matrix, FUN = factoextra::hcut, method = "wss", k.max = 35)
  silo <- factoextra::fviz_nbclust(matrix, FUN = factoextra::hcut, method = "silhouette", k.max = 35)
  stats <- data.frame(k = wss$data$clusters, wss = wss$data$y, silo = silo$data$y)
  stats <- stats %>% mutate(rank_wss = dense_rank(wss), rank_silo = dense_rank(dplyr::desc(silo)))
  stats <- stats %>% mutate(total_rank = rank_wss + rank_silo)
  a <- ggplot(stats, aes(x = k, y = wss)) +
    geom_point(shape = 21, fill = "gray60", color = "black", size = 3, alpha = 0.8) +
    scale_x_discrete(breaks = c(1, seq(5, 30, 5))) +
    colorspace::scale_fill_discrete_sequential("Grays", name = "Statistic") +
    labs(x = "k Partition", y = "Within sum-of-squares") +
    PointTheme(14)
  b <- ggplot(stats, aes(x = k, y = silo)) +
    geom_point(shape = 21, fill = "gray60", color = "black", size = 3, alpha = 0.8) +
    scale_x_discrete(breaks = c(1, seq(5, 30, 5))) +
    colorspace::scale_fill_discrete_sequential("Grays", name = "Statistic") +
    labs(x = "k Partition", y = "Silhouette Coefficient") +
    PointTheme(14)
  heuristics_plot <- cowplot::plot_grid(a, b)
  return(list(plot = heuristics_plot, rangek = which.min(wss$data$y[wss$data$y > min(wss$data$y)+(max(wss$data$y)-min(wss$data$y))*0.05]),
    rankk = which.min(stats$total_rank)))
}

CollateMarkers <- function(clusters, markers, required_nsets = 0.33) {

  grouped_markers <- map(.x = names(clusters), .f = ~ {
    ui_info("Retrieving genes for cluster {.x}")
    cluster_genes <- list()
    for (i in seq_along(1:nrow(clusters[[.x]]))) {
      dataset <- clusters[[.x]][i, "Dataset"]
      cluster <- clusters[[.x]][i, "Orig_Cluster"]
      ui_info("Pulling markers from {dataset}, {cluster}")
      genes <- markers[[dataset]] %>% dplyr::filter(group == cluster)
      cluster_genes[[i]] <- genes
    }
    return(cluster_genes)
  })

  grouped_markers <- lapply(grouped_markers, function(x) {
    return(Reduce(rbind, x))
  })
  names(grouped_markers) <- names(clusters)

  top_markers <- pblapply(names(grouped_markers), function(x) {
    grouped_markers[[x]] <- grouped_markers[[x]] %>%
      group_by(feature) %>%
      mutate(n_sets = n(), median_auc = median(auc), median_logFC = median(logFC))
    threshold <- max(c(round(max(grouped_markers[[x]]$n_sets)*required_nsets), 2))
    top_markers <- grouped_markers[[x]] %>%
      filter(n_sets >= threshold) %>%
      group_by(feature) %>%
      arrange(desc(median_logFC)) %>%
      distinct(feature, .keep_all = TRUE)
    top_markers$partition <- x
    return(top_markers)
  })

  top_markers_df <- Reduce(rbind, top_markers)
  return(top_markers_df)
}


AddScore <- function(
  object,
  features,
  pool = NULL,
  nbin = 24,
  ctrl = 100,
  k = FALSE,
  assay = NULL,
  name = 'Cluster',
  seed = 1,
  search = FALSE,
  ...
) {
  if (!is.null(x = seed)) {
    set.seed(seed = seed)
  }
  assay.old <- DefaultAssay(object = object)
  assay <- assay %||% assay.old
  DefaultAssay(object = object) <- assay
  assay.data <- GetAssayData(object = object)
  rowmeans <- Matrix::rowMeans(assay.data, na.rm = TRUE)
  assay.data <- assay.data[rowmeans > 0, ]
  features.old <- features
  if (k) {
    .NotYetUsed(arg = 'k')
    features <- list()
    for (i in as.numeric(x = names(x = table(object@kmeans.obj[[1]]$cluster)))) {
      features[[i]] <- names(x = which(x = object@kmeans.obj[[1]]$cluster == i))
    }
    cluster.length <- length(x = features)
  } else {
    if (is.null(x = features)) {
      stop("Missing input feature list")
    }
    features <- lapply(
      X = features,
      FUN = function(x) {
        missing.features <- setdiff(x = x, y = rownames(x = assay.data))
        if (length(x = missing.features) > 0) {
          warning(
            "The following features are not present in the object: ",
            paste(missing.features, collapse = ", "),
            ifelse(
              test = search,
              yes = ", attempting to find updated synonyms",
              no = ", not searching for symbol synonyms"
            ),
            call. = FALSE,
            immediate. = TRUE
          )
          if (search) {
            tryCatch(
              expr = {
                updated.features <- UpdateSymbolList(symbols = missing.features, ...)
                names(x = updated.features) <- missing.features
                for (miss in names(x = updated.features)) {
                  index <- which(x == miss)
                  x[index] <- updated.features[miss]
                }
              },
              error = function(...) {
                warning(
                  "Could not reach HGNC's gene names database",
                  call. = FALSE,
                  immediate. = TRUE
                )
              }
            )
            missing.features <- setdiff(x = x, y = rownames(x = assay.data))
            if (length(x = missing.features) > 0) {
              warning(
                "The following features are still not present in the object: ",
                paste(missing.features, collapse = ", "),
                call. = FALSE,
                immediate. = TRUE
              )
            }
          }
        }
        return(intersect(x = x, y = rownames(x = assay.data)))
      }
    )
    cluster.length <- length(x = features)
  }
  if (!all(LengthCheck(values = features))) {
    warning(paste(
      'Could not find enough features in the object from the following feature lists:',
      paste(names(x = which(x = !LengthCheck(values = features)))),
      'Attempting to match case...'
    ))
    features <- lapply(
      X = features.old,
      FUN = CaseMatch,
      match = rownames(x = object)
    )
  }
  if (!all(LengthCheck(values = features))) {
    stop(paste(
      'The following feature lists do not have enough features present in the object:',
      paste(names(x = which(x = !LengthCheck(values = features)))),
      'exiting...'
    ))
  }
  pool <- rownames(x = assay.data)
  data.avg <- Matrix::rowMeans(x = assay.data[pool, ])
  data.avg <- data.avg[order(data.avg)]
  data.cut <- cut_number(x = data.avg + rnorm(n = length(data.avg))/1e30, n = nbin, labels = FALSE, right = FALSE)
  #data.cut <- as.numeric(x = Hmisc::cut2(x = data.avg, m = round(x = length(x = data.avg) / (nbin + 1))))
  names(x = data.cut) <- names(x = data.avg)
  ctrl.use <- vector(mode = "list", length = cluster.length)
  for (i in 1:cluster.length) {
    features.use <- features[[i]]
    for (j in 1:length(x = features.use)) {
      ctrl.use[[i]] <- c(
        ctrl.use[[i]],
        names(x = sample(
          x = data.cut[which(x = data.cut == data.cut[features.use[j]])],
          size = ctrl,
          replace = FALSE
        ))
      )
    }
  }
  ctrl.use <- lapply(X = ctrl.use, FUN = unique)
  ctrl.scores <- matrix(
    data = numeric(length = 1L),
    nrow = length(x = ctrl.use),
    ncol = ncol(x = object)
  )

  all_genes <- c(unique(unlist(ctrl.use)), unique(unlist(features)))
  ui_todo("Scaling data...")
  assay.data <- t(apply(assay.data[all_genes, ], 1, scale, center = FALSE))
  ui_done("Done")
  gc()
  for (i in 1:length(ctrl.use)) {
    features.use <- ctrl.use[[i]]
    ctrl.scores[i, ] <- Matrix::colMeans(x = assay.data[features.use, ])
  }

  features.scores <- matrix(
    data = numeric(length = 1L),
    nrow = cluster.length,
    ncol = ncol(x = object)
  )
  for (i in 1:cluster.length) {
    features.use <- features[[i]]
    data.use <- assay.data[features.use, , drop = FALSE]
    features.scores[i, ] <- Matrix::colMeans(x = data.use)
  }
  assertthat::assert_that(!all(is.nan(features.scores)))
  assertthat::assert_that(!all(is.nan(ctrl.scores)))
  features.scores.use <- features.scores - ctrl.scores
  rownames(x = features.scores.use) <- paste0(name, 1:cluster.length)
  features.scores.use <- as.data.frame(x = t(x = features.scores.use))
  rownames(x = features.scores.use) <- colnames(x = object)
  object[[colnames(x = features.scores.use)]] <- features.scores.use
  CheckGC()
  DefaultAssay(object = object) <- assay.old
  return(object)
}

LengthCheck <- function(values, cutoff = 0) {
  return(vapply(
    X = values,
    FUN = function(x) {
      return(length(x = x) > cutoff)
    },
    FUN.VALUE = logical(1)
  ))
}

CheckGC <- function() {
  if (getOption(x = "Seurat.memsafe")) {
    gc(verbose = FALSE)
  }
}


rowScale = function(x,
  center = TRUE,
  scale = TRUE,
  add_attr = TRUE,
  rows = NULL,
  cols = NULL) {

  if (!is.null(rows) && !is.null(cols)) {
    x <- x[rows, cols, drop = FALSE]
  } else if (!is.null(rows)) {
    x <- x[rows, , drop = FALSE]
  } else if (!is.null(cols)) {
    x <- x[, cols, drop = FALSE]
  }

  ################
  # Get the column means
  ################
  cm = rowMeans(x, na.rm = TRUE)
  ################
  # Get the column sd
  ################
  if (scale) {
    csd = apply(x, )
  } else {
    # just divide by 1 if not
    csd = rep(1, length = length(cm))
  }
  if (!center) {
    # just subtract 0
    cm = rep(0, length = length(cm))
  }
  x = (x - cm) / csd
  if (add_attr) {
    if (center) {
      attr(x, "scaled:center") <- cm
    }
    if (scale) {
      attr(x, "scaled:scale") <- csd
    }
  }
  return(x)
}

TestScores <- function(
  .x,
  iterations,
  cdf,
  scaled_scores
) {

  # pull datasets and program
  source_i <- iterations[.x, 1]
  program_i <- iterations[.x, 2]
  ui_todo("Processing {source_i} and {program_i}...")
  datasets <- as.character(cdf %>% filter(source == source_i) %>% pull(set))

  # compile scores
  model_scores <- scaled_scores %>%
    group_by(set, group) %>%
    filter(geneset == as.character(program_i) & set %in% datasets)

  # define threshold
  mm <- mixtools::normalmixEM(model_scores$score, k = 2)
  min_peak <- which.min(mm$mu)
  threshold <- mm$mu[min_peak] + 2*mm$sigma[min_peak]

  # pull top clusters
  passing_clusters <- model_scores %>%
    group_by(set, group) %>%
    summarize(medscore = mean(score), set_cluster = unique(set_cluster)) %>%
    top_n(2, medscore) %>%
    #filter(medscore > 0) %>%
    pull(set_cluster)

  # filter scores
  fmodel_scores <- model_scores %>% filter(set_cluster %in% passing_clusters)

  # sample stats
  split_model_scores <- split(fmodel_scores$score, fmodel_scores$disease_classification)
  reps = 1000
  sample_size = round(min(table(fmodel_scores$disease_classification))/4)
  sampled_stat <- furrr::future_map_dfr(1:reps, ~ {
    sh <- sample(split_model_scores[["health"]], size = sample_size, replace = TRUE)
    sd <- sample(split_model_scores[["disease"]], size = sample_size, replace = TRUE)
    #return(var(sh)/mean(sh)-var(sd)/mean(sd))
    #return(mean(sd) - mean(sh))
    #return(sd(sd)-sd(sh))
    return(data.frame(qd = quantile(sd, .9), qh = quantile(sh, 0.9), diff = quantile(sd, .9) - quantile(sh, 0.9),
      dthres = sum(sd > threshold)/sample_size, hthres = sum(sh > threshold)/sample_size,
      thresdiff = sum(sd > threshold)/sample_size - sum(sh > threshold)/sample_size))
  }, split_model_scores = split_model_scores, sample_size = sample_size)

  # fit linear model
  ui_todo("Fitting model with {nrow(fmodel_scores)} cells...")
  linear <- nlme::lme(fixed = score ~ disease_var, random =  ~ 1 | batch, data = fmodel_scores)
  linear$data <- NULL
  attr(linear$terms, ".Environment") <- NULL
  linear$modelStruct <- NULL
  #format(object.size(linear), "Mb")
  #saveRDS(linear, "test.rds")
  linear_tidy <- broom::tidy(linear)
  #p_linear <- lmtest::lrtest(linear, 1)[2, 5]
  p_linear <- nlme::anova.lme(linear, Terms = 2)[1, 4]

  # compile results
  kpis <- data.frame(source = source_i)
  kpis$program <- program_i
  kpis$plinear <- p_linear
  #models <- list(lm = linear)

  Clean()
  return(list(stats = sampled_stat, model = linear, model_tidy = linear_tidy, pvalue = kpis))
}

ScoreObject <- function(object, geneset, name, group, ...) {

  geneset <- geneset[geneset %in% rownames(object)]
  count_subset <- GetAssayData(object, slot = "counts")[geneset, ]
  percents <- Matrix::colSums(count_subset > 0)
  percents <- percents/nrow(count_subset)
  object <- AddModuleScore(object, features = list(geneset), name = "Geneset", ...)
  if (grepl("name", colnames(object@meta.data))) {
    ui_info("{name} already present in object, replacing")
    object@meta.data[grep("name", colnames(object@meta.data))] <- NULL
  }
  colnames(object@meta.data)[grep("Geneset1", colnames(object@meta.data))] <- name

  df <- object@meta.data %>% select(!!sym(group), batch, !!sym(name), nCount_RNA)
  df <- df %>% rownames_to_column("cell_barcode")
  df <- df %>% gather(key = "program", value = "score", -!!sym(group), -cell_barcode, -batch, -nCount_RNA)
  df$percent <- percents
  df$set <- unique(object$set)
  return(df)
}

PlotGenesetTrends <- function(dataset, geneset, scores, num_top) {
  tp <- scores[[paste(dataset, geneset, "_")]]
  tpscores <- bind_rows(scores[as.character(cdf$set[cdf$source == dataset])]) %>%
    filter(program == geneset)
  passing_clusters <- tpscores %>%
    group_by(set, merged_leiden) %>%
    summarize(mean = median(score), set_cluster = unique(set_cluster)) %>%
    top_n(num_top, mean) %>%
    select(set_cluster)
  tpscores <- tpscores %>% mutate(plot = ifelse(set_cluster %in% passing_clusters$set_cluster, TRUE, FALSE))

  dt <- diptest::dip.test(tpscores$score)
  dt$p.value
  if (dt$p.value >= 0.5) {
    median_score <- median(tpscores$score)
    mad_score <- mad(tpscores$score)
    threshold <- median_score + 1.5*mad_score
  } else {
    mm <- mixtools::normalmixEM(tpscores$score, k = 2)
    plot(mm, 2)
    min_mm <- which.min(mm$mu)
    threshold <- mm$mu[min_mm] + mm$sigma[min_mm]*2
  }

  sum_scores <- tpscores %>%
    group_by(set, program) %>%
    summarize(n = n(), mean_score = mean(score),
      n_above = sum(score > threshold)/n)

  a <- ggplot(tpscores, aes(x = score, y = category, fill = category)) +
    ggridges::geom_density_ridges(alpha = 0.3, color = "black") +
    geom_vline(xintercept = threshold, size = 1, linetype = "dotted") +
    scale_fill_manual(values = c("Gray60", "#D55E00")) +
    guides(fill = FALSE) +
    scale_y_discrete(expand = expansion(mult = c(0.25, 1), add = c(0, 0.5)),
      labels = snakecase::to_title_case(rev(unique(as.character(tpscores$category))))) +
    annotate(geom = "text", label = glue("{round(sum_scores$n_above[1], 2)}% above threshold"),
      y = 2.5, x = 1.25, fontface = "bold") +
    annotate(geom = "text", label = glue("{round(sum_scores$n_above[2], 2)}% above threshold"),
      y = 1.5, x = 1.25, fontface = "bold") +
    labs(x = "Score", y = "", subtitle = "% cells above global threshold (1.5 MADs)") +
    DensityTheme(14)

  b <- ggplot(tpscores %>% filter(plot == TRUE),
    aes(x = nCount_RNA, y = score , color = category)) +
    geom_point(alpha = 0.5, size = 0.2) +
    geom_density2d(size = 0.5) +
    scale_x_continuous(trans = scales::log10_trans(),
      breaks = scales::trans_breaks("log10", function(x) 10^x),
      labels = scales::trans_format("log10", scales::math_format(10^.x))) +
    scale_color_manual(values = c("Gray60", "#D55E00")) +
    labs(x = "log10(Library size)", y = "Score", subtitle = "Gene set score vs. library size (top 2 clusters)") +
    ViolinTheme() + NoLegend()
  b <- ggExtra::ggMarginal(b, tpscores %>% filter(plot == TRUE),
    x = score, y = nCount_RNA, type = "density", groupFill = TRUE)

  numeric_flvls <- gsub(glue("{dataset}_(health|disease)_"), "", unique(tpscores$set_cluster))
  lvls <- unique(tpscores$set_cluster)[order(as.numeric(numeric_flvls))]
  lvls <- c(lvls[grep("disease", lvls)], lvls[grep("health", lvls)])
  tpscores$set_cluster <- factor(tpscores$set_cluster, levels = rev(lvls))
  lvl_labels <- gsub(glue("{dataset}_(health|disease)_"), "", rev(lvls))
  #lvl_labels <- snakecase::to_title_case(gsub("_", "", lvl_labels))

  c <- ggplot(tpscores,
    aes(x = score, y = set_cluster, fill = category)) +
    #geom_jitter(shape = 21, alpha = 0.25, width = 0.1, size = 0.2) +
    geom_boxplot(alpha = 0.75, width = 0.75, position = position_dodge(width = 0.5), outlier.color = "Gray60", outlier.alpha = 0.25) +
    scale_fill_manual(values = c("Gray60", "#d55e00")) +
    scale_y_discrete(labels = lvl_labels) +
    scale_x_continuous(limits = c(min(tpscores$score),
      ceiling(max(tpscores$score)))) +
    # annotate(geom = "text", label = glue("P = {format.pval(unique(kpis$plogit))}, \u03b2 = {round(unique(kpis$blogit), 2)}"), x = "health",
    #   y = ceiling(max(scores_to_plot$score)), hjust = 0) +
    labs(y = "Cluster", x = "Score", subtitle = c("Cluster-specific expression of gene set")) +
    ViolinTheme(base_size = 14) + NoLegend()

  title <- cowplot::ggdraw() +
    cowplot::draw_label(glue("Marker {gsub('M_', '', geneset)}, {gsub('_', ' ', dataset)}"), fontface = 'bold', x = 0, hjust = 0) +
    theme(plot.margin = margin(0, 0, 0, 24))

  grid <- cowplot::plot_grid(title, cowplot::plot_grid(a, b, c, nrow = 1), nrow = 2, rel_heights = c(0.1, 1))
  SavePlot("plots/{dataset}_{geneset}.png", grid, base_asp = 3, base_height = 6)
}

ScoreObjectScaled <- function(object_ls, set_df, set_source, genesets, names, group) {
  ui_todo("Merging object...")
  sets <- set_df$set[set_df$source == set_source]
  obj <- merge(object_ls[[sets[1]]], object_ls[[sets[2]]], add.cell.ids = c("d", "h"), merge.data = TRUE)
  ui_done("Object merged")

  ui_todo("Preparing gene sets...")
  genesets <- map(names, ~ {
    genesets[[.x]] <- genesets[[.x]][genesets[[.x]] %in% rownames(obj)]
    return(genesets[[.x]])
  }, obj = obj)
  percents <- map_dfr(.x = names, ~ {
    count_subset <- GetAssayData(obj, slot = "counts")[genesets[[.x]], ]
    percents <- Matrix::colSums(count_subset > 0)
    percents <- percents/nrow(count_subset)
    percents <- data.frame(percents = percents, geneset = .x, cell_barcode = colnames(obj))
    return(percents)
  }, obj = obj)
  ui_done("Gene sets prepared")

  assertthat::assert_that(all(lapply(genesets, length) > 0))

  ui_todo("Scoring object...")
  obj <- AddScore(object = obj, features = genesets, nbin = 30, ctrl = 30, name = "Geneset", seed = 1)
  assertthat::assert_that(length(grep("Geneset[[:digit:]]", colnames(obj@meta.data))) == length(genesets))
  colnames(obj@meta.data)[grep("Geneset[[:digit:]]", colnames(obj@meta.data))] <- names
  ui_done("Object scored")

  df <- obj@meta.data[, c("nCount_RNA", "disease_classification", "set", group, "batch", names)]
  df <- df %>% rownames_to_column("cell_barcode")
  df <- df %>% gather(key = "program", value = "score", names)
  df <- cbind(df, percents[, c("percents", "geneset")])
  df$source <- set_source
  df
  colnames(df)[5] <- "group"
  colnames(df)[10] <- "percents_geneset"
  colnames(df)[7] <- "geneset"
  gc()
  return(df)
}

#' @param object Object (SingleCellExperiment) to classify
#' @param indices Named list of indices to use
#' @param threhsold Similarity threshold to call labels
ClusterClassifier <- function(object, indices, threshold = 0.5) {
  ui_info("Classifying object against {names(indices)}...")
  mapped <- scmapCluster(
    projection = object,
    index_list = indices,
    threshold = threshold
  )
  usethis::ui_done("Classification complete")
  #usethis::ui_oops("Number of cells unclassified: {10}")
  return(mapped)
}

ScoreObjectScaled <- function(object_ls, set_df, set_source, genesets, names, group) {
  ui_todo("Merging object...")
  sets <- set_df$set[set_df$source == set_source]
  obj <- merge(object_ls[[sets[1]]], object_ls[[sets[2]]], add.cell.ids = c("d", "h"), merge.data = TRUE)
  ui_done("Object merged")

  ui_todo("Preparing gene sets...")
  genesets <- map(names, ~ {
    genesets[[.x]] <- genesets[[.x]][genesets[[.x]] %in% rownames(obj)]
    return(genesets[[.x]])
  }, obj = obj)
  percents <- map_dfr(.x = names, ~ {
    count_subset <- GetAssayData(obj, slot = "counts")[genesets[[.x]], ]
    percents <- Matrix::colSums(count_subset > 0)
    percents <- percents/nrow(count_subset)
    percents <- data.frame(percents = percents, geneset = .x, cell_barcode = colnames(obj))
    return(percents)
  }, obj = obj)
  ui_done("Gene sets prepared")

  assertthat::assert_that(all(lapply(genesets, length) > 0))

  ui_todo("Scoring object...")
  obj <- AddScore(object = obj, features = genesets, nbin = 30, ctrl = 30, name = "Geneset", seed = 1)
  assertthat::assert_that(length(grep("Geneset[[:digit:]]", colnames(obj@meta.data))) == length(genesets))
  colnames(obj@meta.data)[grep("Geneset[[:digit:]]", colnames(obj@meta.data))] <- names
  ui_done("Object scored")

  df <- obj@meta.data[, c("nCount_RNA", "disease_classification", "set", group, "batch", names)]
  df <- df %>% rownames_to_column("cell_barcode")
  df <- df %>% gather(key = "program", value = "score", names)
  df <- cbind(df, percents[, c("percents", "geneset")])
  df$source <- set_source
  df
  colnames(df)[5] <- "group"
  colnames(df)[10] <- "percents_geneset"
  colnames(df)[7] <- "geneset"
  gc()
  return(df)
}

ScoreLigandActivity <- function(object, object_name, genesets, geneset_names) {

  # score object
  object <- AddModuleScore(object, genesets, nbin = 30, ctrl = 30, name = "Geneset")
  assertthat::assert_that(length(grep("Geneset[[:digit:]]", colnames(object@meta.data))) == length(genesets))
  colnames(object@meta.data)[grep("Geneset[[:digit:]]", colnames(object@meta.data))] <- geneset_names

  colnames(object[[]])
  unique(object$mnp_withclusters)
  max_ident <- max(as.numeric(unique(object$mnp_withclusters)), na.rm = TRUE)
  as.numeric(object$base_clustering[object$mnp_withclusters == "NonMP"]) + max_ident
  object@meta.data$mnp_withclusters[object$mnp_withclusters == "NonMP"] <- as.numeric(object$base_clustering[object$mnp_withclusters == "NonMP"]) + max_ident
  #DimPlot(object, group.by = "mnp_withclusters", reduction = "harmony_umap", label = TRUE)
  #table(object$mnp_withclusters)

  object <- SetIdent(object, value = "mnp_withclusters")
  idents_to_test <- 1:max_ident
  names(idents_to_test) <- as.character(1:max_ident)
  idents_to_test <- as.character(idents_to_test)

  # pull receiver clusters
  sum_scores <- object@meta.data %>% select_if(names(.) %in% geneset_names | grepl("mnp_withclusters", names(.))) %>%
    group_by(mnp_withclusters) %>%
    summarize_all(mean) %>%
    filter(mnp_withclusters %in% idents_to_test) %>%
    gather(key = "GEP", value = "score", -mnp_withclusters) %>%
    group_by(GEP) %>%
    top_n(2, score) %>%
    filter(score > 0)

  # grab activities per program
  programs_run <- unique(sum_scores$GEP)
  specific_activities <- map(programs_run, ~ {

    ident <- sum_scores$mnp_withclusters[sum_scores$GEP == .x]
    expressed_genes_receiver <- get_expressed_genes(ident, object, pct = 0.05)
    background_expressed_genes <- expressed_genes_receiver %>% .[. %in% rownames(ligand_target_matrix)]

    sender_celltypes <- setdiff(Idents(object), ident)
    sender_celltypes <- setdiff(sender_celltypes, names(table(object$mnp_withclusters))[table(object$mnp_withclusters) < 50])
    ui_info("Pulling genes for {.x}")
    list_expressed_genes_sender <-  sender_celltypes %>%
      unique() %>%
      pblapply(get_expressed_genes, object, 0.05)
    expressed_genes_sender <- list_expressed_genes_sender %>% unlist() %>% unique()

    geneset_oi <-  genesets[[.x]] %>% .[. %in% rownames(ligand_target_matrix)]

    ligands <- lr_network %>% pull(from) %>% unique()
    receptors <- lr_network %>% pull(to) %>% unique()

    expressed_ligands <- intersect(ligands, expressed_genes_sender)
    expressed_receptors <- intersect(receptors, expressed_genes_receiver)
    potential_ligands <- lr_network %>% filter(from %in% expressed_ligands & to %in% expressed_receptors) %>% pull(from) %>% unique()

    ui_info("Predicting activities for {.x}")
    ligand_activities <- predict_ligand_activities(geneset = geneset_oi,
      background_expressed_genes = background_expressed_genes,
      ligand_target_matrix = ligand_target_matrix,
      potential_ligands = potential_ligands)
    ligand_activities <-  ligand_activities %>% arrange(-pearson) %>% mutate(rank = rank(desc(pearson)))
    ligand_activities$program <- .x
    ligand_activities$cluster <- list(ident)
    return(ligand_activities)
  })

  sa <- bind_rows(specific_activities)
  sa$set <- object_name
  return(sa)
}