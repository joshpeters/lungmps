# ---
# Description: Miscellaneous functions used in this project.
# Author: Josh Peters
# ---

#' Timestamp
#' Get current timestamp
#' @return timestamp
#' @export
Timestamp <- function() {
  return(format(Sys.time(), '%y%m%d%I%M'))
}

#' Timestamp
#' Get current timestamp
#' @return timestamp
#' @usage GetDate()
#' @export
GetDate <- function() {
  return(format(Sys.Date(), "%y%m%d"))
}

#' Clean
#' gc but invisible
#' @return NULL
#' @export
Clean <- function() {
  return(invisible(gc(verbose = FALSE)))
}

#' Prep environment
#' Common setup for scripts
#' @return list with date and time
#' @export
PrepEnv <- function() {
  library(glue)
  library(tidyverse)
  set.seed(1)
  date <- format(Sys.Date(), "%Y%m%d")
  time <- format(Sys.time(), '%y%m%d%I%M')
  gc()
  return(list(date = date, time = time))
}

#' Convert count matrix to sparse matrix
#' @return dgCMatrix
ConvertToSparseCounts <- function(matrix, gene_column = TRUE) {
  genes <- matrix[, 1, drop = TRUE]
  barcodes <- colnames(matrix)[2:ncol(matrix)]
  matrix <- matrix[, -1]
  matrix <- Matrix::Matrix(as.matrix(matrix), sparse = TRUE)
  colnames(matrix) <- barcodes
  rownames(matrix) <- genes
  return(matrix)
}


#' Load and format the HGNC database
#' @return Data frame with HGNC information
#' @export
LoadHGNC <- function() {
  hgnc <- read_delim("data/hgnc_database.txt", delim = "\t")
  hgnc <- janitor::clean_names(hgnc)
  hgnc_columns <- colnames(hgnc)
  hgnc_approved <- hgnc[, c("approved_symbol", "ensembl_gene_id")]
  hgnc_alias <- hgnc[, c("alias_symbol", "ensembl_gene_id")]
  hgnc_previous <- hgnc[, c("previous_symbol", "ensembl_gene_id")]
  colnames(hgnc_approved) <- c("symbol", "id")
  colnames(hgnc_alias) <- c("symbol", "id")
  colnames(hgnc_previous) <- c("symbol", "id")
  hgnc_total <- rbind(hgnc_approved, hgnc_alias, hgnc_previous)
  return(unique(hgnc_total))
}

#' Write 10x files using DropletUtils
#' @importFrom DropletUtils write10xCounts
#' @export
Write10X <- function(
  object,
  rootpath
) {

  message("Writing files...")
  DropletUtils::write10xCounts(path = rootpath,
    x = object@assays$RNA@counts,
    barcodes = colnames(object@assays$RNA@counts),
    gene.id = rownames(object@assays$RNA@counts),
    gene.symbol = rownames(object@assays$RNA@counts),
    gene.type = "Gene Expression",
    overwrite = TRUE,
    type = "sparse",
    version = "3")

  meta_path <- glue("{rootpath}/meta.tsv")
  write_tsv(object@meta.data, meta_path)
  R.utils::gzip(meta_path)
}

#' Loads list of Seurat objects
#' This function assumes your objects are in a subdirectory called "objects"
#' @param dir.path
#' @param pattern
#' @param verbose
#' @return list of Seurat objects
#' @importFrom future.apply future_lapply
#' @export

LoadObjectList <- function(dir.path, file_pattern, name_pattern, exclude, verbose = TRUE) {

  if (missing(name_pattern)) { name_pattern <- file_pattern }
  object_files <- list.files(path = dir.path, full.names = TRUE)
  object_files_to_load <- object_files[grep(pattern = file_pattern, x = object_files)]
  object_files_to_load <- object_files_to_load[!grepl(pattern = exclude, x = object_files_to_load)]
  object_names <- stringr::str_match(object_files_to_load, pattern = paste0("objects\\/(.*)", name_pattern))[, 2]
  names(object_names) <- object_names
  names(object_files_to_load) <- object_names
  ui_todo("Loading in {length(object_files_to_load)} files")

  object_list <- future.apply::future_lapply(object_names, function(x) {
    if (verbose) usethis::ui_info("\nLoading {x} ({round(file.info(object_files_to_load[x])$size/1e6, 2)} mb)")
    obj <- readRDS(file = object_files_to_load[x])
    if (verbose) usethis::ui_done(("\nCompleted {x}"))
    Clean()
    return(obj)
  })
  names(object_list) <- object_names
  return(object_list)
}

#' Calculate the percentage of a vector above some threshold
#'
#' Seurat 3.1.1
#' @param x Vector of values
#' @param threshold Threshold to use when calculating percentage
#'
#' @return Returns the percentage of `x` values above the given threshold
#' @export
#'
PercentAbove <- function(x, threshold) {
  return(length(x = x[x > threshold]) / length(x = x))
}

#' Convert variable to string
#'
#' This function turns variable names into a string.
#' @param variable variable to convert to a string name
#' @return Returns a string of variable name
#' @usage VarToString(variable)
#' @export
#'
VarToString <- function(variable) {
  return(deparse(substitute(variable)))
}

#' Logarithmic sequence
#'
#' Get a sequence of numbers spaced logarithmically
#' @return sequence
#' @usage LogSeq(1E-6, 1E-1, 5)
#' @export
#'
LogSeq <- function(from, to, length) {
  sequence <- exp(seq(log(from), log(to), length.out = length))
  return(sequence)
}



#' Write message to logfile
#'
#' Writes message to a designated logfile
#' @param message Message for logfile
#' @param logname Name of logfile to open
#' @usage WriteLog(message, logname)
#' @importFrom glue glue
#' @export
WriteLog <- function(
  message,
  logname
) {
  cat(message, file = glue::glue("{logname}.log"), append = TRUE, sep = "\n")
  message(glue::glue("Wrote log to {logname}.log"))
}

#' Print top 5 rows and columns
Check5 <- function(x) {
  print(x[1:5, 1:5])
}


#' Set Future plan to multisession
StartFutureLapply <- function(num.workers = 4, memory.per.worker = 1000) {
  future::plan(strategy = "multisession", workers = num.workers, gc = FALSE)
  options(future.globals.maxSize = memory.per.worker * 1024^2)
}
#' Reset Future plan
StopFutureLapply <- function() {
  future:::ClusterRegistry("stop")
  future::plan(strategy = "sequential")
  ui_done("Multisession stopped.")
}

#' Get colorblind palette from Bang Wong
GetCBPalette <- function(with.gray = TRUE) {
  if (with.gray) return(c("#999999", "#E69F00", "#56B4E9", "#009E73",
    "#F0E442", "#0072B2", "#D55E00", "#CC79A7"))
  if (!with.gray) return(c("#000000", "#E69F00", "#56B4E9", "#009E73",
    "#F0E442", "#0072B2", "#D55E00", "#CC79A7"))
}

#' Display color palette
DisplayPalette <- function(colors) {
  n <- length(colors)
  image(1:n, 1, as.matrix(1:n), col = colors,
    xlab = "", ylab ="", yaxt = "n", bty = "n")
}

ScaleVector <- function(x){
  (x - mean(x, na.rm=TRUE)) / sd(x, na.rm=TRUE)
}

GetMetaReference <- function(object_names) {
  sets <- readxl::read_excel("/Users/jpeters/Projects/lungmps/data/raw/sets.xlsx")
  lung_sets <- sets %>% filter(Tissue == "Lung")
  lung_ids <- lung_sets$ID
  disease_sets <- lung_sets %>% filter(Type != "Healthy") %>% pull(ID)
  healthy_sets <- lung_sets %>% filter(Type == "Healthy") %>% pull(ID)
  names(disease_sets) <- disease_sets
  names(healthy_sets) <- healthy_sets
  dataset_healthy <- data.frame(set = healthy_sets, disease = rep("Control", length(healthy_sets)))
  dataset_disease <- data.frame(set = c("Valenzi_2019", "Morse_2019", "Habermann_2019", "Reyfman_2019",
    "Lambrechts_2018", "Laughney_2020", "Zilionis_2019", "Mayr_2020"),
    disease = c("SS-ILD", "PF", "PF", "PF", "Cancer", "Cancer", "Cancer", "PF"))
  conditions <- stringr::str_match(object_names, "([a-zA-Z]*_[0-9]{4})_(.*)$")
  conditions <- conditions[, c(1,2,3)]
  conditions[12, 3] <- "health"
  conditions <- as.data.frame(conditions)
  colnames(conditions) <- c("ID", "Study", "Condition")
  # conditions$Type <- c("Healthy", "Chronic", "Healthy", "Cancer", "Healthy",
  #   "Cancer", "Healthy", "Healthy", "Chronic", "Healthy",
  #   "Chronic", "Healthy", "Healthy", "Chronic", "Healthy",
  #   "Healthy", "Chronic", "Cancer")
  conditions$Disease <- "Health"
  matched <- match(conditions$Study[conditions$Condition != "health"], dataset_disease$set)
  conditions$Disease[conditions$Condition != "health"] <- as.character(dataset_disease$disease[matched])
  return(conditions)
}

CheckLsNames <- function(list, names) {
  if(is.null(names(list))) {
    ui_oops("List has no names")
    names(list) <- names
  } else {
    ui_done("List has names")
  }
  return(list)
}

RemoveColumnsFromObjList <- function(object.list, object.names, col.names) {
  object.list <- pblapply(object.names, function(y) {
    ui_info("Cleaning up {y}")
    col.names <- paste(col.names, collapse = "|")
    cols_to_remove <- grep(glue("({col.names})"), colnames(object.list[[x]]@meta.data))
    object.list[[x]]@meta.data[, cols_to_remove] <- NULL
    ui_done(("\nCompleted {y}"))
    return(object.list[[x]])
  })
  return(object.list)
}



NormalizeVector <- function(x)
{
  return((x- min(x, na.rm = TRUE)) /(max(x, na.rm = TRUE)-min(x, na.rm = TRUE)))
}

ConvertListToDF <- function(list, prefix = "M") {
  names(list) <- paste0(prefix, seq(1:length(list)))
  df <- lapply(names(list), function(y) {
    df <- data.frame(feature = list[[y]], module = y)
    return(df)
  })
  df <- do.call(rbind, df)
  return(df)
}

# Get density of points in 2 dimensions.
# @param x A numeric vector.
# @param y A numeric vector.
# @param n Create a square n by n grid to compute density.
# @return The density within each square.
GetDensity <- function(x, y, ...) {
  dens <- MASS::kde2d(x, y, ...)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
}

#' Merge sparse matrices
MergeSparseMatrices <- function(...) {
  cnnew <- character()
  rnnew <- character()
  x <- vector()
  i <- numeric()
  j <- numeric()
  for (M in list(...)) {
    cnold <- colnames(M)
    rnold <- rownames(M)
    cnnew <- union(cnnew,cnold)
    rnnew <- union(rnnew,rnold)
    cindnew <- match(cnold,cnnew)
    rindnew <- match(rnold,rnnew)
    ind <- summary(M)
    i <- c(i,rindnew[ind[,1]])
    j <- c(j,cindnew[ind[,2]])
    x <- c(x,ind[,3])
  }
  return(Matrix::sparseMatrix(i=i,j=j,x=x,dims=c(length(rnnew),length(cnnew)),dimnames=list(rnnew,cnnew)))
}

ConvertNamedVecToDF <- function(x) {
  df <- data.frame(name = names(x), value = x)
  return(df)
}

WriteH5 <- function(object, filename) {
  h5createFile(filename)
  h5createGroup(filename, "matrix")
  h5createGroup(filename, "matrix/features")
  matrix <- GetAssayData(object, slot = "counts")
  h5write(as.array(matrix@x), file = filename, name = "matrix/data")
  h5write(as.array(Cells(object)), file = filename, name = "matrix/barcodes")
  h5write(as.array(matrix@Dimnames[[1]]), file = filename, name = "matrix/features/name")
  h5write(as.array(matrix@i), file = filename, name = "matrix/indices")
  h5write(as.array(matrix@p), file = filename, name = "matrix/indptr")
}

#' Determine quantile breaks for colorscheme
#' @param x matrix for quantile calculation
#' @param n number of breaks
QuantileBreaks <- function(x, n = 10) {
  breaks <- quantile(x, probs = seq(0, 1, length.out = n))
  breaks[!duplicated(breaks)]
  return(breaks)
}

MatrixCosSim <- function(x, y){
  mat <- tcrossprod(x, y)
  t1 <- sqrt(apply(x, 1, crossprod))
  t2 <- sqrt(apply(y, 1, crossprod))
  sim <- mat / outer(t1, t2)
  return(sim)
}

L2Norm <- function(x) {
  return(sqrt(sum(x^2)))
}

SelfName <- function(list) {
  names(list) <- list
  return(list)
}

Brillouin <- function(x) {
  N <- sum(x)
  (log(factorial(N)) - sum(log(factorial(x))))/N
}

