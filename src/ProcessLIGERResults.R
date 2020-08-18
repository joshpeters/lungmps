# ---
# Description: Process LIGER results from cluster
# Author: Josh Peters
# ---

# [1] Basic setup ----

# load libraries
rm(list = ls())
library(future)
library(pbapply)
library(Seurat)
source('src/Functions.R')
source('src/Utilities.R')
stamp <- PrepEnv()

# [2] Load files ----

# CRITICAL - set specific directory path
dir_path = ""
inmf_files <- list.files(dir_path)
file_names <- unique(str_match(string = inmf_files, pattern = "(.*)_(2018|2019|2020)_(health|disease)")[, 1, drop = TRUE])
file_names <- SelfName(file_names)

# turn files into list by dataset
inmf_files_byset <- map(file_names, ~ {
  inmf_files[grep(pattern = .x, inmf_files)]
})

# load all files
StartFutureLapply()
all_w <- pblapply(file_names, function(x) {
  set_w <- future.apply::future_lapply(inmf_files_byset[[x]], function(y, dir_path) {
    w <- readRDS(glue("{dir_path}/{y}"))
    w <- w$W
    return(w)
  }, dir_path = dir_path)
  return(set_w)
})
StopFutureLapply()

# save raw loaded files from server
saveRDS(all_w, file = "data/interim/compiled_W.rds")

# [3]  Compile factors for each dataset ----
safe_CompileFactorReps <- purrr::safely(CompileFactorReps)
median_w <- pblapply(file_names, function(x, w) {
  median_w <- safe_CompileFactorReps(w[[x]])
  return(median_w)
}, w = all_w)
ws <- lapply(file_names, function(x, ws) { ws[[x]]$result }, ws = median_w)

# factors = all_w[[file_names[5]]]
# kfactors = 20
# num_reps = 30
# rho = 0.3
# quant_thresh = 0.90
#
# factors <- Reduce(rbind, factors)
# ui_info("{nrow(factors)} factors by {ncol(factors)} genes")
# norm_w <- apply(factors, 1, function(x) x/L2Norm(x))
# norm_w <- t(norm_w)
#
# knn <- rho*num_reps
# nn <- RANN::nn2(data = norm_w, query = norm_w, k = knn)
# avg_dist_to_nns <- apply(nn$nn.dists, 1, function(x) sum(x)/knn)
# ui_info("Removing factors with average distance > {signif(quantile(avg_dist_to_nns, quant_thresh), 2)}")
# fnorm_w <- norm_w[avg_dist_to_nns < quantile(avg_dist_to_nns, quant_thresh), ]
# dist_w <- dist(fnorm_w)
# kfnorm_w <- kmeans(dist_w, centers = kfactors, nstart = 30, iter.max = 30)
# ui_todo("Compiling median factor loadings...")
# median_w <- map(unique(kfnorm_w$cluster), ~ {
#   factors <- as.numeric(names(kfnorm_w$cluster)[kfnorm_w$cluster == .x])
#   subset_w <- fnorm_w[factors, ]
#   median_factors <- apply(subset_w, 2, median)
#   return(median_factors)
# })
# median_w <- Reduce(rbind, median_w)
# median_w <- apply(median_w, 2, function(x) x/sum(x))
# rownames(median_w) <- paste0("F_", rep(1:nrow(median_w)))
# ws[[5]] <- median_w

# save median, compiled W
saveRDS(ws, file = "data/interim/median_W.rds")
