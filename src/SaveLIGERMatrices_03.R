# ---``
# Description: Save LIGER matrices for clusters
# Author: Josh Peters
# ---

# [1] Basic setup ----

# load libraries
rm(list = ls())
library(future)
library(pbapply)
library(Seurat)
source("src/Functions.R")
source("src/Utilities.R")
source("src/Plotting.R")
stamp <- PrepEnv()

# [2] Save liger matrices ----

dir_path <- "data/objects"
file_pattern <- "_mnp_final.rds"
name_pattern <- "_mnp_final.rds"
verbose <- TRUE

object_ls <- LoadObjectList("data/objects", file_pattern = "_mnp_final.rds", verbose = TRUE)
ids <- names(object_ls)
ids <- SelfName(ids)

purrr::walk2(.x = object_ls, .y = ids, .f = ~ {
  split_obj <- SplitObject(.x, split.by = "batch")
  print(sapply(split_obj, ncol))
  split_obj <- split_obj[sapply(split_obj, ncol) > 50]
  matrices <- map(.x = split_obj, .f = GetAssayData, slot = "counts", assay = "RNA")
  saveRDS(matrices,
    file = glue::glue("data/objects/{.y}_mnp_matrices.rds"))
})
