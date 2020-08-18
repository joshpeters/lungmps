# ---
# Description: Check dataset formatting and prepare metadata
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

# load metadata
sets <- readxl::read_excel("data/study_metadata.xlsx")
ids <- sets$ID
disease_sets <- sets %>% filter(Type != "Healthy") %>% pull(ID)
disease_sets <- SelfName(disease_sets)
healthy_sets <- sets %>% filter(Type == "Healthy") %>% pull(ID)
healthy_sets <- SelfName(healthy_sets)

# load objects
ids <- SelfName(ids)
ids <- ids[order(ids)]
suffix <- "base_object.rds"
dir <- "data/objects/"
verbose <- TRUE
load_files <- list.files(path = dir, pattern = suffix, full.names = TRUE)
load_files <- load_files[grep(paste0(ids, collapse = "|"), x = load_files)]
ui_todo("Loading in {length(load_files)} files")
names(load_files) <- ids

StartFutureLapply()
object_ls <- future.apply::future_lapply(ids, function(x) {
  if (verbose) usethis::ui_info("\nLoading {x} ({round(file.info(load_files[x])$size/1e6, 2)} mb)")
  obj <- readRDS(file = load_files[x])
  if (verbose) usethis::ui_done(("\nCompleted {x}"))
  Clean()
  return(obj)
})
StopFutureLapply()
names(object_ls) <- ids
object_names <- sets$ID[sets$Train == TRUE]
object_names <- SelfName(object_names)

# [2] Check metadata for each dataset and compile sample level information

all_colnames <- map(object_ls, ~ {colnames(.x[[]])})
table(unlist(all_colnames))
# batch, batch_frequency, cell_barcode, disease_classification
# condition, donor, sample_id, study

sapply(all_colnames, function(x) {"condition" %in% x})
sapply(all_colnames, function(x) {"donor" %in% x})
sapply(all_colnames, function(x) {"sample_id" %in% x})
sapply(all_colnames, function(x) {"study" %in% x})

# [3] Prepare sample metadata ----
sample_metadata <- map(.x = object_names, .f = ~ {
  meta <- object_ls[[.x]]@meta.data %>% select(study, sample_id, donor, batch, batch_frequency, condition, disease_classification, percent_mit, percent_rib, percent_hsp, percent_hgb) %>%
    group_by(donor, sample_id, batch) %>%
    mutate(percent_mit = mean(percent_mit), percent_rib = mean(percent_rib), percent_hsp = mean(percent_hsp), percent_hgb = mean(percent_hgb)) %>%
    ungroup() %>%
    mutate_if(is.factor, as.character) %>%
    mutate_if(is.logical, as.character) %>%
    unique()
  #distinct(batch, .keep_all = TRUE)
  #   study = unique(study), batch_frequency = unique(batch_frequency), condition = unique(condition), disease_classification = unique(disease_classification)) %>%
  # select(study, donor, sample_id, batch, condition, disease_classification, everything())
  return(meta)
})
b_sample_metadata <- bind_rows(sample_metadata)
sum(b_sample_metadata %>% group_by(study, batch) %>% unique() %>% pull(batch_frequency))
write_csv(b_sample_metadata, path = "data/publication/TableS2.csv")