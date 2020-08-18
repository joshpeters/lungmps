# ---
# Description: Develop consensus variable genes
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

# [2] Identify consensus variable features ----
genes <- lapply(mnp_ls, rownames)
all_genes <- unique(unlist(lapply(mnp_ls, rownames)))

var_feature_ls <- pbapply::pblapply(ids, ExtractVariableFeatures, object_ls = mnp_ls, recalc = FALSE)
var_feature_df <- do.call(rbind, c(var_feature_ls, make.row.names = FALSE))
var_feature_df <- var_feature_df[!is.na(var_feature_df$vst.mean), ]
saveRDS(object = var_feature_df, file =  glue("data/all_genes_meta_df.rds"))

# assign top variable based on variance alone
top_var_genes <- var_feature_df %>%
  group_by(set) %>%
  top_n(n = 4000, wt = vst.variance) %>%
  mutate(vst.variable.raw = TRUE) %>%
  right_join(var_feature_df) %>%
  mutate(vst.variable.raw = replace_na(vst.variable.raw, FALSE))
table(top_var_genes$vst.variable, top_var_genes$vst.variable.raw)

# summarize all genes
sum_var_df <- top_var_genes %>%
  select(feature, name, vst.mean, vst.variance, vst.variance.expected, vst.variance.standardized, vst.variable, vst.variable.raw, set) %>%
  group_by(feature) %>%
  summarize(
    n = n(),
    n_var = sum(vst.variable == TRUE),
    n_var_raw = sum(vst.variable.raw == TRUE),
    median_var_std = median(vst.variance.standardized),
    mean_var_std = mean(vst.variance.standardized),
    mean_var = mean(vst.variance),
    median_var = median(vst.variance),
    mean_mean = mean(vst.mean),
    median_mean = median(vst.mean)) %>%
  distinct(feature, .keep_all = TRUE) %>%
  arrange(desc(n_var))

# filter consensus genes and rank by number of datasets they are variable in
numvar_ranking <- sum_var_df %>%
  filter(n >= 14) %>%
  filter(median_var_std > 0) %>%
  arrange(desc(median_var_std)) %>%
  top_n(n = 3750, wt = n_var)
table(numvar_ranking$n_var)
length(numvar_ranking$n_var)
quantile(numvar_ranking$median_var_std)
saveRDS(object = numvar_ranking$feature, file =  glue("data/consensusvargenes_rankbynvar_n3750_in14minimum.rds"))

# filter consensus genes and rank by number of median std. var
medianvar_ranking <- sum_var_df %>%
  filter(n == 18) %>%
  filter(n_var >= 5) %>%
  filter(median_var_std > 0) %>%
  arrange(desc(median_var_std)) %>%
  top_n(n = 4000, wt = median_var_std)
table(medianvar_ranking$n_var)
length(medianvar_ranking$n_var)
quantile(medianvar_ranking$median_var_std)

# compare rankings
table(medianvar_ranking$feature %in% numvar_ranking$feature)
table(numvar_ranking$feature %in% medianvar_ranking$feature)
setdiff(numvar_ranking$feature, medianvar_ranking$feature)
setdiff(medianvar_ranking$feature, numvar_ranking$feature)

# compare to Seurat procedure
seurat_integration_genes <- SelectIntegrationFeatures(mnp_ls, nfeatures = 4000, verbose = TRUE)
table(seurat_integration_genes %in% numvar_ranking$feature)
setdiff(seurat_integration_genes, numvar_ranking$feature)
saveRDS(object = seurat_integration_genes, file =  glue("data/consensusvargenes_seuratmethod_n4000.rds"))

# save CSV for publication
sum_var_df_save <- sum_var_df
sum_var_df_save$consensus <- FALSE
sum_var_df_save$consensus[sum_var_df_save$feature %in% numvar_ranking$feature] <- TRUE
colnames(sum_var_df_save) <- c("Feature", "# Datasets", "# Variable", "# Variable (raw)", "Median std. var.", "Mean std. var.", "Mean var.", "Median var.", "Mean mean expr.", "Median mean expr.", "Consensus")
write_csv(x = sum_var_df_save, path =  "data/publication/TableS3.csv")

# [3] Blacklist filtering ----
# identify blacklist genes
all_genes <- unique(var_feature_df$feature)
acs <- all_genes[grep(pattern = "(AC|AL|AP|AJ|AF|BX|AD|FP|Z|U)[[:digit::]{5,8}", x = all_genes)]
rn7s <- all_genes[grep(pattern = "RN7SL", x = all_genes)]
rnu <- all_genes[grep(pattern = "RNU[[:digit:]]{1}", x = all_genes)]
mt <- all_genes[grep(pattern = "MT-", x = all_genes)]
ct <- all_genes[grep(pattern = "CT(A|B|D|)-", x = all_genes)]
ig <- all_genes[grep("^IG[H/L/V/K].*$", all_genes)]
tcrb <- all_genes[grep("^TRBV.*$", all_genes)]
tcrg <- all_genes[grep("^TRGV.*$", all_genes)]
tcrc <- all_genes[grep("^TRGC.*$", all_genes)]
tcrj <- all_genes[grep("^TRGJ.*$", all_genes)]
tcra <- all_genes[grep("^TRA[C|V].*$", all_genes)]
tcrs <- c(tcrb, tcrg, tcrc, tcrj, tcra)
rps <- all_genes[grep("^RP([[:digit:]]){1,2}-", all_genes)]
blacklist <- c(acs, rn7s, rnu, mt, ct, ig, tcrs, rps)
length(blacklist)
table(blacklist %in% numvar_ranking$feature)
nonMP_markers <- readRDS("data/nonmp_markers_total.rds")

complete_blacklist <- c(blacklist, nonMP_markers)
table(numvar_ranking$feature %in% blacklist)
numvar_clean <- numvar_ranking[!(numvar_ranking$feature %in% blacklist), ]
saveRDS(object = numvar_clean$feature, file =  glue("data/consensusvargenes_rankbynvar_n3750_in14minimum_clean.rds"))

# [4] Generate plots ----
plot_df <- sum_var_df %>%
  mutate(
    label = ifelse(feature %in% blacklist, "Blacklist",
      ifelse(feature %in% nonMP_markers, "Non-MP Marker",
        ifelse(feature %in% numvar_ranking$feature, "Variable",
          ifelse(n <= 1, "Unique", "None")))))
colnames(plot_df)

unhits <- plot_df %>% filter(label == "None") %>% arrange(desc(mean_var_std))

plot_df$label <- factor(plot_df$label, levels = c("None", "Variable", "Unique", "Blacklist", "Non-MP Marker"))
plot_df <- plot_df[order(plot_df$label), ]
vargenes_plot <- ggplot(plot_df %>% filter(label %in% c("Variable", "None")), aes(x = mean_mean, y = mean_var_std, color = label)) +
  geom_point(alpha = 0.5) +
  scale_x_continuous(trans = "log10", limits = c(1E-5, 150), expand = c(0.05, 0), breaks = scales::trans_breaks("log10", function(x) 10^x),
    labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  scale_y_continuous(trans = "log1p", limits = c(0, 30), expand = c(0.05, 0)) +
  scale_color_manual(values = list(Variable = ggthemes::colorblind_pal()(8)[4], None = "#999999")) +
  annotate(geom = "label", x = 90, y = 5, label = "Variable", color = ggthemes::colorblind_pal()(8)[4], size = 4, fontface = "bold") +
  annotate(geom = "label", x = 5, y = 0.1, label = "Undefined", color = "#999999", size = 4, fontface = "bold") +
  ggrepel::geom_text_repel(
    data = plot_df %>% filter(label %in% c("Variable", "None") &
        feature %in% c("CCL17", "AREG", "CXCL10", "SPP1", "AREG", "CHIT1", "RBP4", "SERPINB2", "CCR7", "FOLR2", "FABP4", "IDO1", "MKI67")),
    mapping = aes(x = mean_mean, y = mean_var_std, label = feature),
    size = 4, color = "black"
  ) +
  guides(color = FALSE) +
  labs(x = "log10(Mean expression)", y = "log1p(Mean standardized variance)") +
  theme_classic(base_size = 14) +
  ggeasy::easy_all_text_color("black") +
  theme(
    axis.title.y = element_text(angle = 90, vjust = 0.5, margin = margin(0, 12, 0, 0), face = "bold"),
    axis.title.x = element_text(margin = margin(12, 0, 0, 0), face = "bold"),
    axis.text.y = element_text(margin = margin(0, 4, 0, 0)),
    #axis.ticks = element_blank(),
    axis.line = element_blank(),
    panel.border = element_rect(size = 1, color = "black", fill = "transparent")
    # legend.text = element_text(size = 10),
    # legend.justification = "top",
    # legend.key.size = unit(1, "line")
  )

plot_df$label <- factor(plot_df$label, levels = c("None", "Variable", "Blacklist", "Unique", "Non-MP Marker"))
plot_df <- plot_df[order(plot_df$label), ]
varnonhits_plot <- ggplot(plot_df %>% filter(label %in% c("Unique", "Blacklist", "Non-MP Marker")),
  aes(x = mean_mean, y = mean_var_std, color = label)) +
  geom_point(alpha = 0.5) +
  scale_x_continuous(trans = "log10", limits = c(1E-5, 150), expand = c(0.05, 0), breaks = scales::trans_breaks("log10", function(x) 10^x),
    labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  scale_y_continuous(trans = "log1p", limits = c(0, 30), expand = c(0.05, 0)) +
  scale_color_manual(values = list(Unique = ggthemes::colorblind_pal()(8)[6],
    Blacklist = ggthemes::colorblind_pal()(8)[1], `Non-MP Marker` = ggthemes::colorblind_pal()(8)[7])) +
  annotate(geom = "label", x = 1, y = 15, label = "Unique", color = ggthemes::colorblind_pal()(8)[6], size = 4, fontface = "bold") +
  annotate(geom = "label", x = 20, y = 4, label = "Non-MP Marker", color = ggthemes::colorblind_pal()(8)[7], size = 4, fontface = "bold") +
  annotate(geom = "label", x = 1E-2, y = 0, label = "Blacklist", color = ggthemes::colorblind_pal()(8)[1], size = 4, fontface = "bold") +
  guides(color = FALSE) +
  labs(x = "log10(Mean expression)", y = "log1p(Mean standardized variance)") +
  theme_classic(base_size = 14) +
  ggeasy::easy_all_text_color("black") +
  theme(
    axis.title.y = element_text(angle = 90, vjust = 0.5, margin = margin(0, 12, 0, 0), face = "bold"),
    axis.title.x = element_text(margin = margin(12, 0, 0, 0), face = "bold"),
    axis.text.y = element_text(margin = margin(0, 4, 0, 0)),
    #axis.ticks = element_blank(),
    axis.line = element_blank(),
    panel.border = element_rect(size = 1, color = "black", fill = "transparent")
    # legend.text = element_text(size = 10),
    # legend.justification = "top",
    # legend.key.size = unit(1, "line")
  )
varnonhits_plot
plot <- cowplot::plot_grid(vargenes_plot, varnonhits_plot, align = c("hv"), nrow = 1)
cowplot::save_plot(plot = plot, filename = "plots/fig1_nonvargenes.png", base_height = 6, base_asp = 2.1)
cowplot::save_plot(plot = plot, filename = "plots/fig1_nonvargenes.svg", base_height = 6, base_asp = 2.1)

features <- numvar_clean$feature
library(enrichR)
dbs <- listEnrichrDbs()
dbs <- c("GO_Biological_Process_2018")
results <- enrichr(features, databases = dbs)
results <- results[[dbs]][results[[dbs]]$Adjusted.P.value < 0.0001, ]
parsed <- stringr::str_match(results$Term, "^.*\\((GO\\:\\d{6,8})\\)$")
term_names <- parsed[, 1]
terms <- parsed[, 2]
terms <- data.frame(names = term_names, terms = terms, padj = results$Adjusted.P.value)
#clipr::write_clip(features, breaks = "\n")
clipr::write_clip(terms[, c(2,3)], breaks = "\n")
write.csv(terms, file = "data/interim/variable_enrichment.csv")
# submit csv online to REViGO
# use similarity of 0.7, SimRel, whole UniProt
# read in new csv here
revigo <- read_csv("data/interim/revigo.csv")

revigo_df <- revigo[revigo$eliminated == 0, ]
revigo_df$neglog10padj <- -revigo_df$`log10 p-value`
revigo_df <- revigo_df %>% arrange(neglog10padj)
revigo_df$description <- as.character(revigo_df$description)

# fix long names
revigo_df$description[grep("DNA damage response", revigo_df$description)] <- "DNA damage response"
revigo_df$description[grep("regulation of macrophage derived", revigo_df$description)] <- "macrophage derived foam cell differentiation"
revigo_df$description[grep("positive regulation of NF-kappaB", revigo_df$description)] <- "positive regulation of NF-kB TF activity"
revigo_df$description[grep("transmembrane receptor protein tyrosine", revigo_df$description)] <- "transmembrane receptor tyrosine kinase signaling"
revigo_df$description[grep("antigen processing and presentation", revigo_df$description)] <- "antigen processing and presentation of peptides"
revigo_df$description[grep("positive regulation of cytosolic calcium ion concentration", revigo_df$description)] <- "positive regulation of cytosolic calcium"

revigo_df$description <- as.factor(revigo_df$description)
revigo_df$description <- factor(revigo_df$description, levels = revigo_df$description)
enrplot <- ggplot(revigo_df, aes(x = description, y = neglog10padj, fill = uniqueness)) +
  geom_segment(aes(x = description,
    xend = description,
    y = min(neglog10padj),
    yend = max(neglog10padj)),
    linetype = "dashed",
    size = 0.2, color = "gray40") +
  geom_point(shape = 21, size = 4) +
  colorspace::scale_fill_continuous_sequential("Grays", name = "Uniqueness") +
  labs(x = "Description", y = "-log10(adj. P value)") +
  theme_classic(base_size = 14) +
  ggeasy::easy_all_text_color("black") +
  theme(
    axis.title.y = element_text(angle = 90, vjust = 0.5, margin = margin(0, 4, 0, 0), face = "bold"),
    axis.title.x = element_text(face = "bold"),
    axis.text.y = element_text(margin = margin(0, 4, 0, 0), size = 10),
    axis.ticks.y = element_blank(),
    axis.line = element_blank(),
    panel.border = element_rect(size =1, color = "black", fill = "transparent"),
    legend.position = "none",
    legend.text = element_text(size = 10),
    legend.justification = "top",
    legend.key.size = unit(0.5, "line")
  ) + coord_flip()
enrplot
cowplot::save_plot(plot = enrplot, filename = "plots/fig1_vargenes_enr.svg", base_height = 6, base_asp = 4/3)
cowplot::save_plot(plot = enrplot, filename = "plots/fig1_vargenes_enr.png", base_height = 6, base_asp = 4/3)

saveRDS(revigo_df, "data/vargenes_revigo.rds")
