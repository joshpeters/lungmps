# ---
# Description: Plotting helper functions used in this project.
# Author: Josh Peters
# ---

PlotUMAP <- function(object, column) {
  umap_df <- Embeddings(object, reduction = "harmony_umap")
  meta_df <-  object[[]] %>% select(cell_barcode, batch, !!sym(column))
  assertthat::assert_that(all.equal(rownames(umap_df), rownames(meta_df)))
  meta_df <- merge(meta_df, umap_df, by = 0)
  meta_df <- meta_df[order(meta_df[, column]), ]
  if (label == TRUE)
    a <- ggplot(meta_df, aes(x = hUMAP_1, y = hUMAP_2, fill = !!sym(column))) +
    geom_point(shape = 21, size = 2, alpha = 0.9) +
    labs(x = "UMAP 1", y = "UMAP 2") +
    UMAPTheme()
  return(a)
}

#' Plots embedding of data
#'
#' Plot a manifold representation of a Seurat object calculated from RunUMAP, RunTSNE, RunPCA
#'
#' @param object Seurat object
#' @param group Grouping variable
#' @param reduction Reduction to use
#' @param label Label the centroids of each group
#' @param legend Include a legend
#' @param pt.size Size of points in embedding
#' @param pt.line Size of outline in points
#' @param pt.alpha Alpha of points in embedding
#' @param base.size Base size of text
#' @param base.pt Base point size
#' @param label.size Size of label if label == TRUE
#' @usage PlotEmbedding(object, group, reduction, label, legend, pt.size, pt.line, pt.alpha,
#'     base.size, base.pt, label.size)
#' @return list including the all variable features, top x% genes, and consensus variable features
#'
#' @importFrom dplyr group_by summarise
#' @importFrom ggrepel geom_label_repel
#' @importFrom rlang .data
#' @export

PlotEmbedding <- function(
  object,
  group = "orig.ident",
  reduction = "umap",
  label = TRUE,
  legend = TRUE,
  pt.size = 1,
  pt.line = 0.95,
  pt.alpha = 0.8,
  base.size = 12,
  base.pt = 1,
  label.size = 4
) {

  embeddings <- as.data.frame(object@reductions[[reduction]]@cell.embeddings)

  stopifnot(group %in% names(object[[]]))

  # prepare plot data
  colnames(embeddings) <- c("UMAP1", "UMAP2")
  embeddings$group <- object[[group]][[1]]
  embeddings %>%
    dplyr::group_by(group) %>%
    summarise(UMAP1 = median(UMAP1), UMAP2 = median(UMAP2)) -> centers

  # generate base plot
  embedding_plot <- ggplot(embeddings, aes(UMAP1, UMAP2)) +
    geom_point(size = pt.size*1.25, color = "black",
      stroke = pt.line, alpha = pt.alpha) +
    geom_point(size = pt.size, color = "white", alpha = 1) +
    geom_point(aes(color = group), size = pt.size, alpha = pt.alpha)

  if (label) {
    embedding_plot <- embedding_plot +
      ggrepel::geom_label_repel(centers, mapping = aes(label = group),
        fontface = "bold", label.size = NA, alpha = 0.8, size = label.size, label.r = 0.25,
        label.padding = 0.2, na.rm = TRUE, seed = 17) +
      ggrepel::geom_label_repel(centers, mapping = aes(label = group),
        fontface = "bold", label.size = NA, alpha = 1, size = label.size, label.r = 0.25,
        label.padding = 0.2, na.rm = TRUE, fill = NA, seed = 17)
  }

  embedding_plot <- embedding_plot + theme(
    # plot.title = element_text(size=base.size*1.5, face="bold", hjust=0,
    #   margin=margin(b = base.pt*5)),
    # plot.subtitle = element_text(size=base.size, color="grey20", hjust=0,
    #   margin = margin(b=base.pt*5)),
    # plot.caption = element_text(size=base.size, color="grey20",
    #   margin=margin(t=base.pt*5)),
    axis.title.y = element_text(size=base.size*1.25,
      margin=margin(r=base.pt*5, l=base.pt*5)),
    axis.title.x = element_text(size=base.size*1.25,
      margin=margin(t=base.pt*5, b=base.pt*5)),
    axis.line = element_line(size = 0.75),
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks = element_blank(),
    panel.background = element_blank(),
    #panel.border = element_rect(size = 1.5, fill = NA, colour = "black"),
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.background = element_blank(),
    plot.margin = unit(rep(base.pt*1, 4), "pt"))

  if (!legend) {
    embedding_plot <- embedding_plot + theme(legend.position = "none")
  } else {
    embedding_plot <- embedding_plot + theme(
      legend.title  =  element_text(size = base.size*1.25, face = "bold", hjust = 0.5),
      legend.position = "right",
      legend.text = element_text(size = base.size*1.25),
      legend.key.size = unit(base.pt*1.25, 'lines'),
      legend.background = element_rect(fill = "transparent", color = NA),
      legend.key = element_rect(fill = "transparent", color = NA)
    )
  }

  return(embedding_plot)
}

PlotAnnotatedClusters <- function(object, group.by = NULL) {
  if (is.null(group.by)) {
    num_current_clusters <- length(unique(object@active.ident))
    group_name <- "current identity"
  } else{
    num_current_clusters <- length(unique(object[[group.by, drop = TRUE]]))
    group_name <- group.by
  }
  colors <- colorspace::qualitative_hcl("Dark3", n = num_current_clusters)
  colors_1 <- colors[c(TRUE, FALSE)]
  colors_2 <- rev(colors[c(FALSE, TRUE)])
  colors_use <- c(rbind(colors_1, colors_2))
  a <- DimPlot(object, reduction = "harmony_umap", group.by = group.by,
    label = TRUE, repel = FALSE, label.size = 6) +
    labs(x = "UMAP1", y = "UMAP2", title = glue::glue("UMAP Embedding clustered by Leiden ({group_name})")) +
    scale_color_manual(values = colors_use, name = "Cluster") +
    theme(
      plot.title = element_text(face = "plain", hjust = 0.5, margin = margin(t = 8, b = 8, unit = "pt")),
      axis.title.x = element_text(margin = margin(t = 8, unit = "pt")),
      axis.title.y = element_text(margin = margin(r = 8, unit = "pt")),
      plot.margin = unit(rep(8, 4), "pt"),
      axis.ticks = element_blank(),
      axis.text = element_blank()
    ) +
    ggeasy::easy_all_text_color("black")
  return(a)
}

PointTheme <- function(base_size) {
  theme_classic(base_size = base_size) +
    ggeasy::easy_all_text_color("black") +
    theme(
      axis.title.x = element_text(face = "bold", margin = margin(12, 0, 0, 0)),
      axis.title.y = element_text(face = "bold", margin = margin(0, 12, 0, 0)),
      axis.line = element_blank(),
      plot.title = element_text(size = base_size, margin = margin(1,1,4,1), face = "bold"),
      plot.subtitle = element_text(size = base_size-2, color = "gray40", margin = margin(1,1,2,1)),
      panel.border = element_rect(size = 1, color = "black", fill = "transparent"),
      panel.grid.major = element_line(size = 0.5, color = "gray90"),
      panel.grid.minor = element_line(size = 0.5, color = "gray90"),
      legend.position = "right",
      legend.justification = "top"
    )
}

DensityTheme <- function(base_size) {
  theme_classic(base_size = base_size) +
    ggeasy::easy_all_text_color("black") +
    theme(
      axis.title.x = element_text(face = "bold", margin = margin(12, 0, 0, 0)),
      axis.title.y = element_text(face = "bold", margin = margin(0, 12, 0, 0)),
      axis.line = element_blank(),
      plot.title = element_text(size = base_size, margin = margin(1,1,2,1)),
      plot.subtitle = element_text(size = base_size-2, color = "gray40", margin = margin(1,1,2,1)),
      panel.border = element_rect(size = 1, color = "black", fill = "transparent"),
      panel.grid.major = element_line(size = 0.5, color = "gray90"),
      panel.grid.minor = element_line(size = 0.5, color = "gray90"),
      legend.position = "right",
      legend.justification = "top"
    )
}

ViolinTheme <- function(base_size = 14) {
  theme_classic(base_size = 14) +
    ggeasy::easy_all_text_color("black") +
    theme(
      axis.title.x = element_text(face = "bold", margin = margin(12, 0, 0, 0)),
      axis.title.y = element_text(face = "bold", margin = margin(0, 12, 0, 0)),
      panel.border = element_rect(size = 1, color = "black", fill = "transparent"),
      plot.title = element_text(size = base_size, margin = margin(1,1,2,1)),
      plot.subtitle = element_text(size = base_size-2, color = "gray40", margin = margin(1,1,2,1)),
      #panel.grid.major = element_line(size = 0.5, color = "gray90"),
      #panel.grid.minor = element_line(size = 0.5, color = "gray90"),
      axis.line = element_blank(),
      legend.position = "right",
      validate = TRUE
    )
}

UMAPTheme <- function(base_size) {
  theme_void(base_size = base_size) +
    theme(
      plot.title = element_text(hjust = 0, face = "plain", margin = margin(0, 0, 12, 0, unit = "pt")),
      plot.subtitle = element_text(hjust = 0, face = "plain", margin = margin(0, 0, 12, 0, unit = "pt")),
      plot.caption = element_text(hjust = 1, vjust = 1, face = "plain", margin = margin(0, 0, 0, 0, unit = "pt")),
      axis.title.x = element_text(hjust = 0, face = "bold", margin = margin(12, 0, 0, 0, unit = "pt")),
      axis.title.y = element_text(hjust = 0, angle = 90, face = "bold", margin = margin(0, 12, 0, 0, unit = "pt")),
      panel.background = element_rect(fill = "transparent", color = "black", size = 1),
      plot.background = element_rect(fill = "transparent", color = "transparent", size = 1),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.background = element_rect(fill = "transparent", color = "transparent"),
      legend.box.background = element_rect(fill = "transparent", color = "transparent"),
      legend.position = "right",
      legend.title = element_text(face = "bold"),
      validate = TRUE
    )
}

TileTheme <- function(...) {
  theme_classic(...) +
    ggeasy::easy_all_text_color("black") +
    theme(
      axis.title.y = element_text(angle = 90, vjust = 0.5, margin = margin(0, 8, 0, 0), face = "bold"),
      axis.title.x = element_text(face = "bold"),
      axis.text.y = element_text(margin = margin(0, 4, 0, 0)),
      axis.ticks = element_blank(),
      axis.line = element_blank(),
      panel.border = element_blank(),
      #panel.border = element_rect(size = 1, color = "black", fill = "transparent"),
      #legend.text = element_text(size = 10),
      legend.justification = "top",
      legend.key.size = unit(1, "line"),
      validate = TRUE
    )
}

LoadColors <- function() {
  cb_colors <- ggthemes::colorblind_pal()(8)
  scales::show_col(cb_colors)
}

RemoveY <- function() {
  theme(axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title.y = element_blank(),
    axis.title.x = element_text(margin = margin(12, 0, 0, 0), face = "bold"),
    axis.line = element_blank(),
    panel.border = element_rect(size = 1, color = "black", fill = "transparent"),
    plot.margin = margin(12, 12, 12, 12))
}

SavePlot <- function(filename, plot, base_asp = 1, base_height = 6) {
  cowplot::save_plot(glue::glue(filename), plot, base_asp = base_asp, base_height = base_height)
}

GeneralTheme <- function(base_size, ...) {
  theme_classic(base_size = base_size, ...) +
    ggeasy::easy_all_text_color("black") +
    theme (
      axis.title.x = element_text(face = "bold", margin = margin(12, 0, 0, 0)),
      axis.title.y = element_text(face = "bold", margin = margin(0, 12, 0, 0), angle = 90, hjust = 0.5, vjust = 0.5),
      axis.line = element_blank(),
      plot.title = element_text(size =  base_size, color = "black", face = "bold", margin = margin(0,0,4,0)),
      plot.subtitle = element_text(size = base_size - 2, color = "black", margin = margin(0,0,4,0)),
      panel.background = element_rect(fill = "transparent", color = "black", size = 0.5),
      plot.background = element_rect(fill = "transparent", color = "transparent", size = 0),
      legend.title = element_text(size = base_size, face = "bold"),
      legend.text = element_text(size = base_size - 2),
      legend.background = element_rect(fill = "transparent", color = "transparent"),
      legend.box.background = element_rect(fill = "transparent", color = "transparent"),
      legend.position = "right",
      legend.justification = "top",
      legend.key.size = unit(1, "line"),
      validate = TRUE
    )
}
