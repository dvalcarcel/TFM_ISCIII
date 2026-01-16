library(Seurat)
library(ggplot2)
library(dplyr)
library(patchwork)
setwd("PATH")


SpatialFeaturePlotBlend <- function(cells_obj, column_1, column_2,
                                    combine = TRUE,
                                    idents = NULL,
                                    show.background = FALSE,
                                    pt.size.factor = 0.1,
                                    images = NULL) {
  # --------------------------
  # Optional subsetting
  # --------------------------
  if (!is.null(idents)) {
    cells_obj <- subset(cells_obj, idents = idents)
  }
  
  # --------------------------
  # Helper functions
  # --------------------------
  as_hex <- function(num) {
    hex_str <- as.character(as.hexmode(num))
    if (nchar(hex_str) == 1) {
      hex_str <- paste0("0", hex_str)
    }
    return(hex_str)
  }
  
  metadata_to_hexadecimal <- function(in_dat) {
    apply(in_dat, 2, function(x) {
      x - min(x) # Make minimum 0
    }) %>%
      apply(2, function(x) {
        round(255 * (x / max(x))) # Scale to [0,255]
      }) %>%
      apply(1, function(x) {
        toupper(paste0("#", as_hex(x[1]), as_hex(x[2]), "00"))
      })
  }
  
  blend_plot_theme <- theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5)
  )
  
  # --------------------------
  # Feature plots
  # --------------------------
  bg_alpha <- ifelse(show.background, 1, 0)
  
  plot_list <- lapply(c(column_1, column_2), function(column) {
    max_color <- ifelse(column == column_1, "#FF0000", "#00FF00")
    SpatialFeaturePlot(
      cells_obj, column,
      image.alpha = bg_alpha,
      pt.size.factor = pt.size.factor,
      images = images
    ) +
      scale_fill_gradient(low = "#000000", high = max_color) +
      ggtitle(column) +
      blend_plot_theme
  })
  
  # --------------------------
  # Blend colors
  # --------------------------
  dat <- FetchData(cells_obj, c(column_1, column_2))
  colors <- as.matrix(dat) %>% metadata_to_hexadecimal()
  
  new_md_column <- paste0(column_1, "_vs_", column_2)
  cells_obj[[new_md_column]] <- colors
  names(colors) <- as.character(colors)
  
  plot_list[[3]] <- SpatialDimPlot(
    cells_obj, new_md_column, cols = colors,
    image.alpha = bg_alpha,
    pt.size.factor = pt.size.factor,
    images = images
  ) +
    ggtitle(paste0(column_1, "_", column_2)) +
    blend_plot_theme
  
  # --------------------------
  # Legend
  # --------------------------
  side_length <- 100
  legend_grid <- expand.grid(
    seq(from = min(dat[, column_1]),
        to = max(dat[, column_1]), length.out = side_length),
    seq(from = min(dat[, column_2]),
        to = max(dat[, column_2]), length.out = side_length)
  )
  colnames(legend_grid) <- c(column_1, column_2)
  legend_colors <- metadata_to_hexadecimal(legend_grid)
  legend_grid$color <- legend_colors
  names(legend_colors) <- legend_colors
  
  legend <- ggplot(legend_grid,
                   aes(x = .data[[column_1]], y = .data[[column_2]],
                       color = color)) +
    geom_point(shape = 15, size = 1.9) +
    scale_color_manual(values = legend_colors) +
    coord_cartesian(expand = FALSE) +
    theme(legend.position = "none", aspect.ratio = 1,
          panel.background = element_blank())
  
  plot_list[[4]] <- wrap_plots(
    ggplot() + theme_void(), legend, ggplot() + theme_void(),
    ncol = 1, heights = c(0.2, 0.6, 0.2)
  )
  
  # --------------------------
  # Return plots
  # --------------------------
  if (combine == FALSE) {
    return(plot_list)
  } else {
    p <- wrap_plots(plot_list, nrow = 1,
                    widths = c(0.28, 0.28, 0.28, 0.16))
    return(p)
  }
}



######## Patient 38 #############

load("Patient_38_RCTD_snRNAseq.rda")


png("Patient_38_SPP1_CD44.png", width = 9000, height = 3000, res = 800)
SpatialFeaturePlotBlend(
  cells_obj = Patient_38_RCTD,
  column_1 = "SPP1",
  column_2 = "CD44",
  idents = c("Myeloid", "T_cells", "Stromal", "tumor", "B_cells", "Endothelial_cells"),  # subset to selected clusters
  show.background = TRUE,           # hide histology background
  pt.size.factor = 9.2,             # smaller point size for VisiumHD
) 
dev.off()


######## Patient 06 #############

load("Patient_06_RCTD_snRNAseq.rda")


png("Patient_06_CXCL12_CXCR4.png", width = 9000, height = 3000, res = 800)
SpatialFeaturePlotBlend(
  cells_obj = Patient_06_RCTD,
  column_1 = "CXCL12",
  column_2 = "CXCR4",
  idents = c("Myeloid", "T_cells", "B_cells"),  # subset to selected clusters
  show.background = TRUE,           # hide histology background
  pt.size.factor = 9.2,             # smaller point size for VisiumHD
) 
dev.off()


######## Patient 61 #############

load("Patient_61_RCTD_snRNAseq.rda")


png("Patient_61_SPP1_CD44.png", width = 9000, height = 3000, res = 800)
SpatialFeaturePlotBlend(
  cells_obj = Patient_61_RCTD,
  column_1 = "SPP1",
  column_2 = "CD44",
  idents = c("Myeloid", "T_cells", "Stromal", "tumor", "B_cells", "Endothelial_cells"),  # subset to selected clusters
  show.background = TRUE,           # hide histology background
  pt.size.factor = 9.2,             # smaller point size for VisiumHD
) 
dev.off()


######## Patient 4 #############

load("Patient_4_26_RCTD_snRNAseq.rda")

Patient_4_RCTD <- subset(
  Patient_4_26_RCTD,
  subset = Patient == "ptx4"
)


png("Patient_4_SPP1_CD44.png", width = 9000, height = 3000, res = 800)
SpatialFeaturePlotBlend(
  cells_obj = Patient_4_RCTD,
  column_1 = "SPP1",
  column_2 = "CD44",
  idents = c("Myeloid", "T_cells", "Stromal", "tumor", "B_cells"),  # subset to selected clusters
  show.background = TRUE,           # hide histology background
  pt.size.factor = 9.2,             # smaller point size for VisiumHD
) 
dev.off()

######## Patient 41 #############

load("Patient_41_67_RCTD_snRNAseq.rda")

Patient_41_RCTD <- subset(
  Patient_41_67_RCTD,
  subset = Patient == "ptx41"
)


png("Patient_41_APP_CD74.png", width = 9000, height = 3000, res = 800)
SpatialFeaturePlotBlend(
  cells_obj = Patient_41_RCTD,
  column_1 = "APP",
  column_2 = "CD74",
  idents = c("Myeloid", "Endothelial_cells", "Stromal", "tumor"),  # subset to selected clusters
  show.background = TRUE,           # hide histology background
  pt.size.factor = 9.2,             # smaller point size for VisiumHD
) 
dev.off()


######## Patient 41 #############

load("Patient_41_67_RCTD_snRNAseq.rda")

Patient_67_RCTD <- subset(
  Patient_41_67_RCTD,
  subset = Patient == "ptx67"
)


png("Patient_67_SPP1_CD44.png", width = 9000, height = 3000, res = 800)
SpatialFeaturePlotBlend(
  cells_obj = Patient_67_RCTD,
  column_1 = "SPP1",
  column_2 = "CD44",
  idents = c("Myeloid", "tumor"),  # subset to selected clusters
  show.background = TRUE,           # hide histology background
  pt.size.factor = 9.2,             # smaller point size for VisiumHD
) 
dev.off()


