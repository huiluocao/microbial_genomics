############################
## Packages
############################
# install.packages(c(
#   "tidyverse", "ape", "ggplot2", "patchwork", "cowplot",
#   "circlize", "png", "grid", "scales"
# ))
# BiocManager::install(c("ggtree", "ggtreeExtra"))

library(tidyverse)
library(ape)
library(ggplot2)
library(patchwork)
library(cowplot)
library(circlize)
library(png)
library(grid)
library(scales)
library(ggtree)
library(ggtreeExtra)
library(ggnewscale)

############################
## Colors approximated from image
############################
year_cols <- c(
  "2021" = "black",
  "Previous" = "grey90"
)

fraction_cols <- c(
  "SW" = "#26B6E3",   # seawater
  "O"  = "#E448E8"    # oyster
)

phylo_cols <- c(
  "Yes" = "red",
  "No"  = "white"
)

clade_cols <- c(
  "V1" = "#35D8D8",
  "V2" = "#F0EA00",
  "V3" = "#8E0F9B",
  "V4" = "#FF8C00",
  "V4 sister" = "#F2C65D",
  "V5" = "#980000",
  "V6" = "#8E9200",
  "V7" = "#BDBDBD",
  "V8" = "#4D9800",
  "V9" = "#E013FF",
  "V10" = "#DCCAA0",
  "V11" = "#FF160F",
  "V12" = "#102CFF",
  "V13" = "#178B8B",
  "Not in clade" = "grey90"
)

rgp_cols <- c(
  "Out of RGP"    = "grey30",
  "Unclassified"  = "#BFD37A",
  "Plasmid"       = "#7ED398",
  "Phage"         = "#7EA6FF",
  "Satellite"     = "#D58BE8",
  "+ Integrase"   = "#C46AE7"
)

############################
## Theme
############################
theme_fig <- function() {
  theme_classic(base_size = 16) +
    theme(
      axis.text = element_text(color = "black"),
      axis.title = element_text(color = "black"),
      legend.title = element_text(face = "bold"),
      legend.text = element_text(size = 11),
      plot.tag = element_text(face = "bold", size = 24),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold")
    )
}

############################
## Helper for p-value bracket
############################
add_p_bracket <- function(p, x1, x2, y, label, h = 0.05) {
  p +
    annotate("segment", x = x1, xend = x2, y = y, yend = y) +
    annotate("segment", x = x1, xend = x1, y = y, yend = y - h) +
    annotate("segment", x = x2, xend = x2, y = y, yend = y - h) +
    annotate("text", x = mean(c(x1, x2)), y = y + h * 0.35, label = label, size = 6)
}

############################
## Panel A: phylogeny + side annotations
############################
make_panel_A <- function(tree_file, meta_file) {
  tr <- read.tree(tree_file)
  meta <- read_csv(meta_file, show_col_types = FALSE) %>%
    mutate(
      year_group = factor(year_group, levels = c("2021", "Previous")),
      fraction = factor(fraction, levels = c("SW", "O")),
      phylopathotype = factor(phylopathotype, levels = c("Yes", "No")),
      clade = factor(clade, levels = names(clade_cols))
    )
  
  p <- ggtree(tr, size = 0.35) %<+% meta +
    theme_tree2() +
    xlim(0, max(ggtree(tr)$data$x, na.rm = TRUE) + 1.8) +
    labs(tag = "A") +
    theme(
      legend.position = "left",
      plot.margin = margin(5, 5, 5, 5)
    )
  
  # Year strip
  p <- p +
    geom_fruit(
      data = meta,
      geom = geom_tile,
      mapping = aes(y = tip, x = 1, fill = year_group),
      width = 0.06,
      offset = 0.02
    ) +
    scale_fill_manual(values = year_cols, name = "Year of isolation") +
    ggnewscale::new_scale_fill()
  
  # Fraction strip
  p <- p +
    geom_fruit(
      data = meta,
      geom = geom_tile,
      mapping = aes(y = tip, x = 1, fill = fraction),
      width = 0.06,
      offset = 0.10
    ) +
    scale_fill_manual(values = fraction_cols, name = "Fraction") +
    ggnewscale::new_scale_fill()
  
  # Phylopathotype strip
  p <- p +
    geom_fruit(
      data = meta,
      geom = geom_tile,
      mapping = aes(y = tip, x = 1, fill = phylopathotype),
      width = 0.06,
      offset = 0.18,
      color = "black",
      size = 0.1
    ) +
    scale_fill_manual(values = phylo_cols, name = "Phylopathotype") +
    ggnewscale::new_scale_fill()
  
  # Clade strip
  p <- p +
    geom_fruit(
      data = meta,
      geom = geom_tile,
      mapping = aes(y = tip, x = 1, fill = clade),
      width = 0.09,
      offset = 0.26
    ) +
    scale_fill_manual(values = clade_cols, name = "Clade") +
    ggnewscale::new_scale_fill()
  
  # Genome size strip
  p <- p +
    geom_fruit(
      data = meta,
      geom = geom_tile,
      mapping = aes(y = tip, x = 1, fill = genome_size_mbp),
      width = 0.09,
      offset = 0.38
    ) +
    scale_fill_gradientn(
      colours = c("#7BBF6A", "#E7E36B", "#F3A35D", "#EA6A63"),
      name = "Genome size\nMbp"
    ) +
    theme(
      legend.box = "vertical",
      legend.position = "left"
    )
  
  p
}

############################
## Panel B: boxplot genome size
############################
make_panel_B <- function(panelB_file) {
  df <- read_csv(panelB_file, show_col_types = FALSE) %>%
    mutate(fraction = factor(fraction, levels = c("SW", "O")))
  
  p <- ggplot(df, aes(fraction, genome_size_mbp, fill = fraction)) +
    geom_boxplot(width = 0.55, outlier.shape = NA, color = "black") +
    geom_jitter(width = 0.12, size = 1.2, alpha = 0.7, color = "black") +
    scale_fill_manual(values = fraction_cols) +
    labs(x = NULL, y = "Genome size (Mbp)", tag = "B") +
    coord_cartesian(ylim = c(4.5, 7.1)) +
    theme_fig() +
    theme(
      legend.position = "none",
      axis.text.x = element_text(angle = 45, hjust = 1, size = 20),
      axis.text.y = element_text(size = 16),
      axis.title.y = element_text(size = 20)
    )
  
  add_p_bracket(p, 1, 2, y = 7.05, label = "P = 0.0004", h = 0.05)
}

############################
## Panel C: three boxplots
############################
make_panel_C_one <- function(df, yvar, ylab, p_label = "P < 0.0001", tag = NULL, ymax = NULL) {
  p <- ggplot(df, aes(fraction, .data[[yvar]], fill = fraction)) +
    geom_boxplot(width = 0.55, outlier.shape = NA, color = "black") +
    geom_jitter(width = 0.12, size = 1.1, alpha = 0.55, color = "grey30") +
    scale_fill_manual(values = fraction_cols) +
    labs(x = NULL, y = ylab, tag = tag) +
    theme_fig() +
    theme(
      legend.position = "none",
      axis.text.x = element_text(angle = 45, hjust = 1, size = 20),
      axis.title.y = element_text(size = 20)
    )
  
  if (!is.null(ymax)) {
    p <- p + coord_cartesian(ylim = c(0, ymax))
    p <- add_p_bracket(p, 1, 2, y = ymax * 0.98, label = p_label, h = ymax * 0.03)
  } else {
    yy <- max(df[[yvar]], na.rm = TRUE)
    p <- add_p_bracket(p, 1, 2, y = yy * 1.03, label = p_label, h = yy * 0.03)
  }
  
  p
}

make_panel_C <- function(panelC_file) {
  df <- read_csv(panelC_file, show_col_types = FALSE) %>%
    mutate(fraction = factor(fraction, levels = c("SW", "O")))
  
  p1 <- make_panel_C_one(df, "plasmid_regions", "Plasmid regions", tag = "C", ymax = 30.5)
  p2 <- make_panel_C_one(df, "phage_regions", "Phage regions", ymax = 17.5)
  p3 <- make_panel_C_one(df, "satellite_regions", "Satellite regions", ymax = 3.1)
  
  plot_grid(p1, p2, p3, nrow = 1, align = "h")
}

############################
## Panel D: circular chord diagram
## Note: this is an approximation of the visual style
############################
make_panel_D_png <- function(nodes_file, links_file, out_png = "panelD.png") {
  nodes <- read_csv(nodes_file, show_col_types = FALSE) %>%
    mutate(
      clade = factor(clade, levels = names(clade_cols)),
      fraction = factor(fraction, levels = c("SW", "O"))
    )
  
  links <- read_csv(links_file, show_col_types = FALSE)
  
  # Order genomes by clade then fraction
  nodes <- nodes %>%
    arrange(clade, fraction, genome)
  
  grid.col <- setNames(clade_cols[as.character(nodes$clade)], nodes$genome)
  link.col <- alpha(clade_cols[links$source_clade], 0.45)
  
  png(out_png, width = 1800, height = 1800, res = 250, bg = "white")
  par(mar = c(1, 1, 1, 1))
  
  circos.clear()
  circos.par(
    start.degree = 110,
    gap.after = c(rep(1, nrow(nodes) - 1), 8),
    track.margin = c(0.005, 0.005),
    cell.padding = c(0, 0, 0, 0)
  )
  
  chordDiagram(
    x = links %>% select(from, to, weight),
    order = nodes$genome,
    grid.col = grid.col,
    col = link.col,
    transparency = 0.25,
    annotationTrack = "grid",
    preAllocateTracks = list(track.height = 0.04)
  )
  
  # Outer mini-track for fraction
  circos.trackPlotRegion(
    track.index = 1,
    ylim = c(0, 1),
    bg.border = NA,
    panel.fun = function(x, y) {
      sec <- get.cell.meta.data("sector.index")
      xlim <- get.cell.meta.data("xlim")
      frac <- nodes$fraction[match(sec, nodes$genome)]
      circos.rect(
        xleft = xlim[1],
        ybottom = 0.72,
        xright = xlim[2],
        ytop = 1,
        col = fraction_cols[as.character(frac)],
        border = NA
      )
    }
  )
  
  title("D", adj = 0.02, cex.main = 2.2, font.main = 2)
  dev.off()
}

load_panel_D_as_gg <- function(png_file = "panelD.png") {
  ggdraw() +
    draw_image(png_file) +
    draw_plot_label("D", x = 0.02, y = 0.98, hjust = 0, vjust = 1,
                    fontface = "bold", size = 24)
}

############################
## Panel E: stacked log-scale boxplots
############################
make_panel_E_one <- function(df, ylab, letters_vec, show_x = TRUE, tag = NULL) {
  rgp_order <- c("Out of RGP", "Unclassified", "Plasmid", "Phage", "Satellite", "+ Integrase")
  
  df <- df %>%
    mutate(rgp_type = factor(rgp_type, levels = rgp_order))
  
  ymax <- max(df$value, na.rm = TRUE)
  
  letter_df <- tibble(
    rgp_type = factor(rgp_order, levels = rgp_order),
    lab = letters_vec,
    y = ymax * 1.35
  )
  
  p <- ggplot(df, aes(rgp_type, value, fill = rgp_type)) +
    geom_boxplot(width = 0.6, outlier.shape = NA, color = "black") +
    geom_jitter(width = 0.15, size = 1, alpha = 0.45, color = "black") +
    geom_text(data = letter_df, aes(rgp_type, y, label = lab),
              inherit.aes = FALSE, size = 6) +
    scale_fill_manual(values = rgp_cols) +
    scale_y_log10() +
    labs(x = if (show_x) "Type of RGP" else NULL, y = ylab, tag = tag) +
    theme_fig() +
    theme(
      legend.position = "none",
      axis.text.x = if (show_x) element_text(angle = 45, hjust = 1) else element_blank(),
      axis.ticks.x = if (show_x) element_line() else element_blank()
    )
  
  p
}

make_panel_E <- function(defense_file, metal_file) {
  defense <- read_csv(defense_file, show_col_types = FALSE)
  metal <- read_csv(metal_file, show_col_types = FALSE)
  
  p_top <- make_panel_E_one(
    defense,
    ylab = expression("Defense systems kbp"^-1),
    letters_vec = c("a", "b", "b", "c", "d", "e"),
    show_x = FALSE,
    tag = "E"
  )
  
  p_bottom <- make_panel_E_one(
    metal,
    ylab = expression("Metal resistance genes kbp"^-1),
    letters_vec = c("a", "b", "c", "d", "e", "e"),
    show_x = TRUE
  )
  
  plot_grid(p_top, p_bottom, ncol = 1, align = "v", rel_heights = c(1, 1))
}

############################
## Build everything
############################
pA <- make_panel_A("tree.nwk", "meta_tree.csv")
pB <- make_panel_B("panelB.csv")
pC <- make_panel_C("panelC.csv")

make_panel_D_png("panelD_nodes.csv", "panelD_links.csv", out_png = "panelD.png")
pD <- load_panel_D_as_gg("panelD.png")

pE <- make_panel_E("panelE_defense.csv", "panelE_metal.csv")

############################
## Assemble full figure
## Layout approximates the image:
## A | (B + D over C) | E
############################
center_top <- plot_grid(pB, pD, nrow = 1, rel_widths = c(0.8, 1.9), align = "h")
center_col <- plot_grid(center_top, pC, ncol = 1, rel_heights = c(1.15, 1.0))

final_fig <- plot_grid(
  pA, center_col, pE,
  nrow = 1,
  rel_widths = c(1.2, 2.2, 1.0),
  align = "h"
)

ggsave("multi_panel_figure.pdf", final_fig, width = 18, height = 9)
ggsave("multi_panel_figure.png", final_fig, width = 18, height = 9, dpi = 300)
