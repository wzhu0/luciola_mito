################################################################################
#
# R Script: Plotting time-calibrated MCC tree from RevBayes output.
#
# Usage: Rscript 05_revbayes/scripts/plot_rb_tree.R Luciola_mito_timetree_MCC.tree
# 
# Authors  : Wenjie Zhu
#
################################################################################

library(RevGadgets)
library(ggplot2)
library(ggtree)
library(ape)

options(ignore.negative.edge = TRUE)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) stop("Please provide a tree file")
tree_file <- args[1]

tree <- readTrees(paths = paste0("05_revbayes/output/", tree_file))

# Patch plotTreeFull to extend the x-axis to 80 Ma with 10 Ma intervals.
f <- RevGadgets:::plotTreeFull
src <- deparse(body(f))

src <- gsub("            max_age <- tree_height",
            "            max_age <- 81", src)
src <- gsub("            max_age <- max\\(minmax, na.rm = TRUE\\)",
            "            max_age <- 81", src)
src <- gsub("interval <- max_age/5",
            "interval <- 10", src)
src <- gsub("dx <- max_age%%interval",
            "dx <- 0", src)
src <- gsub("xlim = c\\(-max\\(minmax, na.rm = TRUE\\), tree_height/2\\)",
            "xlim = c(-81, tree_height/2)", src)
src <- gsub("        xline <- pretty\\(c\\(0, max_age\\)\\)\\[pretty\\(c\\(0, max_age\\)\\) < ",
            "        xline <- seq(0, 81, by = 10)[seq(0, 81, by = 10) < ", src)
src <- gsub('box_col <- "white"', 'box_col <- "grey97"', src)

body(f) <- parse(text = paste(src, collapse = "\n"))[[1]]
assignInNamespace("plotTreeFull", f, ns = "RevGadgets")

# Rescale branch length and 95% HPD age intervals to million years
tree[[1]][[1]]@phylo$edge.length <- tree[[1]][[1]]@phylo$edge.length / 1000
tree[[1]][[1]]@data$age_0.95_HPD <- lapply(
  tree[[1]][[1]]@data$age_0.95_HPD,
  function(x) if (is.null(x) || all(is.na(x))) x else x / 1000
)

# Extract population labels from tip names
tip_labels <- tree[[1]][[1]]@phylo$tip.label
pop_prefix <- sub("[-_].*", "", tip_labels)
# Merge AlVe and MoDo into one colour group
pop_prefix[pop_prefix %in% c("AlVe", "MoDo")] <- "AlVeMoDo"

tip_pop <- data.frame(
  label      = tip_labels,
  population = pop_prefix,
  stringsAsFactors = FALSE
)

color_table <- read.table("utils/color_code.txt", header = TRUE, sep = "\t",
                          stringsAsFactors = FALSE)
color_table$HEX <- paste0("#", gsub("^#", "", color_table$HEX))
pop_colors <- setNames(color_table$HEX, color_table$POP)

pop_colors["AlVeMoDo"] <- pop_colors["AlVe"]

pop_list   <- split(tip_pop$label, tip_pop$population)
tree_phylo <- tree[[1]][[1]]@phylo

n_tips_all <- length(tree_phylo$tip.label)
n_nodes    <- tree_phylo$Nnode
node_pop   <- rep(NA_character_, n_tips_all + n_nodes)
node_pop[1:n_tips_all] <- tip_pop$population[match(tree_phylo$tip.label, tip_pop$label)]

for (nd in (n_tips_all + 1):(n_tips_all + n_nodes)) {
  desc_tips <- which(tree_phylo$tip.label %in%
                       extract.clade(tree_phylo, nd)$tip.label)
  pops      <- unique(node_pop[desc_tips])
  node_pop[nd] <- if (length(pops) == 1) pops else "mixed"
}

# Plot the tree with original edges set to fully transparent.
p <- plotTree(tree = tree,
              timeline                     = TRUE,
              tip_labels_remove_underscore = FALSE,
              node_age_bars                = TRUE,
              line_width                   = 0.8,
              tip_labels                   = FALSE,
              branch_color                 = "#00000000")  # fully transparent

p$data$population <- node_pop[match(p$data$node, seq_along(node_pop))]

edge_df            <- as.data.frame(p$data)
edge_df$population <- node_pop[match(edge_df$node, seq_along(node_pop))]
parent_coords      <- edge_df[match(edge_df$parent, edge_df$node), c("x", "y")]
edge_df$xparent    <- parent_coords$x
edge_df$yparent    <- parent_coords$y
edge_df            <- edge_df[!is.na(edge_df$xparent), ]

h_seg <- geom_segment(data        = edge_df,
                      aes(x = xparent, y = y, xend = x, yend = y,
                          color = population),
                      linewidth   = 0.8,
                      inherit.aes = FALSE)

v_seg <- geom_segment(data        = edge_df,
                      aes(x = xparent, y = yparent, xend = xparent, yend = y,
                          color = population),
                      linewidth   = 0.8,
                      inherit.aes = FALSE)

layer_classes <- sapply(p$layers, function(l) class(l$geom)[1])
first_seg_idx <- which(layer_classes == "GeomSegment")[1]

p$layers <- c(
  p$layers[seq_len(first_seg_idx - 1)],
  list(h_seg, v_seg),
  p$layers[seq(first_seg_idx, length(p$layers))]
)

p <- p +
  scale_color_manual(values   = c(pop_colors, "mixed" = "#000000"),
                     na.value = "#000000",
                     guide    = "none")

p$coordinates$limits$x <- c(-81, 20)
p$scales$scales[[2]]$name <- "Age (Ma)                                   "

p <- p +
  geom_tiplab(aes(color = population),
              size = 1.8, fontface = "bold", geom = "text",
              offset = 0.5) +
  geom_nodelab(
    aes(label = ifelse(!isTip & !is.na(posterior) & as.numeric(posterior) < 0.99,
                       formatC(as.numeric(posterior), digits = 2, format = "f"),
                       "")),
    size = 1.8, fontface = "bold", nudge_x = -2, nudge_y = 1.0
  ) +
  theme(plot.margin = margin(t = -8, r = -10, b = 10, l = 10, unit = "mm"))

out_file <- paste0("05_revbayes/plots/", tree_file, ".pdf")
ggsave(out_file, p, width = 8.27, height = 11.69)