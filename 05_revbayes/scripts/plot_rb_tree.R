library(RevGadgets)
library(ggtree)
library(ggplot2)
library(dplyr)
library(tidytree)
library(ape)

options(ignore.negative.edge = TRUE)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) stop("Please provide a tree file")
tree_file <- args[1]

# tree_file <- "Luciola_mito_timetree_MAP_run_2.tree"
tree <- readTrees(paths = paste0("05_revbayes/output/", tree_file))

# Extract population prefix (everything before first - or _)
tip_labels <- tree[[1]][[1]]@phylo$tip.label
pop_prefix <- sub("[-_].*", "", tip_labels)

tip_pop <- data.frame(
  label      = tip_labels,
  population = pop_prefix,
  stringsAsFactors = FALSE
)

# Read color codes
color_table <- read.table("utils/color_code.txt", header = TRUE, sep = "\t",
                          stringsAsFactors = FALSE)
color_table$HEX <- paste0("#", gsub("^#", "", color_table$HEX))
pop_colors <- setNames(color_table$HEX, color_table$POP)

# Build population groupings and tag tips
pop_list   <- split(tip_pop$label, tip_pop$population)
tree_phylo <- tree[[1]][[1]]@phylo

n_tips_all <- length(tree_phylo$tip.label)
n_nodes    <- tree_phylo$Nnode
node_pop   <- rep(NA_character_, n_tips_all + n_nodes)
node_pop[1:n_tips_all] <- tip_pop$population[match(tree_phylo$tip.label, tip_pop$label)]

for (nd in (n_tips_all + 1):(n_tips_all + n_nodes)) {
  desc_tips <- which(tree_phylo$tip.label %in%
                       extract.clade(tree_phylo, nd)$tip.label)
  pops <- unique(node_pop[desc_tips])
  node_pop[nd] <- if (length(pops) == 1) pops else "mixed"
}

tree_grouped <- groupOTU(tree_phylo, pop_list, "population")

# X-axis breaks in Ma (tree is in kyr)
x_max        <- max(ape::node.depth.edgelength(tree_phylo))
breaks_kyr   <- seq(0, ceiling(x_max / 10000) * 10000, by = 10000)
break_pos    <- x_max - breaks_kyr
break_labels <- breaks_kyr / 1000

# Plot base timetree
p <- ggtree(tree_grouped, aes(color = population),
            ladderize = TRUE,
            linewidth = 1) +
  scale_color_manual(values = c(pop_colors, "0" = "#000000", "mixed" = "#000000")) +
  geom_tiplab(aes(color = population),
              size     = 3,
              fontface = "bold",
              align    = FALSE,
              offset   = 100) +
  theme_tree2() +
  theme(
    legend.position = "none",
    plot.margin     = margin(20, 20, 20, 20),
    axis.text.x     = element_text(size = 12, face = "bold"),
    axis.title.x    = element_text(size = 14, face = "bold")
  ) +
  scale_x_continuous(
    name   = "Age (Ma)",
    breaks = break_pos,
    labels = break_labels,
    expand = expansion(mult = c(0.02, 0.18))
  )

# Count descendant tips per node to distinguish deep vs shallow internal nodes
n_tips     <- length(tree_phylo$tip.label)
desc_count <- sapply(1:(n_tips + tree_phylo$Nnode), function(nd) {
  if (nd <= n_tips) return(1L)
  length(extract.clade(tree_phylo, nd)$tip.label)
})
node_desc  <- data.frame(node = 1:(n_tips + tree_phylo$Nnode), n_desc = desc_count)
tree_anno  <- left_join(as_tibble(tree[[1]][[1]]), node_desc, by = "node")

# Deep internal nodes (>40 descendant tips): prominent bold black labels
p <- p %<+% tree_anno +
  geom_nodelab(
    aes(label = ifelse(!isTip & !is.na(posterior) & posterior < 0.95 & n_desc > 40,
                       round(posterior, 2),
                       NA_character_)),
    size     = 3,
    fontface = "bold",
    nudge_x  = 300,
    nudge_y  =  0,
    na.rm    = TRUE,
    color    = "black"
  ) +
  # Shallow internal nodes (<=40 descendant tips): subtle grey labels
  geom_nodelab(
    aes(label = ifelse(!isTip & !is.na(posterior) & posterior < 0.95 & n_desc <= 40,
                       round(posterior, 2),
                       NA_character_)),
    size     = 2,
    fontface = "plain",
    nudge_x  = -500,
    nudge_y  =  0,
    na.rm    = TRUE,
    color    = "grey40"
  )

out_file <- paste0("05_revbayes/output/", tree_file, "_plot.pdf")
ggsave(out_file, p, width = 16, height = 24)
