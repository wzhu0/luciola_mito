################################################################################
# Preprocess REvBayes unrooted MCC tree for mPTP:
#   1. Reroot using outgroup
#   2. Convert to Newick format
#   3a. Write full tree (all samples)
#   3b. Drop GenBank tips, write field-only tree
#
# Usage: Rscript 06_mptp/scripts/01_preprocess_tree.R <tree_file>
#
# Example:
#   Rscript 06_mptp/scripts/01_preprocess_tree.R 05_revbayes/output/Luciola_mito_unrooted_non-clock_MCMC150000_MCC.tree 
#
################################################################################

library(ape)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) stop("Usage: Rscript 01_preprocess_tree.R <tree_file>")

tree_file     <- args[1]
outgroup_file <- "utils/outgroup_list.txt"
out_dir       <- "06_mptp/data"

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

basename <- sub("\\.tree$", "", basename(tree_file))
out_full  <- file.path(out_dir, paste0(basename, "_rooted.nwk"))
out_field <- file.path(out_dir, paste0(basename, "_rooted_field_only.nwk"))

# read tree
cat("Reading tree:", tree_file, "\n")
tree <- read.nexus(tree_file)
cat("Tips:", length(tree$tip.label), "\n")

# read outgroups
outgroups <- readLines(outgroup_file)
outgroups <- outgroups[outgroups != ""]

# only keep outgroups that exist in tree
outgroups_present <- intersect(outgroups, tree$tip.label)
cat("Outgroups found in tree:", length(outgroups_present), "/", length(outgroups), "\n")
if (length(outgroups_present) == 0) stop("No outgroups found in tree")

# reroot
cat("Rerooting...\n")
outgroup_node <- getMRCA(tree, outgroups_present)
tree_rooted <- root(tree, node = outgroup_node, resolve.root = TRUE)

# remove root edge length (mPTP cannot handle root branch lengths)
tree_rooted$root.edge <- NULL

# write full rooted newick
write.tree(tree_rooted, file = out_full)
cat("Full rooted tree written:", out_full, "\n")

# drop GenBank tips (identified by absence of hyphen in name)
genbank_tips <- tree_rooted$tip.label[!grepl("-", tree_rooted$tip.label)]
cat("Pruning", length(genbank_tips), "GenBank tips:\n")
cat(paste(" ", genbank_tips, collapse = "\n"), "\n")

tree_field <- drop.tip(tree_rooted, genbank_tips)
write.tree(tree_field, file = out_field)
cat("Field-only rooted tree written:", out_field, "\n")