#!/usr/bin/env Rscript
# =====================================================
# supermatrix_tree.R
#
# Assemble the concatenated supermatrix, write codon-position partitions for IQ-TREE, and
# (after IQ-TREE has been run) root, plot, and export the maximum-likelihood tree.
#
#   concat_loci(aln_list)              -> list(super = <DNAStringSet>, partitions = <data.frame>)
#   write_codon_partitions(part, file)  write a NEXUS charset block split by codon position
#   root_and_plot(treefile, ...)       -> rooted phylo; writes a PDF with branch-support labels
#   export_tree(rt, newick, support)    write the rooted tree as newick + a branch-support table
#
# IMPORTANT: this file intentionally does NOT wrap IQ-TREE. The ML search is run as a documented,
# commented bash step in 2_build_danaini_supermatrix_tree.Rmd so the heavy inference stays
# explicit and reproducible outside the knit.
#
# Provenance:
#   Written by N. Weksler for the consolidated Danaini tree-building workflow. Follows the
#   concatenation / codon-partition / Heliconius-rooting scheme used in Marianna Karageorgi's
#   original Danaini phylogeny pipeline; export_tree() (rooted newick + branch-support table) and
#   the on-figure support labels are new here.
# =====================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(Biostrings)
  library(ape)
  library(ggtree)
})

# Concatenate a named list of per-locus alignments (DNAStringSets) into one supermatrix.
# Species missing a locus get that block filled with gaps (treated as missing by IQ-TREE).
concat_loci <- function(aln_list) {
  aln_list <- aln_list[!vapply(aln_list, is.null, logical(1))]
  all_species <- sort(unique(unlist(lapply(aln_list, names))))
  superseq <- setNames(rep("", length(all_species)), all_species)
  parts <- data.frame(); pos <- 1
  for (nm in names(aln_list)) {
    aln <- aln_list[[nm]]; w <- unique(width(aln)); stopifnot(length(w) == 1)
    block <- setNames(rep(strrep("-", w), length(all_species)), all_species)
    block[names(aln)] <- as.character(aln)
    superseq <- paste0(superseq, block)   # NOTE: paste0 drops names (order preserved)
    parts <- rbind(parts, data.frame(locus = nm, start = pos, end = pos + w - 1))
    pos <- pos + w
  }
  super <- DNAStringSet(superseq); names(super) <- all_species   # restore names paste0 dropped
  list(super = super, partitions = parts)
}

# Write a NEXUS charset block that splits every locus into its three codon positions, so
# IQ-TREE can fit a separate model to each (with MFP+MERGE).
write_codon_partitions <- function(part, outfile) {
  con <- file(outfile, open = "wt"); on.exit(close(con))
  writeLines(c("#nexus", "begin sets;"), con)
  for (i in seq_len(nrow(part))) {
    s <- part$start[i]; e <- part$end[i]
    if (((e - s + 1) %% 3) != 0) stop("block not divisible by 3: ", part$locus[i])
    for (pos in 1:3)
      writeLines(paste0("  charset ", part$locus[i], "_pos", pos, " = ",
                        s + (pos - 1), "-", e, "\\3;"), con)
  }
  writeLines("end;", con)
}

# Root the ML tree on the chosen outgroup species and plot it, with the ultrafast-bootstrap
# branch-support values shown at the internal nodes. Returns the rooted tree.
root_and_plot <- function(treefile, out_pdf, title, root_taxon = ROOT_TAXON) {
  tr <- read.tree(treefile)
  rt <- if (root_taxon %in% tr$tip.label)
    root(tr, outgroup = root_taxon, resolve.root = TRUE) else { warning("root species absent"); tr }
  p <- ggtree(rt) +
    geom_tiplab(size = 2) +
    geom_nodelab(size = 1.8, hjust = -0.2, color = "#2c7fb8") +  # UFBoot support at each node
    ggtitle(title) + theme_tree2() + xlim(0, NA)
  ggsave(out_pdf, p, width = 10, height = 14, limitsize = FALSE)
  rt
}

# Export the rooted tree as newick, and write a table of internal-node branch-support values
# (the node labels carried over from IQ-TREE's ultrafast bootstrap).
export_tree <- function(rt, out_newick, out_support_csv) {
  write.tree(rt, file = out_newick)
  n_tip <- length(rt$tip.label)
  lab <- rt$node.label
  if (is.null(lab) || !length(lab)) lab <- rep(NA_character_, rt$Nnode)
  support <- tibble(
    node    = (n_tip + 1L):(n_tip + rt$Nnode),
    label   = lab,
    support = suppressWarnings(as.numeric(lab))
  )
  write_csv(support, out_support_csv)
  invisible(support)
}
