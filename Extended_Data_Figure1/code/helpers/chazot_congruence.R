#!/usr/bin/env Rscript
# =====================================================
# chazot_congruence.R
#
# Compare our ML tree against the Chazot et al. 2021 reference tree and draw a congruence figure.
#
#   load_chazot(nexus_path)                 -> list(ref = <phylo>, ref_sp = <Genus_species vec>)
#   congruence_vs_chazot(tree, chazot, ...) -> list(plot, pct_recovered, n_shared, rf, rf_norm,
#                                                    node_df, our_pruned, ref_pruned)
#
# The figure: our tree pruned to the species shared with Chazot, tips ordered to match Chazot
# (node rotations only, via ape::rotateConstr — topology unchanged), internal nodes coloured by
# whether that exact clade is recovered in Chazot and sized by bootstrap support.
#
# Provenance:
#   Adapted by N. Weksler from Marianna Karageorgi's congruence_ge1_vs_chazot.R and the
#   congruence_plot() style in 5_validation_viz.R (Extended_data_fig1/code). Packaged here as
#   reusable functions that return the plot and the congruence statistics.
# =====================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(ape)
  library(ggtree)
  library(phangorn)   # RF.dist
})
options(ignore.negative.edge = TRUE)

COL_OK  <- "#1b9e77"   # colorblind-safe Dark2 pair (matches 5_validation_viz)
COL_BAD <- "#d95f02"

# A Chazot tip label like "Danaus_plexippus_specimenCode" -> "Danaus_plexippus" (or NA).
parse_sp <- function(lbl) {
  toks <- strsplit(lbl, "_")[[1]]
  gi <- which(grepl("^[A-Z][a-z]{2,}$", toks))
  if (!length(gi)) return(NA_character_)
  i <- gi[1]
  if (i < length(toks) && grepl("^[a-z]+$", toks[i + 1]))
    paste(toks[i], toks[i + 1], sep = "_") else NA_character_
}

# Read the Chazot NEXUS tree and parse a Genus_species label per tip.
load_chazot <- function(nexus_path) {
  ref <- read.nexus(nexus_path)
  if (inherits(ref, "multiPhylo")) ref <- ref[[1]]
  ref_sp <- vapply(ref$tip.label, parse_sp, character(1))
  list(ref = ref, ref_sp = ref_sp)
}

# tip labels of `phy` in plotting order (bottom -> top)
.ggtree_tip_order <- function(phy) {
  d <- ggtree::fortify(phy); d <- d[d$isTip, ]; d$label[order(d$y)]
}

#' Congruence of `tree` with the Chazot reference.
#' @param tree       our phylo (will be rooted on root_taxon if present)
#' @param chazot     output of load_chazot()
#' @param title      plot title
#' @param root_taxon outgroup species to root on
congruence_vs_chazot <- function(tree, chazot, title = "Congruence with Chazot et al. 2021",
                                 root_taxon = "Heliconius_erato") {
  ref <- chazot$ref; ref_sp <- chazot$ref_sp
  if (root_taxon %in% tree$tip.label)
    tree <- root(tree, outgroup = root_taxon, resolve.root = TRUE)

  shared <- intersect(tree$tip.label, unique(na.omit(ref_sp)))
  if (length(shared) < 3) stop("Fewer than 3 species shared with Chazot — cannot compare.")

  # prune Chazot to the shared species (relabel specimen tips to Genus_species)
  idx <- match(shared, ref_sp); idx <- idx[!is.na(idx)]
  rp <- keep.tip(ref, ref$tip.label[idx])
  rp$tip.label <- ref_sp[match(rp$tip.label, ref$tip.label)]
  mt <- keep.tip(tree, shared); rp <- keep.tip(rp, shared)

  # order our tips to match Chazot (rotations only), reflect at draw time if still inverted
  ref_order <- .ggtree_tip_order(ladderize(rp))
  rank_of   <- setNames(seq_along(ref_order), ref_order)
  mt <- rotateConstr(mt, mt$tip.label[order(rank_of[mt$tip.label])])
  drawn <- .ggtree_tip_order(mt)
  rho   <- suppressWarnings(cor(match(drawn, drawn), rank_of[drawn], method = "spearman"))
  reflect <- !is.na(rho) && rho < 0

  # per-node congruence (is this clade recovered in Chazot?) + bootstrap support.
  # Compare unrooted bipartitions directly instead of via prop.clades(): against a
  # single reference tree prop.clades() returns NA (never 0) for an absent clade,
  # so it cannot tell a genuine conflict from the trivial, root-adjacent split that
  # an arbitrary rooting introduces. For each internal node we take its descendant
  # tip set; the split is "recovered" if that set (or its complement) is a clade in
  # Chazot. Trivial splits (one side is a single tip or the whole tree) carry no
  # topological information: they are left NA — uncoloured and excluded from pct —
  # so neither the rooting choice nor the root taxon is miscounted as a conflict.
  all_tips <- sort(mt$tip.label)
  pp       <- prop.part(rp); labs <- attr(pp, "labels")
  ref_keys <- new.env(parent = emptyenv())
  for (ci in pp) {
    s <- sort(labs[ci]); comp <- setdiff(all_tips, s)
    if (length(s))    assign(paste(s,    collapse = "\r"), TRUE, ref_keys)
    if (length(comp)) assign(paste(comp, collapse = "\r"), TRUE, ref_keys)
  }
  nodes  <- (length(mt$tip.label) + 1):(length(mt$tip.label) + mt$Nnode)
  status <- vapply(nodes, function(nd) {
    s <- sort(extract.clade(mt, nd)$tip.label); comp <- setdiff(all_tips, s)
    if (length(s) <= 1 || length(comp) <= 1) return(NA_character_)      # trivial split
    if (exists(paste(s, collapse = "\r"), envir = ref_keys, inherits = FALSE))
      "recovered in Chazot" else "conflicting"
  }, character(1))
  node_df <- tibble(
    node      = nodes,
    congruent = status,
    support   = suppressWarnings(as.numeric(mt$node.label))
  )
  pct <- round(100 * mean(status == "recovered in Chazot", na.rm = TRUE))

  p <- ggtree(mt) %<+% node_df +
    geom_tiplab(size = 2.4, offset = 0.002) +
    geom_nodepoint(aes(colour = congruent, size = support)) +
    scale_colour_manual(values = c("recovered in Chazot" = COL_OK, "conflicting" = COL_BAD),
                        name = NULL, na.translate = FALSE) +
    scale_size_continuous(range = c(1, 4), name = "bootstrap", limits = c(0, 100)) +
    ggtitle(sprintf("%s (n = %d shared species)", title, length(shared))) +
    theme_tree2() + theme(legend.position = "right") +
    xlim(0, max(node.depth.edgelength(mt)) * 1.35)
  if (reflect) p <- p + scale_y_reverse()

  list(plot = p, pct_recovered = pct, n_shared = length(shared),
       rf = RF.dist(mt, rp), rf_norm = RF.dist(mt, rp, normalize = TRUE),
       node_df = node_df, our_pruned = mt, ref_pruned = rp)
}
