# Figure 1 — Host use and ATPα ancestral states across Danaini

Code and data for Figure 1: ancestral reconstruction of host use (none / milkweed
clade / native-range *Asclepias*) along the Danaini butterfly phylogeny, together
with the maximum-likelihood ancestral reconstruction of the ATPα functional
genotype across the same tree.

Two notebooks in `code/`, run with the working directory at `code/` (all paths are
relative to that folder), built on prepared inputs carried over from Extended Data
Figure 1 (so Figure 1, Extended Data Fig. 1, and Extended Data Fig. 2 all rest on
the same curated host-use data and the same butterfly tree):

- **`code/1_host_use_ancestral_state_reconstruction.Rmd`** — a three-state host-use
  ancestral state reconstruction (none / milkweed clade / native *Asclepias*),
  fitted with a symmetric (`SYM`) Mk model selected by AIC. The full-tree figure
  carries the three-state node-pie reconstruction, two host-use tip bars, and an
  inner **ATPα-state tip strip** (the curated per-species genotype, with unsequenced
  taxa shown as "Unknown"), so host use and the ATPα state can be read side by side.

- **`code/2_atpalpha_ASR.Rmd`** — ancestral state reconstruction of the ATPα
  functional states across the Danaini, plus the sensitivity analyses. Each taxon's
  state is the three-residue code at ATPα positions 1 / 9 / 12 (VSH, VSN, LSN, QAN,
  LAN). Sections:
  1. **Whole-tree reconstruction** — the reconstruction on the full Danaini tree,
     unsampled taxa entered as missing data (`phangorn`, LG); produces a per-node
     state table (code + support + descendant coverage) and a whole-tree figure in
     which low-coverage / low-support nodes are visually downweighted.
  2. **Sequenced-subset reconstruction (IQ-TREE)** — LG+Γ, marginal (empirical-Bayes),
     on the tree pruned to taxa with sequence data, rooted on the outgroups
     *Heliconius erato* and *Mechanitis polymnia*; cross-checked against `phangorn`
     and against the Section 1 whole-tree reconstruction.
  3. **Bootstrap** — the subset reconstruction repeated across the 1000
     ultrafast-bootstrap trees of the Danaini backbone (phylogenetic uncertainty).
  4. **Jackknife** — each sequenced taxon dropped in turn (taxon sampling).

The ATPase ouabain-resistance (IC₅₀) analysis lives in `../Figure2/`; the
phylogenetic-signal (Fritz & Purvis D-statistic) analysis lives in
`../Extended_Data_Figure2/`.

## Directory layout

```
Figure1/
├── README.md
├── code/
│   ├── 1_host_use_ancestral_state_reconstruction.Rmd
│   └── 2_atpalpha_ASR.Rmd
├── input_data/
│   ├── Danaini_in_HOSTS.newick          # Danaini species tree (from Extended Data Fig. 1)
│   ├── danaini_reference_order.newick   # curated tip order shared with Ext. Fig. 1 / Figure 2
│   ├── danaini_host_use_states.csv      # per-species host-use states + ATPα (from Ext. Fig. 1 notebook 5)
│   └── alignment.fasta                  # ATPα functional-motif alignment (incl. 4 congener/subspecies-sourced tips)
├── output_data/
│   ├── host_use_milkweed_asclepias_ASR_node_probabilities.csv   # (notebook 1) 3-state ASR node probabilities
│   ├── atpalpha_ASR_node_posteriors.csv              # (notebook 2) IQ-TREE per-node posteriors
│   ├── atpalpha_ASR_bootstrap_robustness.csv         # (notebook 2) per-node state consistency across boot trees
│   ├── atpalpha_ASR_jackknife_stability.csv          # (notebook 2) per-node leave-one-out stability
│   ├── atpalpha_ASR_coverage.csv                     # (notebook 2, §1) whole-tree node-state table (code + support + coverage)
│   ├── atpalpha_ASR_full_tree_posteriors.csv         # (notebook 2) full-tree per-node posteriors
│   └── iqtree_asr/                                   # (notebook 2) IQ-TREE run files (.state, .treefile, …)
└── visuals/
    ├── host_use_milkweed_asclepias_ASR.png          # (notebook 1) full-tree 3-state reconstruction + ATPα strip
    ├── atpalpha_ASR_whole_tree.png                  # (notebook 2, §1) whole-tree reconstruction (size=coverage, opacity=support)
    ├── atpalpha_ASR_tree.png                        # (notebook 2, §2) IQ-TREE subset ancestral states
    ├── atpalpha_ASR_bootstrap.png                   # (notebook 2, §3) bootstrap state consistency
    └── atpalpha_ASR_jackknife.png                   # (notebook 2, §4) jackknife state stability
```

## Inputs (prepared in Extended Data Figure 1)

The three host-use input files are copied from `../Extended_Data_Figure1/`:

- **`Danaini_in_HOSTS.newick`** — the filtered Danaini ML tree, including the
  outgroups *Mechanitis polymnia* and *Heliconius erato*.
- **`danaini_host_use_states.csv`** — one host-use state per species (the
  native-range *Asclepias* flag, the milkweed-clade flag, the collapsed state
  coding, and the ATPα genotype).
- **`danaini_reference_order.newick`** — the curated tip order (shared with
  Extended Data Fig. 1 and Figure 2); notebook 1 rotates the tree to this order so
  the figures match the orientation used elsewhere.
- **`alignment.fasta`** — the ATPα functional-motif alignment used by notebook 2.

A tree tip with no row in the host-use table (here *Anetia thirza*) has no
documented host record and is coded as no milkweed / *Asclepias* use.

## Model note (notebook 1)

The three-state host-use reconstruction uses a **symmetric (`SYM`)** Mk model, chosen
by AIC over equal-rates (`ER`) and all-rates-different (`ARD`). `ER` is a markedly
worse fit and, because its equal rates imply a uniform equilibrium, inflates the rare
`asclepias` state at deep nodes (root `asclepias` ≈ 0.29 under `ER` vs ≈ 0.05 under
`SYM`); the notebook prints the AIC comparison.

## External dependencies (notebook 2)

**Notebook 2 requires IQ-TREE and the Danaini bootstrap trees.** IQ-TREE is called
via `system2()`; set the `IQTREE_BIN` environment variable to your binary, otherwise
it defaults to `/opt/anaconda3/envs/danaini-phylo/bin/iqtree` and falls back to
`iqtree`/`iqtree2`/`iqtree3` on `PATH`. The bootstrap section reads
`../../Extended_Data_Figure1/output_data/danaini_tree/supermatrix_ge1/supermatrix_ge1.ufboot`.

## Running

```r
setwd("<path>/Figure1/code")

# Notebook 1 — host-use reconstruction + ATPα-strip figure:
rmarkdown::render("1_host_use_ancestral_state_reconstruction.Rmd")

# Notebook 2 — ATPα ancestral states + sensitivity (needs IQ-TREE on PATH/IQTREE_BIN):
rmarkdown::render("2_atpalpha_ASR.Rmd")
```

## R package dependencies

- **Notebook 1:** `tidyverse`, `ape`, `phytools`, `ggtree`, `ggnewscale`
- **Notebook 2:** `ape`, `phangorn`, `seqinr`, `tidyverse` (plus the IQ-TREE binary)
