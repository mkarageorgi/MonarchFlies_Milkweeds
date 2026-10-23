# Figure 1 — Host-use ancestral state reconstruction across Danaini

Code and data for Figure 1: reconstruction of ancestral host use (none / milkweed
clade / native-range *Asclepias*) along the Danaini butterfly phylogeny.

A single notebook, built on prepared inputs carried over from Extended Data Figure 1
(so Figure 1, Extended Data Fig. 1, and Extended Data Fig. 2 all rest on the same
curated host-use data and the same butterfly tree):

- **`code/1_host_use_ancestral_state_reconstruction.Rmd`** — a three-state host-use
  ancestral state reconstruction (none / milkweed clade / native *Asclepias*),
  fitted with a symmetric (`SYM`) Mk model selected by AIC, with the full-tree
  figure.

The ATPα-focused figure reuses this reconstruction and lives in `../Figure2/`. The
related phylogenetic-signal (Fritz & Purvis D-statistic) analysis lives in
`../Extended_Data_Figure2/`.

The notebook runs with the working directory at `code/`; all paths are relative to
that folder.

## Directory layout

```
Figure1/
├── README.md
├── code/
│   └── 1_host_use_ancestral_state_reconstruction.Rmd
├── input_data/
│   ├── Danaini_in_HOSTS.newick          # Danaini species tree (from Extended Data Fig. 1)
│   ├── danaini_reference_order.newick   # curated tip order shared with Ext. Fig. 1 / Figure 2
│   └── danaini_host_use_states.csv      # per-species host-use states + ATPα (from Ext. Fig. 1 notebook 5)
├── output_data/
│   └── host_use_milkweed_asclepias_ASR_node_probabilities.csv   # 3-state ASR node probabilities
└── visuals/
    └── host_use_milkweed_asclepias_ASR.png          # full-tree 3-state reconstruction
```

## Inputs (prepared in Extended Data Figure 1)

All three inputs are copied from `../Extended_Data_Figure1/`:

- **`Danaini_in_HOSTS.newick`** — the filtered Danaini ML tree (56 tips, including
  the outgroups *Mechanitis polymnia* and *Heliconius erato*).
- **`danaini_host_use_states.csv`** — one host-use state per species (the
  native-range *Asclepias* flag, the milkweed-clade flag, the collapsed state
  coding, and the ATPα genotype).
- **`danaini_reference_order.newick`** — the curated tip order (shared with
  Extended Data Fig. 1 and Figure 2); the notebook rotates the tree to this order so
  the figures match the orientation used elsewhere.

A tree tip with no row in the host-use table (here *Anetia thirza*) has no
documented host record and is coded as no milkweed / *Asclepias* use.

## Model note

The three-state reconstruction uses a **symmetric (`SYM`)** Mk model, chosen by AIC
over equal-rates (`ER`) and all-rates-different (`ARD`). `ER` is a markedly worse
fit and, because its equal rates imply a uniform equilibrium, inflates the rare
`asclepias` state at deep nodes (root `asclepias` ≈ 0.29 under `ER` vs ≈ 0.05 under
`SYM`); the notebook prints the AIC comparison.

## Running

```r
setwd("<path>/Figure1/code")
rmarkdown::render("1_host_use_ancestral_state_reconstruction.Rmd")
```

R package dependencies: `tidyverse`, `ape`, `phytools`, `ggtree`, `ggnewscale`.
