# Figure 2 — Na⁺/K⁺-ATPase ouabain resistance (IC₅₀)

Code and data for Figure 2. A single notebook in `code/`, run with the working
directory at `code/` (all paths are relative to that folder):

- **`1_butterfly_vs_fly_IC50s.Rmd`** — ouabain IC₅₀ dose–response analysis comparing
  Na⁺/K⁺-ATPase resistance of milkweed butterflies (Petschenka et al. 2013) with
  engineered *Drosophila melanogaster* ATPα knock-in lines (Karageorgi et al. 2019).
  Fits one four-parameter log-logistic curve per biological preparation, then
  compares IC₅₀ values by ATPα state. Produces the supplementary IC₅₀ table and the
  IC₅₀-by-state figure.

## Directory layout

```
Figure2/
├── README.md
├── code/
│   └── 1_butterfly_vs_fly_IC50s.Rmd
├── input_data/
│   ├── Petschenka2013_NaKATPase_raw_absorbance.csv   # raw colorimetric assay data (Dryad)
│   └── Karageorgi2019_Fig3a_NaKATPase_activity.csv   # Fig. 3a source data (engineered Drosophila)
├── output_data/
│   └── Supplementary_Table_ATPase_IC50_replicates.xlsx  # per-preparation IC50 estimates
└── visuals/
    └── IC50_by_ATPase_state_means.png                # IC50 by ATPα state
```

## Running

```r
setwd("<path>/Figure2/code")
rmarkdown::render("1_butterfly_vs_fly_IC50s.Rmd")
```

## R package dependencies

`tidyverse`, `drc`, `openxlsx`, `purrr`, `tidyr`
