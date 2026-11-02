# Tidy CKMR notebooks

Start with [fit_delfi_A.Rmd](fit_delfi_A.Rmd), a static R Markdown companion to
Mark's Acme Archipelago example. Preparation is shown inline, with
inspectable tables, figures, and mathematics. The only analysis function is the
likelihood used repeatedly by the optimizer.

Open the Rmd in RStudio and run its chunks in order or click **Knit**. Required
packages are `dplyr` (1.1.0 or newer), `tidyr`, `ggplot2`, `tibble`, `rmarkdown`,
and `knitr`. Rendering HTML also requires Pandoc, included with RStudio.

Change `data_path` in the first data-loading chunk, or edit the `data_path` YAML
parameter. The default is `../.can2024/datasets/samp_delfi_A.rda`, relative to
the Rmd's folder. In this checkout the repository link is actually named
`.can24`; the notebook explicitly falls back to that spelling when the unchanged
default is unavailable. Absolute paths work too. The simulated data are not
copied into this folder; participants need their own authorized copy of Mark's
course repository.

To render from the CKMRshiny project root:

```r
rmarkdown::render(
  "tidy_rmd_examples/fit_delfi_A.Rmd",
  output_dir = "tidy_rmd_examples/output"
)
```

Here `output_dir` is relative to the console's working directory (the project
root), while relative data paths inside chunks resolve from the Rmd's folder.
To override the data location explicitly:

```r
rmarkdown::render(
  "tidy_rmd_examples/fit_delfi_A.Rmd",
  params = list(data_path = "/absolute/path/to/can2024/datasets/samp_delfi_A.rda",
                check_mark_code = FALSE),
  output_dir = "tidy_rmd_examples/output"
)
```

The final optional section compares preparation and likelihood values with
Mark's original functions, and checks the fitted parameter estimates. It uses
`mvbutils`, `offarray`, and `atease` only when available and when
`check_mark_code` is true. It sources `r-scripts/ckmr_funs.r` from the repository
containing the chosen dataset; it does not run Mark's interactive driver script.

Generated HTML files are kept in `output/`, which is ignored by Git.
