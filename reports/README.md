# reports

| File | Content |
|---|---|
| `genotype_methods.md` | genotyping, sample and SNP QC of the GUSTO trio arrays |
| `Placenta_POE_PFAS_methods.Rmd` | methods write-up for the POE and PFAS interaction analyses |
| `all_PFAS_report.Rmd` | per-PFAS interaction results report |
| `PFAS_POE_SRS_NCsim_annotated.Rmd` | simulation study of the planned PFAS -> imprinted isoform -> SRS mediation analysis |
| `figures/MANIFEST.tsv` | every figure file the two figure-based reports need, with its source |
| `stage_figures.sh` | copies those figures into `figures/methods/` and `figures/interaction/` |
| `render_reports.R` | re-knits the two figure-based reports into `docs/reports/` |

The rendered HTML of each Rmd is in `docs/reports/` and is linked from the
project site.

## Re-knitting the figure-based reports

Both Rmds only arrange existing PNGs (cowplot + magick); they run no analysis.
Each takes a `fig_dir` parameter, by default a folder under `reports/figures/`:

| Report | `fig_dir` default | Files | Produced by |
|---|---|---|---|
| `all_PFAS_report.Rmd` | `figures/interaction` | 128: `{gene,tx}_diffPOE_<cond>_BetaInt{,_pvals,_pQQ}.png`, `{gene,tx}_forest_<cond>.png` for the 16 `PFAS_{mat,child}` conditions | stage 08 `interaction_report.R`, under the same names in `FIG_DIR/interaction/` |
| `Placenta_POE_PFAS_methods.Rmd` | `figures/methods` | 22: `figure2` ... `figure16`, plus `figure7b, 8b, 10b, 12b, 13b, 14b, 15b` | pre-assembled panels kept in the local project folder `Reorganize/`; see the table below |

`figures/MANIFEST.tsv` lists, for every file, its path in the local project
folder (`local_source`) and the pipeline figure it corresponds to
(`pipeline_source`, relative to `FIG_DIR`).

### 1. Stage the figures

From the local project folder (exact panels the current HTML was built from):

```bash
reports/stage_figures.sh --local /path/to/Placenta_POE
```

or from pipeline output (`FIG_DIR` from `config/paths.sh`, or given explicitly):

```bash
reports/stage_figures.sh --pipeline                 # uses $FIG_DIR
reports/stage_figures.sh --pipeline /path/to/figures
reports/stage_figures.sh --local /path/to/Placenta_POE   # then fill figure5.png
```

`--check` lists anything still missing; existing files are kept unless
`--force` is given.

### 2. Render

```bash
Rscript reports/render_reports.R            # both, HTML, into docs/reports/
Rscript reports/render_reports.R all_PFAS   # or: methods
```

`all_PFAS_report.Rmd` can also read stage 08 output in place, without staging:

```bash
REPORT_INTERACTION_DIR=$FIG_DIR/interaction Rscript reports/render_reports.R all_PFAS
```

In RStudio, use *Knit with Parameters* and set `fig_dir`. Note that the
methods Rmd lists `pdf_document` first, so a plain *Knit* builds the PDF (needs
LaTeX); `render_reports.R` builds the HTML used on the site.

R packages: `rmarkdown`, `knitr`, `cowplot`, `magick`, `ggplot2`.

### Methods panels and their pipeline equivalents

| Panel | Content | Pipeline figure (`FIG_DIR/...`) |
|---|---|---|
| figure2, figure3 | genotype PCA: reported ethnicity, KNN-inferred ancestry | `PCA_infants_reported_ethnicity.png`, `PCA_infants_inferred_ancestry.png` (restyled) |
| figure4 | heterozygous variants per transcript | `txomic_variants_per_tx.png` (no summary annotations) |
| figure5 | heterozygous variants per consensus transcript | **none**: local copy only |
| figure6, figure7 | within-sample biased genes per sample: counts, proportions | `significant_POE_totals.png`, `significant_POE_props.png` (restyled) |
| figure7b | majority-allele ratio vs abs(beta) | `beta_by_majority_allele.png` (restyled) |
| figure8, 8b, 9 | gene-level baseline POE: p-values, QQ, volcano | `genePOE_Bpdist.png`, `genePOE_pQQ.png`, `genePOE.png` |
| figure10, 10b, 11 | transcript-level baseline POE | `isoPOE_Bpdist.png`, `isoPOE_pQQ.png`, `isoPOE.png` |
| figure12, 12b, 13, 13b | gene x maternal PFBS: p-values, QQ, volcano, forest | `interaction/gene_diffPOE_PFBS_mat_BetaInt{_pvals,_pQQ,}.png`, `interaction/gene_forest_PFBS_mat.png` |
| figure14, 14b, 15, 15b | transcript x maternal PFBS | `interaction/tx_...PFBS_mat...` (same pattern) |
| figure16 | DAG | `PFAS_POE_DAG.png` (stage 10, `DAG.R`) |

"Restyled" means the pipeline draws the same quantity with different titles or
colours, so a report built from `--pipeline` will not look identical to the
current HTML.

## Simulation report

`PFAS_POE_SRS_NCsim_annotated.Rmd` runs its own simulation (long; it caches
`results_list.RData` in the working directory). Stage 10 is the cluster
version of the same simulation.
