# placenta-poe-pfas

Parent-of-origin effects (POE) on placental gene and isoform expression, and
their modification by prenatal PFAS exposure, in GUSTO mother-father-infant
trios.

The pipeline goes from trio SNP-array genotypes and placental RNA-seq to
cohort-level POE estimates per gene and transcript and to PFAS x POE
interaction effects:

1. trio genotypes -> GRCh38, ancestry assignment, pedigree-aware phasing (SHAPEIT5)
2. SEAD imputation of the phased infants
3. allele-specific quantification against a long-read placenta transcriptome (Ornaments)
4. maternal / paternal assignment of infant haplotypes from the parents' genotypes
5. negative-binomial mixed models (glmmTMB) for baseline POE and for allele x PFAS / GDM interactions
6. power analysis and a simulation of the planned mediation analysis

**Start with [`WORKFLOW.md`](WORKFLOW.md)**: stage diagram, inputs, how to run,
methods per stage, outputs and caveats. The project site (`docs/`, built from
`book/`) collects the workflow and the rendered analysis reports.

## Quick start

```bash
# 1. edit config/paths.sh (data, references, software environments)
# 2. look at the plan, then submit
./run_pipeline.sh all
./run_pipeline.sh --submit all
./run_pipeline.sh status
```

The pipeline runs on an LSF cluster. It pauses after phasing for the external
SEAD imputation step (`WORKFLOW.md`, section 4.1) and resumes from completion
receipts when the same command is run again.

## Layout

| Path | Content |
|---|---|
| `run_pipeline.sh` | LSF submitter with dependencies, receipts, pilot mode |
| `config/` | paths, parameters and LSF resources |
| `00_reference/` ... `10_mediation_sim/` | one folder per stage |
| `R/` | shared R library (models, tiers, plots) |
| `bin/`, `tools/` | job helpers, resource report |
| `reports/` | Rmd sources and genotyping methods |
| `book/`, `docs/` | project site source and rendered site |

## Data

No data or participant identifiers are included. GUSTO data are available
under the cohort's data access procedures.

## Funding

NIH R21 (NIH RePORTER project 11437481).
