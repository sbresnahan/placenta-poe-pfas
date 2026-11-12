# Genotyping and genotype QC (GUSTO trios)

## DNA and arrays

Maternal DNA was extracted from blood collected at mid-gestation, infant DNA
from umbilical cord tissue or blood, and paternal DNA from buccal swabs. DNA
extraction for mothers and infants was as described previously [1-3]. Isohelix
DNA buccal swabs stored at -80 C were equilibrated to room temperature for 1 h;
1 mL of ATL lysis buffer (Qiagen) was added and incubated at room temperature
for 30 min, then 100 uL of proteinase K (Qiagen) was added and incubated at
60 C with shaking for 30-60 min. DNA was extracted from the lysates with
QIAsymphony DNA kits according to the manufacturer's instructions.

Paternal buccal-swab DNA (n = 732), maternal buffy-coat DNA collected at 26-28
weeks of gestation (n = 1079) and umbilical-cord DNA from the children
(n = 1073) were genotyped on the Illumina OmniExpressExome-8 v1.2 array
(~1 million SNPs) [3]. Hybridisation and scanning were performed by Expression
Analysis, Inc. (Morrisville, North Carolina). Data were processed with the
GenomeStudio Genotyping Module (Illumina): genotypes were called by GenCall
(GenTrain clustering, Bayesian calling), and genotypes with a GenCall score
< 0.15 were set to missing.

## Array-level QC (before this pipeline)

* Principal components analysis confirmed self-reported ethnicity.
* Samples with call rate < 97%, cryptic relatedness, or sex / ethnicity
  discrepancies were excluded.
* Alleles are reported on hg19.
* SNP QC was done separately in the Chinese, Malay and Indian subgroups: a SNP
  was dropped if call rate < 0.95, MAF < 0.05 or Hardy-Weinberg P < 1e-6.

The resulting PLINK files (one each for infants, mothers and fathers) and the
array manifest are the genotype inputs of stage 01.

## Processing in this pipeline

| Step | What it does |
|---|---|
| `01_genotype_prep/prep_ped.R` | trio pedigree from the three `.fam` files |
| `01_genotype_prep/plink2vcf.lsf` | merge, TOP strand -> hg19 forward strand (`bcftools +fixref`), REF check |
| `01_genotype_prep/liftover.lsf` | hg19 -> GRCh38, biallelic SNPs with GRCh38-matching REF |
| `01_genotype_prep/check_mendel.lsf` | Mendelian error rates |
| `02_ancestry_phasing/infer_ancestry.R` | PCA + 3-NN ancestry for infants without reported ethnicity |
| `02_ancestry_phasing/split_by_ancestry.lsf` | per-group filters: missingness <= 0.05, MAF >= 0.05, HWE p >= 1e-6 |
| `02_ancestry_phasing/phase.lsf` | SHAPEIT5 pedigree-aware phasing per group and chromosome |
| `03_imputation_SEAD/` | SEAD imputation of phased infants; biallelic, MAF > 0.01, Rsq >= 0.3 |

See `WORKFLOW.md` for details.

## References

1. Chen L, Pan H, Tuan TA, et al. Brain-derived neurotrophic factor (BDNF)
   Val66Met polymorphism influences the association of the methylome with
   maternal anxiety and neonatal brain volumes. *Dev Psychopathol* 2015;
   27(1):137-150.
2. Lin X, Teh AL, Chen L, et al. Choice of surrogate tissue influences neonatal
   EWAS findings. *BMC Med* 2017; 15(1):211.
3. Sampathkumar A, Tan KM, Chen L, et al. Genetic link determining the
   maternal-fetal circulation of vitamin D. (manuscript)
