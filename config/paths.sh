# =============================================================================
# config/paths.sh -- every path, environment and file name used by the pipeline
#
# Sourced by run_pipeline.sh, by every job script (via bin/job_env.sh) and by
# R (via R/setup.R, which imports these variables into Sys.getenv()).
# Each value can be overridden by exporting the variable before running, or
# by editing the default here.
# =============================================================================

# ---- cluster / software -----------------------------------------------------
: "${LAB_DIR:=/rsrch5/home/epi/bhattacharya_lab}"
: "${SCRATCH_DIR:=/rsrch5/home/epi/stbresnahan/scratch/Placenta_POE}"
: "${CONDA_EXE_PATH:=/risapps/rhel8/miniforge3/24.5.0-0/bin/conda}"
: "${BCFTOOLS_ENV:=bcftools-1.16}"                 # conda env with bcftools, tabix, bgzip
: "${SAMTOOLS_ENV:=samtools-1.16.1}"               # conda env with samtools
: "${LIFTOVER_BCFTOOLS:=${LAB_DIR}/software/bcftools_1.16_build/bcftools/bcftools}"  # build with the +liftover plugin
: "${PLINK_MODULE:=plink/1.90-beta}"
: "${SHAPEIT5_VENV:=${LAB_DIR}/software/shapeit5}"
: "${ORNAMENTS_VENV:=${LAB_DIR}/software/ornaments}"
: "${PYTHON_ENV:=${BCFTOOLS_ENV}}"                 # conda env with python3 + pysam (stage 05)
: "${RSTUDIO_SIF:=/risapps/singularity/repo/RStudio/4.3.1/rstudio_4.3.1.sif}"
: "${SING_BIND:=/rsrch5}"
: "${R_LIB_DIR:=/home/stbresnahan/R/ubuntu/4.3.1}"

# ---- run bookkeeping ---------------------------------------------------------
: "${RUN_DIR:=${SCRATCH_DIR}/run}"                 # receipts/ and plans/ live here
: "${LOG_DIR:=${SCRATCH_DIR}/logs}"
: "${FIG_DIR:=${SCRATCH_DIR}/figures}"             # all figures written by R drivers
: "${WORK_DIR:=${SCRATCH_DIR}/work}"               # tables written by R drivers

# ---- references ---------------------------------------------------------------
: "${REF_DIR:=${LAB_DIR}/data/GenomicReferences}"
: "${HG19_FASTA:=${REF_DIR}/genome/hg19_chr.fasta}"                       # chr-prefixed hg19
: "${HG38_FASTA:=${REF_DIR}/genome/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna}"
: "${HG19_TO_HG38_CHAIN:=${REF_DIR}/liftover/hg19ToHg38.over.chain}"
: "${GENETIC_MAP_DIR:=${REF_DIR}/genome/shapeit4_maps}"                  # chr<N>.b38.gmap.gz
: "${DBSNP_HG19_VCF:=${LAB_DIR}/data/dbSNP/hg19/common_all_20180423.vcf.gz}"
: "${DBSNP_HG38_VCF:=${LAB_DIR}/data/dbSNP/common_all_20180418.vcf.gz}"

# ---- long-read placenta transcriptome (ESPRESSO + SQANTI3) -----------------
: "${ASSEMBLY_DIR:=${LAB_DIR}/data/Placenta_scRNAseq/salmon}"
: "${TX_FASTA:=${ASSEMBLY_DIR}/ESPRESSO_corrected_SC_filtered.fasta}"
: "${TX_GTF:=${ASSEMBLY_DIR}/ESPRESSO_corrected_SC_filtered.gtf}"
: "${TX_GTF_EXONS:=${ASSEMBLY_DIR}/ESPRESSO_corrected_SC_filtered_withExons.gtf}"
: "${TX_CLASSIF:=${ASSEMBLY_DIR}/ESPRESSO_corrected_SC_filtered_classification.txt}"

# ---- GUSTO inputs -------------------------------------------------------------
: "${GUSTO_DIR:=${LAB_DIR}/data/mapqtl/GUSTO}"
: "${TRIOS_DIR:=${GUSTO_DIR}/genotypes/trios}"
: "${PLINK_PREFIXES:=infant/Infant1073 mother/Mother1079 father/Father732}"   # relative to TRIOS_DIR
: "${ARRAY_MANIFEST:=${TRIOS_DIR}/SNP_Table_8v1_2.txt}"
: "${INFANT_META_CSV:=${TRIOS_DIR}/n200_infant_metadata.csv}"        # SubjectID, maternal ethnicity
: "${FASTQ_DIR:=${GUSTO_DIR}/rnaseq/fastq}"
: "${FASTQ_ID_MAP:=${TRIOS_DIR}/infant/infant_samples_n200_fastqIDs.txt}"  # <SubjectID> <fastq id>
: "${RNASEQ_COVARS_CSV:=${LAB_DIR}/data/Restore08072025/Placenta_LRRNAseq/STB/infrv_GUSTO/20220823-Full_200_RNAseq_covars_v3.csv}"
: "${PFAS_MAT_CSV:=${SCRATCH_DIR}/inputs/GUSTO_PFAS_maternalBlood.csv}"
: "${PFAS_CHILD_CSV:=${SCRATCH_DIR}/inputs/GUSTO_PFAS_childBlood.csv}"
: "${PFAS_MAT_COHORT_XLSX:=${GUSTO_DIR}/covariates/mat_pfas_imputed.xlsx}"
: "${PFAS_CORD_COHORT_XLSX:=${GUSTO_DIR}/covariates/cord_pfas_imputed.xlsx}"

# ---- genotype intermediates (stages 01-03) ---------------------------------
: "${TRIO_PED:=${TRIOS_DIR}/GUSTO_trios.ped}"                         # FID IID PID MID SEX PHENO
: "${TRIO_FAM_SHAPEIT5:=${TRIOS_DIR}/GUSTO_trios_shapeit5.fam}"       # kid father mother
: "${TRIO_VCF_HG19:=${TRIOS_DIR}/merged/trios_hg19.vcf.gz}"
: "${TRIO_VCF_HG38:=${TRIOS_DIR}/merged/trios_hg38.vcf.gz}"           # biallelic, used downstream
: "${ANCESTRY_DIR:=${TRIOS_DIR}/sample_lists_by_ancestry}"            # samples_<group>.txt
: "${BY_ANCESTRY_DIR:=${TRIOS_DIR}/trios_by_ancestry}"
: "${SEAD_UPLOAD_DIR:=${TRIOS_DIR}/infant/phased_by_ancestry}"        # <group>/<group>_chr<N>.vcf.gz
: "${SEAD_DOWNLOAD_DIR:=${TRIOS_DIR}/infant/imputed_by_ancestry}"     # <group>_SEAD/chr<N>.SEAD.GRCh38.dose.vcf.gz
: "${IMPUTED_DIR:=${TRIOS_DIR}/infant/imputed_SEAD}"                  # chr<N>.biallelic.filtered.vcf.gz
: "${ANCESTRY_GROUPS:=chinese indian malay}"
: "${AUTOSOMES:=1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22}"

# ---- allele-specific expression (stages 04-06) -----------------------------
: "${PROJECT_DIR:=${LAB_DIR}/data/Placenta_POE}"
: "${ORN_DIR:=${PROJECT_DIR}/ornaments}"
: "${TXVCF_DIR:=${ORN_DIR}/tx_vcf_SEAD}"
: "${INFANT_SAMPLES:=${TRIOS_DIR}/infant/infant_samples_rnaseq.txt}" # VCF IDs of RNA-seq infants
: "${TXVCF_MERGED:=${TXVCF_DIR}/txomic.merged.infants.vcf.gz}"              # bgzipped + indexed
: "${HAP_ORIGIN:=${TXVCF_DIR}/hap_origin_matrix.tsv}"
: "${QUANT_DIR:=${ORN_DIR}/quant}"                                   # <sample>/ per infant
: "${TXOME_DIR:=${ORN_DIR}/txomes}"                                  # personalised transcriptomes + indexes
: "${KEEP_TXOME_INDEX:=0}"                                           # 1 = keep each ornaments index after quant
: "${POE_GENE_COUNTS:=${ORN_DIR}/POE_gene_SNP_counts.tsv}"
: "${POE_TX_COUNTS:=${ORN_DIR}/POE_tx_SNP_counts.tsv}"

# ---- models (stages 07-09) ---------------------------------------------------
: "${MODEL_DIR:=${ORN_DIR}/models}"
: "${MODEL_INPUTS:=${MODEL_DIR}/model_inputs.rds}"                    # matrices + metadata
: "${BATCH_DIR:=${MODEL_DIR}/batches}"                                # per-batch inputs and results
: "${N_BATCHES_GENE:=80}"
: "${N_BATCHES_TX:=150}"
: "${POE_CONDITIONS:=PFBA_mat PFOA_mat PFNA_mat PFDA_mat PFUnDA_mat PFBS_mat PFHxS_mat PFOS_mat PFBA_child PFOA_child PFNA_child PFDA_child PFUnDA_child PFBS_child PFHxS_child PFOS_child GDM}"
: "${N_BOOT:=100}"
: "${IMPRINT_DB:=${WORK_DIR}/geneimprint_db_ensembl.txt}"
: "${PLACENTA_DB:=${WORK_DIR}/placenta_db_ensembl.txt}"
: "${EXCLUDE_SAMPLES_FILE:=}"                                         # optional, one sample ID per line; removed from all models
: "${DOSE_RESPONSE_CONDITIONS:=PFBS_mat}"                             # stage 08 dose-response plots
: "${DOSE_RESPONSE_EXCLUDE:=HSPA9 DIS3L2}"                            # gene symbols left out of gene-level dose-response plots

# ---- simulation (stage 10) ---------------------------------------------------
: "${SIM_DIR:=${PROJECT_DIR}/PFAS_simulation}"
: "${SIM_N_RUNS:=100}"
