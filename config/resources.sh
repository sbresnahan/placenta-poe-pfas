# =============================================================================
# config/resources.sh -- LSF resources per step
#
# RES_<step>   : bsub flags (wall time, cores, memory in GB, queue)
# ARR_<step>   : array concurrency limit (the %K in "name[1-N]%K"); 0 = no limit
# Job scripts read the number of cores LSF actually granted from
# LSB_DJOB_NUMPROC, so changing "-n" here is the only edit needed to scale a step.
# Tune the Ornaments values after a pilot: `run_pipeline.sh 04_ornaments --pilot 2`
# then `tools/lsf_usage_report.sh 04_ornaments`.
# =============================================================================

: "${LSF_QUEUE_SHORT:=short}"     # <= 3 h
: "${LSF_QUEUE_MEDIUM:=medium}"   # <= 24 h
: "${LSF_QUEUE_LONG:=long}"

res() { echo "-W $1 -n $2 -M $3 -R rusage[mem=$3] -R span[hosts=1] -q $4"; }

# ---- 00 reference ----
RES_00_reference_index=$(res 1:00 2 16 "$LSF_QUEUE_SHORT")
RES_00_imprint_db=$(res 1:00 1 8 "$LSF_QUEUE_SHORT")

# ---- 01 genotype preparation ----
RES_01_pedigree=$(res 0:30 1 4 "$LSF_QUEUE_SHORT")
RES_01_plink2vcf=$(res 2:00 4 16 "$LSF_QUEUE_SHORT")
RES_01_liftover=$(res 3:00 4 16 "$LSF_QUEUE_SHORT")
RES_01_check_mendel=$(res 1:00 4 16 "$LSF_QUEUE_SHORT")
RES_01_manifest_qc=$(res 3:00 12 32 "$LSF_QUEUE_SHORT")

# ---- 02 ancestry + phasing ----
RES_02_infer_ancestry=$(res 2:00 12 32 "$LSF_QUEUE_SHORT")
RES_02_split=$(res 2:00 8 16 "$LSF_QUEUE_SHORT")
RES_02_phase=$(res 3:00 8 16 "$LSF_QUEUE_SHORT")
RES_02_subset=$(res 2:00 8 16 "$LSF_QUEUE_SHORT")
ARR_02_phase=0

# ---- 03 SEAD imputation (filter + merge) ----
RES_03_filter=$(res 1:00 4 8 "$LSF_QUEUE_SHORT")
RES_03_merge=$(res 2:00 4 16 "$LSF_QUEUE_SHORT")
RES_03_qc=$(res 2:00 8 32 "$LSF_QUEUE_SHORT")
ARR_03_filter=0

# ---- 04 transcriptomic VCF + Ornaments ----
RES_04_prep_tx=$(res 1:00 1 8 "$LSF_QUEUE_SHORT")
RES_04_g2t=$(res 3:00 1 16 "$LSF_QUEUE_SHORT")
RES_04_g2t_merge=$(res 3:00 2 32 "$LSF_QUEUE_SHORT")
RES_04_ornaments=$(res 12:00 4 16 "$LSF_QUEUE_MEDIUM")
ARR_04_ornaments=50

# ---- 05 parent-of-origin assignment ----
RES_05_assign_poe=$(res 2:00 22 32 "$LSF_QUEUE_SHORT")

# ---- 06 feature matrices ----
RES_06_summarize=$(res 3:00 16 64 "$LSF_QUEUE_SHORT")

# ---- 07 baseline POE / 08 PFAS interaction ----
RES_07_inputs=$(res 1:00 1 32 "$LSF_QUEUE_SHORT")
RES_07_within_sample=$(res 2:00 16 48 "$LSF_QUEUE_SHORT")
RES_07_fit=$(res 1:00 8 16 "$LSF_QUEUE_SHORT")
RES_07_combine=$(res 0:30 1 16 "$LSF_QUEUE_SHORT")
RES_07_report=$(res 1:00 1 32 "$LSF_QUEUE_SHORT")
RES_08_fit=$(res 3:00 8 16 "$LSF_QUEUE_SHORT")
RES_08_combine=$(res 1:00 1 32 "$LSF_QUEUE_SHORT")
RES_08_report=$(res 2:00 1 48 "$LSF_QUEUE_SHORT")
ARR_07_fit=0
ARR_08_fit=0

# ---- 09 power ----
RES_09_power=$(res 3:00 16 32 "$LSF_QUEUE_SHORT")

# ---- 10 mediation simulation ----
RES_10_sim=$(res 12:00 12 32 "$LSF_QUEUE_MEDIUM")
RES_10_aggregate=$(res 0:30 1 8 "$LSF_QUEUE_SHORT")
ARR_10_sim=0
