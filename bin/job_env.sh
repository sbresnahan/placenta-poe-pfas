# =============================================================================
# bin/job_env.sh -- common environment for every job script
#
# Usage (first lines of a job script):
#   source "${REPO:?REPO not set; submit through run_pipeline.sh}/bin/job_env.sh"
#
# Provides
#   NCPU              cores granted by LSF (LSB_DJOB_NUMPROC; 1 outside LSF)
#   IDX               array index (LSB_JOBINDEX; 1 outside an array)
#   use_bcftools      activate the conda env with bcftools/tabix/bgzip
#   use_conda ENV     activate any conda env
#   rscript ARGS...   Rscript inside the RStudio Singularity image
#   nth_word N LIST   N-th (1-based) word of a space-separated list
#   receipt STEP [I]  mark STEP (array element I) complete
#   run_step_r FILE   rscript "$REPO/FILE" with the repo on the R path
# =============================================================================
set -euo pipefail

set -a
# shellcheck source=../config/paths.sh
source "${REPO}/config/paths.sh"
set +a

export NCPU="${LSB_DJOB_NUMPROC:-${NCPU:-1}}"
export IDX="${LSB_JOBINDEX:-${IDX:-1}}"
[[ "${IDX}" == "0" ]] && IDX=1
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export REPO

mkdir -p "${RUN_DIR}/receipts" "${LOG_DIR}" "${FIG_DIR}" "${WORK_DIR}"

log() { printf '[%s] %s\n' "$(date '+%F %T')" "$*" >&2; }

nth_word() { local n=$1; shift; echo "$*" | tr -s ' \t\n' '\n' | sed -n "${n}p"; }
n_words()  { echo "$*" | wc -w; }

_conda_hook_loaded=0
use_conda() {
  if [[ ${_conda_hook_loaded} -eq 0 ]]; then
    [[ -f /etc/profile.d/modules.sh ]] && source /etc/profile.d/modules.sh
    set +u; eval "$("${CONDA_EXE_PATH}" shell.bash hook)"; set -u
    _conda_hook_loaded=1
  fi
  set +u; conda activate "$1"; set -u
}
use_bcftools() { use_conda "${BCFTOOLS_ENV}"; }
use_venv()     { set +u; source "$1/bin/activate"; set -u; }
use_module()   { [[ -f /etc/profile.d/modules.sh ]] && source /etc/profile.d/modules.sh; module load "$@"; }

rscript() {
  singularity exec --bind "${SING_BIND}" "${RSTUDIO_SIF}" Rscript "$@"
}
run_step_r() { local f=$1; shift; rscript "${REPO}/${f}" "$@"; }

receipt() {
  local step=$1 i=${2:-1}
  date '+%F %T' > "${RUN_DIR}/receipts/${step}.${i}.ok"
  log "done: ${step} [${i}]"
}

# group/chromosome helpers for (group x autosome) arrays
group_of() { nth_word $(( ($1 - 1) / $(n_words ${AUTOSOMES}) + 1 )) ${ANCESTRY_GROUPS}; }
chr_of()   { nth_word $(( ($1 - 1) % $(n_words ${AUTOSOMES}) + 1 )) ${AUTOSOMES}; }

log "job ${LSB_JOBID:-local}[${IDX}] on $(hostname) with ${NCPU} cores"
