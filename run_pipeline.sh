#!/usr/bin/env bash
# =============================================================================
# run_pipeline.sh -- submit the placenta POE x PFAS pipeline to LSF
#
#   ./run_pipeline.sh [options] TARGET [TARGET ...]
#   ./run_pipeline.sh status
#   ./run_pipeline.sh list
#
# TARGET   a step name (e.g. 02_phase), a stage number (e.g. 02),
#          `all` (every main-path step), `qc` (optional QC steps) or `sim`
#          (stage 10 simulation).
#
# Options
#   --submit      submit the jobs (default is a dry run that prints the bsub plan)
#   --force       ignore completion receipts and re-run every selected element
#   --pilot N     submit at most N elements of each selected array step;
#                 steps downstream of a piloted step are held back
#   --lenient     downstream steps start when upstream jobs *end* (ended())
#                 instead of only when every element succeeded (done())
#   --tag NAME    job-name prefix (default: poe). Use a different tag to run
#                 two copies side by side.
#
# How it works
#   * Every step is one LSF job or one LSF job array. Dependencies are LSF
#     job dependencies (-w "done(<tag>_<step>)"), so a single call queues the
#     whole chain and nothing polls.
#   * Each finished job/element writes $RUN_DIR/receipts/<step>.<index>.ok.
#     Steps whose receipts are all present are skipped; array steps resubmit
#     only the missing elements. Re-running the same command therefore resumes.
#   * Stage 03 needs the SEAD imputation results, which come from an external
#     server. Until all chr<N>.SEAD.GRCh38.dose.vcf.gz files exist, `all`
#     submits everything up to the SEAD upload files and stops there.
# =============================================================================
set -euo pipefail

REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
export REPO
set -a
source "${REPO}/config/paths.sh"
source "${REPO}/config/resources.sh"
set +a

SUBMIT=0; FORCE=0; PILOT=0; LENIENT=0; TAG=poe
TARGETS=()
while [[ $# -gt 0 ]]; do
  case "$1" in
    --submit)  SUBMIT=1 ;;
    --force)   FORCE=1 ;;
    --pilot)   PILOT="$2"; shift ;;
    --lenient) LENIENT=1 ;;
    --tag)     TAG="$2"; shift ;;
    -h|--help) sed -n '2,33p' "$0"; exit 0 ;;
    *)         TARGETS+=("$1") ;;
  esac
  shift
done
[[ ${#TARGETS[@]} -gt 0 ]] || { sed -n '2,33p' "$0"; exit 1; }

# -----------------------------------------------------------------------------
# Step registry
#   name | job script | array kind | dependencies | flags
#   array kinds: - (single job), groups, group_chr, chr, samples,
#                levels, batches, sim
#   flags: main (part of `all`), qc, sim, gate_sead
# -----------------------------------------------------------------------------
REGISTRY=$(cat <<'EOF'
00_reference_index | 00_reference/reference_index.lsf        | -         | -                                   | main
00_imprint_db      | 00_reference/imprint_db.lsf             | -         | -                                   | main
01_pedigree        | 01_genotype_prep/pedigree.lsf           | -         | -                                   | main
01_plink2vcf       | 01_genotype_prep/plink2vcf.lsf          | -         | 00_reference_index                  | main
01_liftover        | 01_genotype_prep/liftover.lsf          | -         | 01_plink2vcf                        | main
01_check_mendel    | 01_genotype_prep/check_mendel.lsf       | -         | 01_liftover 01_pedigree       | main
01_manifest_qc     | 01_genotype_prep/qc/manifest_qc.lsf     | -         | 00_reference_index                  | qc
02_infer_ancestry  | 02_ancestry_phasing/infer_ancestry.lsf  | -         | 01_liftover 01_pedigree       | main
02_split           | 02_ancestry_phasing/split_by_ancestry.lsf | groups  | 02_infer_ancestry                   | main
02_phase           | 02_ancestry_phasing/phase.lsf           | group_chr | 02_split 01_pedigree                | main
02_subset          | 02_ancestry_phasing/subset_by_role.lsf  | groups    | 02_phase                            | main
03_filter          | 03_imputation_SEAD/filter_sead.lsf      | group_chr | -                                   | main gate_sead
03_merge           | 03_imputation_SEAD/merge_sead.lsf       | chr       | 03_filter                           | main gate_sead
03_qc              | 03_imputation_SEAD/imputation_qc.lsf    | -         | 03_merge                            | qc
04_prep_tx         | 04_ornaments_ase/prep_transcriptome.lsf | -         | 02_infer_ancestry                   | main
04_g2t             | 04_ornaments_ase/g2t_vcf.lsf            | chr       | 03_merge 04_prep_tx                 | main
04_g2t_merge       | 04_ornaments_ase/g2t_merge.lsf          | -         | 04_g2t                              | main
04_ornaments       | 04_ornaments_ase/ornaments.lsf          | samples   | 04_g2t_merge                        | main
04_txomic_stats    | 04_ornaments_ase/txomic_stats.lsf       | -         | 04_g2t_merge                        | qc
05_assign_poe      | 05_parental_origin/assign_poe.lsf       | -         | 03_merge 01_liftover 01_pedigree    | main
06_summarize       | 06_feature_matrices/summarize.lsf       | levels    | 04_ornaments 05_assign_poe 04_prep_tx | main
07_inputs          | 07_poe_baseline/model_inputs.lsf        | -         | 06_summarize                        | main
07_within_sample   | 07_poe_baseline/within_sample.lsf       | -         | 07_inputs 00_imprint_db             | main
07_fit             | 07_poe_baseline/fit_baseline.lsf        | batches   | 07_inputs                           | main
07_combine         | 07_poe_baseline/combine_baseline.lsf    | -         | 07_fit                              | main
07_report          | 07_poe_baseline/baseline_report.lsf     | -         | 07_combine 00_imprint_db            | main
08_fit             | 08_pfas_interaction/fit_interaction.lsf | batches   | 07_inputs                           | main
08_combine         | 08_pfas_interaction/combine_interaction.lsf | -     | 08_fit                              | main
08_report          | 08_pfas_interaction/interaction_report.lsf | -      | 08_combine 00_imprint_db            | main
09_power           | 09_power/power.lsf                      | -         | 07_combine                          | main
10_sim             | 10_mediation_sim/sim.lsf                | sim       | -                                   | sim
10_aggregate       | 10_mediation_sim/aggregate.lsf          | -         | 10_sim                              | sim
EOF
)

field() { echo "$1" | awk -F'|' -v k="$2" '{gsub(/^[ \t]+|[ \t]+$/, "", $k); print $k}'; }
ALL_STEPS=$(echo "${REGISTRY}" | awk -F'|' '{gsub(/[ \t]/, "", $1); print $1}')
row_of()  { echo "${REGISTRY}" | awk -F'|' -v s="$1" '{n=$1; gsub(/[ \t]/, "", n)} n==s'; }

nwords() { echo "$*" | wc -w; }

array_size() {
  case "$1" in
    -)         echo 0 ;;
    groups)    nwords ${ANCESTRY_GROUPS} ;;
    group_chr) echo $(( $(nwords ${ANCESTRY_GROUPS}) * $(nwords ${AUTOSOMES}) )) ;;
    chr)       nwords ${AUTOSOMES} ;;
    levels)    echo 2 ;;
    batches)   echo $(( N_BATCHES_GENE + N_BATCHES_TX )) ;;
    sim)       echo "${SIM_N_RUNS}" ;;
    samples)
      # One element per RNA-seq infant. The list is written by 04_prep_tx; before it
      # exists, size the array from the metadata table (elements past the end exit 0).
      if [[ -s "${INFANT_SAMPLES}" ]]; then wc -l < "${INFANT_SAMPLES}"
      elif [[ -s "${INFANT_META_CSV}" ]]; then echo $(( $(wc -l < "${INFANT_META_CSV}") - 1 ))
      else echo 200; fi ;;
  esac
}

missing_indices() {   # step size [force] -> space-separated indices without a receipt
  local step=$1 n=$2 force=${3:-0} i out=()
  local m=$(( n > 0 ? n : 1 ))
  for (( i = 1; i <= m; i++ )); do
    if [[ ${force} -eq 1 || ! -e "${RUN_DIR}/receipts/${step}.${i}.ok" ]]; then out+=("$i"); fi
  done
  echo "${out[*]:-}"
}

compress_ranges() {   # "1 2 3 5 7 8" -> "1-3,5,7-8"
  awk '{ s=$1; p=$1; out="";
         for (i=2; i<=NF; i++) { if ($i==p+1) { p=$i; continue }
           out=out (out?",":"") (s==p ? s : s"-"p); s=$i; p=$i }
         out=out (out?",":"") (s==p ? s : s"-"p); print out }' <<< "$1"
}

sead_ready() {
  local g c
  for g in ${ANCESTRY_GROUPS}; do for c in ${AUTOSOMES}; do
    [[ -s "${SEAD_DOWNLOAD_DIR}/${g}_SEAD/chr${c}.SEAD.GRCh38.dose.vcf.gz" ]] || return 1
  done; done
}

active_jobs() {      # running/pending jobs with this name (0 if bjobs unavailable)
  command -v bjobs > /dev/null 2>&1 || { echo 0; return; }
  bjobs -J "$1" -noheader 2> /dev/null | grep -cE ' (RUN|PEND|PSUSP|USUSP|SSUSP) ' || true
}

# ---- status / list ------------------------------------------------------------
if [[ "${TARGETS[0]}" == "list" ]]; then
  echo "${REGISTRY}" | awk -F'|' '{printf "%-20s %-10s %-40s %s\n", $1, $3, $4, $5}' | tr -s ' '
  exit 0
fi
if [[ "${TARGETS[0]}" == "status" ]]; then
  printf '%-20s %10s %8s\n' STEP DONE ACTIVE
  for s in ${ALL_STEPS}; do
    r=$(row_of "$s"); n=$(array_size "$(field "$r" 3)"); m=$(( n > 0 ? n : 1 ))
    d=$(compgen -G "${RUN_DIR}/receipts/${s}.*.ok" | wc -l || true)
    printf '%-20s %10s %8s\n' "$s" "${d}/${m}" "$(active_jobs "${TAG}_${s}")"
  done
  exit 0
fi

# ---- expand targets -----------------------------------------------------------
SELECTED=()
for t in "${TARGETS[@]}"; do
  for s in ${ALL_STEPS}; do
    flags=$(field "$(row_of "$s")" 5)
    case "$t" in
      all) [[ " ${flags} " == *" main "* ]] && SELECTED+=("$s") ;;
      qc)  [[ " ${flags} " == *" qc "*   ]] && SELECTED+=("$s") ;;
      sim) [[ " ${flags} " == *" sim "*  ]] && SELECTED+=("$s") ;;
      [0-9][0-9]) [[ "$s" == "${t}_"* ]] && SELECTED+=("$s") ;;
      *)   [[ "$s" == "$t" ]] && SELECTED+=("$s") ;;
    esac
  done
done
[[ ${#SELECTED[@]} -gt 0 ]] || { echo "No step matches: ${TARGETS[*]} (see ./run_pipeline.sh list)" >&2; exit 1; }
is_selected() { local x; for x in "${SELECTED[@]}"; do [[ "$x" == "$1" ]] && return 0; done; return 1; }

# ---- build the plan -------------------------------------------------------------
mkdir -p "${RUN_DIR}/plans" "${RUN_DIR}/receipts" "${LOG_DIR}"
PLAN="${RUN_DIR}/plans/plan_$(date +%Y%m%d_%H%M%S).sh"
declare -A STATE      # step -> queued | pilot | complete | blocked
SEAD_OK=1; sead_ready || SEAD_OK=0
DEP_FN=done; [[ ${LENIENT} -eq 1 ]] && DEP_FN=ended

{
  echo "#!/usr/bin/env bash"
  echo "# plan written $(date '+%F %T') by run_pipeline.sh ${TARGETS[*]}"
  echo "set -euo pipefail"
} > "${PLAN}"

for s in ${ALL_STEPS}; do
  is_selected "$s" || continue
  row=$(row_of "$s")
  script=$(field "$row" 2); kind=$(field "$row" 3); deps=$(field "$row" 4); flags=$(field "$row" 5)
  [[ "${deps}" == "-" ]] && deps=""

  # gate: external SEAD imputation must be downloaded
  if [[ " ${flags} " == *" gate_sead "* && ${SEAD_OK} -eq 0 ]]; then
    STATE[$s]=blocked
    echo "# ${s}: waiting for SEAD imputation results in ${SEAD_DOWNLOAD_DIR}" >> "${PLAN}"
    continue
  fi

  # dependencies
  wexpr=(); blocked=0; why=""
  for d in ${deps}; do
    if [[ "${STATE[$d]:-}" == "blocked" ]]; then blocked=1; why="upstream ${d} is blocked"; break; fi
    if [[ "${STATE[$d]:-}" == "pilot" ]];   then blocked=1; why="${d} runs as a pilot only"; break; fi
    if [[ "${STATE[$d]:-}" == "queued" ]]; then wexpr+=("${DEP_FN}(${TAG}_${d})"); continue; fi
    if [[ "${STATE[$d]:-}" == "complete" ]]; then continue; fi
    # dependency not part of this submission: it must already be complete
    drow=$(row_of "$d"); dn=$(array_size "$(field "${drow}" 3)")
    if [[ -n "$(missing_indices "$d" "$dn" 0)" ]]; then
      echo "ERROR: ${s} needs ${d}, which is incomplete and not selected. Add it to the targets." >&2
      exit 1
    fi
  done
  if [[ ${blocked} -eq 1 ]]; then
    STATE[$s]=blocked; echo "# ${s}: blocked (${why})" >> "${PLAN}"; continue
  fi

  n=$(array_size "${kind}")
  todo=$(missing_indices "$s" "$n" "${FORCE}")
  if [[ -z "${todo}" ]]; then
    STATE[$s]=complete; echo "# ${s}: complete (receipts present)" >> "${PLAN}"; continue
  fi
  if [[ "$(active_jobs "${TAG}_${s}")" != "0" ]]; then
    echo "ERROR: jobs named ${TAG}_${s} are already queued or running. Wait, bkill them, or use --tag." >&2
    exit 1
  fi
  piloted=0
  if [[ ${PILOT} -gt 0 && ${n} -gt 0 && $(nwords ${todo}) -gt ${PILOT} ]]; then
    todo=$(echo ${todo} | tr ' ' '\n' | head -n "${PILOT}" | tr '\n' ' '); piloted=1
  fi

  resvar="RES_${s}"; arrvar="ARR_${s}"
  res="${!resvar:?no RES_${s} in config/resources.sh}"
  if [[ ${n} -gt 0 ]]; then
    throttle="${!arrvar:-0}"; th=""; [[ "${throttle}" != "0" ]] && th="%${throttle}"
    jname="${TAG}_${s}[$(compress_ranges "${todo}")]${th}"
    logs="-o ${LOG_DIR}/${s}.%J.%I.out -e ${LOG_DIR}/${s}.%J.%I.err"
  else
    jname="${TAG}_${s}"
    logs="-o ${LOG_DIR}/${s}.%J.out -e ${LOG_DIR}/${s}.%J.err"
  fi
  dep=""
  if [[ ${#wexpr[@]} -gt 0 ]]; then
    dep="-w \"$(printf '%s && ' "${wexpr[@]}" | sed 's/ && $//')\""
  fi
  echo "bsub -J \"${jname}\" ${res} ${logs} -env \"all,REPO=${REPO}\" ${dep} < \"${REPO}/${script}\"" >> "${PLAN}"
  if [[ ${piloted} -eq 1 ]]; then STATE[$s]=pilot; else STATE[$s]=queued; fi
done

if [[ ${SEAD_OK} -eq 0 ]] && grep -q 'waiting for SEAD' "${PLAN}"; then
  cat >> "${PLAN}" <<EOF
# ---------------------------------------------------------------------------
# NEXT (manual): upload ${SEAD_UPLOAD_DIR}/<group>/<group>_chr<N>.vcf.gz to the
# SEAD imputation server (one job per ancestry group), download the results with
#   03_imputation_SEAD/retrieve_SEAD.sh <group> <url_list.txt>
# and run this command again; finished steps are skipped.
# ---------------------------------------------------------------------------
EOF
fi

cat "${PLAN}"
if [[ ${SUBMIT} -eq 1 ]]; then
  bash "${PLAN}"
  echo "Submitted. Plan saved to ${PLAN}. Check progress with ./run_pipeline.sh status"
else
  echo
  echo "Dry run. Plan saved to ${PLAN}. Re-run with --submit to queue these jobs."
fi
