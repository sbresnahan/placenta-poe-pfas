#!/bin/bash
# =============================================================================
# tools/lsf_usage_report.sh -- peak memory and run time of finished jobs
#
# Usage: tools/lsf_usage_report.sh <step-name-prefix>     e.g. 04_ornaments
# Reads the LSF job reports in LOG_DIR/<step>*.out and prints one line per job
# plus the maximum, to set RES_<step> in config/resources.sh after a pilot.
# =============================================================================
set -euo pipefail
REPO=${REPO:-$(cd "$(dirname "$0")/.." && pwd)}
source "${REPO}/config/paths.sh"
STEP=${1:?usage: lsf_usage_report.sh <step-name-prefix>}

shopt -s nullglob
files=("${LOG_DIR}/${STEP}"*.out)
[[ ${#files[@]} -gt 0 ]] || { echo "no logs matching ${LOG_DIR}/${STEP}*.out" >&2; exit 1; }

printf '%-45s %-10s %12s %12s %10s\n' log status max_mem_MB run_time_s cpu_time_s
for f in "${files[@]}"; do
  awk -v f="$(basename "$f")" '
    /^Successfully completed/ { st = "ok" }
    /^Exited with/            { st = "failed" }
    /Max Memory :/            { mem = $4; if ($5 ~ /GB/) mem *= 1024 }
    /Run time :/              { rt = $4 }
    /CPU time :/              { cpu = $4 }
    END { printf "%-45s %-10s %12.0f %12.0f %10.0f\n", f, st, mem, rt, cpu }' "$f"
done | sort -k3,3n | tee /dev/stderr \
  | awk '{ if ($3 > m) m = $3; if ($4 > t) t = $4 }
         END { printf "\nmax memory %.0f MB (%.1f GB), max run time %.0f s (%.1f h)\n", m, m/1024, t, t/3600 }'
