#!/usr/bin/env bash
# =============================================================================
# reports/stage_figures.sh -- copy the figure files the report Rmds need into
# reports/figures/{methods,interaction}/, following reports/figures/MANIFEST.tsv.
#
# Usage:
#   reports/stage_figures.sh --local    <Placenta_POE folder>   # exact panels from the local project tree
#   reports/stage_figures.sh --pipeline [FIG_DIR]               # regenerate from pipeline output
#   reports/stage_figures.sh --check                            # list what is still missing
#
# --local    uses the "local_source" column (e.g. Reorganize/figure2.png,
#            Reorganize/all_PFAS/gene_forest_PFBS_mat.png). These are the files
#            the rendered reports in docs/reports/ were built from.
# --pipeline uses the "pipeline_source" column, relative to FIG_DIR (default:
#            $FIG_DIR, else the value in config/paths.sh). Methods panels from the
#            pipeline are restyled, and figure5.png has no pipeline source; any
#            file already present in reports/figures/ is kept (use --local for those).
# Existing files are overwritten only with --force.
# =============================================================================
set -euo pipefail
here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repo="$(dirname "$here")"
manifest="$here/figures/MANIFEST.tsv"
dest="$here/figures"

mode="${1:-}"; src="${2:-}"; force=0
for a in "$@"; do [[ "$a" == "--force" ]] && force=1; done
[[ "$src" == "--force" ]] && src=""

case "$mode" in
  --local)
    [[ -n "$src" && -d "$src" ]] || { echo "usage: $0 --local <Placenta_POE folder> [--force]" >&2; exit 2; }
    col=3 ;;
  --pipeline)
    if [[ -z "$src" ]]; then
      src="$(bash -c "source '$repo/config/paths.sh' >/dev/null 2>&1; echo \"\$FIG_DIR\"")"
    fi
    [[ -d "$src" ]] || { echo "FIG_DIR not found: '$src'" >&2; exit 2; }
    col=4 ;;
  --check) col=0 ;;
  *) sed -n '3,19p' "$0" | sed 's/^# \{0,1\}//'; exit 2 ;;
esac

copied=0; kept=0; missing=()
while IFS=$'\t' read -r report file local_src pipe_src note; do
  [[ "$report" == "report" ]] && continue
  out="$dest/$file"
  if [[ $col -eq 0 ]]; then
    [[ -f "$out" ]] || missing+=("$file")
    continue
  fi
  rel="$local_src"; [[ $col -eq 4 ]] && rel="$pipe_src"
  if [[ -f "$out" && $force -eq 0 ]]; then kept=$((kept+1)); continue; fi
  [[ "$rel" == "-" ]] && rel=""
  if [[ -z "$rel" || ! -f "$src/$rel" ]]; then
    missing+=("$file  <-  ${rel:-<no source>}"); continue
  fi
  mkdir -p "$(dirname "$out")"
  cp -p "$src/$rel" "$out"; copied=$((copied+1))
done < "$manifest"

[[ $col -ne 0 ]] && echo "copied: $copied   already present (kept): $kept"
if (( ${#missing[@]} )); then
  echo "missing: ${#missing[@]}"; printf '  %s\n' "${missing[@]}"; exit 1
fi
echo "all $(($(wc -l < "$manifest") - 1)) report figures present in reports/figures/"
