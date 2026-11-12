#!/bin/bash
# =============================================================================
# Stage 03 | retrieve_SEAD.sh      (run by hand on a login/transfer node)
# Download one ancestry group's results from the SEAD imputation server.
#
# Usage:
#   REPO=/path/to/placenta-poe-pfas bash retrieve_SEAD.sh <group> <url_list.txt>
#     group        : one of ANCESTRY_GROUPS (chinese | indian | malay)
#     url_list.txt : one download URL per line, copied from the SEAD job page.
#                    The URLs carry job-specific access tokens: keep the list
#                    outside the repository (config/sead_urls_*.txt is gitignored).
#
# Out: SEAD_DOWNLOAD_DIR/<group>_SEAD/chr<N>.SEAD.GRCh38.dose.vcf.gz (+ .info, .md5)
# Files must carry exactly these names; rename them if the server's differ.
# =============================================================================
set -euo pipefail
REPO=${REPO:-$(cd "$(dirname "$0")/.." && pwd)}
source "${REPO}/config/paths.sh"

GROUP=${1:?usage: retrieve_SEAD.sh <group> <url_list.txt>}
URLS=$(readlink -f "${2:?usage: retrieve_SEAD.sh <group> <url_list.txt>}")
[[ " ${ANCESTRY_GROUPS} " == *" ${GROUP} "* ]] || { echo "unknown group: ${GROUP}" >&2; exit 1; }

OUT="${SEAD_DOWNLOAD_DIR}/${GROUP}_SEAD"
mkdir -p "${OUT}"
cd "${OUT}"

# resumable download; --content-disposition keeps server-side file names
wget -c --content-disposition -i "${URLS}"

# verify checksums where .md5 files were provided
shopt -s nullglob
for m in *.md5; do md5sum -c "${m}"; done

missing=0
for c in ${AUTOSOMES}; do
  [[ -s "chr${c}.SEAD.GRCh38.dose.vcf.gz" ]] || { echo "missing: chr${c}.SEAD.GRCh38.dose.vcf.gz" >&2; missing=1; }
done
[[ ${missing} -eq 0 ]] && echo "${GROUP}: all dose VCFs present in ${OUT}"
