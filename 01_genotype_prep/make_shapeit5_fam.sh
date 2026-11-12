#!/bin/bash
# =============================================================================
# Stage 01 | make_shapeit5_fam.sh
# Convert the PLINK-style trio pedigree into the SHAPEIT5 --pedigree format:
# one line per offspring with >= 1 genotyped parent, columns
#   <offspring> <father> <mother>, "NA" for a missing parent.
# IDs are written as they appear in the VCF (FID_IID, PLINK's VCF naming);
# set ID_MODE=iid if the VCF uses bare IIDs.
# In : TRIO_PED     Out: TRIO_FAM_SHAPEIT5
# =============================================================================
set -euo pipefail
if [[ -z "${TRIO_PED:-}" ]]; then set -a; source "$(dirname "$0")/../config/paths.sh"; set +a; fi
: "${ID_MODE:=fid_iid}"

awk -v mode="${ID_MODE}" 'BEGIN{OFS="\t"}
  $0 ~ /^#/ || NF < 4 { next }
  {
    fid=$1; iid=$2; pid=$3; mid=$4
    if (pid == "0" && mid == "0") next                  # founders: no line
    if (mode == "fid_iid") {
      kid = fid "_" iid
      fa  = (pid == "0") ? "NA" : fid "_" pid
      mo  = (mid == "0") ? "NA" : fid "_" mid
    } else {
      kid = iid
      fa  = (pid == "0") ? "NA" : pid
      mo  = (mid == "0") ? "NA" : mid
    }
    print kid, fa, mo
  }' "${TRIO_PED}" > "${TRIO_FAM_SHAPEIT5}"

echo "wrote ${TRIO_FAM_SHAPEIT5}: $(wc -l < "${TRIO_FAM_SHAPEIT5}") offspring;" \
     "$(awk '$2=="NA"' "${TRIO_FAM_SHAPEIT5}" | wc -l) without father," \
     "$(awk '$3=="NA"' "${TRIO_FAM_SHAPEIT5}" | wc -l) without mother" >&2
