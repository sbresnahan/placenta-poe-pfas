#!/usr/bin/env python3
# =============================================================================
# Stage 05 | assign_POE.py
# Assign each phased infant haplotype (left "L" / right "R" of the GT field)
# to its maternal or paternal origin, per chromosome, from the parents'
# unphased genotypes.
#
# Per chromosome:
#   1. parental panel  = biallelic SNPs genotyped in >= (100 - MISSING_RATE_PCT)%
#                        of all parents in the unphased trio VCF
#   2. informative     = panel sites where this infant's mother and father are
#                        homozygous for different alleles (REF/ALT must match
#                        between the trio VCF and the imputed infant VCF)
#   3. vote            = at each informative site where the infant is a phased
#                        heterozygote, the L allele votes "maternal" if it equals
#                        the mother's allele, otherwise "paternal"; R takes the
#                        opposite. Majority wins; ties are "undefined".
#
# Genotypes are streamed with `bcftools query` (one pass over each VCF per
# chromosome); chromosomes run in parallel.
#
# Usage:
#   assign_POE.py --unphased_vcf TRIO_VCF_HG38 --ped TRIO_PED \
#                 --phased_dir IMPUTED_DIR --out HAP_ORIGIN --workers N
# Requires: python >= 3.8, bcftools on PATH
# Out: TSV, one row per chromosome (chr1..chr22), columns <FID>_<IID>_L, <FID>_<IID>_R
#      (infant VCF sample ID + haplotype),
#      values maternal | paternal | undefined
#      plus <out>.votes.tsv: per (chr, infant) informative sites used and votes
# =============================================================================
import argparse
import os
import subprocess
import sys
import tempfile
import time
from concurrent.futures import ProcessPoolExecutor

MISSING_RATE_PCT = 10
PHASED_PATTERN = "{chrom}.biallelic.filtered.vcf.gz"


def log(msg):
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", file=sys.stderr, flush=True)


def parse_ped(path):
    """infant IID -> (FID, father IID or None, mother IID or None)."""
    ped = {}
    with open(path) as fh:
        for line in fh:
            f = line.split()
            if len(f) < 4:
                continue
            fid, iid, pid, mid = f[:4]
            ped[iid] = (fid, None if pid == "0" else pid, None if mid == "0" else mid)
    return ped


def vcf_id(fid, iid):
    return f"{fid}_{iid}"


def vcf_samples(path):
    out = subprocess.run(["bcftools", "query", "-l", path], check=True, capture_output=True, text=True)
    return out.stdout.split()


def query(path, chrom, samples, targets=None):
    """Yield (pos, ref, alt, [GT, ...]) for one chromosome and a sample subset."""
    with tempfile.NamedTemporaryFile("w", suffix=".txt", delete=False) as fh:
        fh.write("\n".join(samples) + "\n")
        sfile = fh.name
    cmd = ["bcftools", "query", "-r", chrom, "-S", sfile, "-f", "%POS\t%REF\t%ALT[\t%GT]\n"]
    if targets:
        cmd[4:4] = ["-T", targets]
    try:
        with subprocess.Popen(cmd + [path], stdout=subprocess.PIPE, text=True, bufsize=1 << 20) as p:
            for line in p.stdout:
                f = line.rstrip("\n").split("\t")
                yield int(f[0]), f[1], f[2], f[3:]
        if p.returncode:
            raise RuntimeError(f"bcftools query failed ({p.returncode}): {' '.join(cmd)} {path}")
    finally:
        os.unlink(sfile)


def hom_allele(gt):
    """0 or 1 for a homozygous REF/ALT call, else None."""
    if gt in ("0/0", "0|0"):
        return 0
    if gt in ("1/1", "1|1"):
        return 1
    return None


def called(gt):
    return len(gt) == 3 and gt[0] != "." and gt[2] != "."


def decide(m, p):
    if m > p:
        return "maternal"
    if p > m:
        return "paternal"
    return "undefined"


def run_chrom(chrom, unphased_vcf, phased_path, parents, trios, tmpdir):
    """trios: list of (infant IID, infant VCF id, mother VCF id, father VCF id)."""
    t0 = time.time()
    if not os.path.exists(phased_path):
        log(f"{chrom}: missing {phased_path}")
        return chrom, {}, 0

    # ---- pass 1: parental panel + informative sites per trio ----
    pidx = {s: i for i, s in enumerate(parents)}
    trio_cols = [(k, pidx[m], pidx[d]) for k, (_, _, m, d) in enumerate(trios)]
    max_missing = MISSING_RATE_PCT / 100.0
    n_par = len(parents)
    informative = {}     # pos -> (ref, alt, [(trio k, mother allele), ...])
    n_panel = 0
    for pos, ref, alt, gts in query(unphased_vcf, chrom, parents):
        if len(ref) != 1 or len(alt) != 1:
            continue
        present = sum(1 for g in gts if called(g))
        if 1 - present / n_par >= max_missing:
            continue
        n_panel += 1
        hits = []
        for k, mi, di in trio_cols:
            mh = hom_allele(gts[mi])
            if mh is None:
                continue
            dh = hom_allele(gts[di])
            if dh is None or dh == mh:
                continue
            hits.append((k, mh))
        if hits:
            informative[pos] = (ref, alt, hits)

    if not informative:
        log(f"{chrom}: panel {n_panel:,} sites, none informative")
        return chrom, {}, n_panel

    # ---- pass 2: phased infant genotypes at informative sites ----
    targets = os.path.join(tmpdir, f"{chrom}.targets.tsv")
    with open(targets, "w") as fh:
        for pos in sorted(informative):
            fh.write(f"{chrom}\t{pos}\n")
    infants = [t[1] for t in trios]
    mat = [0] * len(trios)   # votes for "L is maternal"
    pat = [0] * len(trios)   # votes for "L is paternal"
    for pos, ref, alt, gts in query(phased_path, chrom, infants, targets):
        site = informative.get(pos)
        if site is None or site[0] != ref or site[1] != alt:
            continue
        for k, mh in site[2]:
            g = gts[k]
            if len(g) != 3 or g[1] != "|" or g[0] == g[2] or "." in g:
                continue                       # phased heterozygotes only
            if int(g[0]) == mh:
                mat[k] += 1
            else:
                pat[k] += 1

    res = {trios[k][0]: (decide(mat[k], pat[k]), decide(pat[k], mat[k]), mat[k], pat[k])
           for k in range(len(trios))}
    n_def = sum(1 for v in res.values() if v[0] != "undefined")
    log(f"{chrom}: panel {n_panel:,}, informative {len(informative):,}, "
        f"assigned {n_def}/{len(trios)} infants ({time.time() - t0:.0f}s)")
    return chrom, res, n_panel


def main():
    ap = argparse.ArgumentParser(description=__doc__ or "Assign parental origin of phased infant haplotypes.")
    ap.add_argument("--unphased_vcf", required=True, help="unphased trio VCF (bgzipped + indexed)")
    ap.add_argument("--ped", required=True, help="pedigree: FID IID PID MID ... (0 = missing)")
    ap.add_argument("--phased_dir", required=True, help="directory with chr<N>.biallelic.filtered.vcf.gz")
    ap.add_argument("--out", required=True, help="output TSV")
    ap.add_argument("--workers", type=int, default=1)
    ap.add_argument("--chroms", default=",".join(f"chr{i}" for i in range(1, 23)))
    a = ap.parse_args()
    chroms = a.chroms.split(",")

    ped = parse_ped(a.ped)
    if not ped:
        sys.exit("no rows in pedigree")
    unphased_ids = set(vcf_samples(a.unphased_vcf))
    first = next((os.path.join(a.phased_dir, PHASED_PATTERN.format(chrom=c)) for c in chroms
                  if os.path.exists(os.path.join(a.phased_dir, PHASED_PATTERN.format(chrom=c)))), None)
    if first is None:
        sys.exit(f"no phased VCFs in {a.phased_dir}")
    phased_ids = set(vcf_samples(first))

    trios, parents = [], set()
    for iid in sorted(ped):
        fid, dad, mom = ped[iid]
        if not dad or not mom:
            continue
        inf, m, d = vcf_id(fid, iid), vcf_id(fid, mom), vcf_id(fid, dad)
        if inf in phased_ids and m in unphased_ids and d in unphased_ids:
            trios.append((iid, inf, m, d))
    for fid, dad, mom in ped.values():   # panel missingness is over all genotyped parents
        for p in (dad, mom):
            if p and vcf_id(fid, p) in unphased_ids:
                parents.add(vcf_id(fid, p))
    parents = sorted(parents)
    log(f"pedigree infants {len(ped)}, complete trios with imputed infant {len(trios)}, "
        f"parents {len(parents)}, workers {a.workers}")

    results = {}
    with tempfile.TemporaryDirectory() as tmp, ProcessPoolExecutor(max_workers=max(1, a.workers)) as ex:
        futs = [ex.submit(run_chrom, c, a.unphased_vcf,
                          os.path.join(a.phased_dir, PHASED_PATTERN.format(chrom=c)),
                          parents, trios, tmp) for c in chroms]
        for f in futs:
            chrom, res, _ = f.result()
            results[chrom] = res

    infants = sorted(ped)
    col = {i: vcf_id(ped[i][0], i) for i in infants}
    os.makedirs(os.path.dirname(os.path.abspath(a.out)), exist_ok=True)
    with open(a.out, "w") as out, open(a.out + ".votes.tsv", "w") as vout:
        out.write("chr\t" + "\t".join(f"{col[i]}_{s}" for i in infants for s in ("L", "R")) + "\n")
        vout.write("chr\tinfant\tL_maternal_votes\tL_paternal_votes\n")
        for c in chroms:
            row = [c]
            for i in infants:
                L, R, m, p = results[c].get(i, ("undefined", "undefined", 0, 0))
                row += [L, R]
                if i in results[c]:
                    vout.write(f"{c}\t{col[i]}\t{m}\t{p}\n")
            out.write("\t".join(row) + "\n")
    log(f"wrote {a.out}")


if __name__ == "__main__":
    main()
