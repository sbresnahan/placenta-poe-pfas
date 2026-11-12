#!/usr/bin/env python3
# =============================================================================
# Stage 01 (QC) | convert_manifest_to_vcf.py
# Align the Illumina array manifest (SNP_Table) to the hg19 reference: for each
# SNP try the alleles as given, swapped, strand-flipped and flipped+swapped and
# record which correction made REF match the genome.
# Usage: convert_manifest_to_vcf.py <manifest_subset.txt> <out.vcf>
# Env  : HG19_FASTA (indexed), NCPU (worker processes)
# =============================================================================
import csv
import pysam
from tqdm import tqdm
from concurrent.futures import ProcessPoolExecutor
import multiprocessing
import sys
import os

# --- Configuration ---
# Indexed hg19 FASTA (chr-prefixed or not; both names are tried).
REF_FASTA = os.environ["HG19_FASTA"]


# --- Global Variables ---
# Reverse complement dictionary for flipping alleles
REV_COMP = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N'}


def init_worker(fasta_path):
    """
    Initializer for each worker process.
    This function creates a new, independent pysam.FastaFile handle for each
    worker, preventing file handle corruption from concurrent access.
    """
    global ref
    try:
        # Each worker process gets its own 'ref' object.
        ref = pysam.FastaFile(fasta_path)
    except FileNotFoundError:
        print(f"WORKER ERROR: Reference FASTA file or its index (.fai) not found at {fasta_path}", file=sys.stderr)
        # Exit the worker process if the reference can't be opened.
        sys.exit(1)


def reverse_complement(base):
    """Calculates the reverse complement of a single DNA base."""
    return REV_COMP.get(base, "N")

def process_row(row):
    """
    Processes a single row from the input manifest to align it to the reference genome.
    The logic ignores any strand information and attempts all possibilities for each SNP:
    1. Parse the chromosome, position, and alleles from the row.
    2. Fetch the reference base from the reference FASTA file using the process-local handle.
    3. Try to match the manifest alleles (no change or swapped).
    4. If no match, flip the alleles and try to match again (flipped or flipped+swapped).
    5. Report the final REF/ALT alleles and the correction applied in the INFO field.
    """
    # --- 1. Parse and Validate Input ---
    try:
        chrom_in = row["Chr"].strip()
        pos = int(row["Position"])
        rsid = row["Name"].strip()
        # The 'SNP' column is expected in '[A/G]' format
        allele1, allele2 = row["SNP"].strip("[]").split("/")
        allele1 = allele1.strip().upper()
        allele2 = allele2.strip().upper()
    except (ValueError, KeyError, AttributeError):
        # Return an error if parsing fails for any reason
        return [
            row.get("Chr", "."), row.get("Position", "."), row.get("Name", "."), ".", ".", ".", ".",
            "Correction=error_parsing_row"
        ]

    # --- 2. Fetch Reference Base ---
    chrom = chrom_in
    try:
        # pysam is 0-based, VCF is 1-based. fetch(pos-1, pos) gets the single base at `pos`.
        ref_base = ref.fetch(chrom, pos - 1, pos).strip().upper()
    except ValueError:
        # If that fails, it's often due to a 'chr' prefix mismatch.
        # We try to fix this automatically by adding/removing 'chr'.
        fixed_chrom = 'chr' + chrom if not chrom.startswith('chr') else chrom.replace('chr', '')
        try:
            ref_base = ref.fetch(fixed_chrom, pos - 1, pos).strip().upper()
            chrom = fixed_chrom # Update chrom name if the fix worked
        except ValueError:
            # If both attempts fail, we can't resolve the reference base.
            return [chrom_in, pos, rsid, allele1, allele2, ".", ".", "Correction=error_ref_fetch"]

    # --- 3. Attempt to Resolve Alleles ---

    # Step 3a: Try manifest alleles
    if allele1 == ref_base:
        return [chrom, pos, rsid, allele1, allele2, ".", ".", "Correction=none"]
    if allele2 == ref_base:
        return [chrom, pos, rsid, allele2, allele1, ".", ".", "Correction=swapped"]

    # Step 3b: If no match, flip both alleles and try again
    f1 = reverse_complement(allele1)
    f2 = reverse_complement(allele2)

    if f1 == ref_base:
        return [chrom, pos, rsid, f1, f2, ".", ".", "Correction=flipped"]
    if f2 == ref_base:
        return [chrom, pos, rsid, f2, f1, ".", ".", "Correction=flipped+swapped"]

    # --- 4. Final Fallback ---
    # If none of the above combinations match the reference, mark as unresolved.
    return [
        chrom, pos, rsid, ref_base, f"{allele1},{allele2}", ".", ".",
        "Correction=unresolved_mismatch"
    ]


def read_manifest(path):
    """Reads the tab-delimited input file into a list of dictionaries."""
    try:
        with open(path, 'r', encoding='utf-8') as f:
            reader = csv.DictReader(f, delimiter='\t')
            return list(reader)
    except FileNotFoundError:
        print(f"FATAL: Input file not found at {path}", file=sys.stderr)
        sys.exit(1)


def write_vcf_header(outf, ref_path):
    """Writes a standard VCF header, including contigs from the reference."""
    # We need a temporary FastaFile handle here just to write the header contigs.
    # This doesn't interfere with the worker processes.
    with pysam.FastaFile(ref_path) as header_ref:
        outf.write("##fileformat=VCFv4.2\n")
        outf.write(f"##reference=file://{ref_path}\n")
        # Add contig headers from the reference file
        for contig in header_ref.references:
            outf.write(f"##contig=<ID={contig},length={header_ref.get_reference_length(contig)}>\n")
    outf.write('##INFO=<ID=Correction,Number=1,Type=String,Description="Describes the correction applied to the manifest alleles to match the reference. Can be none, swapped, flipped, flipped+swapped, or an error status.">\n')
    outf.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")


def main():
    """Main function to orchestrate reading, processing, and writing."""
    if len(sys.argv) != 3:
        print(f"Usage: python {sys.argv[0]} <input_manifest.txt> <output.vcf>", file=sys.stderr)
        sys.exit(1)

    input_file = sys.argv[1]
    output_file = sys.argv[2]

    print(f"Reading manifest: {input_file}")
    rows = read_manifest(input_file)
    total = len(rows)
    if total == 0:
        print("Input file is empty. Exiting.")
        return

    num_workers = max(1, int(os.environ.get("NCPU", "1")))
    print(f"Processing {total} SNPs using {num_workers} cores...")
    
    with open(output_file, "w", newline='') as out_f:
        writer = csv.writer(out_f, delimiter='\t')
        write_vcf_header(out_f, REF_FASTA)
        
        # Use an initializer to create a new FastaFile handle in each worker process
        with ProcessPoolExecutor(
            max_workers=num_workers,
            initializer=init_worker,
            initargs=(REF_FASTA,)
        ) as executor:
            # Use executor.map for parallel processing and tqdm for a progress bar
            results = tqdm(executor.map(process_row, rows, chunksize=5000), total=total, desc="Aligning SNPs to reference")
            writer.writerows(results)

    print(f"\nProcessing complete. Output written to {output_file}")


if __name__ == "__main__":
    main()

