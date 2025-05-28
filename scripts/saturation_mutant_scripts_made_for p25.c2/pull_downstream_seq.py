#!/usr/bin/ipython

# the goal of this script is to take the barcode characterization files and pull sequence 1,000bp downstream from the location of the insertion. This script is used for designing primers for each of the insertion mutants.

#this had to be run on all 37 plates with the script run_pull_downstream_seq.sh




from Bio import SeqIO
from Bio.Seq import Seq
import subprocess
import glob
import os 
import gzip
import sys
# pysam failed to install repeatedly


plate = sys.argv[1]  # e.g. 'plate01'


# === Config ===
base_dir = "/uufs/chpc.utah.edu/common/home/karasov-group1/lab_members/esymeonidi/p25c2_characterization_library/Oct24/data"
genome="/uufs/chpc.utah.edu/common/home/karasov-group1/lab_members/esymeonidi/p25c2_characterization_library/reference_p25c2/p25.c2_correct_ref_gen.fasta"
#bam_file="/uufs/chpc.utah.edu/common/home/karasov-group1/lab_members/esymeonidi/p25c2_characterization_library/Oct24/data/plate01/plate01_CKDL240032762-1A_22GKHFLT4_L3_2.fq.gz.bam"

#barcode_file="/uufs/chpc.utah.edu/common/home/karasov-group1/lab_members/esymeonidi/p25c2_characterization_library/Oct24/data/plate01/trimmed_final_plate01_R1.fq.gz"

plate_dir = os.path.join(base_dir, plate)

# Match BAM using wildcard since it contains flowcell info
bam_candidates = glob.glob(os.path.join(plate_dir, f"{plate}_*.bam"))

# Construct FASTQ filename using plate name
fq_file = os.path.join(plate_dir, f"trimmed_final_{plate}_R1.fq.gz")

if not bam_candidates or not os.path.exists(fq_file):
    print(f"⚠️ Skipping {plate}: missing BAM or FASTQ")
    sys.exit(1)

bam_file = bam_candidates[0]



output_dir = "/uufs/chpc.utah.edu/common/home/karasov-group1/saturation_library/barcode_pos_map"
os.makedirs(output_dir, exist_ok=True)  # make sure it exists

output_fasta = os.path.join(output_dir, f"{plate}.fasta")

plate_dir = os.path.join(base_dir, plate)
bam_candidates = glob.glob(os.path.join(plate_dir, "*.bam"))
fq_candidates = glob.glob(os.path.join(plate_dir, "*.fastq")) + \
                glob.glob(os.path.join(plate_dir, "*.fastq.gz")) + \
                glob.glob(os.path.join(plate_dir, "*.fq.gz"))

if not bam_candidates or not fq_candidates:
    print(f"⚠️ Skipping {plate}: missing files")
    sys.exit(0)

bam_file = bam_candidates[0]
fq_file = fq_file #fq_candidates[0]



# Path to your BAM file
# output_tsv = "/uufs/chpc.utah.edu/common/home/karasov-group1/saturation_library/barcode_pos_map/bam_read_positions.tsv"






# === Load genome ===
genome = SeqIO.to_dict(SeqIO.parse(genome, "fasta"))


# === Barcode loader ===
# Load barcodes
def read_fastq_barcodes(fq_file):
    open_func = gzip.open if fq_file.endswith(".gz") else open
    barcodes = {}
    with open_func(fq_file, "rt") as f:
        while True:
            name = f.readline().strip()
            seq = f.readline().strip()
            f.readline()
            f.readline()
            if not name or not seq:
                break
            read_id = name[1:].split()[0]
            barcodes[read_id] = seq
    return barcodes

barcodes = read_fastq_barcodes(fq_file)
seen_regions = set()

with subprocess.Popen(
    ["samtools", "view", bam_file],
    stdout=subprocess.PIPE,
    text=True
) as proc, open(output_fasta, "w") as out_f:

    for line in proc.stdout:
        fields = line.strip().split("\t")
        if len(fields) < 5:
            continue

        original_read_id = fields[0]
        flag = int(fields[1])
        ref = fields[2]
        pos = int(fields[3]) - 1

        if ref not in genome:
            continue

        barcode = barcodes.get(original_read_id)
        if not barcode:
            continue

        strand = "-" if (flag & 16) else "+"

        if strand == "+":
            start = pos
            end = pos + 1000
        else:
            end = pos + 1
            start = max(0, end - 1000)

        if start < 0 or end > len(genome[ref]):
            continue

        region_key = (ref, start, end, strand)
        if region_key in seen_regions:
            continue
        seen_regions.add(region_key)

        seq = genome[ref].seq[start:end]

        header = f"{barcode}_{ref}:{start+1}-{end}({strand})_{plate}"
        out_f.write(f">{header}\n{seq}\n")

print(f"✅ Wrote FASTA for {plate} to {output_fasta}")



# the trimmed final file is the barcode associated with each of the reads. The key is to match the reads in the bam with each of the barcodes according to the id.
barcode_label="/uufs/chpc.utah.edu/common/home/karasov-group1/lab_members/esymeonidi/p25c2_characterization_library/Oct24/data/plate01/trimmed_final_plate01_R1.fq.gz"