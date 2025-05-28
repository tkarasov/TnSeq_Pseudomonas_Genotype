#!/usr/bin/ipython

# the goal of this script is to take the output from saturation_ordered_mutant_bc.py and pull_downstream_seq.py into an output dataframe that has the barcode the gene and 1000bp downstream from the insertion site.



from Bio import SeqIO




# === Config ===
base_dir = "/uufs/chpc.utah.edu/common/home/karasov-group1/saturation_library/barcode_pos_map"
txt_file= base_dir+"/barcode_identity.txt"
fasta_file=base_dir+"/all_plates_combined.fasta"
output_file = base_dir+"/barcode_identity_annotated_output.txt"



# Step 1: Map barcodes to sequences from FASTA
barcode_to_seq = {}
for record in SeqIO.parse(fasta_file, "fasta"):
    barcode = record.id.split("_")[0]
    barcode_to_seq[barcode] = str(record.seq)

# Step 2: Annotate and clean
with open(txt_file, "r") as infile, open(output_file, "w") as out:
    for line in infile:
        parts = line.strip().split("\t")
        if len(parts) < 3:
            continue  # skip malformed lines

        # Clean all fields of "ID=" if present
        parts = [field.replace("ID=", "") for field in parts]

        barcode = parts[0]
        sequence = barcode_to_seq.get(barcode, "NA")

        # Insert sequence at column 4 (index 3)
        new_parts = parts[:3] + [sequence] + parts[3:]
        out.write("\t".join(new_parts) + "\n")

print(f"✅ Cleaned and annotated file written to: {output_file}")