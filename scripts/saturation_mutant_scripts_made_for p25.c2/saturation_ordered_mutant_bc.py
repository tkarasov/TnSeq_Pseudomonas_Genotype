#!/usr/bin/ipython

#this script is meant to take Effie's top barcode IDS

# step 1 
# build a dictionary for ever position and every plate.
from Bio import Seq
import os
from collections import defaultdict
import pickle
import csv



plate_dir="/uufs/chpc.utah.edu/common/home/karasov-group1/lab_members/esymeonidi/novogene_dec24/01.RawData/analysis/final_bc_files/plate_files"
pos_dir="/uufs/chpc.utah.edu/common/home/karasov-group1/lab_members/esymeonidi/novogene_dec24/01.RawData/analysis/final_bc_files/position_files"

#cd /uufs/chpc.utah.edu/common/home/karasov-group1/saturation_library/barcode_pos_map/

def extract_sequences_by_prefix(filepath):
    prefix_sequences = defaultdict(list)

    with open(filepath, 'r', encoding='utf-8') as file:
        for line in file:
            line = line.strip()
            if not line:
                continue
            parts = line.split()
            if len(parts) != 2:
                continue

            full_id, sequence = parts
            id_prefix = full_id.split(':')[0]         # e.g., pl21ID014
            clean_prefix = id_prefix.split('ID')[0].split('id')[0]   # e.g., pl21

            prefix_sequences[clean_prefix].append(sequence)

    return prefix_sequences






def load_and_combine_cpks_recursive(root_dir):
    combined = defaultdict(list)

    for dirpath, _, filenames in os.walk(root_dir):
        for filename in filenames:
            if filename.endswith('.cpk'):
                filepath = os.path.join(dirpath, filename)
                print(f"Loading: {filepath}")
                with open(filepath, 'rb') as f:
                    try:
                        current_dict = pickle.load(f)
                        for key, value_list in current_dict.items():
                            combined[key].append(value_list)
                    except Exception as e:
                        print(f"❌ Failed to load {filepath}: {e}")
    
    return dict(combined)

# write a function to relabel the barcode with the gene annotation.
def add_function(inverted_dict, cpk_combined):
	for rec in inverted_dict.keys():
		try:
			functions=('').join(cpk_combined[rec])
			inverted_dict[rec].append(functions)
		except KeyError:
			pass
	return inverted_dict


# Set your directory path here
directory_path = plate_dir

# This will hold all sequences grouped by prefix (e.g., pl21, pl22, etc.)
pl_data = defaultdict(list)

# Read all files in the directory for plates
for filename in os.listdir(directory_path):
    filepath = os.path.join(directory_path, filename)
    if os.path.isfile(filepath):
        file_sequences = extract_sequences_by_prefix(filepath)
        for prefix, sequences in file_sequences.items():
            pl_data[prefix].extend(sequences)

# Convert to regular dict (optional)
pl_data = dict(pl_data)


#now do the same for the position data
directory_path = pos_dir
# This will hold all sequences grouped by prefix (e.g., pl21, pl22, etc.)
pos_data = defaultdict(list)

# Read all files in the directory for plates
for filename in os.listdir(directory_path):
    filepath = os.path.join(directory_path, filename)
    if os.path.isfile(filepath):
        file_sequences = extract_sequences_by_prefix(filepath)
        for prefix, sequences in file_sequences.items():
            pos_data[prefix].extend(sequences)

# Convert to regular dict (optional)
pos_data = dict(pos_data)

# Now I have two dictionaries pl_data and pos_data and want to find overlap.
plate_pos_list=defaultdict(list)

has_dup=defaultdict(list)
# now find overlap
for pos in pos_data.keys():
	for plate in pl_data.keys():
		#find overlap
		recs=[rec for rec in pos_data[pos] if rec in pl_data[plate]]
		if len(recs)>1:
			has_dup[plate, pos] = recs
		plate_pos_list[pos,plate]=[rec for rec in pos_data[pos] if rec in pl_data[plate]]


from collections import Counter

# Aggregate all sequences from all lists into one big list
all_sequences = []
for seq_list in plate_pos_list.values():
    all_sequences.extend(seq_list)

# Count how many times each unique sequence appears
sequence_counts = Counter(all_sequences)
#78 barcodes are observed in more than one pairing.

#first let's get rid of the worst offenders
bad_list=['GTTCGAGTCCGTAAGGGGTT',
 'GTGGGCTGCCTGCGCTCCAG',
 'AGGTCGGGTCGCGGGGGTAT',
 'ATGCATTACGTGGAATGTCT',
 'TATAGGCATTCATTTGCATG',
 'ATCGTTTGGGCCCCGTGTAA',
 'TTTCCTTCATAAGGAAGGCG',
 'TGGATGGGAGAACTTTGAGA']
 #bad_list=[rec for rec in sequence_counts if sequence_counts[rec]>3]

#let's remove these barcodes from all of the samples
to_remove = bad_list  # Convert to set for faster lookup

for key in plate_pos_list:
    plate_pos_list[key] = [seq for seq in plate_pos_list[key] if seq not in to_remove]


# now invert the dictionary so we have the key is the barcode and the plate,pos is the value

from collections import defaultdict

# Invert the dictionary
inverted = defaultdict(list)

for key, sequences in plate_pos_list.items():
    for seq in sequences:
        inverted[seq].append(key)

# Convert to regular dict if you want
inverted = dict(inverted)

#now i want to write inverted dictionary to file
with open('/uufs/chpc.utah.edu/common/home/karasov-group1/saturation_library/barcode_pos_map/bc_pl_position_04032025.txt', 'w') as f:
    for key, value_list in inverted.items():
        f.write(f"{key}\t")
        for tup in value_list:
            f.write(f"  {tup[1]},{tup[0]}\t")
        f.write("\n")




# Now let's match the barcodes to their gene. The gene characteization libraries are here:

charact_dict='/uufs/chpc.utah.edu/common/home/karasov-group1/lab_members/esymeonidi/p25c2_characterization_library/Oct24/data'


combined_dict = load_and_combine_cpks_recursive('/uufs/chpc.utah.edu/common/home/karasov-group1/lab_members/esymeonidi/p25c2_characterization_library/Oct24/data')

#I need to curate the combined dict
subset = {k: v for k, v in combined_dict.items() if isinstance(k, Seq)}
for key in subset:
    subset[key] = ('').join(subset[key]).split(";")[0]



relabel=add_function(inverted, subset)

#reorder
# def flip_and_clean(d):
#    for key in d:
#        new = [(b.split(';')[1], a) for [(a, b),c] in d[key]]
#    return d

# finished=flip_and_clean(relabel)
genes=[rec[-1] for rec in relabel.values()]



with open('/uufs/chpc.utah.edu/common/home/karasov-group1/saturation_library/barcode_pos_map/barcode_identity.txt', 'w') as f:
    for key, value in relabel.items():
        # Flatten the value list (tuples -> items, strings -> just as-is)
        flat = []
        for item in value:
            if isinstance(item, tuple):
                flat.extend(item)
            else:
                flat.append(item)
        line = '\t'.join([key] + flat)
        f.write(line + '\n')


