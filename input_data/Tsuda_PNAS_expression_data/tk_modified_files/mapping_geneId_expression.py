/Users/talia/opt/anaconda3/bin/ipython
#this file maps the gene IDs from the Tsuda (which were genbank ids) to the protein ids from Aubrey (which were Refseq protein ids). Talia did this on 7/22/2024. Seems to be imperfect for putting in the RS numbers but the mapping between PSPTO and WP seems to be OK for now.

expression=[line.strip().split('\t') for line in open("/Users/talia/Documents/GitHub/TnSeq_Pseudomonas_Genotype/input_data/Tsuda_PNAS_expression_data/tk_modified_files/pnas.1800529115_tk_modified.sd06.txt").readlines()]
	
gene_mappings=[line.strip().split('\t') for line in open("/Users/talia/Documents/GitHub/TnSeq_Pseudomonas_Genotype/input_data/Tsuda_PNAS_expression_data/tk_modified_files/GCF_000007805.1_ASM780v1_feature_table_tlk.txt").readlines()]


#for every line in the expression file I want to pull out the corresponding id in the gene mappings. 

for line in expression[2:]:
	PSPTO=line[0]
	mapped=[rec for rec in gene_mappings if rec[-1]==PSPTO]
	if(len(mapped)>=1):
		print("yes")
		WP=mapped[0][10]
		try:
			PSPTO_RS=mapped[0][16]
			line.append(WP)
			line.append(PSPTO_RS)
		except:
			IndexError()

#Write this mapping to file
handle="/Users/talia/Documents/GitHub/TnSeq_Pseudomonas_Genotype/input_data/Tsuda_PNAS_expression_data/tk_modified_files/expression_gene_mapping.txt"

with open(handle, "w", newline='') as file:
	writer=csv.writer(file, delimiter="\t")
	for line in expression:
		writer.writerow(line)

#handle.close()
