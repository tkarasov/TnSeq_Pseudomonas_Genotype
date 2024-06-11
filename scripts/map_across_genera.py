'''the goal of this script is to take the mappings of panX run across genera (for example: within gammaproteobacteria) and to connect the panX output results between two species.
'''

import os,sys,time,csv,glob, pickle, cPickle
from collections import defaultdict
from Bio import GenBank
from Bio import SeqIO


def convert_fasta_genbank(path_to_faa_dir):
	for file in os.listdir(path_to_faa_dir):
		if 'faa' in file:
			temp = list(SeqIO.parse(open(file), 'fasta'))
			keep_all=[]
			for rec in temp:
				new = SeqRecord(Seq(str(rec.seq), IUPAC.protein))
				new.id = rec.id.split(":")[1] 
				new.name = rec.name.split(":")[1]
				new.description = rec.description

				#rec.id=rec.id.split(":")[1]
				#rec.name = rec.name.split(":")[1]
				keep_all.append(new)
			file_gbk = file.strip('.faa')+'.gbk'
			SeqIO.write(keep_all, file_gbk, "genbank")




def load_sorted_clusters(path_to_pangenome_dir):
    '''
    load gene clusters and sort 1st by abundance and then by clusterID
    '''
    geneClusterPath='%s%s'%(path_to_pangenome_dir,'protein_faa/diamond_matches/')
    geneCluster_dt=load_pickle(geneClusterPath+'allclusters_postprocessed.cpk')
    from operator import itemgetter
    # sort by decreasing abundance (-v[0], minus to achieve decreasing)
    # followed by increasing strain count
    return sorted(geneCluster_dt.iteritems(),
                key=lambda (k,v): (-itemgetter(0)(v),itemgetter(2)(v)), reverse=False)
    #return sorted(geneCluster_dt.iteritems(),
    #            key=lambda (k,v): (-itemgetter(0)(v),itemgetter(2)(v)), reverse=False)


def map_redone(g_map):
	redone={}
	for rec in g_map.keys():
		temp=rec.split(":")[1]
		new_key=g_map[rec][0]
		redone[new_key]=temp
	return redone

def edit_distance(seq1, seq2):
	count = sum(1 for a, b in zip(seq1, seq2) if a != b)
	return count


def map_locus_tags():



path_to_pangenome_dir='/ebio/abt6_projects9/tnseq/wdata/panX-tnseq/data/Gammaproteobacteria/'
path_to_fitness_tables='/ebio/abt6_projects9/tnseq/data/fitness_datasets/fitness_tables'
path_to_faa_dir='/ebio/abt6_projects9/tnseq/data/fitness_datasets/tnseq_ref_genomes'
clusters = path_to_pangenome_dir + "/allclusters_final.tsv"

#cd /ebio/abt6_projects9/tnseq/data/across_genera

#list of every ortholog, the genes assigned to that orthology group and the gene number
sorted_genelist = load_sorted_clusters(path_to_pangenome_dir)

#SNP mappings
SNP_matrix=list(SeqIO.parse(open("/ebio/abt6_projects9/tnseq/wdata/panX-tnseq/data/Gammaproteobacteria/geneCluster/SNP_whole_matrix.aln"), 'fasta'))
SNP_dist=pd.DataFrame(columns=[species.name for species in SNP_matrix], index=[species.name for species in SNP_matrix])
for species in SNP_matrix:
	for spec in SNP_matrix:
		count=edit_distance(species.seq, spec.seq)
		SNP_dist[species.name][spec.name]=float(count)/len(spec.seq)

SNP_dist.to_csv("distance_matrix_across_strains.txt")



def red_conditions(fit_table):
	#reduce condition names
	cond=fit_table.columns
	cond_dict = {}
	for rec in cond:
		if 'set' in rec:
			new = rec.partition(' ')[2]
			cond_dict[rec] = new
	return cond_dict

def matching_conditions(words, cond_dict):
	for word in words:
		is_match_value = [s for s in cond_dict.values() if word in s]#[rec for rec in cond_dict if word in cond_dict[rec]]
		is_match_key =  [key for key in cond_dict if cond_dict[key] in is_match_value]
		if len(is_match_key)>0:
			return is_match_key[0]
		else:
			print "No match"



def rs_genecluster_map(sorted_genelist):
	gene_dir = {}
	for rec in sorted_genelist:
		gc = rec[0]
		for val in rec[1][1]:
			strain = val.split("|")[0]
			gene = val.split("|")[1]
			if 'RS' in gene:
				rgene = gene.split("RS")[1]
			elif "PP" in gene:
				rgene = gene.split("_")[1]
			else:
				rgene = gene 
			if "SHEW" in rgene:
				print rgene
			gene_dir[(strain, rgene.lstrip("0"))]= gc
			#gene_dir[strain, gene]= gc
	return gene_dir

def rename_gc(fit_table, str_fit_table, gene_dir):
	genes = [str(rec) for rec in fit_table.index]
	i=0
	z=0
	for gene in genes:
		if str(gene).isdigit():
			gene_a = str(gene).lstrip('0')
			print gene_a
		if "_" in gene:
			gene_a = gene.split("_")[1].lstrip('0')
			gene_a = str(re.sub("[^0-9]", "", gene_a).lstrip('0'))
		else:
			gene_a = str(re.sub("[^0-9]", "", gene.lstrip('0')))
		try:
			gc = gene_dir[(str_fit_table, str(re.sub("[^0-9]", "", gene_a)))]
			genes[i] = gc
			#print "YES"
			#print gc
			z=z+1
		except KeyError:
			poop.append(gene_a)
		i=i+1
			#print "PRB"
			genes[i] = gene_a
			'''
		try:
			gc = gene_dir[(str_fit_table, gene)]
			genes[i] = gc
		except KeyError:
			print "NO"
			'''
		i=i+1
	fit_table.index = genes
	return fit_table


def remove_duplicated(pd_dataframe):
	#take genes that are duplicated (index and returns only first)
	dups = pd_dataframe.index.duplicated()#(keep="first")
	non_dups = []
	for rec in dups:
		if rec == True:
			non_dups.append(False)
		if rec == False:
			non_dups.append(True)
	keep=pd_dataframe.ix[non_dups]
	return keep

def cond_map(cond, cond_dict, my_strains):
	words=cond
	combined_cond = []
	for cond_dict in all_conds:
		combined_cond.append(matching_conditions(words, cond_dict))
	non_duplicated= remove_duplicated(my_strains[0])
	keep_mat=pd.DataFrame(index = non_duplicated.index, columns = combined_cond)
	for strain in my_strains:
		for cond in strain.columns: 
			if cond in combined_cond:
				non_duplicated= remove_duplicated(strain)
				non_dup2 = [ind for ind in non_duplicated.index if ind in keep_mat.index]
				keep_mat[cond][non_dup2] = non_duplicated[cond][non_dup2]
	return keep_mat




#load SNP datasets
NC_002947=pd.read_table("NC_002947_fit_organism_Putida.txt", index_col=2) #worked
NC_004347=pd.read_table("NC_004347_fit_organism_MR1.txt", index_col=3) #doesn't work
NC_008577=pd.read_table("NC_008577_fit_organism_ANA3.txt", index_col=3) #mismatch between labeling
NC_008700=pd.read_table("NC_008700_fit_organism_SB2B.txt", index_col=3) #no
NC_009092=pd.read_table("NC_017506_fit_organism_Marino.txt", index_col=3) #not really
NC_019936=pd.read_table("NC_019936_fit_organism_psRCH2.txt", index_col=3) #no
NZ_CP007637=pd.read_table("NZ_CP007637_fit_organism_WCS417.txt", index_col=3) #yes
NZ_CP009273=pd.read_table("NZ_CP009273_fit_organism_Keio.txt", index_col=3) #no
NZ_CP012830=pd.read_table("NZ_CP012830_fit_organism_pseudo3_N2E3.txt", index_col=3) #yes
NZ_CP015225=pd.read_table("NZ_CP015225_fit_organism_pseudo6_N2E2.txt", index_col=2) #yes
NZ_JIBD01000001=pd.read_table("NZ_JIBD01000001_fit_organism_Dyella79.txt", index_col=2) #no
NZ_KB893650=pd.read_table("NZ_KB893650_fit_organism_Kang.txt", index_col=2) #no
NZ_LKBJ01000002=pd.read_table("NZ_LKBJ01000002_fit_organism_pseudo13_GW456_L13.txt", index_col=2) #no


all_strains=[NC_002947,NC_004347,NC_008577,NC_008700,NC_009092,NC_019936,NZ_CP007637,NZ_CP009273,NZ_CP012830,NZ_CP015225,NZ_JIBD01000001,NZ_KB893650,NZ_LKBJ01000002]
all_strains_str=['NC_002947','NC_004347','NC_008577','NC_008700','NC_009092','NC_019936','NZ_CP007637','NZ_CP009273','NZ_CP012830','NZ_CP015225','NZ_JIBD01000001','NZ_KB893650','NZ_LKBJ01000002']

#rename to geneclusters
gene_dir = rs_genecluster_map(sorted_genelist)
for i in range(0,len(all_strains)):
	all_strains[i] = rename_gc(all_strains[i], all_strains_str[i], gene_dir)

#find shared conditions
#shared=[all_conds[0][x] for x in all_conds[0] if all_conds[0][x] in all_conds[1].values() and in all_conds[2].values() ]

#that doesn't work. Let's choose a few shared conditions:
cond1=["D-Glucose", "D-Glucose (C)"]
cond2 = ['D-Fructose']
cond3 = ['L-Glutamine']
cond4 = ['Tetracycline'] #take the highest value of Tetracycline
cond5 = ['Chloramphenicol'] #take the highest values of Chloramphenicol

all_conds = []
joint_values = red_conditions(all_strains[0]).values()
all_conds.append(red_conditions(all_strains[0]))



for rec in sorted_genelist:
    for val in rec[1][1]:
        if val.split("|")[0]=='NC_008577':
            gene = val.split("|")[1]
        if "RS" not in gene:
        	print "PROBLEM"
        #if "RS" in gene:
        #	print gene.split("RS")[1]
        



'''
for rec in [all_strains[6], all_strains[8], all_strains[9]]:
	all_conds.append(red_conditions(rec))
	joint_values = set(joint_values) & set(red_conditions(rec).values())
'''

my_strains = [all_strains[0], all_strains[6], all_strains[8]]
Glucose = cond_map(cond1, cond_dict, my_strains).astype('float64')
Fructose = cond_map(cond2, cond_dict, my_strains).astype('float64')
Glutamine = cond_map(cond3, cond_dict, my_strains).astype('float64')
Chlor = cond_map(cond5, cond_dict, my_strains).astype('float64')
Tetracycline = cond_map(cond4, cond_dict, my_strains).astype('float64')

Glucose.to_csv("Glucose_mapping.csv")
Tetracycline.to_csv("Tetracycline_mapping.csv")
Fructose.to_csv()"Tetracycline_mapping.csv"

'''
#output co-map of desired genomes, output only one gene from ortholog group
genus2="NC_020063"
genus1="NZ_CP012830"
keep=[]
for line in sorted_genelist:
	isin=[s for s in line[1][1] if genus1 in s or genus2 in s]
	if len(isin)>0:
		print isin
		isin1=[s for s in isin if genus1 in s]
		isin2=[s for s in isin if genus2 in s]
		if len(isin2)>0 and len(isin1)>0:
			keep.append([isin1[0], isin2[0]])


#between the Pfl and K12 1421 overlapping ortholog groups (very small!!)
pickle.dump(keep, open("NC_020063_NZ_CP011567_ortholog_map.cpk", "w"))



#Now load the gene data for the two genera
genus1='/ebio/abt6_projects9/tnseq/data/combined_datasets/Pseudmonadales/FW300_N1B4/'
genus2='/ebio/abt6_projects9/tnseq/data/combined_datasets/Enterobacteriales/Ecoli_BW25113/'


g1_map=pickle.load(open(genus1+"tnseq_panx_mapped_dict.cpk"))
g2_map=pickle.load(open(genus2+"tnseq_panx_mapped_dict.cpk"))
together=pickle.load(open("NC_020063_NZ_CP011567_ortholog_map.cpk"))
g1_redone=map_redone(g1_map)
g2_redone=map_redone(g2_map)
together_redone=[[g1_redone[line[1]], g2_redone[line[0]]] for line in together]


#Build table with gene information across genera
g1_fit=pickle.load(open(genus1+"tnseq_panx_fitness_data.cpk"))
g2_fit=pickle.load(open(genus2+"tnseq_panx_fitness_data.cpk"))
'''


