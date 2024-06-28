#!/usr/bin/env python
'''the goal of this script is to take the pan genome for a given bacterial group, and to assign diversity statistics to each ortholog group. Looking good with full_tag_pd output file on 6/2/2024'''

import sys
import os
import pickle
import ete3
import pandas as pd
from Bio import *
from Bio import Phylo
#from newick import read

# now build object for every gene in tnseq genome


class tnseq_gene:

    def __init__(self, tnid, pgid, pggc, pg_num, ttree, gene_ev, gene_div):
        self.tnid = tnid
        self.pgid = pgid
        self.pggc = pggc
        self.pg_num = pg_num
        # self.pident = pident
        # self.len_gene = len_gene
        # self.gc_content = gc_content
        # self.gc_bias =
        self.ttree = ttree
        self.gene_ev = gene_ev
        self.gene_div = gene_div


def getGC(seq):
    return(str(int(round((sum([1.0 for nucl in seq if nucl in ['G', 'C']]) / len(seq)) * 100))))


def fit_object(file_name):
    fit = pandas.read_table(file_name, index_col=2)
    # Time spent in phylogenetic tree
    tnseq_stats = {}
    for key, value in tnseq_panx.iteritems():
        try:
            tnid = key.split(":")[1]
        except IndexError:
            print("ERROR")
        tnid = tnseq_panx[key][0].split("_")[2]
        pgid = value[0]
        pggc = value[1]
        pgnum = value[2]
        pident = gene_events[pgnum]
        ttree = ts_tree[pgnum]
        gene_div = genediv[pggc]
        gene_ev = gene_events[pgnum]
        tnseq_stats[tnid] = tnseq_gene(tnid, pgid, pggc, pgnum, ttree, gene_ev, gene_div)
    fit.index = [str(rec) for rec in fit.index]
    fit['ttree'] = [tnseq_stats[str(rec)].ttree for rec in fit.index]
    fit['gene_ev'] = [tnseq_stats[rec].gene_ev for rec in fit.index]
    fit['gene_div'] = [tnseq_stats[rec].gene_div for rec in fit.index]
    with open("/ebio/abt6_projects9/tnseq/tnseq_function/data/" + file_name + "_tnseq_panx_fitness_data.cpk", 'wb') as handle:
        pickle.dump(fit, handle)


# path_to_pangenome_dir = '/ebio/ag-neher/share/users/wding/panX-refseq/data/Pseudomonadales'#sys.argv[1]
path_to_pangenome_dir="/Users/talia/Library/CloudStorage/GoogleDrive-tkarasov@gmail.com/My Drive/Utah_Professorship/projects/Tnseq/panX_data_the50/MySpecies_50_subset/"
clusters = path_to_pangenome_dir + "/allclusters_final.tsv"
# tnseq_genome = list(SeqIO.parse(os.popen("ls | grep _cds_from_genomic.fna").read().strip(), "fasta"))

# tnseq_panx = pickle.load(open("tnseq_panx_mapped_dict.cpk", 'rb')) '''I think this is the deprecated mapping'''
#gene_mapping = [line.strip().split(", ") for line in open("/ebio/abt6_projects9/tnseq/tnseq_function/fitness_tables/gammaproteobacteria_mappings.txt").readlines()]

# note that the gene_mapping file starts at 1 whereas the panX files start at 0. I am subtracting 1 but must check that this mapping is correct.
#gene_map_dict = {(int(line[0]) - 1): line[1] for line in gene_mapping}
#Load gene identities
geneID = pickle.load(open(path_to_pangenome_dir + "/geneID_to_description.cpk", 'rb'))
gene_ID_seq_mapping = pickle.load(open(path_to_pangenome_dir +"geneID_to_geneSeqID.cpk", 'rb'))


#allclusters is the gene orthology mappings
allclusters = [line.strip().split() for line in open(path_to_pangenome_dir + "allclusters_final.tsv", 'r').readlines()]
# I want to keep the listing for p25.c2 (if there is one) for each line
p25_c2 = {}
i=0
for line in allclusters:
    record="NA"
    for rec in line:
        if 'p25.C2' in str(rec):
            record = rec.split("|")[1]
            p25_c2[i]=record
    if record == "NA":
        #p25_c2[i]="NA"
        p25_c2[i]=line[0].split("|")[1]
    i=i+1




# Nucleotide diversity of gene. For some reason this is substantially longer (~3000 genes that is the)
genediv = pickle.load(open(path_to_pangenome_dir + "/geneCluster/gene_diversity.cpk", 'rb'))

# Build Gene Cluster mappings. This is mapping the GC number to the gene name in p25.c2. We need to iterate through the whole folder of vis/geneCluster, look in every file ending in aa_aln.fa.gz  and find the genes in the different strains associated with that gene cluster

GC_dict = {}
# DACKMO:GC
for filename in os.listdir(path_to_pangenome_dir+"/vis/geneCluster/"):
    if filename=="strain_tree.nwk":
        pass
    elif filename.endswith(".nwk"):
        GC=filename.strip(".nwk")
        tree = Phylo.read(path_to_pangenome_dir+"/vis/geneCluster/"+filename, "newick")
        for leaf in tree.get_terminals(): 
            GC_dict[leaf.name.split("|")[1].split("-")[0]] = GC
    else:
        pass

# Now connect the gene name to the divergene
GC_name_divergence = {}
for key,value in GC_dict.items():
    GC_name_divergence[key]=genediv[value]





# GC content -- Not done yet.

# Deletion/Duplication events pickle.load(open(path_to_pangenome_dir+"/geneCluster/dt_geneEvents.cpk"))
uncode_gene_events = pickle.load(open(path_to_pangenome_dir + "/geneCluster/dt_geneEvents.cpk", 'rb'))
#gene_events = {gene_map_dict[k]: v for k, v in uncode_gene_events.items()}

# time_spent in tree


uncode_ts_tree = pickle.load(open("/Users/talia/Documents/GitHub/TnSeq_Pseudomonas_Genotype/output_data/pan_genome/branch_length.cpk", 'rb')) #was branch_gene before
gene_map_dict = {p25_c2[item]:uncode_ts_tree[item] for item in uncode_ts_tree.keys()}
ts_tree = gene_map_dict #{gene_map_dict[k]: v for k, v in uncode_ts_tree.items()}

tnseq_stats = {}
for number,gene in p25_c2.items():
    if gene.startswith("DAKCFMEO_"):
        #right now we are only adding the orthologs that have a p25.c2 gene in them
        print(gene)
        tnid = gene
        pgid = 'NULL'
        pggc = 'NULL'
        pgnum = 'NULL'
        pident = uncode_gene_events[number]
        ttree = ts_tree[gene]
        gene_div = float(GC_name_divergence[gene])
        gene_ev = float(uncode_gene_events[number])
        tnseq_stats[tnid] = tnseq_gene(tnid, pgid, pggc, pgnum, ttree, gene_ev, gene_div)

# Now attach the gene object information to the large fitness output file. Everything in one place!

full_tag_pd = pd.DataFrame(index=[gene for gene in tnseq_stats], columns = ["Time in tree", "Genetic_Diversity", "Num_gene_events", "GC_content"])
#read_csv('/Users/talia/Documents/GitHub/TnSeq_Pseudomonas_Genotype/output_data/pan_genome/full_tag_pd.csv', sep=",", index_col=1)
not_measured = {}


# This takes a few hours interactively
for gene in tnseq_stats:
    print(gene)
    try:
        full_tag_pd["Time in tree"].loc[gene] = tnseq_stats[gene].ttree
        full_tag_pd["Genetic_Diversity"].loc[gene] = tnseq_stats[gene].gene_div
        full_tag_pd["Num_gene_events"].loc[gene] = tnseq_stats[gene].gene_ev
        full_tag_pd["GC_content"].loc[gene] = tnseq_stats[gene].pggc
    except KeyError:
        not_measured[gene] = tnseq_stats[gene]


full_tag_pd.to_csv("/Users/talia/Documents/GitHub/TnSeq_Pseudomonas_Genotype/output_data/pan_genome/full_tag_pd.csv")


