#!/usr/bin/env python
'''the goal of this script is to take the pan genome for a given bacterial group, and to assign diversity statistics to each ortholog group'''

import sys
import os
import pickle
import ete3
import pandas as pd
from Bio import SeqIO

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
path_to_pangenome_dir='/Users/talia/Library/CloudStorage/GoogleDrive-tkarasov@gmail.com/My Drive/Utah_Professorship/projects/Tnseq/panX_data_the50/MySpecies_50_subset/'
clusters = path_to_pangenome_dir + "/allclusters_final.tsv"
# tnseq_genome = list(SeqIO.parse(os.popen("ls | grep _cds_from_genomic.fna").read().strip(), "fasta"))

# tnseq_panx = pickle.load(open("tnseq_panx_mapped_dict.cpk", 'rb')) '''I think this is the deprecated mapping'''
gene_mapping = [line.strip().split(", ") for line in open("/ebio/abt6_projects9/tnseq/tnseq_function/fitness_tables/gammaproteobacteria_mappings.txt").readlines()]

# note that the gene_mapping file starts at 1 whereas the panX files start at 0. I am subtracting 1 but must check that this mapping is correct.
gene_map_dict = {(int(line[0]) - 1): line[1] for line in gene_mapping}

# Nucleotide diversity of gene
genediv = pickle.load(open(path_to_pangenome_dir + "/geneClusterx/gene_diversity.cpk", 'rb'))

# GC content

# Deletion/Duplication events pickle.load(open(path_to_pangenome_dir+"/geneCluster/dt_geneEvents.cpk"))
uncode_gene_events = pickle.load(open(path_to_pangenome_dir + "/geneClusterx/dt_geneEvents.cpk", 'rb'))
gene_events = {gene_map_dict[k]: v for k, v in uncode_gene_events.items()}

# time_spent in tree
uncode_ts_tree = pickle.load(open("/Users/talia/Documents/GitHub/TnSeq_Pseudomonas_Genotype/output_data/pan_genome/branch_gene.cpk", 'rb'))
ts_tree = {gene_map_dict[k]: v for k, v in uncode_ts_tree.items()}

tnseq_stats = {}
for gene in gene_map_dict.values():
    tnid = gene
    pgid = 'NULL'
    pggc = 'NULL'
    pgnum = 'NULL'
    pident = gene_events[gene]
    ttree = ts_tree[gene]
    gene_div = float(genediv[gene])
    gene_ev = float(gene_events[gene])
    tnseq_stats[tnid] = tnseq_gene(tnid, pgid, pggc, pgnum, ttree, gene_ev, gene_div)

# Now attach the gene object information to the large fitness output file. Everything in one place!

full_tag_pd = pd.read_csv('/Users/talia/Documents/GitHub/TnSeq_Pseudomonas_Genotype/output_data/pan_genome/full_tag_pd.csv', sep=",", index_col=1)
not_measured = {}

full_tag_pd = full_tag_pd.append(pandas.Series(name="Time in tree"))
full_tag_pd = full_tag_pd.append(pandas.Series(name="Genetic_Diversity"))
full_tag_pd = full_tag_pd.append(pandas.Series(name="Num_gene_events"))
full_tag_pd = full_tag_pd.append(pandas.Series(name="GC_content"))

# This takes a few hours interactively
for gene in tnseq_stats:
    print(gene)
    try:
        full_tag_pd[gene].loc["Time in tree"] = tnseq_stats[gene].ttree
        full_tag_pd[gene].loc["Genetic_Diversity"] = tnseq_stats[gene].gene_div
        full_tag_pd[gene].loc["Num_gene_events"] = tnseq_stats[gene].gene_ev
        full_tag_pd[gene].loc["GC_content"] = tnseq_stats[gene].pggc
    except KeyError:
        not_measured[gene] = tnseq_stats[gene]


full_tag_pd.to_csv("/Users/talia/Documents/GitHub/TnSeq_Pseudomonas_Genotype/output_data/pan_genome/full_tag_pd_with_fitness.csv")


'''
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

# Fitness info
# fit = pandas.read_table("fit_organism_pseudo1_N1B4.txt", index_col =2)
fit = pandas.read_table("NZ_CP009273_fit_organism_Keio.txt", index_col=2)

'''


'''
fit.index = [str(rec) for rec in fit.index]
fit['ttree'] = [tnseq_stats[str(rec)].ttree for rec in fit.index]
fit['gene_ev'] = [tnseq_stats[rec].gene_ev for rec in fit.index]
fit['gene_div'] = [tnseq_stats[rec].gene_div for rec in fit.index]


with open("/ebio/abt6_projects9/tnseq/tnseq_function/data/Keio_tnseq_panx_fitness_data.cpk", 'wb') as handle:
    pickle.dump(fit, handle)
'''
