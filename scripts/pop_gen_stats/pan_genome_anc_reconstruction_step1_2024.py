#!/usr/bin/python3

'''the goal of this script is to do ancestral state reconstruction of the loss and gain of genes along the provided phylogenetic tree. TaliaK repurposed script in 6/2024'''
import os
import sys
import copy
import numpy as np
from collections import defaultdict
sys.path.append('./')
from treetime import *
#from treetime.treetime import treeanc as ta
#import treetime.treeanc as ta
from treetime.gtr import GTR
#from treetime import io: didn't work on 3.27.2019 but didn't seem necesseary
from treetime import seq_utils
from Bio import Phylo, AlignIO
import pickle
#from sf_miscellaneous import write_json, write_pickle
#from sf_geneCluster_align_makeTree import load_sorted_clusters


# path_to_pangenome_dir='/ebio/ag-neher/share/users/wding/panX-refseq/data/Pseudomonadales'#sys.argv[1]
path_to_pangenome_dir = '/Users/talia/Documents/GitHub/TnSeq_Pseudomonas_Genotype/input_data/panX_data_the50/MySpecies_50_subset'
# path_to_pangenome_dir='/ebio/ag-neher/share/users/wding/panX-refseq/data/Enterobacteriales/'
#t = path_to_pangenome_dir + '/vis/strain_tree.nwk'


def infer_gene_gain_loss(path, rates=[1.0, 1.0]):
    '''code nabbed and edited from panX'''

    # initialize GTR model with default parameters
    mu = np.sum(rates)
    gene_pi = np.array(rates) / mu
    gain_loss_model = GTR.custom(pi=gene_pi, mu=mu,
                                 W=np.ones((2, 2)),
                                 alphabet=np.array(['0', '1']))
    # add "unknown" state to profile
    gain_loss_model.profile_map['-'] = np.ones(2)
    #root_dir = os.path.dirname(os.path.realpath(__file__))

    # define file names for pseudo alignment of presence/absence patterns as in 001001010110
    # path_to_pangenome_dir='/ebio/ag-neher/share/users/wding/panX-refseq/data/Enterobacteriales/'
    # path_to_pangenome_dir='/ebio/ag-neher/share/users/wding/panX-refseq/data/Pseudomonadales'#sys.argv[1]
    nwk = path_to_pangenome_dir + "/vis/strain_tree.nwk"
    fasta = path_to_pangenome_dir + "/geneClusterx/genePresence.aln"

    # instantiate treetime with custom GTR
    t = ta.TreeAnc(nwk, gtr=gain_loss_model, verbose=2)
    # fix leaves names since Bio.Phylo interprets numeric leaf names as confidence
    for leaf in t.tree.get_terminals():
        if leaf.name is None:
            leaf.name = str(leaf.confidence)
    t.aln = fasta
    t.tree.root.branch_length = 0.0001
    t.reconstruct_anc(method='ml')
    for n in t.tree.find_clades():
        n.genepresence = n.sequence

    return t


def num_events(t):
    mutation_dict = {}
    for n in t.tree.find_clades():
        mut = n.mutations
        for mutation in mut:
            try:
                mutation_dict[mutation[1]] += 1
            except KeyError:
                mutation_dict[mutation[1]] = 1
    handle = open("/Users/talia/Documents/GitHub/TnSeq_Pseudomonas_Genotype/output_data/pan_genome/pres_abs_events_gene.cpk", "wb")
    pickle.dump(mutation_dict, handle)
    handle.close()


def get_parent(tree, child_clade):
    node_path = tree.get_path(child_clade)
    return node_path[-2]


def calc_branch_length_in_tree(t):
    '''takes output from running ancestral reconstruction in treetime, iterates through tree and if a gene present in node is present in parent, adds sum of branch length. Basically measure the amount of time (units in branch length) that gene has spent in the tree'''
    # do a breadth first search. Starting at the MRCA of tree, going to each child, store the branch length of parent to that child. Then go through every gene,
    # initialize
    parent = t.tree.root.sequence
    tracking_gene = {}
    len_dict = {}
    for i in range(0, len(t.tree.root.sequence)):
        tracking_gene[i] = t.tree.root.sequence[i]
        len_dict[i] = 0

    for n in t.tree.find_clades():
        try:
            p = get_parent(t.tree, n)
            pn = get_parent(t.tree, n).sequence
        except IndexError:
            continue
        for i in range(0, len(t.tree.root.sequence)):
            tracking_gene[i] = p.sequence[i]
        branch = n.branch_length
        events = n.mutations
        # change in tracking_gene only those locations that have a putative mutation
        for mutation in events:
            print(tracking_gene[mutation[1]], mutation[2])
            tracking_gene[mutation[1]] = mutation[2]
        for gene in tracking_gene.keys():
            if tracking_gene[gene] == '1':
                len_dict[gene] = len_dict[gene] + branch

        # dump out branch lengths
        handle = open("/Users/talia/Documents/GitHub/TnSeq_Pseudomonas_Genotype/output_data/pan_genome/", "wb")
        pickle.dump(len_dict, handle)
        handle.close()


t = infer_gene_gain_loss(path_to_pangenome_dir, rates=[1.0, 1.0])
calc_branch_length_in_tree(t)  # (path_to_pangenome_dir)
num_events(t)
