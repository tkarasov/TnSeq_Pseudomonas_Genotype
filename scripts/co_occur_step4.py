#!/usr/bin/env python

'''the goal of this script is to look at the co-occurence of genes and ask whether this predicts similar functions'''
import os,sys,copy;import numpy as np
from collections import defaultdict
sys.path.append('./')
from Bio import Phylo, AlignIO
import pickle
import seaborn as sns
import pandas as pd
from joblib import Parallel, delayed
from Bio import SeqIO


def distance_genes(gene1, gene2):
	if gene1.scaffoldId == gene2.scaffoldId:
		b1 = (float(gene1.begin) + float(gene1.end))/2
		b2 = (float(gene2.begin) + float(gene2.end))/2
		s1 = float(gene1.begin) - float(gene1.end)
		s2 = float(gene1.begin) - float(gene1.end)
		if gene1.strand == gene2.strand:
			strand = True
		else:
			strand = False
		if s1 > 0 and s2 >0:
			direc = True
		elif s1 < 0 and s2 < 0:
			direc = True
		else:
			direc = False
		distance = abs(b1 - b2)
		return [distance, strand, direc]

	else:
		return ["NA", "NA", "NA"]


def co_occur(pa_cor_keep, fit_keep, thresh):
	'''calculates the fitness correlation between genes that show a given threshold of p/a polymorphism'''
	co_ocur_cor = {}
	for rec in pa_cor_keep.columns:
		fit_keep_av = np.mean(fit_keep[rec])
		enough = [val for val in pa_cor_keep.index if rec!=val and pa_cor_keep[rec][val]> thresh]
		fit_keep_spec = np.mean(fit_keep[rec].ix[[val for val in enough]])
		co_ocur_cor[rec]= [fit_keep_av, fit_keep_spec]
		print rec

	return co_ocur_cor

def empirical_dist(pd_series, thresh):
	pd_sort=pd_series.sort_values(ascending=False)
	num_include=int(thresh*len(pd_series))
	pd_keep=pd_sort[0:num_include]
	return pd_keep

def add_attributes(fit, together):
	'''need to add on other attributes of genes in consideration'''
	together['gene_div'] = (np.array(fit['gene_div'][together.index].astype('float'))+np.array(fit['gene_div'][together['gene']].astype('float')))/2
	together['gene_ev'] = (np.array(fit['gene_ev'][together.index].astype('float'))+np.array(fit['gene_ev'][together['gene']].astype('float')))/2
	together['ttree'] = (np.array(fit['ttree'][together.index].astype('float'))+np.array(fit['ttree'][together['gene']].astype('float')))/2
	return together
fit = pickle.load(open("tnseq_panx_fitness_data.cpk"))
fit_reduced = fit[fit.columns[3:142]]
fit_red_fin = empirical_dist(fit_reduced.mean(axis=1), 0.2)
fit_reduced1 = fit_reduced.ix[fit_red_fin.index]
#choose tails of variance


#correlation between genes
fit_cor = fit_reduced1.T.corr()

#distance between genes
location = pd.DataFrame([line.strip().split('\t') for line in open("organism_pseudo1_N1B4_genes.tab").readlines()])
location.columns = location.iloc[0]
location.index = location[location.columns[0]]
loc_new = location.drop(location.index[0])
loc_new2 = loc_new.drop(loc_new.columns[0], axis = 1)


#now map co-occurence
path_to_pangenome_dir='/ebio/ag-neher/share/users/wding/panX-refseq/data/Pseudomonadales'#sys.argv[1]
fasta=list(SeqIO.parse(open(path_to_pangenome_dir+"/geneCluster/genePresence.aln"), 'fasta'))
pa_df = pd.DataFrame(columns = range(0, len(fasta[0].seq)), index = [rec.name for rec in fasta])
for line in fasta:
	pa_df.ix[line.name] = [val for val in line.seq]


panx_gene_identified = pickle.load(open("tnseq_panx_mapped_dict.cpk")) # gene number to name mapping!!!
num_name_map = {}
for name in panx_gene_identified.keys():
	num = panx_gene_identified[name][2]
	new = name.split(':')[1]
	num_name_map[num] = new

#there are about 900 genes that map to the same ortholog cluster...we must get rid of them...but also ask whether they have the same behavior


tnseq_loci = list(set([line[2] for line in panx_gene_identified.values()]))
pa_tnseq = pa_df[tnseq_loci]
pa_tnseq.columns = [num_name_map[rec] for rec in pa_tnseq.columns]
pa_tnseq = pa_tnseq.convert_objects(convert_numeric = True)

#pa_cor is the correlation matrix for abundances 
pa_cor = pa_tnseq.corr()


#compare fit_cor and pa_cor
pa_cor_keep = pa_cor[[rec for rec in fit_cor.index if rec in pa_cor.index]].ix[[rec for rec in fit_cor.index if rec in pa_cor.index]]
pa_cor_keep = pa_cor_keep[pa_cor_keep.dropna(how='all').index].ix[pa_cor_keep.dropna(how='all').index]
fit_keep = fit_cor[pa_cor_keep.columns].ix[pa_cor_keep.index]


#now map distance between every pair of genes of subset
'''pairwise_dist = pd.DataFrame(columns = pa_cor_keep.index, index = pa_cor_keep.index)

#this loop takes hours. Figure out how to make faster
for rec in pa_cor_keep.index:
    for val in pa_cor_keep.index:
        pairwise_dist[rec][val] = distance_genes(loc_new2.ix[rec], loc_new2.ix[val])
        print rec, val

with open("tnseq_pairwise_dist.cpk", 'wb') as handle:
    pickle.dump(pairwise_dist, handle)

#all_comb = [[loc_new2.ix[rec], loc_new2.ix[val]] for rec in loc_new2.index for val in loc_new2.index ]

#results = Parallel(n_jobs=-1, verbose=verbosity_level, backend="threading")(map(delayed(distance_gene), all_comb))
'''
pairwise_dist = pickle.load(open("tnseq_pairwise_dist.cpk"))
pairwise = pairwise_dist[fit_keep.columns].ix[fit_keep.columns]

#now for linear modeling of relationships


pa_melt = pd.melt(pa_cor_keep)
fit_melt = pd.melt(fit_keep)
pairwise_melt = pd.melt(pairwise)



together = pd.DataFrame({'gene':pairwise_melt['variable'], 'dist':[line[0] for line in pairwise_melt['value']], 'fit':fit_melt['value'], 'pa':pa_melt['value']})
together.index=list(fit_keep.index)*len(fit_keep.index)

handle=open("pan_genome_general_stats.txt", 'w')
together.to_csv(handle, sep=",")
handle.close()

x = numpy.array(fit['ttree'])
y = numpy.array(fit['set5IT028 Minimal media pH9'])
#plt.scatter(x, y)
#ys = lowess(x, y)
m9 = sns.regplot(x, y, fit_reg=False, scatter_kws={"color": "black"})) #, x_bins = 10, x_ci = 'sd')
m9.set_xlim(0,)
m9_fig = m9.get_figure()
m9_fig.savefig("m9_fitness.pdf")

