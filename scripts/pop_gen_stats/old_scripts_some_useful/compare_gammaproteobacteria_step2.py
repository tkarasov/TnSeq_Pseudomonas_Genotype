
'''the goal of this script is to take the gammaproteobacteria output from panx
and to look at the relationship between divergence and the correlation of the
phenotype. Should be used to process the output of script that processes blast output or actually for the panX output overall'''
import json
import os
import pandas as pd
from Bio import SeqIO
import pickle
import numpy as np
from scipy import stats
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt


def prepare_tag_file():
    # read in the files of orthology assignment
    with open("/ebio/abt6_projects9/tnseq/wdata/panX-tnseq/data/Gammaproteobacteria/vis/geneCluster.json") as f:
        annot = json.load(f)

    # get all homologs associated with a gene
    annot_new = []

    for gene in annot:
        stuff = [str(gene['geneId']), gene['msa'], gene['GName']]
        stuff.extend(gene['locus'].split(" "))
        annot_new.append(stuff)

    # write out mappings between geneID and all the corresponding loci
    with open("/ebio/abt6_projects9/tnseq/data/pan_genome_datasets/gammaproteobacteria_mappings.txt", 'w') as myfile:
        for item in annot_new:
            myfile.write((', ').join((item)) + "\n")

    # now build a dictionary with every gene tag and its corresponding mapping
    tag_cluster = {}
    for line in annot_new:
        for rec in line[3:]:
            tag_cluster["-".join(rec.split("_", 2)[2:3])] = line[1]

    # build a dictionary for all locus_tag to set mapping. This takes more than an hour. I originally had a problem because the set names were used more than once
    file_names = os.listdir("/ebio/abt6_projects9/tnseq/data/fitness_datasets/fitness_tables/")
    full_tag_dict = {}
    for file in file_names:
        if "fit" in file and "amended" in file:
            full_tag_dict = process_tag_file(file, full_tag_dict, tag_cluster)
    #(pd.DataFrame.from_dict(data=full_tag_dict, orient='index').to_csv('/ebio/abt6_projects9/tnseq/data/fitness_datasets/fitness_tables/full_tag_dict_file.csv', header=False))
    with open('/ebio/abt6_projects9/tnseq/tnseq_function/data/full_tag_dict_file_273_2019.cpk', "wb") as file:
        file.write(pickle.dumps(full_tag_dict))


def process_tag_file(file, tag_dict, tag_cluster):
    # this function processes the tag_file to add to the dictionary
    file_in = pd.read_csv(file)
    file_in.index = file_in['locusID']
    organism = list(file_in['orgID'])[0]
    for locus in file_in['locusID']:
        for condition in file_in.columns:
            if 'set' in condition:
                try:
                    if type(file_in[condition][locus]) == numpy.float64:
                        # there is a problem here. Some loci have pandas series while others floats. I mustgo back to determine why
                        # tag_dict[(condition.split(" ")[0], tag_cluster[locus])] = file_in[condition][locus] ***this was changed 3/26/2019. I think this is what caused bad overlap
                        new_cond = condition.split()[0]
                        tag_dict[(new_cond, organism, tag_cluster[locus])] = file_in[condition][locus]
                        # issue with N515DRAFT in addition ot several others
                        print(locus)
                        print(1)
                except KeyError:
                    print(locus)
                    print(0)
                    pass
    return tag_dict


def choose_condition(conditions_table, condition, full_tag_dict):
    focal_condition = [rec for rec in conditions_table if condition in rec[-1]]

    # keep only gammaproteobacteria under consideration
    keep_condition = [rec for rec in focal_condition if rec[0] in [line[1] for line in genome_id]]

    # now build dataframe with set as index and locus as column
    full_tag_subset = {}
    for k, v in full_tag_dict.items():
        if [k[1], k[0]] in [line[0:2] for line in keep_condition]:
            full_tag_subset[k] = v
    full_tag_subset = {k: v for k, v in full_tag_dict.items() if [k[1], k[0]] in [line[0:2] for line in keep_condition]}
    condition_df = pd.DataFrame(columns=set([line[2] for line in full_tag_subset.keys()]), index=set([line[0] + "_" + line[1] for line in full_tag_subset.keys()]))
    for k, v, x in full_tag_subset.keys():
        print(k, v, x)
        condition_df[x][k + v] = full_tag_subset[(k, v, x)].astype(float)
    return condition_df.T.astype(float)


def pd_fill_diagonal(df_matrix, value=0):
    mat = df_matrix.values
    n = mat.shape[0]
    mat[range(n), range(n)] = value
    return pd.DataFrame(mat)


def calculate_divergence_genomes():
    # this takes the SNP matrix from panx and calculates a divergence matrix
    align = list(SeqIO.parse("/ebio/abt6_projects9/tnseq/data/fitness_datasets/fitness_tables/SNP_whole_matrix.aln", "fasta"))
    pairwise_divergence = pd.DataFrame(columns=[rec.name for rec in align], index=[rec.name for rec in align])
    for s1 in align:
        for s2 in align:
            n1 = s1.name
            n2 = s2.name
            count = sum(1 for a, b in zip(s1.seq, s2.seq) if a != b)
            dist = float(count) / len(s1.seq)
            pairwise_divergence[n1][n2] = dist
    return pairwise_divergence


def pairwise_condition(pairwise_divergence, condition_corr):
    genomes = [conditions_dict[rec][0] for rec in condition_corr.corr()]  # [conditions_dict[rec][0].strip(".gbk") for rec in condition_corr.corr()]
    present = map(lambda x: x in pairwise_divergence.columns, genomes)
    keep_conditions = list(filter(lambda x: conditions_dict[x][0].strip(".gbk") in pairwise_divergence.columns, condition_corr.corr()))
    keep_mat = condition_corr.corr()[keep_conditions].loc[keep_conditions]
    keep_pairwise = pairwise_divergence[[conditions_dict[rec][0].strip(".gbk") for rec in keep_mat]].loc[[conditions_dict[rec][0].strip(".gbk") for rec in keep_mat]]
    # assign NA to diagonals
    keep_mat = pd_fill_diagonal(keep_mat, numpy.nan)
    keep_pairwise = pd_fill_diagonal(keep_pairwise, numpy.nan)
    return(keep_pairwise, keep_mat)


# Assuming we have previously run prepare_tag_file(), we just need to read in the output otherwise: prepare_tag_file()
full_tag_dict = pickle.load(open('/ebio/abt6_projects9/tnseq/tnseq_function/data/full_tag_dict_file_273_2019.cpk', 'rb'))

# generate pandas dataframe of the tag_dict. This takes a long time!
#my_index = set([(k[0], k[1]) for k in list(full_tag_dict.keys())])
#my_col = set([k[2] for k in list(full_tag_dict.keys())])
#full_tag_pd = pd.DataFrame(index=my_index, columns=my_col)
# for rec in list(full_tag_dict.items()):
#    ind = rec[0][0:2]
#    col = rec[0][2]
#   full_tag_pd[col].loc[ind] = rec[1]

# full_tag_pd.to_csv("/ebio/abt6_projects9/tnseq/tnseq_function/data/full_tag_pd.csv")

#'/ebio/abt6_projects9/tnseq/data/fitness_datasets/fitness_tables/full_tag_dict_file.cpk', 'rb'))

genome_id = [line.strip().split('\t') for line in open("/ebio/abt6_projects9/tnseq/tnseq_function/data/genome_list.txt").readlines()]
conditions = [line.strip().split('\t') for line in open("/ebio/abt6_projects9/tnseq/tnseq_function/fitness_tables/all_conditions_amended.txt")]
conditions_dict = {}
id_dict = {line[1]: line[0] for line in genome_id}
for line in conditions:
    if len(line) > 1:
        try:
            print(line[1] + "_" + id_dict[line[0]])
            conditions_dict[line[1] + "_" + line[0]] = [id_dict[line[0]], line[-1]]
        except KeyError:
            conditions_dict[line[1]] = [line[0], line[-1]]


tetc = choose_condition(conditions, "Tetracycline", full_tag_dict)
spec = choose_condition(conditions, "Spectinomycin", full_tag_dict)
gentamicin = choose_condition(conditions, "Gentamicin", full_tag_dict)
nalidixic = choose_condition(conditions, "Nalidixic", full_tag_dict)
Doxycycline = choose_condition(conditions, "Doxycycline", full_tag_dict)
Bacitracin = choose_condition(conditions, "Bacitracin", full_tag_dict)


pairwise_divergence = calculate_divergence_genomes()

fig = plt.figure(figsize=(16, 12))
#f, (ax1, ax2, ax3) = plt.subplots(3, 1, sharex='col')
ax1 = plt.subplot(231)
x, y = pairwise_condition(pairwise_divergence, tetc)
ax1.scatter(x, y)
x_flat = np.array(x).flatten().astype('float32')
y_flat = np.array(y).flatten().astype('float32')
mask = ~np.isnan(x_flat) & ~np.isnan(y_flat)
slope, intercept, r_value, p_value, std_err = stats.linregress(x_flat[mask], y_flat[mask])
line = slope * x_flat[mask] + intercept
ax1.plot(x_flat[mask], line, 'r', label='y={:.2f}x+{:.2f}'.format(slope, intercept))
ax1.set_xlabel("% Sequence Divergence")
ax1.set_ylabel("Pearson Correlation Coefficient (R)")
ax1.set_title("Tetracycline")


plt.subplot(232, sharex=ax1, sharey=ax1)
x, y = pairwise_condition(pairwise_divergence, gentamicin)
plt.scatter(x, y)
x_flat = np.array(x).flatten().astype('float32')
y_flat = np.array(y).flatten().astype('float32')
mask = ~np.isnan(x_flat) & ~np.isnan(y_flat)
slope, intercept, r_value, p_value, std_err = stats.linregress(x_flat[mask], y_flat[mask])
line = slope * x_flat[mask] + intercept
plt.plot(x_flat[mask], line, 'r', label='y={:.2f}x+{:.2f}'.format(slope, intercept))
#ax2.xlabel("% Sequence Divergence")
#ax2.ylabel("Pearson Correlation Coefficient (R)")
plt.title("Gentamicin")

plt.subplot(2, 3, 3, sharex=ax1, sharey=ax1)
fin_pair, spec_pair = pairwise_condition(pairwise_divergence, spec)
plt.scatter(fin_pair, spec_pair)
x_flat = np.array(x).flatten().astype('float32')
y_flat = np.array(y).flatten().astype('float32')
mask = ~np.isnan(x_flat) & ~np.isnan(y_flat)
slope, intercept, r_value, p_value, std_err = stats.linregress(x_flat[mask], y_flat[mask])
line = slope * x_flat[mask] + intercept
plt.plot(x_flat[mask], line, 'r', label='y={:.2f}x+{:.2f}'.format(slope, intercept))
plt.title("Spectinomycin")

plt.subplot(2, 3, 4)
x, y = pairwise_condition(pairwise_divergence, Doxycycline)
plt.scatter(x, y)
x_flat = np.array(x).flatten().astype('float32')
y_flat = np.array(y).flatten().astype('float32')
mask = ~np.isnan(x_flat) & ~np.isnan(y_flat)
slope, intercept, r_value, p_value, std_err = stats.linregress(x_flat[mask], y_flat[mask])
line = slope * x_flat[mask] + intercept
plt.plot(x_flat[mask], line, 'r', label='y={:.2f}x+{:.2f}'.format(slope, intercept))
plt.title("Doxycycline")

plt.subplot(2, 3, 5)
x, y = pairwise_condition(pairwise_divergence, nalidixic)
plt.scatter(x, y)
x_flat = np.array(x).flatten().astype('float32')
y_flat = np.array(y).flatten().astype('float32')
mask = ~np.isnan(x_flat) & ~np.isnan(y_flat)
slope, intercept, r_value, p_value, std_err = stats.linregress(x_flat[mask], y_flat[mask])
line = slope * x_flat[mask] + intercept
plt.plot(x_flat[mask], line, 'r', label='y={:.2f}x+{:.2f}'.format(slope, intercept))
plt.title("Nalidixic acid")

plt.subplot(2, 3, 6)
x, y = pairwise_condition(pairwise_divergence, Bacitracin)
plt.scatter(x, y)
x_flat = np.array(x).flatten().astype('float32')
y_flat = np.array(y).flatten().astype('float32')
mask = ~np.isnan(x_flat) & ~np.isnan(y_flat)
slope, intercept, r_value, p_value, std_err = stats.linregress(x_flat[mask], y_flat[mask])
line = slope * x_flat[mask] + intercept
plt.plot(x_flat[mask], line, 'r', label='y={:.2f}x+{:.2f}'.format(slope, intercept))
plt.title("Bacitracin")

plt.tight_layout()
plt.savefig("antibiotic.pdf")

Lserine = choose_condition(conditions, "L-Serine", full_tag_dict)
DGlucose = choose_condition(conditions, "D-Glucose", full_tag_dict)
Sucrose = choose_condition(conditions, "Sucrose", full_tag_dict)
LLeucine = choose_condition(conditions, "L-Leucine", full_tag_dict)
DAlanine = choose_condition(conditions, "D-Alanine", full_tag_dict)
LAlanine = choose_condition(conditions, "L-Alanine", full_tag_dict)

fig = plt.figure(figsize=(16, 12))
ax1 = plt.subplot(231)
x, y = pairwise_condition(pairwise_divergence, Lserine)
ax1.scatter(x, y)
ax1.set_xlabel("% Sequence Divergence")
ax1.set_ylabel("Pearson Correlation Coefficient (R)")
x_flat = np.array(x).flatten().astype('float32')
y_flat = np.array(y).flatten().astype('float32')
mask = ~np.isnan(x_flat) & ~np.isnan(y_flat)
slope, intercept, r_value, p_value, std_err = stats.linregress(x_flat[mask], y_flat[mask])
line = slope * x_flat[mask] + intercept
ax1.plot(x_flat[mask], line, 'r', label='y={:.2f}x+{:.2f}'.format(slope, intercept))
ax1.set_title("L-serine")

plt.subplot(232, sharex=ax1, sharey=ax1)
x, y = pairwise_condition(pairwise_divergence, DGlucose)
plt.scatter(x, y)
x_flat = np.array(x).flatten().astype('float32')
y_flat = np.array(y).flatten().astype('float32')
mask = ~np.isnan(x_flat) & ~np.isnan(y_flat)
slope, intercept, r_value, p_value, std_err = stats.linregress(x_flat[mask], y_flat[mask])
line = slope * x_flat[mask] + intercept
plt.plot(x_flat[mask], line, 'r', label='y={:.2f}x+{:.2f}'.format(slope, intercept))
plt.title("D-Glucose")

plt.subplot(2, 3, 3, sharex=ax1, sharey=ax1)
x, y = pairwise_condition(pairwise_divergence, Sucrose)
plt.scatter(x, y)
x_flat = np.array(x).flatten().astype('float32')
y_flat = np.array(y).flatten().astype('float32')
mask = ~np.isnan(x_flat) & ~np.isnan(y_flat)
slope, intercept, r_value, p_value, std_err = stats.linregress(x_flat[mask], y_flat[mask])
line = slope * x_flat[mask] + intercept
plt.plot(x_flat[mask], line, 'r', label='y={:.2f}x+{:.2f}'.format(slope, intercept))
plt.title("Sucrose")

plt.subplot(2, 3, 4, sharex=ax1, sharey=ax1)
x, y = pairwise_condition(pairwise_divergence, LLeucine)
plt.scatter(x, y)
x_flat = np.array(x).flatten().astype('float32')
y_flat = np.array(y).flatten().astype('float32')
mask = ~np.isnan(x_flat) & ~np.isnan(y_flat)
slope, intercept, r_value, p_value, std_err = stats.linregress(x_flat[mask], y_flat[mask])
line = slope * x_flat[mask] + intercept
plt.plot(x_flat[mask], line, 'r', label='y={:.2f}x+{:.2f}'.format(slope, intercept))plt.title("L-Leucine")

plt.subplot(2, 3, 5, sharex=ax1, sharey=ax1)
x, y = pairwise_condition(pairwise_divergence, DAlanine)
plt.scatter(x, y)
x_flat = np.array(x).flatten().astype('float32')
y_flat = np.array(y).flatten().astype('float32')
mask = ~np.isnan(x_flat) & ~np.isnan(y_flat)
slope, intercept, r_value, p_value, std_err = stats.linregress(x_flat[mask], y_flat[mask])
line = slope * x_flat[mask] + intercept
plt.plot(x_flat[mask], line, 'r', label='y={:.2f}x+{:.2f}'.format(slope, intercept))
plt.title("D-Alanine")

plt.subplot(2, 3, 6, sharex=ax1, sharey=ax1)
x, y = pairwise_condition(pairwise_divergence, LAlanine)
plt.scatter(x, y)
x_flat = np.array(x).flatten().astype('float32')
y_flat = np.array(y).flatten().astype('float32')
mask = ~np.isnan(x_flat) & ~np.isnan(y_flat)
slope, intercept, r_value, p_value, std_err = stats.linregress(x_flat[mask], y_flat[mask])
line = slope * x_flat[mask] + intercept
plt.plot(x_flat[mask], line, 'r', label='y={:.2f}x+{:.2f}'.format(slope, intercept))
plt.title("L-Alanine")

plt.tight_layout()
plt.savefig("/ebio/abt6_projects9/tnseq/tnseq_function/figs/nutrient.pdf")
