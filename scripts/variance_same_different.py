from matplotlib import pyplot as plt
import numpy as np
import pandas as pd
import pickle
import math
'''The goal of this script is to calculate an amended correlation coefficient for the correlation of gene fitness effects between environments. In this case the denominator is the multiplied covariances of the replicate environments. Only environments that have replicates are considered'''

conditions = pd.read_csv("/ebio/abt6_projects9/tnseq/tnseq_function/fitness_tables/all_conditions_amended.txt", sep="\t")

# this is a dictionary with the set and gene cluster as the key and the fitness measurement as a value
#full_tag_dict = pickle.load(open('/ebio/abt6_projects9/tnseq/data/fitness_datasets/fitness_tables/full_tag_dict_file.cpk', 'rb'))

# this is a pandas dataframe with the set (experiment) as the columns and the gene cluster (output from panx) as the row
full_tag_pd = pd.read_csv('//ebio/abt6_projects9/tnseq/tnseq_function/data/full_tag_pd.csv', sep=",")

full_name = full_tag_pd['Unnamed: 0'] + "_" + full_tag_pd['Unnamed: 1']
full_tag_pd.index = full_name
full_tag_new = full_tag_pd.drop(['Unnamed: 0'], axis=1).drop(['Unnamed: 1'], axis=1).transpose()
# mean centered and divided by standard deviation
full_tag_pd_centered = full_tag_new.apply(lambda x: x - x.mean()) / full_tag_new.std()


def give_shared_cond(conditions):
    '''Find experiments in which the same organism was measured in the same condition more than once'''
    shared = {}
    for strain in set(conditions['Organism']):
        for condition in set(conditions['Condition']):
            subset = []
            subset = conditions.loc[(conditions['Organism'] == strain) & (conditions['Condition'] == condition)]
            if len(subset) > 1:
                shared["_".join((strain + "_" + str(condition)).split(" "))] = [rec + "_" + strain for rec in list(subset.Name)]
                #shared["_".join((strain + "_" + str(condition)).split(" "))] = [(strain, rec) for rec in list(subset.Name)]
                #shared["_".join((strain + "__" + str(condition)).split(" "))] = list(subset.Name)
            else:
                pass
    return shared


def give_shared_descr(conditions):
    '''Find experiments in which the same organism was measured in the same condition more than once'''
    shared = {}
    for strain in set(conditions['Organism']):
        for condition in set(conditions['Condition']):
            subset = []
            subset = conditions.loc[(conditions['Organism'] == strain) & (conditions['Description'] == condition)]
            if len(subset) > 1:
                shared["_".join((strain + "__" + str(condition)).split(" "))] = [(strain, rec) for rec in list(subset.Name)]
                #shared["_".join((strain + "__" + str(condition)).split(" "))] = list(subset.Name)
            else:
                pass
    return shared


def calc_corr_within(shared, full_tag_pd_centered):
    '''Calculate covariance between fitness of same gene in "same condition"'''
    # full covariance matrix
    cov_pd = full_tag_pd_centered.astype(float).cov()
    cov_dict = {}
    group_dict = {}
    # first find relevant subsets
    for key, group in shared.items():
        temp = []
        i = 0
        z = 1
        while z < len(group):
            try:
                temp.append(cov_pd[group[i]].loc[group[z]])
                print("success")
            except KeyError:
                print("Error in group")
                # errror here for ~200 paired conditions
                # print(group[i])
            i = i + 1
            z = z + 1
        if len(temp) > 0:
            cov_dict[key] = np.max(temp)
            for val in group:
                group_dict[val] = np.max(temp)
    group_dict = dict([rec for rec in group_dict.items() if not math.isnan(rec[1])])
    return group_dict


def correct_large_matrix(group_dict, full_tag_pd_centered):
    '''this function takes the max covariance observed within a condition replicate group, and calculates cov(xa, xb)/sqrt(rho_a^2*rho_b^2). I am unhappy with this matrix however, as there are many values that exceed unity. This likely means that there are conditions that have better correlations than the replicate condition itself. What do we do with this?'''
    cov_pd = full_tag_pd_centered.astype(float).cov()
    cov_normalized = pd.DataFrame(columns=[rec for rec in group_dict.keys() if rec in cov_pd.index], index=[rec for rec in group_dict.keys() if rec in cov_pd.index])
    for col in cov_normalized.columns:
        col_cor = "NaN"
        ind_cor = "NaN"
        yes = "no"
        for ind in cov_normalized.index:
            if col in list(group_dict.keys()):
                col_cor = group_dict[col]
                yes = "yes"
            else:
                col_cor = "NaN"
                yes = "no"
                cov_normalized[col].loc[ind] = "NaN"
                continue

            if ind in list(group_dict.keys()):
                ind_cor = group_dict[ind]
                yes = "yes"
                denom = np.sqrt(ind_cor**2 * col_cor**2)
                cov_normalized[col].loc[ind] = cov_pd[col].loc[ind] / denom
            else:
                ind_cor = "NaN"
                yes = "no"
                cov_normalized[col].loc[ind] = "NaN"
                continue
    return cov_normalized



# 1424 conditions have some repetition. But....only about 800 of them seem to be recognized in the creation of group_dict. Need to figure out why this is. Seems something wrong with Dylella and only Pseudomonas showing up.
shared = give_shared_cond(conditions)
shared_desc = give_shared_descr(conditions)
group_dict = calc_corr_within(shared, full_tag_pd_centered)
# group_dict_desc = calc_corr_within(shared_desc, full_tag_pd_centered) too few conditions
cov_normalized = correct_large_matrix(group_dict, full_tag_pd_centered)
# cov_normalized_desc = correct_large_matrix(group_dict_desc, full_tag_pd_centered) too few conditions
cov_normalized.to_csv("/ebio/abt6_projects9/tnseq/tnseq_function/fitness_tables/cov_normalized_tag_dict_file.csv")
