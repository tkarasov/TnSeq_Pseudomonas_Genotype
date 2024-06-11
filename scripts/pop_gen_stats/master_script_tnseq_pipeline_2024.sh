#!/usr/bin/env python



#performs anc reconstruction to infer events. Need output from panX
#3.28.2019 works fine
python pan_genome_anc_reconstruction_step1.py

#This next script takes output from panX and the fitness tables. It builds a large dictionary of (locus_id, condition):fitness_value mappings
#3.28.2019 Works fine. Builds full_tag_dict and can make figures on divergence vs. fitness correlation
python compare_gammaproteobacteria_step2.py

#creates tnseq gene object with attributes. 
python pan_genome_stats_step3.py

#statistical analysis with fitness data and gene attributes. Need fitness + panX
python tnseq_fitness_class_pan_genome_stats_step4.py

#compare the covariance of replicates to related environments. Working on 29.3.2019 but the convariance matrix can exceed 1.m
python variance_same_different.py



#Below is deprecated
#first script maps between tnseq and panX results. I don't see how to use this script if not using blast to assign.
#python map_tnseq_pangenome.py

#performs anc reconstruction to infer events
#python pan_genome_anc_reconstruction_step3.py

#creates tnseq gene object with a ton of attributes
#python pan_genome_stats_step2.py



#co-occurence analysis
#python co_occur_step4.py