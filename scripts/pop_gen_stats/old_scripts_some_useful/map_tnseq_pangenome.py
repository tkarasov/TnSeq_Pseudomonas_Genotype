#!/usr/bin/env python
import os,sys,time,csv,glob, pickle, cPickle
from collections import defaultdict
from Bio import GenBank/Users/tkarasov/work_main/abt6_projects9/tnseq/tnseq_function/code/map_tnseq_pangenome.py
from Bio import SeqIO



'''The goal of this script is to take the output from a panX run and the gene data from a tnseq experiment and to assign the genes in the tnseq experiment to ortholog categories in the panX run
input 1: genome file of the tnseq strain (fasta format)
input 2: Output from pan
'''
#path_to_pangenome_dir='/ebio/ag-neher/share/users/wding/panX-refseq/data/Pseudomonadales/'#sys.argv[1]
path_to_pangenome_dir='/ebio/ag-neher/share/users/wding/panX-refseq/data/Enterobacteriales/'#sys.argv[1]

#Make a blastdb of all of the genes in the pan genome

def fetch_present(genome_name):
    ## https://www.ncbi.nlm.nih.gov/genome/doc/ftpfaq/#allcomplete
    ## download RefSeq data for all complete bacterial genomes?
    ## ftp://ftp.ncbi.nlm.nih.gov/genomes/genbank/bacteria/assembly_summary.txt
    ## ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria/assembly_summary.txt
    ## download refseq
    #this downloads a genome of interest specified genome_name
    os.system('wget -c ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria/assembly_summary.txt')
    genome_rec=genome_rec=os.popen('less assembly_summary.txt | grep %s'% genome_name).read().split('\t')[19]
    os.system('wget -c %s/**_translated_cds.faa.gz'%genome_rec) #cds_from_genomic.fna.gz'%genome_rec)
    os.system('gunzip *_translated_cds.faa.gz')
    #also get fna file
    os.system('wget -c  %s/**_cds_from_genomic.fna.gz'%genome_rec)
    os.system('gunzip *_cds_from_genomic.fna.gz'%genome_rec)
    
    

def make_combined_fasta(path_to_pangenome_dir):
    #makes file with representative member of pan genome
    filenames=[line.strip() for line in  os.popen('ls %s/nucleotide_fna/*.fna' % (path_to_pangenome_dir)).readlines()]
    with open("pan_genome_fna.fasta", 'w') as outfile:
        for fname in filenames:
             with open(fname, 'r') as readfile:
                 outfile.write(readfile.read() + "\n\n")

def make_blast_db():
    os.system('/usr/bin/makeblastdb -dbtype nucl -in pan_genome_fna.fasta -out pan_genome.db')

def blast_tnseq_panxdb(tnseq_genome):
    #tnseq_genome = os.popen(input_faa).readlines()[0].strip()
    #os.system('/usr/bin/tblastn -db /ebio/abt6_projects9/tnseq/data/combined_datasets/Pseudmonadales/pan_genome.db -query %s -out pan_genome_blast.txt -outfmt "7 qseqid sseqid pident length qlen slen"  -max_target_seqs 1' %(tnseq_genome))
    os.system('/usr/bin/tblastn -db /ebio/abt6_projects9/tnseq/data/combined_datasets/Enterobacteriales/pan_genome.db -query %s -out pan_genome_blast.txt -outfmt "7 qseqid sseqid pident length qlen slen"  -max_target_seqs 1' %(tnseq_genome))

def rename_blast_results():
    '''parse blast results to output locus_tag'''
    locus_tag = {}
    blast_results = [line.strip().split() for line in open("pan_genome_blast.txt")]
    for line in blast_results: 
        if line[0:2]==['#', 'Query:']:
            locus_tag[line[2]] = line[3].strip(']').strip('[').strip("locus_tag=")
    return locus_tag


def filter_blast_results(identity_cutoff, len_cutoff, locus_tag):
    '''filter results based off of desired cutoffs'''
    blast_results=[line.strip().split() for line in open("pan_genome_blast.txt") if line[0]!="#"]
    pass_filter={}
    for line in blast_results:
        query = line[0] #locus_tag[line[0]]
        subject = line[1]
        pident = float(line[2])
        if pident >= identity_cutoff:
            len_align = float(line[3])/float(line[4])
            if len_align >= len_cutoff:
                pass_filter[query] = [subject, pident, len_align]
            else:
                pass_filter[query] = "NA"
        else:
            pass_filter[query]="NA"

    return pass_filter
    

def load_pickle(filename):
    f = open(filename,"rb")
    p = cPickle.load(f)
    f.close()
    return(p)

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




#fetch genome of interest
#fetch_present(genome_name)

#Make combined fasta file for every genome in pangenome and blastdb
make_combined_fasta(path_to_pangenome_dir)
make_blast_db()#"pan_genome_fna.fasta")

#cd /ebio/abt6_projects9/tnseq/data/combined_datasets/Pseudmonadales/FW300_N1B4
#tnseq_genome = "organism_pseudo1_N1B4.faa"

#cd /ebio/abt6_projects9/tnseq/data/combined_datasets/Enterobacteriales/Ecoli_BW25113
tnseq_genome = "organism_Keio.faa"#GCF_000482265.1_EC_K12_MG1655_Broad_SNP_protein.faa"




#blast tnseq genome against ortholog group library
blast_tnseq_panxdb(tnseq_genome)

#filter blast results based on identity criteria:
locus_tag = rename_blast_results()
pass_filter=filter_blast_results(80, 0.8, locus_tag)

#output file mapping tnseq gene (first column) to ortholog group assignment
with open("/ebio/abt6_projects9/tnseq/tnseq_function/data/tnseq_pangenome_map.cpk", 'wb') as handle:
    pickle.dump(pass_filter, handle)

#list of every ortholog, the genes assigned to that orthology group and the gene number
sorted_genelist = load_sorted_clusters(path_to_pangenome_dir)

#map gene_cluster to number
gc_num = {}
for ind, (clusterID, gene) in enumerate(sorted_genelist):
    gc_num[clusterID] = ind


panx_gene_identifier = {}
for key,value in pass_filter.iteritems():
    gene = value[0]
    cluster = [rec[0] for rec in sorted_genelist if gene in rec[1][1]]
    if len(cluster) > 0: #for ind, (clusterID, gene) in enumerate(sorted_genelist):
        gc_numbering = gc_num[cluster[0]]
        panx_gene_identifier[key] = [gene, cluster[0], gc_numbering]



#output dictionary that has every tnseq result and the corresponding gene cluster
with open("/ebio/abt6_projects9/tnseq/tnseq_function/data/tnseq_panx_mapped_dict.cpk", 'wb') as handle:
    pickle.dump(panx_gene_identifier, handle)


