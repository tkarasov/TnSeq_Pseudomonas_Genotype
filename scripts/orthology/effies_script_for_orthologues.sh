#convert_p25_to_dict.ipynb = jupyter for p25 empty dictionary
#making_dictionaries_blastp.ipynb = jupyter for making p25_DC3000, p25_p7 dictionaries, final dictionary with p25(key) and corresponding values and convertin gthe final dictionary into a pandas dataframe
#
#  script_for_aubrey.sh
#  
#
#  Created by Effie Symeonidi on 11/27/20.
#
blast=/Users/efthymiasymeonidi/Softwares/ncbi-blast-2.10.1+/bin/
#all original references are in /karasov-group1/tnseq/data/ref_database/references_sequences/
#the indexed faa files for this can be found in the folder /uufs/chpc.utah.edu/common/home/karasov-group1/tnseq/data/esymeonidi/references:
#for making the database for each pairwise comparison
#for DC3000
/Users/efthymiasymeonidi/Softwares/ncbi-blast-2.10.1+/bin/makeblastdb -in ./DC3000_GCF_000007805.1/protein.faa -dbtype prot
#similar for the others
#
#out 7 gives also comment lines while out 6 is plain tabular
#p25_c2 vs DC3000
$blast/blastp -query ../p25_c2/annotation/plate25.C2.annotation.faa -db ../DC3000_GCF_000007805.1/protein.faa -evalue 1e-6 -outfmt '6 qseqid sseqid length qlen slen qstart qend sstart send evalue pident mismatch ppos bitscore' -out p25_c2_vs_DC3000_out6.txt
#p25_c2 vs p7_g9
$blast/blastp -query ../p25_c2/annotation/plate25.C2.annotation.faa -db ../p7_g9/annotation/plate7.G9.annotation.faa -evalue 1e-6 -outfmt '6 qseqid sseqid length qlen slen qstart qend sstart send evalue pident mismatch ppos bitscore' -out p25_c2_vs_p7_g9_out6.txt
#
#70% identiy of query over 70% length similarity of query
#output all hits and not only the best ones so we can maintiain the grouping of the different orthologues
#all vs all to identify true orthologues
#filter function
#make p25.c2 dictionary for all genes (this dictionary will have keys but no values)
#output in file
#filter
#use original dictionary to fill the values
#
#(https://www.biostars.org/p/187230/) A bit score is another prominent statistical indicator used in addition to the Evalue in a BLAST output. The bit score measures sequence similarity independent of query sequence length and database size and is normalized based on the rawpairwise alignment score. The bit score (S) is determined by the following formula: S = (λ × S − lnK)/ ln2 where λ is the Gumble distribution constant, S is the raw alignment score, and K is a constant associated with the scoring matrix used. Clearly, the bit score (S) is linearly related to the rawalignment score (S). Thus, the higher the bit score, the more highly significant the match is. The bit score provides a constant statistical indicator for searching different databases of different sizes or for searching the same database at different times as the database enlarges.

#fitlering
#adding the percentage of length alignment as extra column (15)
awk -f div.awk p25_c2_vs_DC3000_out6_nohead.txt > p25_c2_vs_DC3000_out6_nohead_add15.txt
#div.awk
BEGIN { OFS="\t"}

{
    if ($4 != 0) {
        quotient = $3 / $4
        $15 = sprintf("%.3f", quotient)
    }
    else {
        $15 = "UND"
    }
        print
}
# filtering parameters
awk '{if($11>=60 && $13>=70 && $15 > 0.7)print$0}' p25_c2_vs_DC3000_out6_nohead_add15.txt > p25_c2_vs_DC3000_out6_nohead_add15_filter.txt
#
#
#keep only the top hit of each p25 protein
#sorting outputs for best hits based on  (query > bit score >e value > pident)
sort -k1,1 -k14,14gr -k10,10gr -k11,11gr p25_c2_vs_DC3000_out6_nohead_add15_filter.txt > p25_c2_vs_DC3000_out6_nohead_add15_filter_sorted.txt
#grabing the top hit for each pair
#grep needs the capital F in order to look for full string names
for next in $(cut -f1 p25_c2_vs_DC3000_out6_nohead_add15_filter_sorted.txt| sort -u)
do grep -F -m 1 "$next" p25_c2_vs_DC3000_out6_nohead_add15_filter_sorted.txt
done > p25_c2_vs_DC3000_out6_nohead_add15_filter_sorted_tophit.txt
#
#do the same for column 2 (dc3000 entries)
p25_c2_vs_DC3000_out6_nohead_add15_filter_sorted_tophit_2.txt  
#
#
#get all the p25 proteins in order to make a dictionary with these proteins are keys and empty values:
#isolate the name sof -25 proteisn:
cat plate25.C2.annotation.faa | grep ">" > p25_proteins.txt
#and grab only the name sof the preoteins (whatever is between the > and the first space)
awk -F '[> ]' '{print $2}' p25_proteins.txt > p25_proteins_names.txt

#grep is a useful commant if you want to subset information form a file. It can recognise a pattern in the whole document or in a specific column or it can even use a list of elements that the user wants to subset. It can also work in reverse, for example keep everything but the information tha tthe user provides.

######
#for subsetting column 1 (p25 protein) and column2 (Dc protein)
cat p25_c2_vs_DC3000_out6_nohead_add15_filter_sorted_tophit.txt | cut -f1 > p25_c2_vs_DC3000_out6_col1_query.txt
cat p25_c2_vs_DC3000_out6_nohead_add15_filter_sorted_tophit.txt | cut -f2 > p25_c2_vs_DC3000_out6_col2_sub.txt
#
#
