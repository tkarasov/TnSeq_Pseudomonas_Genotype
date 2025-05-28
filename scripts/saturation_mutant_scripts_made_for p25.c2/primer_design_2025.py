from Bio import SeqIO
from Bio.Seq import Seq

# Parameters
primer_len = 20
min_tm = 55
max_tm = 65

def calculate_tm(seq):
    seq = seq.upper()
    a = seq.count('A')
    t = seq.count('T')
    g = seq.count('G')
    c = seq.count('C')
    return 2 * (a + t) + 4 * (g + c)

def calculate_gc(seq):
    seq = seq.upper()
    gc = seq.count('G') + seq.count('C')
    return gc / len(seq) * 100

def design_primer(seq):
    for i in range(len(seq) - primer_len + 1):
        primer = seq[i:i + primer_len]
        tm = calculate_tm(primer)
        if min_tm <= tm <= max_tm:
            gc = calculate_gc(primer)
            return primer, tm, gc
    return None, None, None

def reverse_complement(seq):
    return str(Seq(seq).reverse_complement())

def load_genome(fasta_path):
    record = SeqIO.read(fasta_path, "fasta")
    return str(record.seq).upper()

def design_primers_from_gbk(genome_fasta, gbk_file, output_file):
    genome_seq = load_genome(genome_fasta)

    with open(output_file, 'w') as out:
        out.write("gene_id\tstart_primer\tstart_tm\tstart_gc\tend_primer\tend_tm\tend_gc\n")

        for record in SeqIO.parse(gbk_file, "genbank"):
            for feature in record.features:
                if feature.type == "gene":
                    start = int(feature.location.start)
                    end = int(feature.location.end)
                    strand = feature.location.strand
                    gene_id = feature.qualifiers.get("gene", ["unknown_gene"])[0]

                    gene_seq = genome_seq[start:end]

                    if strand == 1:
                        start_seq = gene_seq[:primer_len + 10]
                        end_seq = gene_seq[-(primer_len + 10):]
                        start_primer, start_tm, start_gc = design_primer(start_seq)
                        end_primer, end_tm, end_gc = design_primer(end_seq)
                    else:
                        gene_seq_rc = reverse_complement(gene_seq)
                        start_seq = gene_seq_rc[:primer_len + 10]
                        end_seq = gene_seq_rc[-(primer_len + 10):]
                        start_primer, start_tm, start_gc = design_primer(start_seq)
                        end_primer, end_tm, end_gc = design_primer(end_seq)

                    out.write(f"{gene_id}\t"
                              f"{start_primer or 'NA'}\t{start_tm or 'NA'}\t{start_gc or 'NA'}\t"
                              f"{end_primer or 'NA'}\t{end_tm or 'NA'}\t{end_gc or 'NA'}\n")

# Example usage
cd /uufs/chpc.utah.edu/common/home/karasov-group1/saturation_library/barcode_pos_map/p25.C2.gbk
gbk_file="/uufs/chpc.utah.edu/common/home/karasov-group1/saturation_library/barcode_pos_map/"
design_primers_from_gbk("p25.C2.fna", "p25.C2.gbk", "gene_primers.tsv")
