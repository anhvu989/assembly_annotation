from Bio import SeqIO
from Bio.Seq import Seq
import os

gff_path = 'gff'
fasta_path = 'fasta'
kmer_file = 'annotation7/gene_hits.txt'
outpath = 'cds'

def extract_CDS_ID(kmer_file):
    id_querry = []
    with open(kmer_file, 'r') as genefile:
        next(genefile)
        id_querry = [line.strip().split('\t')[0] for line in genefile]
    return id_querry

def extract_CDS(gff_file_path, fasta_file_path,cds_id_of_interest,output_path):
    # Initialize variables to store the extracted CDS coordinates
    os.makedirs(output_path, exist_ok=True)
    cds_start = None
    cds_end = None

    # Open and parse the GFF file to find the coordinates of the CDS of interest
    with open(gff_file_path, 'r') as gff_file:
        for line in gff_file:
            # Skip comment lines and header lines
            if line.startswith('#'):
                continue

            # Split the line into fields
            fields = line.strip().split('\t')

            # Check if the feature is a CDS
            if fields[2] == 'CDS':
                # Extract the attributes column
                attributes = fields[8]

                # Split the attributes into key-value pairs
                attribute_pairs = attributes.split(';')

                # Initialize variables to store ID and Parent attributes
                id_attribute = None
                parent_attribute = None

                # Extract ID and Parent attributes
                for pair in attribute_pairs:
                    key, value = pair.strip().split('=')
                    if key == 'ID':
                        id_attribute = value
                    elif key == 'Parent':
                        parent_attribute = value

                # Check if the CDS ID matches the one of interest
                if id_attribute == cds_id_of_interest:
                    cds_start = int(fields[3])  # Start position
                    cds_end = int(fields[4])    # End position
                    break

    # Check if CDS coordinates were found
    if cds_start is not None and cds_end is not None:
        # Open and parse the FASTA file
        fasta_sequences = SeqIO.to_dict(SeqIO.parse(fasta_file_path, 'fasta'))
        # Extract the nucleotide sequence of the CDS
        cds_sequence = fasta_sequences[fields[0]].seq[cds_start - 1 : cds_end]
        out_filename = os.path.join(output_path, f'{cds_id_of_interest}.fasta')
        with open(out_filename, 'w') as cds_writer:
            print(f'>{cds_id_of_interest}', file = cds_writer)
            print(cds_sequence, file = cds_writer)

querry = extract_CDS_ID(kmer_file)
for id_num in querry:
    sample = id_num.split("__")[0]
    gff_file_path = os.path.join(gff_path, f'{sample}.gff')
    fasta_file_path = os.path.join(fasta_path, f'{sample}.contigs.fasta')
    extract_CDS(gff_file_path,fasta_file_path,id_num,outpath)