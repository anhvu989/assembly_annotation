from Bio import SeqIO
from Bio.Seq import Seq
import os

gff_path = 'gff'
kmer_file = 'annotation7/gene_hits.txt'

def extract_CDS_ID(kmer_file):
    id_querry = []
    with open(kmer_file, 'r') as genefile:
        next(genefile)
        id_querry = [line.strip().split('\t')[0] for line in genefile]
    return id_querry

def extract_attribute(gff_file_path,cds_id_of_interest):
    # Initialize variables to store the extracted CDS coordinates
    ec_number = None
    GO = None
    uniprot = None
    gene = None
    product = None
    cog = None

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

                # Check if the CDS ID matches the one of interest
                if id_attribute == cds_id_of_interest:
                    for pair in attribute_pairs:
                        key, value = pair.strip().split('=')
                        if key =='eC_number':
                            ec_number = value
                        elif key == 'gene':
                            gene = value
                        elif key == 'inference':
                            uniprot = value.split(',')[-1].split(':')[-1]
                        elif key == 'product':
                            product = value
                        elif key == 'db_xref':
                            cog = value
                    break
    return f'{cds_id_of_interest}\t{ec_number}\t{gene}\t{uniprot}\t{cog}\t{product}'


sum_file = 'summary_attribute.txt'
querry = extract_CDS_ID(kmer_file)
for protein in querry:
    sample = protein.split('__')[0]
    gff_file_path = os.path.join(gff_path, f'{sample}.gff')
    result = extract_attribute(gff_file_path, protein)
    with open(sum_file,'a') as outfile:
        print(result, file=outfile)   
