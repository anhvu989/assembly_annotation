import sys
from Bio import SeqIO # conda install -c conda-forge biopython

def fasta_stats(file):
#create a dictionary with key is name of contig and value is a dict with of contig's length and contig's coverage
    
    contigs_gt500bp = []
    contigs_gt0bp = []

    for contig in SeqIO.parse(file, "fasta"):
        if len(contig) >= 500:
            contigs_gt500bp.append(len(contig))
        contigs_gt0bp.append(len(contig))
        
    contigs_gt500bp.sort(reverse=True)
    contigs_gt0bp.sort(reverse=True)
    
    n50 = 0
    n90 = 0
    n50_value = None
    n90_value = None
    for length in contigs_gt0bp:
        n50 += length
        n90 += length
        if n50_value == None and n50 >= sum(contigs_gt0bp) / 2:
            n50_value = length
            continue
        if n90_value == None and n90 >= sum(contigs_gt0bp) * 0.9:
            n90_value = length
            continue

    return {
        "length_gt0bp": sum(contigs_gt0bp),
        "length_gt500bp": sum(contigs_gt500bp),
        "num_gt0bp" : len(contigs_gt0bp),
        "num_gt500bp" : len(contigs_gt500bp),
        "largest_contig": max(contigs_gt0bp),
        "smallest_contig": min(contigs_gt0bp),
        "avg_contig_size" : sum(contigs_gt0bp)/len(contigs_gt0bp),
        "n50": n50_value,
        "n90" : n90_value
    }

if __name__ == "__main__":
    file = sys.argv[1]
    result = fasta_stats(file)  
    print(f"Total length of contigs >0bp: {result['length_gt0bp']}")
    print(f"Total length of contigs >500bp: {result['length_gt500bp']}")
    print(f"Number of contigs >0bp: {result['num_gt0bp']}")
    print(f"Number of contigs >500bp: {result['num_gt500bp']}")
    print(f"Largest contig: {result['largest_contig']}")
    print(f"Smallest contig: {result['smallest_contig']}")
    print(f"Average contig size: {result['avg_contig_size']}")
    print(f"N50: {result['n50']}")
    print(f"N90: {result['n90']}")