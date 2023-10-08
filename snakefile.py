import os
import glob

# load file containing names of samples with 70-99% of L.crispatus reads
list_samples_file = "/pasteur/zeus/projets/p02/metasig/gitlab/genome_comparison/list_samples.txt"

with open(list_samples_file, 'r') as file:
    SAMPLES = [line.strip() for line in file]
READS = ["1","2"]

#check whether sample's reads are in rawreads folder
new_samples = []
odd_samples_name = []
for sam in SAMPLES:
    sam = sam[6:]
    alist = []
    if int(sam) in range(476, 571):
        new_sam = str(int(sam) + 1000)
        sam_dir = os.path.join("rawreads",f'*{new_sam}*')
    else:
        sam_dir = os.path.join("rawreads",f'*{sam}*')
    alist = glob.glob(sam_dir)
    if len(alist) == 0:
        sam = sam[1:]
        sam_dir = os.path.join("rawreads",f'*_{sam}_*')
        alist = glob.glob(sam_dir)
    if len(alist) >= 2:
        if len(sam) == 3:
            odd_samples_name.append(sam)
        new_samples.append("A4138_" + sam.zfill(4))

SAMPLES = new_samples

#create directory for each sample's results
for sam in SAMPLES:
    sample_dir = os.path.join("results", sam)
    subdirectories = ["bowtie2", "spades","prokka"]  # List of subdirectory names
    os.makedirs(sample_dir, exist_ok=True)
    for subdir in subdirectories:
        subdir_path = os.path.join(sample_dir, subdir)
        os.makedirs(subdir_path, exist_ok=True)

def getrawreads(sample):
    sample = sample[6:]
    if sample[1:] in odd_samples_name:
        sample = sample[1:]
    raw_read_file = []
    pattern = []
    if int(sample) in range(476, 571):
        new_sample = str(int(sample) + 1000)
        pattern.append(f'A4138_*{new_sample}*')
    elif int(sample) in range (1476, 1571):
        pattern1 = f'??_*{sample}*'
        pattern2 = f'???_*{sample}*'    
        pattern.append(pattern1)
        pattern.append(pattern2)
    else:
        pattern.append(f'*_{sample}_*')
    for pat in pattern:
        sample_path = os.path.join("rawreads",pat)
        raw_read_file.extend(glob.glob(sample_path))
    if len(raw_read_file) == 2:
        return raw_read_file
    elif len(raw_read_file) >2:
        with open('unknown.txt', 'a') as outfile:
            print(sample, file = outfile)
            return raw_read_file

rule all:
    input:
        "pipeline_end.txt"

rule aligned_via_bowtie2: 
    input:
        rawreads = lambda wildcards: getrawreads(wildcards.sample)
    params:
        btindex = "mapping/reference/concat_ref.btindex",
        log = "results/{sample}/bowtie2/{sample}_bowtie2.log"
    output:
        samfile = "results/{sample}/bowtie2/{sample}.sam"
    shell:
        """
        bowtie2 -p 4 \
        -x {params.btindex} -1 {input.rawreads[0]} -2 {input.rawreads[1]} \
        -S {output.samfile} 2> {params.log} 
        """

rule filtered_bam:
    input: "results/{sample}/bowtie2/{sample}.sam"
    output: "results/{sample}/bowtie2/{sample}.filtered.bam"
    shell:
    # exclude reads where neither pair is mapped with -F 12 (UNMAP, MUNMAP)
    # -f include, -F exclude
    # use 'samtools flags 12' to check. Can use any number
        """
        samtools view -SbF 12 {input} | \
        samtools sort -o {output}
        """

rule extract_mapped_reads:
    input: "results/{sample}/bowtie2/{sample}.filtered.bam"
    output:
        read1 = "results/{sample}/bowtie2/{sample}.aligned.1.fastq.gz",
        read2 = "results/{sample}/bowtie2/{sample}.aligned.2.fastq.gz"
    shell:
    ### samtools collate – shuffles and groups reads together by their names
    ### -u : Write uncompressed BAM output
    ### -O : Output to stdout.
    ### samtools fastq - converts a SAM/BAM/CRAM file to FASTA or FASTQ
    ### -n By default, either '/1' or '/2' is added to the end of read names where the corresponding READ1 or READ2 FLAG bit is set. 
    ### Using -n causes read names to be left as they are.
    ### -s write singleton reads to FILE
    ### -0 Write reads where the READ1 and READ2 FLAG bits set are either both set or both unset to FILE instead of outputting them.
        """
        samtools collate -Ou {input} | \
        samtools fastq -1 {output.read1} -2 {output.read2} -0 /dev/null -s /dev/null -n
        """

rule assembly_via_spades:
    input: 
        read1 = rules.extract_mapped_reads.output.read1,
        read2 = rules.extract_mapped_reads.output.read2
    output: 
        contigs = "results/{sample}/spades/{sample}.contigs.fasta",
        sample_dir = directory("results/{sample}/spades/")
    threads:
        8
    shell:
        """
        spades.py --careful --cov-cutoff auto --threads {threads} \
        -1 {input.read1} -2 {input.read2} -o {output.sample_dir} &&
        mv "results/{wildcards.sample}/spades/contigs.fasta" "results/{wildcards.sample}/spades/{wildcards.sample}.contigs.fasta"
        """

rule annotation_via_prokka:
    input:
        contigs = "results/{sample}/spades/{sample}.contigs.fasta"
    output:
        result_file = "results/{sample}/prokka/{sample}.gff",
        result_dir = directory("results/{sample}/prokka/")
    shell:
        """
        prokka --force --kingdom Bacteria --prefix {wildcards.sample} --locustag {wildcards.sample}_ --outdir {output.result_dir} {input.contigs}
        """

rule run_completed:
    input:
        assemblies= expand("results/{sample}/prokka/{sample}.gff", sample=SAMPLES)
    output:
        end_report = "pipeline_end.txt"
    shell:
        """
        touch {output.end_report}
        """
