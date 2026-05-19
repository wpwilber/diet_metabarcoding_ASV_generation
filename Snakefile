# Define primers:
ITS1_F = "GATATCCGTTGCCGAGAGTC"
ITS1_R = "CCGAAGGCGTCAAGGAACAC"
trnL_F = "GGGCAATCCTGAGCCAA"
trnL_R = "CCATTGAGTCTCTGCACCTATC"

# Auto-discover sample names from fastq directory
SAMPLES = glob_wildcards("fastq/{sample}_R1_001.fastq.gz").sample

############################################
# FINAL TARGET
############################################

rule all:
    input:
        expand("trim_clean_qc/demux/{sample}_{amp}_{R}.fastq.gz", sample=SAMPLES, amp=["ITS1","trnL"], R=["R1","R2"])

############################################
# DEMULTIPLEXING
############################################
# Demultiplexing in this case refers to separating fastq files by the target amplicon that was sequenced. This is accomplished through pattern matching of the ITS1 or trnL forward primer defined at the head of this file. This is processed in paired end mode, such that R2 is retained when the forward primer is recognized in the paired R1. Up to 10% error rate in matching the defined forward primer is allowed (defined by -e). Read pairs that do not contain either forward primer are separated into "unnasigned" files. Forward primer matching is unanchored so that reads that do not start precisely on the forward primer are still retained. 

rule demux_by_primer:
    conda: "envs/cutadapt.yaml"
    input:
        r1 = "fastq/{sample}_R1_001.fastq.gz",
        r2 = "fastq/{sample}_R2_001.fastq.gz"
    output:
        ITS1_r1 = "trim_clean_qc/demux/{sample}_ITS1_R1.fastq.gz",
        ITS1_r2 = "trim_clean_qc/demux/{sample}_ITS1_R2.fastq.gz",
        trnL_r1 = "trim_clean_qc/demux/{sample}_trnL_R1.fastq.gz",
        trnL_r2 = "trim_clean_qc/demux/{sample}_trnL_R2.fastq.gz",
        un_r1   = "trim_clean_qc/demux/{sample}_unassigned_R1.fastq.gz",
        un_r2   = "trim_clean_qc/demux/{sample}_unassigned_R2.fastq.gz"
    threads: 8
    shell:
        r"""
        cutadapt \
            -j {threads} \
            -e 0.10 \
            --action=none \
            \
            -g ITS1={ITS1_F} \
            -g trnL={trnL_F} \
            \
            --pair-filter=any \
            \
            --untrimmed-output        {output.un_r1} \
            --untrimmed-paired-output {output.un_r2} \
            \
            -o trim_clean_qc/demux/{wildcards.sample}_{{name}}_R1.fastq.gz \
            -p trim_clean_qc/demux/{wildcards.sample}_{{name}}_R2.fastq.gz \
            {input.r1} {input.r2}
        """

