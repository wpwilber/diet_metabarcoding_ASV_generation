# Define primers:
ITS1_F = "GATATCCGTTGCCGAGAGTC"
ITS1_R = "CCGAAGGCGTCAAGGAACAC"
trnL_F = "GGGCAATCCTGAGCCAA"
trnL_R = "CCATTGAGTCTCTGCACCTATC"
# Reverse complements:
ITS1_F_RC = "GACTCTCGGCAACGGATATC"
ITS1_R_RC = "GTGTTCCTTGACGCCTTCGG"
trnL_F_RC = "TTGGCTCAGGATTGCCC"
trnL_R_RC = "GATAGGTGCAGAGACTCAATGG"
# Auto-discover sample names from fastq directory
SAMPLES = glob_wildcards("fastq/{sample}_R1_001.fastq.gz").sample

############################################
# FINAL TARGET
############################################

rule all:
    input:
        expand("trim_clean_qc/trimmed/{sample}_{amp}_{R}.primertrim.fastq.gz", sample=SAMPLES, amp=["ITS1","trnL"], R=["R1","R2"]),
        expand("trim_clean_qc/trimmed_reports/{sample}_{amp}_cutadapt.txt", sample=SAMPLES, amp=["ITS1","trnL"]),
        expand("trim_clean_qc/qc/{sample}_{amp}_{R}.primertrim_fastqc.html", sample=SAMPLES, amp=["ITS1","trnL"], R=["R1","R2"]),
        expand("trim_clean_qc/qc/{sample}_{amp}_{R}.primertrim_fastqc.zip", sample=SAMPLES, amp=["ITS1","trnL"], R=["R1","R2"])
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

############################################
# TRIMMING
############################################
# This rule trims primers from the forward and reverse reads using identical error allowance as was permitted by demuxing (10% error for each primer and unanchored primer locations). Because the target amplicons are short, targets are sequenced beyond the length of the amplicon. To isolate the amplicon, I trim from the forward primer to the reverse complement of the reverse primer on R1, and from the reverse primer to the reverse complement of the forward primer on R2. Reads that do not contain both primers at the allowed error rate are discarded. Strict N filtering is applied so that variable regions with ambiguous bases are filtered from the data set. A modest minimum length is also applied that helps to remove a few outlying erroneous reads that were identified during quality assessment. Filtering occurs after trimming so that sequences are not discarded on the basis of low quality tails and real biological variation in target amplicon length is preserved. Quality filtering is not handled here and is instead handled by DADA2 during ASV inference as this information may improve DADA2's results. Filtering and trimming results for each sample can be found in trim_clean_qc/trimmed_reports. 
rule trim_primers_cutadapt:
    conda: "envs/cutadapt.yaml"
    input:
        r1 = "trim_clean_qc/demux/{sample}_{amp}_R1.fastq.gz",
        r2 = "trim_clean_qc/demux/{sample}_{amp}_R2.fastq.gz"
    output:
        r1_trim = "trim_clean_qc/trimmed/{sample}_{amp}_R1.primertrim.fastq.gz",
        r2_trim = "trim_clean_qc/trimmed/{sample}_{amp}_R2.primertrim.fastq.gz",
        report = "trim_clean_qc/trimmed_reports/{sample}_{amp}_cutadapt.txt"
    threads: 4
    params:
        F = lambda wc: {"ITS1": ITS1_F, "trnL": trnL_F}[wc.amp],
        R = lambda wc: {"ITS1": ITS1_R, "trnL": trnL_R}[wc.amp],
        F_RC = lambda wc: {"ITS1": ITS1_F_RC, "trnL": trnL_F_RC}[wc.amp],
        R_RC = lambda wc: {"ITS1": ITS1_R_RC, "trnL": trnL_R_RC}[wc.amp],
        min_len = lambda wc: {"ITS1": 50, "trnL": 10}[wc.amp]
    shell:
        r"""
        mkdir -p trim_clean_qc/trimmed trim_clean_qc/trimmed_reports

        cutadapt \
            -j {threads} \
            -e 0.10 \
            -g "{params.F}...{params.R_RC}" \
            -G "{params.R}...{params.F_RC}" \
            --discard-untrimmed \
            --max-n 0 \
            --minimum-length {params.min_len} \
            -o {output.r1_trim} \
            -p {output.r2_trim} \
            {input.r1} {input.r2} \
            > {output.report} 2>&1
        """

############################################
# QUALITY REPORTING
############################################
# This rule runs fastqc on trimmed sequences for troubleshooting before delivering sequences to DADA2.
rule fastqc_final:
    conda: "envs/fastqc.yaml"
    input:
        r1 = "trim_clean_qc/trimmed/{sample}_{amp}_R1.primertrim.fastq.gz",
        r2 = "trim_clean_qc/trimmed/{sample}_{amp}_R2.primertrim.fastq.gz"
    output:
        "trim_clean_qc/qc/{sample}_{amp}_R1.primertrim_fastqc.html",
        "trim_clean_qc/qc/{sample}_{amp}_R1.primertrim_fastqc.zip",
        "trim_clean_qc/qc/{sample}_{amp}_R2.primertrim_fastqc.html",
        "trim_clean_qc/qc/{sample}_{amp}_R2.primertrim_fastqc.zip"
    threads: 2
    shell:
        r"""
        fastqc --threads {threads} --outdir trim_clean_qc/qc {input.r1} {input.r2}
        """

############################################
# ASV inference
############################################
# The following rules run DADA2 ASV inference on trimmed and filtered fastq files described above. For detailed comments on the DADA2 workflow, reference the script in scripts/run_dada2_{amp}.R. 
rule dada2_trnl:
    conda: "envs/DADA2.yaml"
    input:
        expand("trim_clean_qc/trimmed/{sample}_trnL_R1.primertrim.fastq.gz", sample=SAMPLES),
        expand("trim_clean_qc/trimmed/{sample}_trnL_R2.primertrim.fastq.gz", sample=SAMPLES)
    output:
        tsv = "dada2/trnL_seqtab_all.tsv",
        fasta = "dada2/trnL_ASVs.fasta"
    threads: 8
    shell:
        r"""
        mkdir -p dada2
        Rscript scripts/run_dada2_trnL.R
        """

# The following rule runs DADA2 sample inference on ITS1 samples.

rule dada2_its1:
    conda: "envs/DADA2.yaml"
    input:
        expand("trim_clean_qc/trimmed/{sample}_ITS1_R1.primertrim.fastq.gz", sample=SAMPLES),
        expand("trim_clean_qc/trimmed/{sample}_ITS1_R2.primertrim.fastq.gz", sample=SAMPLES)
    output:
        tsv = "dada2/ITS1_seqtab_all.tsv",
        fasta = "dada2/ITS1_ASVs.fasta"
    threads: 8
    shell:
        r"""     
        mkdir -p dada2
        Rscript scripts/run_dada2_ITS1.R
        """

