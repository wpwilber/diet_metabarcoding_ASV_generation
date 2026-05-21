This is a component of the final pipeline for handling data from the diet metabarcoding 24480 project. Unlike earlier versions of this pipeline, this project only handles preprocessing of sequenced data through ASV generation without taxonomic assignment.

Reproducing this analysis requires Conda and a Unix environment.

This project is designed for HPC deployment on the Notre Dame CRC using the Univa Grid Engine job scheduler. This can be adapted for use on SLURM or other schedulers through modification to the snakemake profile found under /profiles/hpc/config.yaml. It will also execute on compatible local systems through use of the local profile. 

How to execute this pipeline:

1: Install and activate the contained snakemake environment.

```
conda env create -f envs/snakemake-env.yaml
conda activate snakemake-env
```

2: Add raw fastq files to the fastq/ directory.

3: Execute the snakefile.

```
snakemake --profile profiles/hpc --jobs 20
```
The job number specified when running on hpc is dependent on local etiquette. The number of simultaneous jobs can scale as high as the number of fasta files you are processing. 

Contact Willi Wilber (wwilber@nd.edu) for support in reproducing this analysis. 
