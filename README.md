This is a component of the final pipeline for handling data from the diet metabarcoding 24480 project. Unlike earlier versions of this pipeline, this project only handles preprocessing of sequenced data through ASV generation without taxonomic assignment.

Reproducing this analysis requires Conda and a Unix environment.

This project is designed for HPC deployment on the Notre Dame CRC using the Univa Grid Engine job scheduler. This can be adapted for use on SLURM or other schedulers through modification to the snakemake profile found under /profiles/hpc/config.yaml. It will also execute on compatible local systems through use of the local profile. 

How to execute this pipeline:

1: Install and activate the contained snakemake environment.

```
conda env create -f envs/snakemake-env.yaml
conda activate snakemake-env
```

2: Add raw fastq files to the fastq/ directory organized into subdirectories by sequencing run number. This pipeline was designed for a sequencing project of 16 runs that were processed individually. This allows for run specific error modeling in DADA2. 

Example directory structure:
```
mkdir fastq/Run_01 fastq/Run_02
```

3: Execute the snakefile.

```
snakemake --profile profiles/hpc --jobs 20 --config project=Run_name
```
Run_name is the name of the sequencing run specific subdirectory. Each run is processed individually. 

For example, given the example directory structure above you would first run:
```
snakemake --profile profiles/hpc --jobs 20 --config project=Run_01
```
The job number specified when running on HPC is dependent on local etiquette. The number of simultaneous jobs can scale as high as the number of fasta files you are processing. 

4: Merge ASV tables and remove chimeras for each marker

```
qsub -N merge_chim_ITS1 scripts/merge_chim.sh ITS1
qsub -N merge_chim_trnL scripts/merge_chim.sh trnL
```
Revise the subission script as necessary to specify your own local dada2 environment and HPC credentials. A dada2 environment can be created from the included yaml:
```
conda env create -f envs/dada2.yaml
```

If running locally you can run the R scripts directly

```
Rscript scripts/merge_runs_remove_chimera.R ITS1
Rscript scripts/merge_runs_remove_chimera.R trnL
```

This finishing step can only be completed after ASVs have been generated for each sequencing run.

Contact Willi Wilber (wwilber@nd.edu) for support in reproducing this analysis. 
