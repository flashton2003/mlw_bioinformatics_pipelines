# mlw_bioinformatics_pipelines
This repository contains bioinformatics pipelines for sequencing data from both nanopore and illumina platforms.

## Running the nanopore script
* Navigate to the test directory via
```
cd /home/ubuntu/data/belson/test_seq_service/nanopore
```

* First activate the snakemake env via
```
conda activate snakemake
```
* Run the nanopore pipeline via

```
snakemake --snakefile nanopore.smk --cores 1
```



