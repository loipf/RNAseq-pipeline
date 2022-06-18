# RNAseq pipeline

a RNAseq quantification pipeline from `.fastq` files to a gene counts matrix using kallisto for coding and non-coding RNA.


---
### set up pipeline


before running, you have to set up the attached Docker image (will take ~30 min):
```sh
docker build -t rnaseq-pipeline https://raw.githubusercontent.com/loipf/RNAseq-pipeline/master/docker/Dockerfile
```

now either replace the Docker container hash (last output line from previous build command) in `nextflow.config` or run nextflow with the `-with-docker rnaseq-pipeline` argument.


---
### run quantification pipeline

no pre-processing or quality improvement is performed and must be done by the user! (check file `kallisto_aligned_reads_qc.csv` for `p_pseudoaligned` >70% and `DESeq2_size_factor` around 0.6-1.4). sums up transcripts to gene_symbols (without haplotype and scaffold genes).

it can be run locally with downloaded github-repo and edited `nextflow.config` file with:
```sh
nextflow run main.nf
```

or

```sh
nextflow run loipf/RNAseq-pipeline -r main --project_dir /path/to/folder --reads_dir /path/to/samples --ensembl_release 101 --num_threads 10 -with-docker rnsaseq-pipeline
```
for this execution to work properly, you have to be in the current project directory.


optional extendable with:
```sh
-resume
-with-report report_RNAseq-pipeline
-with-timeline timeline_RNAseq-pipeline
-w work_dir
```

by default, all output will be saved into the `data` folder of the current directory.
best to run with a new clear folder structure as not all new results do overwrite old ones.

check quality reports in `data/quality_reports` to exclude problematic samples.


additional, an 3' and 5' adapter sequence (file) can be specified with the nextflow arguments `--adapter_3_seq_file [sequence|file.fasta]` and `--adapter_5_seq_file [sequence|file.fasta]` or in the `main.nf` file. if a file is provided, it must be structured like the following example:
```
> adapter_3_batch_01
AANTGG
> adapter_3_batch_02
GATCGG
```





