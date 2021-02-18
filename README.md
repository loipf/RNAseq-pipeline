# RNAseq pipeline

a RNAseq quantification pipeline from `.fasta` files to a gene counts matrix using kallisto for coding and non-coding RNA.


---
### set up pipeline


before running, you have to set up the attached Docker image (will take ~30 min):
```sh
docker build -t rnaseq-pipeline https://raw.githubusercontent.com/loipf/RNAseq-pipeline/master/docker/Dockerfile
```

now either replace the Docker container hash (last output line from previous build command) in `nextflow.config` or run nextflow with the `-with-docker rnaseq-pipeline` argument.


---
### run quantification pipeline

no preprocessing or quality fitering is performed and need to be done by the user! (check file `kallisto_aligned_reads_qc.csv` for `p_pseudoaligned` >70% and `DESeq2_size_factor` around 0.7-1.3, maybe keep only genes on main chromosoms)

it can be run locally with downloaded github-repo and edited `nextflow.config` file with:
```sh
nextflow run main.nf
```

or

```sh
nextflow run loipf/RNAseq-pipeline --project_dir /path/to/folder --reads_dir /path/to/samples --ensembl_release 101 --num_threads 10 -with-docker rnsaseq-pipeline
```
for this execution to work properly, you have to be in the current project directory.


optional extendable with:
```sh
-resume
-with-report report_RNAseq-pipeline
-with-timeline timeline_RNAseq-pipeline
-w work_dir
```

for newer nextflow version you may have to specify the used version, like:
```sh
nextflow run loipf/RNAseq-pipeline -r main ...
```


by default, all output will be saved into the `data` folder of the current directory.
best to run with a new clear folder structure as not all new results do overwrite old ones.

check quality reports in `data/quality_reports` to exclude problematic samples.







