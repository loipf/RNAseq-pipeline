# RNAseq pipeline

a RNAseq quantification pipeline from `.fasta` files to a gene counts matrix using kallisto.


---
### set up pipeline


before running, you have to set up the attached Docker image which includes OpenCRAVAT and 2 others (all will take quite some time and overall need >50 GB of memory, you can add/remove OpenCRAVAT modules according to your needs in the Dockerfile, currently adapted for breast cancer):
```sh
docker build -t dnavc-pipeline https://raw.githubusercontent.com/loipf/DNAvc-pipeline/master/docker/Dockerfile
docker pull google/deepvariant:1.1.0
docker pull quay.io/mlin/glnexus:v1.2.7
```

now either replace the Docker container hash (last output line from previous build command) in `nextflow.config` or run nextflow with the `-with-docker dnavc-pipeline` argument.


---
### run mapping pipeline

no preprocessing or quality fitering is performed and need to be done by the user! (check file `kallisto_aligned_reads_qc.csv` for `p_pseudoaligned` >70% and `DESeq2_size_factor` around 0.7-1.3)

it can be run locally with downloaded github-repo and edited `nextflow.config` file with:
```sh
nextflow run main.nf
```

or

```sh
nextflow run loipf/DNAvc-pipeline --project_dir /path/to/folder --reads_mapped_dir /path/to/samples --num_threads 10 -with-docker dnavc-pipeline
```
for this execution to work properly, you have to be in the current project directory.


optional extendable with:
```sh
-resume
-with-report report_DNAseq-pipeline
-with-timeline timeline_DNAseq-pipeline
-w work_dir
```


by default, all output will be saved into the `data` folder of the current directory.
best to run with a new clear folder structure as not all new results do overwrite old ones.







