

// need params.data_dir declared in main.nf, not implemented in DSL2 yet
params.data_dir	= "$launchDir/data"


process CREATE_KALLISTO_INDEX { 
	publishDir "$params.data_dir/kallisto_index", mode: "copy"

	input:
		val ensembl_release

	output:
		path "kallisto_transcripts.idx", emit: kallisto_index
		path "Homo_sapiens.GRCh38.cdna_ncrna.fa.gz", emit: raw_transcripts

	shell:
	'''
	curl --remote-name ftp://ftp.ensembl.org/pub/release-!{ensembl_release}/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz
	curl --remote-name ftp://ftp.ensembl.org/pub/release-!{ensembl_release}/fasta/homo_sapiens/ncrna/Homo_sapiens.GRCh38.ncrna.fa.gz

	gunzip Homo_sapiens.GRCh38.cdna.all.fa.gz
	gunzip Homo_sapiens.GRCh38.ncrna.fa.gz

	### combine coding and uncoding transcripts
	cat Homo_sapiens.GRCh38.cdna.all.fa Homo_sapiens.GRCh38.ncrna.fa > Homo_sapiens.GRCh38.cdna_ncrna.fa
	gzip -v Homo_sapiens.GRCh38.cdna_ncrna.fa

	kallisto index -k 31 -i kallisto_transcripts.idx Homo_sapiens.GRCh38.cdna_ncrna.fa.gz
	'''
}



process CREATE_T2G_list { 
	publishDir "$params.data_dir", mode: "copy"

	input:
		path raw_transcripts

	output:
		path "transcript_to_gene_list.csv", emit: t2g_list

	"""
	gunzip -c $raw_transcripts > raw_transcripts.fa

	perl -lne 'print \$1 while /^>(.+gene_symbol:\\S+)/gi' raw_transcripts.fa | sed 's/\\s/,/g' > transcript_to_gene_list.csv
	sed  -i '1i transcript_id,transcript_biotype,chromosome,gene_id,gene_biotype,transcript_biotype,gene_symbol' transcript_to_gene_list.csv  ### add header
	sed -i 's/chromosome://g; s/gene://g; s/gene_biotype://g; s/transcript_biotype://g; s/gene_symbol://g;' transcript_to_gene_list.csv ### remove pre-strings

	"""
}


process PREPROCESS_READS { 
	tag "$sample_id"
	publishDir "$params.data_dir/reads_prepro", pattern:"*cutadapt_output.txt", mode: "copy", saveAs: { filename -> "${sample_id}/$filename" }
	stageInMode = 'copy'   // avoids permission denied error

	input:
		tuple val(sample_id), path(reads) 
		val num_threads

	output:
		tuple val(sample_id), path("${sample_id}_prepro_1.fastq.gz"), path("${sample_id}_prepro_2.fastq.gz"), emit: reads_prepro
		path "${sample_id}_cutadapt_output.txt", emit: cutadapt

	shell:
	'''

	cutadapt --cores=!{num_threads} --max-n 0.1 --pair-filter=any --minimum-length 10 -o !{sample_id}_prepro_1.fastq.gz -p !{sample_id}_prepro_2.fastq.gz !{reads} > !{sample_id}_cutadapt_output.txt

	'''
}



// not possible to run dynamically fastqc with same name
process FASTQC_READS_RAW { 
	tag "$sample_id"
	publishDir "$params.data_dir/reads_raw", mode: "copy", overwrite: false, saveAs: { filename -> "${sample_id}/$filename" }
	stageInMode = 'copy'   // avoids permission denied error

	input:
		tuple val(sample_id), path(reads) 
		val num_threads

	output:
		path "*.zip", emit: reports
		path "*.html"

	shell:
	'''
	/home/stefan/tools/FastQC/fastqc -t !{num_threads} --noextract !{reads}
	'''
}



process FASTQC_READS_PREPRO { 
	tag "$sample_id"
	publishDir "$params.data_dir/reads_prepro", mode: "copy", overwrite: false, saveAs: { filename -> "${sample_id}/$filename" }

	input:
		tuple val(sample_id), path(reads_prepro) 
		val num_threads

	output:
		path "*.zip", emit: reports
		path "*.html"

	shell:
	'''
	/home/stefan/tools/FastQC/fastqc -t !{num_threads} --noextract !{reads_prepro}
	'''
}


process QUANT_KALLISTO { 
	tag "$sample_id"
	publishDir "$params.data_dir/reads_quant", mode: 'copy', saveAs: { filename -> "${sample_id}/$filename" }

	input:
		tuple val(sample_id), path(reads_prepro) 
		val num_threads
		path kallisto_index

	output:
		path "*", emit: all


	shell:
	'''

	kallisto quant -i !{kallisto_index} -o . -t !{num_threads} !{reads_prepro}

	'''
}




process DEEPTOOLS_ANALYSIS { 
	publishDir "$params.data_dir/reads_mapped/_deepTools", mode: 'copy'

	input:
		path reads_mapped
		path reads_mapped_index
		val num_threads

	output:
		path "*", emit:all


	shell:
	'''
	multiBamSummary bins -p !{num_threads} --smartLabels --bamfiles !{reads_mapped} -o multiBamSummary.npz

	plotCorrelation --corData multiBamSummary.npz --corMethod spearman --whatToPlot heatmap --outFileCorMatrix plotCorrelation_matrix.tsv

	plotPCA --corData multiBamSummary.npz --outFileNameData plotPCA_matrix.tsv

	plotCoverage -p !{num_threads} --ignoreDuplicates --smartLabels --bamfiles !{reads_mapped} --outRawCounts plotCoverage_rawCounts_woDuplicates.tsv > plotCoverage_output.tsv

	bamPEFragmentSize -p !{num_threads} --smartLabels --bamfiles !{reads_mapped} --table bamPEFragment_table.tsv --outRawFragmentLengths bamPEFragment_rawLength.tsv

	estimateReadFiltering -p !{num_threads} --smartLabels --bamfiles !{reads_mapped} > estimateReadFiltering_output.tsv

	'''
}






process MULTIQC_RAW { 
	publishDir "$params.data_dir/quality_reports", mode: "copy"

	input:
		path stat_files

	output:
		path "*"

	shell:
	'''
	multiqc -f -o reads_raw .
	'''
}


process MULTIQC_PREPRO { 
	publishDir "$params.data_dir/quality_reports", mode: "copy"

	input:
		path stat_files

	output:
		path "*"

	shell:
	'''
	multiqc -f -o reads_prepro .
	'''
}


process MULTIQC_MAPPED { 
	publishDir "$params.data_dir/quality_reports", mode: "copy"

	input:
		path stat_files

	output:
		path "*"

	shell:
	'''
	multiqc -f -o reads_mapped .
	'''
}


















