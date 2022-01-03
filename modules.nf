
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



process CREATE_T2G_LIST { 
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


// removes duplicate transcripts before merging into gene names
process RM_DUPLICATE_TRANSCRIPTS { 
	publishDir "$params.data_dir/kallisto_index", pattern:"Homo*", mode: "copy"

	input:
		path raw_transcripts 

	output:
		path "Homo_sapiens.GRCh38.cdna_ncrna_oneline_unique.txt", emit: trans_oneline_unique
		path "kallisto_removal_info.txt", emit: removal_info

	shell:
	'''
	gunzip -c !{raw_transcripts} | cut -f1 -d" " | awk '/^>/ {printf("\\n%s\\n",$0);next; } { printf("%s",$0);} END {printf("\\n");}' | tail -n +2 | sed 'N;s/\\n/ /' > Homo_sapiens.GRCh38.cdna_ncrna_oneline.txt
	
	remove_duplicate_transcripts.R Homo_sapiens.GRCh38.cdna_ncrna_oneline.txt
	'''
}



process PREPROCESS_READS { 
	tag "$sample_id"
	publishDir "$params.data_dir/reads_prepro", pattern:"*cutadapt_output.txt", mode: "copy", saveAs: { filename -> "${sample_id}/$filename" }
	stageInMode = 'copy'   // avoids permission denied error
	cache false

	input:
		tuple val(sample_id), path(reads) 
		val num_threads
		path adapter_3_seq_file
		path adapter_5_seq_file

	output:
		tuple val(sample_id), path("${sample_id}_prepro_1.fastq.gz"), path("${sample_id}_prepro_2.fastq.gz"), emit: reads_prepro
		path "${sample_id}_cutadapt_output.txt", emit: cutadapt

	shell:
	'''
	
	### combine multiple seq files in the same sample directory with same direction together
	### theoretically sorting not needed but to be ordered
	reads_sorted=$(echo !{reads} | xargs -n1 | sort | xargs)
	reads_sorted_1=$(find $reads_sorted -name "*_1.fq.gz" -o -name "*_1.fastq.gz")
	reads_sorted_2=$(find $reads_sorted -name "*_2.fq.gz" -o -name "*_2.fastq.gz")
	cat $reads_sorted_1 > raw_reads_connected_1.fastq.gz
	cat $reads_sorted_2 > raw_reads_connected_2.fastq.gz

	if [[ (!{adapter_3_seq_file} == "NO_FILE") || (!{adapter_5_seq_file} == "NO_FILE2") ]]
	then
		### run with adapter sequences
		
		### adapter 3' input as string or file
		if [[ !{adapter_3_seq_file} == *".fasta"* ]]
		then
	  		ADAPTER_3=file:!{adapter_3_seq_file}
		else
			ADAPTER_3="!{adapter_3_seq_file}"
		fi
		### adapter 5' input as string or file
		if [[ !{adapter_5_seq_file} == *".fasta"* ]]
		then
	  		ADAPTER_5=file:!{adapter_5_seq_file}
		else
			ADAPTER_5="!{adapter_5_seq_file}"
		fi
		
		cutadapt --cores=!{num_threads} --max-n 0.1 --discard-trimmed --pair-filter=any --minimum-length 10 -a $ADAPTER_3 -A $ADAPTER_5 -o !{sample_id}_prepro_1.fastq.gz -p !{sample_id}_prepro_2.fastq.gz raw_reads_connected_1.fastq.gz raw_reads_connected_2.fastq.gz > !{sample_id}_cutadapt_output.txt
		
	else
		### run without adapter sequence
		cutadapt --cores=!{num_threads} --max-n 0.1 --discard-trimmed --pair-filter=any --minimum-length 10 -o !{sample_id}_prepro_1.fastq.gz -p !{sample_id}_prepro_2.fastq.gz raw_reads_connected_1.fastq.gz raw_reads_connected_2.fastq.gz > !{sample_id}_cutadapt_output.txt
	fi
	
	'''
}



// not possible to run dynamically fastqc with same name
process FASTQC_READS_RAW { 
	tag "$sample_id"
	publishDir "$params.data_dir/reads_raw", mode: "copy", overwrite: false, saveAs: { filename -> "${sample_id}/$filename" }
	stageInMode = 'copy'   // avoids permission denied error
	cache false

	input:
		tuple val(sample_id), path(reads) 
		val num_threads

	output:
		path "*.zip", emit: reports
		path "*.html"

	shell:
	'''
	fastqc -t !{num_threads} --noextract !{reads}
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
	fastqc -t !{num_threads} --noextract !{reads_prepro}
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
		path "${sample_id}_abundance.h5", emit: kallisto_abundance
		path "${sample_id}_run_info.json", emit: kallisto_json
		path "${sample_id}_kallisto_output.txt", emit: kallisto_output

	shell:
	'''

	kallisto quant -i !{kallisto_index} -o . -t !{num_threads} !{reads_prepro} &> kallisto_output.txt

	for f in *; do mv -- "$f" "!{sample_id}_${f%}"; done   # rename all with sample_id prefix

	'''
}



// removes duplicate transcripts before merging into gene names
process CREATE_KALLISTO_QC_TABLE { 
	input:
		path kallisto_json

	output:
		path "kallisto_aligned_reads_qc.csv", emit: kallisto_qc_table

	shell:
	'''
	create_kallisto_qc_table.R
	'''
}


process CREATE_GENE_MATRIX { 
	publishDir "$params.data_dir", mode: "copy", saveAs: { filename -> filename.endsWith(".rds") ? "reads_quant/$filename" : filename }

	input:
		path kallisto_qc_table
		path removal_info
		path trans_oneline_unique
		path t2g_list
		path kallisto_abundance
		
	output:
		path "kallisto_gene_counts.csv", emit: gene_matrix
		path "kallisto_gene_counts_norm_sf_vst.csv", emit: gene_matrix_vst
		path "kallisto_aligned_reads_qc.csv", emit: kallisto_qc_table
		path "all_kallisto_abundance_obj.rds"
		path "kallisto_removal_info.txt"
		path "transcript_to_gene_list.csv"

	shell:
	'''
	create_kallisto_gene_matrix.R !{kallisto_qc_table} !{removal_info} !{trans_oneline_unique} !{t2g_list} 
	'''
}




process CREATE_GENE_COUNT_PLOTS { 
	publishDir "$params.data_dir/quality_reports", mode: "copy"

	input:
		path kallisto_qc_table
		path gene_matrix
		path gene_matrix_vst
		
	output:
		path "gene_count_analysis_plots.html"
	
	shell:
	'''
	#!/usr/bin/env Rscript

	curr_dir = getwd()
	rmarkdown_file_path = system("which create_gene_count_plots.R", intern = TRUE)
	
	params_list = list(curr_dir = curr_dir, reads_qc_table_file='!{kallisto_qc_table}', gene_matrix_file='!{gene_matrix}',gene_matrix_file_vst='!{gene_matrix_vst}')

	rmarkdown::render(rmarkdown_file_path, output_file=file.path(curr_dir,'gene_count_analysis_plots.html'), params= params_list )

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


process MULTIQC_QUANT { 
	publishDir "$params.data_dir/quality_reports", mode: "copy"

	input:
		path stat_files

	output:
		path "*"

	shell:
	'''
	multiqc -f -o reads_quant .
	'''
}


















