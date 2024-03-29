
/* 
 * RNA-SEQ PIPELINE 
 * for paired-end reads
 */


/* 
 * import modules 
 */

nextflow.enable.dsl=2

include { 
	CREATE_KALLISTO_INDEX;
	CREATE_T2G_LIST;
	RM_DUPLICATE_TRANSCRIPTS;
	PREPROCESS_READS;
	FASTQC_READS_RAW;
	FASTQC_READS_PREPRO;
	QUANT_KALLISTO;
	CREATE_KALLISTO_QC_TABLE;
	CREATE_GENE_MATRIX;
	CREATE_GENE_COUNT_PLOTS;
	MULTIQC_RAW;
	MULTIQC_PREPRO;
	MULTIQC_QUANT;
} from './modules.nf' 




/*
 * default parameters
 */ 

params.dev_samples = -1

params.project_dir	= "$projectDir"
params.reads_dir	= "$params.project_dir/data/reads_raw"

params.reads		= "$params.reads_dir/*/*{1,2}.{fastq,fq,fastq.gz,fq.gz}"
params.data_dir	= "$params.project_dir/data"

// either give both or none
params.adapter_3_seq_file	= file("NO_FILE")
params.adapter_5_seq_file	= file("NO_FILE2")


/*
 * other parameters
 */

params.num_threads	= 5
params.ensembl_release	= 101
params.include_ncrna	= true  // false
params.nextflow_stageInMode	= "symlink" // "copy" avoids permission denied error




log.info """\
RNA-SEQ PIPELINE
===================================================
reads			: $params.reads
data_dir		: $params.data_dir
ensembl_version	: $params.ensembl_release
adapter_3_seq_file	: $params.adapter_3_seq_file
adapter_5_seq_file	: $params.adapter_5_seq_file
===================================================

"""



/* 
 * main pipeline logic
 */
workflow {
	// extension from .fromFilePairs( params.reads ) to deal with multiple sequencing files in the same folder which get combined in module
	channel_reads = Channel
			.fromPath(params.reads)
			.map{ files -> tuple(files.getParent().getName(), files) }
			.groupTuple()
			.ifEmpty { error "cannot find any reads matching: ${params.reads}" }
			.take( params.dev_samples )  // only consider a few files for debugging

	CREATE_KALLISTO_INDEX(params.ensembl_release) 
	CREATE_T2G_LIST(CREATE_KALLISTO_INDEX.out.raw_transcripts)
	RM_DUPLICATE_TRANSCRIPTS(CREATE_KALLISTO_INDEX.out.raw_transcripts)

	PREPROCESS_READS(channel_reads, params.num_threads, params.adapter_3_seq_file, params.adapter_5_seq_file)
	channel_reads_prepro = PREPROCESS_READS.out.reads_prepro.map{ it -> tuple(it[0], tuple(it[1], it[2])) }
	FASTQC_READS_RAW(channel_reads, params.num_threads)
	FASTQC_READS_PREPRO(channel_reads_prepro, params.num_threads)
	
	QUANT_KALLISTO(channel_reads_prepro, params.num_threads, CREATE_KALLISTO_INDEX.out.kallisto_index)
	CREATE_KALLISTO_QC_TABLE(QUANT_KALLISTO.out.kallisto_json.collect())
	CREATE_GENE_MATRIX(CREATE_KALLISTO_QC_TABLE.out.kallisto_qc_table, RM_DUPLICATE_TRANSCRIPTS.out.removal_info, RM_DUPLICATE_TRANSCRIPTS.out.trans_oneline_unique, CREATE_T2G_LIST.out.t2g_list, QUANT_KALLISTO.out.kallisto_abundance.collect() )

	CREATE_GENE_COUNT_PLOTS(CREATE_GENE_MATRIX.out.kallisto_qc_table, CREATE_GENE_MATRIX.out.gene_matrix, CREATE_GENE_MATRIX.out.gene_matrix_vst)

	MULTIQC_RAW(FASTQC_READS_RAW.out.reports.collect() )
	MULTIQC_PREPRO(FASTQC_READS_PREPRO.out.reports.concat(PREPROCESS_READS.out.cutadapt).collect() )
	MULTIQC_QUANT(QUANT_KALLISTO.out.kallisto_output.collect())

}




workflow.onComplete { 
	println ( workflow.success ? "\ndone! check the quality reports in --> $params.data_dir/quality_reports\n" : "oops .. something went wrong" ) } 

















