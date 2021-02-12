
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
	CREATE_T2G_list;
	RM_DUPLICATE_TRANSCRIPTS;
	PREPROCESS_READS;
	FASTQC_READS_RAW;
	FASTQC_READS_PREPRO;
	QUANT_KALLISTO;
	CREATE_KALLISTO_QC_TABLE;
	MULTIQC_RAW;
	MULTIQC_PREPRO;
	MULTIQC_QUANT;
} from './modules.nf' 




/*
 * default parameters
 */ 

params.dev_samples = 2

params.project_dir	= "$projectDir"
params.reads_dir	= "$params.project_dir/data/reads_raw"

params.reads		= "$params.reads_dir/*/*_{1,2}.{fastq,fq}.gz"
params.data_dir		= "$params.project_dir/data"
params.scripts_dir	= "$params.project_dir/scripts"


/*
 * other parameters
 */

params.num_threads		= 2
params.ensembl_release	= 101




log.info """\
RNA-SEQ PIPELINE
===================================================
reads			: $params.reads
data_dir		: $params.data_dir
ensembl_version	: $params.ensembl_release

===================================================

"""



/* 
 * main pipeline logic
 */
workflow {
	channel_reads = Channel
			.fromFilePairs( params.reads )
			.ifEmpty { error "cannot find any reads matching: ${params.reads}" }
			.take( params.dev_samples )  // only consider a few files for debugging


	CREATE_KALLISTO_INDEX(params.ensembl_release) 
	CREATE_T2G_list(CREATE_KALLISTO_INDEX.out.raw_transcripts)
	RM_DUPLICATE_TRANSCRIPTS(CREATE_KALLISTO_INDEX.out.raw_transcripts)

	PREPROCESS_READS(channel_reads, params.num_threads)
	channel_reads_prepro = PREPROCESS_READS.out.reads_prepro.map{ it -> tuple(it[0], tuple(it[1], it[2])) }
	FASTQC_READS_RAW(channel_reads, params.num_threads)
	FASTQC_READS_PREPRO(channel_reads_prepro, params.num_threads)
	
	QUANT_KALLISTO(channel_reads_prepro, params.num_threads, CREATE_KALLISTO_INDEX.out.kallisto_index)
	CREATE_KALLISTO_QC_TABLE(QUANT_KALLISTO.out.kallisto_json.collect())



	MULTIQC_RAW(FASTQC_READS_RAW.out.reports.collect() )
	MULTIQC_PREPRO(FASTQC_READS_PREPRO.out.reports.concat(PREPROCESS_READS.out.cutadapt).collect() )
	MULTIQC_QUANT(QUANT_KALLISTO.out.kallisto_output.collect())

}




workflow.onComplete { 
	println ( workflow.success ? "\ndone! check the quality reports in --> $params.data_dir/quality_reports\n" : "oops .. something went wrong" ) } 

















