#!/usr/bin/env Rscript

############################################
### formats single kallisto entries to count matrix

pacman::p_load(data.table, tximport, rhdf5, dplyr, tidyr)

args = commandArgs(trailingOnly=TRUE)
READS_QC_TABLE_FILE = args[1]
PROCESSING_INFO_FILE = args[2]
UNIQUE_TRANS_FILE = args[3]
TRANS2GENE_FILE = args[4]

# READS_QC_TABLE_FILE = "kallisto_aligned_reads_qc.csv"
# PROCESSING_INFO_FILE = "kallisto_removal_info.txt"
# UNIQUE_TRANS_FILE = "Homo_sapiens.GRCh38.cdna_ncrna_oneline_unique.txt"
# TRANS2GENE_FILE = "transcript_to_gene_list.csv"

KALLISTO_FILES = list.files(path=".", pattern="*.h5", full.names=T)


############################################
### read in all necessary tables

tr2g_matrix = fread(TRANS2GENE_FILE)
colnames(tr2g_matrix)[2] = "transcript_type"

unique_trans_table = data.frame(fread(UNIQUE_TRANS_FILE, header = F), stringsAsFactors = FALSE)
colnames(unique_trans_table) = "transcript_id"

reads_qc_table =  data.frame(fread(READS_QC_TABLE_FILE, header = T), stringsAsFactors = FALSE, row.names = 1)

############################################
### read in all kallisto files and create gene counts

processed_samples = sapply(KALLISTO_FILES, function(file_name) { strsplit(basename(file_name),"_")[[1]][1] } )
kallisto_abundance_obj <- tximport(KALLISTO_FILES, type = "kallisto", txOut = TRUE,  importer = tximport:::read_kallisto_h5)

### add sample names to colnames
kallisto_abundance_obj = sapply(kallisto_abundance_obj, function(x) {
  if(is.matrix(x)) {
    colnames(x) = processed_samples
    x
  } else { x }
}, simplify=FALSE, USE.NAMES = T)


### output raw file
saveRDS(kallisto_abundance_obj, file="all_kallisto_abundance_obj.rds")


############################################
### create gene count matrix and filter it

sink(PROCESSING_INFO_FILE)
cat("\n\n### remove 100% sequence identity transcripts to avoid gene mapping bias\n")

### consider only unique transcripts to avoid double gene mapping bias
merged_unique_trans_gene = merge(unique_trans_table, tr2g_matrix[,c("transcript_id","gene_id")])
identical_transcripts_removed = setdiff(tr2g_matrix[["gene_id"]],merged_unique_trans_gene[["gene_id"]])
cat("genes with identical transcripts removed: ", length(identical_transcripts_removed),"/",unique(length(tr2g_matrix[["gene_id"]])),"\n")

kallisto_gene_obj = summarizeToGene(kallisto_abundance_obj, merged_unique_trans_gene[,c("transcript_id","gene_id")])
gene_matrix = kallisto_gene_obj[["counts"]]
cat("gene count matrix dimension (genes x samples): ", dim(gene_matrix), "\n")

### remove duplicates
gene_matrix = gene_matrix[!duplicated(rownames(gene_matrix)),]
cat("gene count matrix with genes with same version removed: ", dim(gene_matrix), "\n")

# ### only keep genes located on the chromosomes without scaffold  
# tr2g_matrix_split = tr2g_matrix %>% separate(chromosome, sep=":", into=c("chr_scaffold","chr_genome","chr_num","chr_start","chr_end","chr_direction"), fill="left")
# CHROMOSOME_NAMES = c(as.character(1:22), "X","Y","MT")
# genes_on_chromosomes = subset(tr2g_matrix_split, is.na(chr_scaffold) & chr_num %in% CHROMOSOME_NAMES)[["gene_id"]]
# gene_matrix = gene_matrix[rownames(gene_matrix) %in% unique(genes_on_chromosomes),]
# cat("gene count matrix with only genes located on main chromosomes: ", dim(gene_matrix),"\n")

### sort genes to newest version first, remove .version and remove duplicated gene names while keeping newest
gene_sorted_order = order(nchar(rownames(gene_matrix)), rownames(gene_matrix), decreasing = TRUE)   # to correctly handle .1 .2 .10
gene_matrix = gene_matrix[gene_sorted_order,]   # keep newest
gene_names_wo_version <- gsub('\\.[0-9]*$', "", rownames(gene_matrix))
rownames(gene_matrix) = gene_names_wo_version
gene_matrix = gene_matrix[!duplicated(rownames(gene_matrix)),]
gene_matrix = gene_matrix[order(rownames(gene_matrix), decreasing = F),]  ### sort increasing again
cat("gene count matrix without duplicated gene names: ", dim(gene_matrix),"\n")

cat("\n\n### Rsession info\n")
sessionInfo()   ### to output r-package versions
sink()

  
### output gene count matrix
fwrite(data.frame(gene_matrix), "kallisto_gene_counts.csv", quote=F, row.names=T)  # genes x samples 


############################################
### DESeq2 normalization - problem with package installation in docker
dummy_colData = data.frame(rep(1, ncol(gene_matrix)), row.names=colnames(gene_matrix))
dds = DESeqDataSetFromMatrix(round(gene_matrix), dummy_colData, design = ~1)  # no covariates included
dds <- estimateSizeFactors(dds)
dds <- estimateDispersions(dds)

# ### size factor normalized output
# gene_matrix_norm = counts(dds, normalized=T)
# fwrite(data.frame(gene_matrix_norm), "kallisto_gene_counts_norm_sf.csv", quote=F, row.names=T)

### vst or log transformation ? not necessary - should filter first ?
### https://seqqc.wordpress.com/2015/02/16/should-you-transform-rna-seq-data-log-vst-voom/
gene_matrix_vst = assay(varianceStabilizingTransformation(dds))
fwrite(data.frame(gene_matrix_vst), "kallisto_gene_counts_norm_sf_vst.csv", quote=F, row.names=T)




# 
# ### only output estimated size factors to see filtering
# ### equals DESeq2::estimateSizeFactorsForMatrix but problem install package on docker
# estimateSizeFactorsForMatrix <- function (counts, locfunc = stats::median, geoMeans, controlGenes) {
#   if (missing(geoMeans)) {
#     incomingGeoMeans <- FALSE
#     loggeomeans <- rowMeans(log(counts))
#   }
#   else {
#     incomingGeoMeans <- TRUE
#     if (length(geoMeans) != nrow(counts)) {
#       stop("geoMeans should be as long as the number of rows of counts")
#     }
#     loggeomeans <- log(geoMeans)
#   }
#   if (all(is.infinite(loggeomeans))) {
#     stop("every gene contains at least one zero, cannot compute log geometric means")
#   }
#   sf <- if (missing(controlGenes)) {
#     apply(counts, 2, function(cnts) {
#       exp(locfunc((log(cnts) - loggeomeans)[is.finite(loggeomeans) & 
#                                               cnts > 0]))
#     })
#   }
#   else {
#     if (!(is.numeric(controlGenes) | is.logical(controlGenes))) {
#       stop("controlGenes should be either a numeric or logical vector")
#     }
#     loggeomeansSub <- loggeomeans[controlGenes]
#     apply(counts[controlGenes, , drop = FALSE], 2, function(cnts) {
#       exp(locfunc((log(cnts) - loggeomeansSub)[is.finite(loggeomeansSub) & 
#                                                  cnts > 0]))
#     })
#   }
#   if (incomingGeoMeans) {
#     sf <- sf/exp(mean(log(sf)))
#   }
#   sf
# }
# 
# 
# # gene_matrix = gene_matrix[!rowSums(gene_matrix) > 10,]   ### low level genes
# size_factors = estimateSizeFactorsForMatrix(gene_matrix)
# sf_df = data.frame(size_factors)
# 

############################################
### append size factors to kallisto read info
merge_rownames_df <- function(df1, df2, ...) {
  merge_df = merge(df1, df2, by="row.names", ...)
  data.frame(merge_df, row.names = 1, check.names = F )
}

sf_df = data.frame(sizeFactors(dds))
colnames(sf_df) = c("DESeq2_size_factor")
reads_qc_table[["DESeq2_size_factor"]] <- NULL   ### fixes double read in bug with nextflow
reads_qc_table_merged = merge_rownames_df(reads_qc_table, sf_df)
reads_qc_table_merged[["sample_id"]] = rownames(reads_qc_table_merged)
fwrite(reads_qc_table_merged, "kallisto_aligned_reads_qc.csv", quote=F, row.names=T)







