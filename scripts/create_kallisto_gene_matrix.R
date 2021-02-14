
############################################
### formats single kallisto entries to count matrix

library(data.table)
library(tximport)
library(DESeq2)
library(rhdf5)
library(dplyr)
library(tidyr)


args = commandArgs(trailingOnly=TRUE)

READS_QC_TABLE_FILE = args[1]
PROCESSING_INFO_FILE = args[2]
UNIQUE_TRANS_FILE = args[3] 
TRANS2GENE_FILE = args[4]

KALLISTO_FILES = list.files(path=".", pattern="*.h5", full.names=T)


# KALLISTO_FILES = list.files(path="/home/stefan/Documents/umcg/RNAseq-pipeline/data/reads_quant", pattern="*.h5", recursive = T, full.names = T)
# READS_QC_TABLE_FILE = "/home/stefan/Documents/umcg/RNAseq-pipeline/data/kallisto_aligned_reads_qc.csv"
# UNIQUE_TRANS_FILE = "/home/stefan/Documents/umcg/RNAseq-pipeline/data/kallisto_index/Homo_sapiens.GRCh38.cdna_ncrna_oneline_unique.txt"
# TRANS2GENE_FILE = '/home/stefan/Documents/umcg/RNAseq-pipeline/data/transcript_to_gene_list.csv'


############################################
### read in all necessary tables

tr2g_matrix = fread(TRANS2GENE_FILE)
colnames(tr2g_matrix)[2] = "transcript_type"

unique_trans_table = data.frame(fread(UNIQUE_TRANS_FILE, header = F), stringsAsFactors = FALSE)
colnames(unique_trans_table) = "transcript_id"

reads_qc_table =  data.frame(fread(READS_QC_TABLE_FILE, header = T), stringsAsFactors = FALSE)


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
saveRDS(kallisto_abundance_obj, file="kallisto_abundance_obj.rds")


############################################
### create gene count matrix and filter it

sink(PROCESSING_INFO_FILE)
cat("\n\n### remove 100% sequence identity transcripts to avoid gene mapping bias\n")


### consider only unique transcripts to avoid double gene mapping bias
merged_unique_trans_gene = merge(unique_trans_table, tr2g_matrix[,c("transcript_id","gene_id")])
identical_transcripts_removed = setdiff(tr2g_matrix$gene_id,merged_unique_trans_gene$gene_id)
cat("genes with identical transcripts removed: ", length(identical_transcripts_removed),"/",unique(length(tr2g_matrix$gene_id)),"\n")

kallisto_gene_obj = summarizeToGene(kallisto_abundance_obj, merged_unique_trans_gene[,c("transcript_id","gene_id")])
gene_matrix = kallisto_gene_obj$counts
cat("gene count matrix dimension (genes x samples): ", dim(gene_matrix), "\n")

### remove duplicates
gene_matrix = gene_matrix[!duplicated(rownames(gene_matrix)),]
cat("gene count matrix with genes with same version removed: ", dim(gene_matrix), "\n")

# ### only keep genes located on the chromosomes without scaffold  
# tr2g_matrix_split = tr2g_matrix %>% separate(chromosome, sep=":", into=c("chr_scaffold","chr_genome","chr_num","chr_start","chr_end","chr_direction"), fill="left")
# CHROMOSOME_NAMES = c(as.character(1:22), "X","Y","MT")
# genes_on_chromosomes = subset(tr2g_matrix_split, is.na(chr_scaffold) & chr_num %in% CHROMOSOME_NAMES)$gene_id
# gene_matrix = gene_matrix[rownames(gene_matrix) %in% unique(genes_on_chromosomes),]
# cat("gene count matrix with only genes located on main chromosomes: ", dim(gene_matrix),"\n")

### sort genes to newest version first, remove .version and remove duplicated gene names while keeping newest
gene_sorted_order = order(nchar(rownames(gene_matrix)), rownames(gene_matrix), decreasing = TRUE)   # to correctly handle .1 .2 .10
gene_matrix = gene_matrix[gene_sorted_order,]   # keep newest
gene_names_wo_version <- gsub("\\.[0-9]*$", "", rownames(gene_matrix))
rownames(gene_matrix) = gene_names_wo_version
gene_matrix = gene_matrix[!duplicated(rownames(gene_matrix)),]
gene_matrix = gene_matrix[order(rownames(gene_matrix), decreasing = F),]  ### sort increasing again
cat("gene count matrix without duplicated gene names: ", dim(gene_matrix),"\n")

sink()

### output gene count matrix
fwrite(data.frame(gene_matrix), "kallisto_gene_counts.csv", quote=F, row.names=T) 



############################################
### DESeq2 normalization
dummy_colData = data.frame(rep(1, ncol(gene_matrix)), row.names=colnames(gene_matrix))
dds = DESeqDataSetFromMatrix(round(gene_matrix), dummy_colData, design = ~1)  # no covariates included
dds <- estimateSizeFactors(dds)
dds <- estimateDispersions(dds)

# ### size factor normalized output
# gene_matrix_norm = counts(dds, normalized=T)
# fwrite(data.frame(gene_matrix_norm), "kallisto_gene_counts_norm_sf.csv", quote=F, row.names=T) 

gene_matrix_vst = assay(varianceStabilizingTransformation(dds))
fwrite(data.frame(gene_matrix_vst), "kallisto_gene_counts_norm_sf_vst.csv", quote=F, row.names=T)


### append size factors to kallisto read info
sf_df = data.frame(sizeFactors(dds))
colnames(sf_df) = c("DESeq2_size_factor")
sf_df$sample_id = rownames(sf_df)
reads_qc_table_merged = merge(reads_qc_table, sf_df, by="sample_id")
rownames(reads_qc_table_merged) = reads_qc_table_merged$sample_id
fwrite(reads_qc_table_merged, "kallisto_aligned_reads_qc.csv", quote=F, row.names=T)




# 
# 
# ### only >70% mapping, complete and ILLUMINA
# # sa_passed_filter = subset(sample_annotation, p_pseudoaligned_over_70==TRUE & failed_complete_download==FALSE & instrument_platform=="ILLUMINA")
# # gene_matrix_filter = gene_matrix %>% select(any_of(sa_passed_filter$sample_id))
# 
# # gene_matrix_filter = gene_matrix %>% select(any_of(trans_matrix_samples))  ### select only the one from tpm file
# # print(">70% mapping, complete and ILLUMINA remaining (selected from TPM-file): ")
# # print(dim(gene_matrix_filter))
# 
# 
# ### remove genes with 0 reads detected
# gene_matrix_filter = gene_matrix
# gene_matrix_filter = gene_matrix_filter[rowSums(gene_matrix_filter)>0,]
# print("removed genes with all reads=0: ")
# print(dim(gene_matrix_filter))
# 
# 
# ### remove duplicated gene names with version
# # gene_matrix_filter = unique(gene_matrix_filter)
# gene_matrix_filter = gene_matrix_filter[!duplicated(as.matrix(gene_matrix_filter)),]
# print("duplicated gene names with same version removed: ")
# print(dim(gene_matrix_filter))
# 
# 
# ### remove genes which have identical transcripts - double mapping bias
# # merged_unique_trans_gene = merge(unique_trans_table, tr2g_matrix[,c("transcript_id","gene_id")])
# # gene_matrix_filter = gene_matrix_filter[rownames(gene_matrix_filter) %in% unique(merged_unique_trans_gene$gene_id),]
# # print("genes with identical transcripts removed: ")
# # print(dim(gene_matrix_filter))
# 
# 
# ### only keep genes located on the chromosomes without scaffold  
# # tr2g_matrix_split = tr2g_matrix %>% separate(chromosome, sep=":", into=c("chr_scaffold","chr_genome","chr_num","chr_start","chr_end","chr_direction"), fill="left")
# # CHROMOSOME_NAMES = c(as.character(1:22), "X","Y","MT")  
# # genes_on_chromosomes = subset(tr2g_matrix_split, is.na(chr_scaffold) & chr_num %in% CHROMOSOME_NAMES)$gene_id
# # gene_matrix_filter = gene_matrix_filter[rownames(gene_matrix_filter) %in% unique(genes_on_chromosomes),]
# # print("genes which are not located on chromosomes: ")
# # print(dim(gene_matrix_filter))
# 
# 
# ### sort genes to newest version first, remove .version and remove duplicated gene names while keeping newest
# gene_sorted_order = order(nchar(rownames(gene_matrix_filter)), rownames(gene_matrix_filter), decreasing = TRUE)   # to correctly handle .1 .2 .10
# gene_matrix_filter = gene_matrix_filter[gene_sorted_order,]
# 
# gene_names_wo_version <- gsub("\\.[0-9]*$", "", rownames(gene_matrix_filter))
# rownames(gene_matrix_filter) = gene_names_wo_version
# 
# gene_matrix_filter = gene_matrix_filter[!duplicated(as.matrix(gene_matrix_filter)),]
# # gene_matrix_filter = unique(gene_matrix_filter)
# print("duplicated gene names removed: ")
# print(dim(gene_matrix_filter))
# 
# 
# ### remove genes which are only expressed by one sample - std dev problem  - TODO CHECK AGAIN ?
# # apply(gene_matrix_filter, 1, sd)
# expressed_samples_per_gene = rowSums(gene_matrix_filter!=0)
# 
# # ### check hist to check how many to remove
# # par(mfrow=c(1,2))
# # hist(expressed_samples_per_gene, breaks=100, main="samples expressed per gene")
# # hist(expressed_samples_per_gene[expressed_samples_per_gene<500], breaks= 100, xlim=c(0,500), main="detailed <500")
# # hist(expressed_samples_per_gene[expressed_samples_per_gene<100], breaks= 20, xlim=c(0,100), main="detailed <100")
# # sum(expressed_samples_per_gene<50)
# 
# expressed_samples_per_gene_min_quant = quantile(expressed_samples_per_gene, c(0.01) )
# print("removed everything under 1% expressed")
# print(expressed_samples_per_gene_min_quant)
# 
# gene_matrix_filter = gene_matrix_filter[expressed_samples_per_gene>=expressed_samples_per_gene_min_quant,]   # minimum 2 samples have expressed this gene - one enough ?
# print("removed genes with only one sample expressed: ")
# print(dim(gene_matrix_filter))
# 
# 
# 
# 
# ### remove highly correlated "identical" samples
# # cor_samples = cor(gene_matrix_filter)
# # cor_samples_failed <- findCorrelation(cor_samples, cutoff=0.9999, names=TRUE)
# # gene_matrix_filter = gene_matrix_filter[!colnames(gene_matrix_filter) %in% cor_samples_failed]
# # print("highly correlated samples removed: ")
# # print(dim(gene_matrix_filter))
# 
# 
# 
# ############################################
# ### DESeq2 normalization
# dummy_colData = data.frame(rep(1, ncol(gene_matrix_filter)), row.names=colnames(gene_matrix_filter))
# dds = DESeqDataSetFromMatrix(round(gene_matrix_filter), dummy_colData, design = ~1)  # no covariates included
# dds <- estimateSizeFactors(dds)
# dds <- estimateDispersions(dds)
# gene_matrix_filter_norm = counts(dds, normalized=T)
# print("finished DESeq2 size factors")
# 
# fwrite(data.frame(gene_matrix_filter_norm), file.path(WORKING_DIR, "processed/rnaseq_counts_filtered_samples_deseq2.csv"), quote=F, row.names=T) 
# 
# gene_matrix_filter_vst = assay(varianceStabilizingTransformation(dds))
# print("finished vst normalization")
# 
# 
# ############################################
# ### export
# 
# fwrite(data.frame(gene_matrix_filter_vst), file.path(WORKING_DIR, "processed/rnaseq_counts_filtered_samples_vst_ensembl.csv"), quote=F, row.names=T) 
# 
# 
# 
# 
# ############################################
# ### ENSG to entrez ids
# 
# hgnc_table = data.frame(fread(HGNC_TABLE_FILE), stringsAsFactors = FALSE)
# 
# ### if HGCNC curated id is missing, take the one currently taken from ensembl/NCBI
# hgnc_table_imputed = hgnc_table %>%
#   mutate(Ensembl.gene.ID = ifelse(Ensembl.gene.ID == "", Ensembl.ID.supplied.by.Ensembl., Ensembl.gene.ID)) %>%
#   mutate(NCBI.Gene.ID = ifelse(is.na(NCBI.Gene.ID)| NCBI.Gene.ID=="", NCBI.Gene.ID.supplied.by.NCBI., NCBI.Gene.ID))
# 
# print("amount of ENSG ids not found in hgnc table: ")
# print(length(intersect(rownames(gene_matrix_filter), hgnc_table_imputed$Ensembl.gene.ID)))
# print(length(intersect(rownames(gene_matrix_filter), hgnc_table$Ensembl.gene.ID)))
# 
# ### create common table
# tr2g_matrix_imputed = tr2g_matrix[,c("gene_id","gene_biotype")]
# tr2g_matrix_imputed = tr2g_matrix_imputed[order(nchar(tr2g_matrix_imputed$gene_id), tr2g_matrix_imputed$gene_id, decreasing = TRUE),]   # to correctly handle .1 .2 .10
# tr2g_matrix_imputed$gene_id = gsub("\\.[0-9]*$", "", tr2g_matrix_imputed$gene_id)
# tr2g_matrix_imputed = data.frame(tr2g_matrix_imputed[!duplicated(tr2g_matrix_imputed$gene_id),])
# rownames(tr2g_matrix_imputed) = tr2g_matrix_imputed$gene_id 
# tr2g_matrix_imputed = tr2g_matrix_imputed[rownames(gene_matrix_filter),]
# 
# gene_hgnc_merged = merge(tr2g_matrix_imputed, hgnc_table_imputed[,c("HGNC.ID","Approved.symbol","Approved.name","Alias.symbols","NCBI.Gene.ID","Ensembl.gene.ID")],
#                          by.x="gene_id", by.y="Ensembl.gene.ID", all.x=TRUE)
# colnames(gene_hgnc_merged) = c("ensembl_id","gene_biotype","HGNC_id","HGNC_symbol","HGNC_name","alias_symbol","NCBI_id")
# 
# fwrite(gene_hgnc_merged, file.path(WORKING_DIR, "processed/gene_id_table_ensembl_ncbi.csv"), na="", quote=T, row.names=F) 
# 
# 
# 
# ### create sub-output with only entrez ids
# 
# ### check potential removed ids
# not_found_matrix = subset(gene_hgnc_merged, is.na(NCBI_id) | NCBI_id =="")
# print("gene_biotype removed: ")
# table(not_found_matrix$gene_biotype)
# 
# ncbi_id_table = subset(gene_hgnc_merged, !is.na(NCBI_id) | NCBI_id =="" )
# ncbi_gene_matrix_filter_vst = gene_matrix_filter_vst[ncbi_id_table$ensembl_id,]
# rownames(ncbi_gene_matrix_filter_vst) = ncbi_id_table$NCBI_id
# 
# print("ncbi_id table: ")
# dim(ncbi_gene_matrix_filter_vst)
# 
# fwrite(data.frame(ncbi_gene_matrix_filter_vst), file.path(WORKING_DIR, "processed/rnaseq_counts_filtered_samples_vst_ncbi.csv"), quote=F, row.names=T) 
# 
# 
# 
# 
# 
# 
# 
# 
# ### check log2 against vst for genes
# # library(ggplot2)
# # library(cowplot)
# # vst_matrix = data.frame(fread(file.path(WORKING_DIR,"processed/rnaseq_counts_filtered_samples_vst_ensembl.csv"), header=T), row.names=1, check.names = FALSE)
# # plot_matrix = t(vst_matrix)
# # count_matrix = t(data.frame(fread(file.path(WORKING_DIR,"processed/rnaseq_counts_filtered_samples_deseq2.csv"), header=T), row.names=1, check.names = FALSE))
# # 
# # pre(vst_matrix)
# # 
# # rowMeans(vst_matrix[1:50,])
# # 
# # gene_num = 15987
# # p1 = qplot(plot_matrix[,gene_num], geom="histogram") 
# # p2 = qplot(log2(count_matrix[,gene_num]+1), geom="histogram") 
# # plot_grid(p1, p2, labels = "AUTO")
# # 
# # vst_matrix <- NULL
# # plot_matrix <- NULL
# # count_matrix <- NULL
# 
# 
# 
# 
# 
# 
# 
# 
# ### TODO OUTPUT RAW GENE COUNTS AND VST_TRANSFORMED
# 
# 
# 

