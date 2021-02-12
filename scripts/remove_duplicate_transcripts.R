
############################################
### remove transcripts with 100% sequence identiy

# setwd("/home/stefan/Documents/umcg/RNAseq-pipeline/data/kallisto_index")

library(data.table)

args = commandArgs(trailingOnly=TRUE)

oneline_filepath = args[1]

############################################
### read in
# trans_table <- fread("Homo_sapiens.GRCh38.cdna_ncrna_oneline.txt", header=F)
trans_table <- fread(oneline_filepath, header=F)

num_identical_transcripts = nrow(trans_table[,.N,by=V2][N!=1])

### pick first one of duplicates 
trans_table_short = trans_table[, .SD[1,], by = V2]
fwrite(trans_table_short[, c(2,1)], "Homo_sapiens.GRCh38.cdna_ncrna_oneline_unique.txt", sep="\t", col.names = F) 


### write numbers to output
sink('transcript_removal_info.txt')

cat("### remove 100% sequence identity transcripts to avoid gene mapping bias\n")
cat("number of unique identical transcripts found: ",num_identical_transcripts,"\n")
cat("transcripts before: ",nrow(trans_table),"\n")
cat("transcripts after removal: ", nrow(trans_table_short),"\n")

sink()














