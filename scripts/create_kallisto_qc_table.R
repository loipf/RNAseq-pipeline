############################################
### check mapping percentage reads, creates table with filtered samples for low amount of reads aligned

library(data.table)
library(rjson)


PER_READS_ALIGNED = 70   ### percentage quality threshold

############################################
### read all jsons with kallisto stats

json_file_list = list.files(".", pattern="*.json", full.names = T)

get_json_sample_list <- function(json_file_sample) {
  sample_id = strsplit(basename(json_file_sample), "_")[[1]][1]
  sample_json = fromJSON(file = json_file_sample)
  return(c(sample_id,sample_json$n_processed, sample_json$n_pseudoaligned, sample_json$p_pseudoaligned, sample_json$n_unique, sample_json$p_unique))
}

qc_sample_list = lapply(json_file_list, get_json_sample_list) ### could be parallelized but not necessary
qc_sample_df = data.frame(do.call(rbind, qc_sample_list), stringsAsFactors = F)
colnames(qc_sample_df) = c("sample_id", "kallisto_n_processed", "kallisto_n_pseudoaligned", "kallisto_p_pseudoaligned", "kallisto_n_unique", "kallisto_p_unique")

qc_sample_df[[paste0("p_pseudoaligned_over_", PER_READS_ALIGNED)]] = qc_sample_df$kallisto_p_pseudoaligned >= PER_READS_ALIGNED

write.csv(qc_sample_df, "kallisto_aligned_reads_qc.csv", quote=FALSE, row.names = FALSE)




