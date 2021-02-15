#!/usr/bin/env Rscript

#' ---
#' title: "gene count performance analysis"
#' author: "loipf"
#'
#' output: html_document
#' code_folding: hide
#' ---

pacman::p_load(knitr, data.table, plotly, ggplot2, pheatmap, vsn, hexbin, LSD, cowplot)


# ### DOES NOT WORK WITH NEXTFLOW# args = commandArgs(trailingOnly=TRUE)
# READS_QC_TABLE_FILE = args[1]
# GENE_MATRIX_FILE = args[2]
# GENE_MATRIX_FILE_VST = args[3]

READS_QC_TABLE_FILE = "/home/stefan/Documents/umcg/RNAseq-pipeline/data/kallisto_aligned_reads_qc.csv"
GENE_MATRIX_FILE = "/home/stefan/Documents/umcg/RNAseq-pipeline/data/kallisto_gene_counts.csv"
GENE_MATRIX_FILE_VST = "/home/stefan/Documents/umcg/RNAseq-pipeline/data/kallisto_gene_counts_norm_sf_vst.csv"

READS_QC_TABLE_FILE = "kallisto_aligned_reads_qc.csv"
GENE_MATRIX_FILE = "kallisto_gene_counts.csv"
GENE_MATRIX_FILE_VST = "kallisto_gene_counts_norm_sf_vst.csv"

### read in everything
reads_qc_table =  data.frame(fread(READS_QC_TABLE_FILE, header = T), stringsAsFactors = FALSE, row.names = 1)
gene_matrix =  data.frame(fread(GENE_MATRIX_FILE, header = T), stringsAsFactors = FALSE, row.names = 1)  ### genes x samples
gene_matrix_vst =  data.frame(fread(GENE_MATRIX_FILE_VST, header = T), stringsAsFactors = FALSE, row.names = 1)
gene_matrix_log2 = log2(gene_matrix+1)  ### vst has sf norm, log2 not


############################################
### functions

center_rowwise <- function(x, na.rm=T) {
  x_center = x - rowMeans(x, na.rm=na.rm)
  return(x_center)
}

plot_heatscatter <- function(x, y, ...) {
  plot_df = data.frame(x=as.vector(x), y=as.vector(y))
  plot_df = na.omit(plot_df)

  LSD::heatscatter(plot_df$x, plot_df$y, cor=TRUE, colpal="standard", ...)
  abline(0,1, col="orange")
}

plot_density_genes <- function(gene_matrix, main="", xlab="values") {
  plot_df = stack(gene_matrix)
  colnames(plot_df)[2] = "sample_id"
  p <- ggplot(plot_df, aes(x=values)) + labs(x=xlab, title=main) +
    stat_density(aes(color=sample_id),position="identity",geom="line") +
    theme_bw() + theme(legend.position='none')
  p
}

plot_density_genes_cumulative <- function(gene_matrix, main="", xlab="values") {
  plot_df_list = sapply(colnames(gene_matrix), function(sample_id) {
    x_value = sort(gene_matrix[,sample_id])
    cumul_value = cumsum(x_value)/sum(x_value)
    data.frame(x=x_value, y=cumul_value, sample_id=sample_id)
  }, simplify = F, USE.NAMES = F)
  plot_df = do.call("rbind", plot_df_list)

  p = ggplot(plot_df, aes(x=x, y=y, color=sample_id)) + geom_line() +
    theme_bw() + labs(y="density", x=xlab, title=main) + theme(legend.position='none')
  p
}



#' ### kallisto performance
mean_kallisto_p_pseudoaligned = mean(reads_qc_table$kallisto_p_pseudoaligned)
p = ggplot(reads_qc_table, aes(x = kallisto_p_pseudoaligned, label=sample_id)) +
  geom_rug(sides="b") + geom_line(stat="density") +
  theme_bw() + labs(x="kallisto pseudoaligned [%]", title="kallisto mapping performance") +
  geom_vline(aes(xintercept = mean_kallisto_p_pseudoaligned),col='orange',size=0.5) +
  geom_text(aes(label=paste0("mean: ",round(mean_kallisto_p_pseudoaligned,4)),y=0,x=ifelse(mean_kallisto_p_pseudoaligned > 50, mean_kallisto_p_pseudoaligned-10, mean_kallisto_p_pseudoaligned+10)), vjust=-1,col='orange',size=4) +
  xlim(0,100)
ggplotly(p)


#' ### DESeq2 size factors (should be around 0.5-1.5)
p = ggplot(reads_qc_table, aes(x = DESeq2_size_factor, label=sample_id, label2=kallisto_p_pseudoaligned)) +
  geom_rug(sides="b") + geom_line(stat="density") +
  theme_bw() + labs(x="size factors", title="DESeq2 size factor estimation") +
  geom_vline(aes(xintercept = 1),col='orange',size=0.5)
ggplotly(p)


#' ### sample correlation
COLOR_PALETTE <- c("#2166AC","#67A9CF","#D1E5F0","#F7F7F7","#FDDBC7","#EF8A62","#B2182B") ## rev(brewer.pal(n = 7, name = "RdBu"))
gene_matrix_vst_center = center_rowwise(gene_matrix_vst)
corr_matrix = cor(gene_matrix_vst_center)
#+ fig.width=10, fig.height=8
pheatmap(corr_matrix, main="sample correlation gene_counts",
         breaks = seq(-1,1, length.out = 100), na_col = "grey", color=colorRampPalette(COLOR_PALETTE)(100) )


#' ### PCA
pca_obj = prcomp(gene_matrix_vst, center=T, scale.=T)
pca_df = data.frame(pca_obj$rotation, sample_id = colnames(gene_matrix_vst))
pca_explained_var = signif(pca_obj$sdev^2/sum(pca_obj$sdev^2)*100,3)

p = ggplot(pca_df, aes(x=PC1, y=PC2, label=sample_id)) + geom_point() + theme_bw() +
  labs(x=paste0("PC1 [",pca_explained_var[1],"% variance]"), y=paste0("PC2 [",pca_explained_var[2],"% variance]"), title="PCA samples")
ggplotly(p)


#' ### gene count distribution
boxplot(gene_matrix_log2, main="gene count distribution", ylab="log2(expression+1)")
boxplot(gene_matrix_vst, main="gene count distribution [sf+vst]", ylab="normalized expression")


#' ### compare first 2 samples
plot_heatscatter(gene_matrix_log2[,1], gene_matrix_log2[,2], main="gene count sample comparison [log2(expression+1)]",
                 xlab=colnames(gene_matrix_log2)[1], ylab=colnames(gene_matrix_log2)[2])
plot_heatscatter(gene_matrix_vst[,1], gene_matrix_vst[,2], main="gene count sample comparison [sf + vst]",
                 xlab=colnames(gene_matrix_log2)[1], ylab=colnames(gene_matrix_log2)[2])


#' ### raw to vst heteroscedasticity
#+ fig.width=10, fig.height=5
p1 = meanSdPlot(as.matrix(gene_matrix), plot=F)$gg + ggtitle("raw") + theme_bw()
p2 = meanSdPlot(as.matrix(gene_matrix_log2), plot=F)$gg + ggtitle("log2 normalized") + theme_bw()
p3 = meanSdPlot(as.matrix(gene_matrix_vst), plot=F)$gg + ggtitle("sf+vst normalized") + theme_bw()
p_grid = plot_grid(p1,p2,p3, labels="AUTO", ncol=3)
title <- ggdraw() + draw_label("gene counts mean-sd", fontface='bold')
plot_grid(title, p_grid, ncol=1, rel_heights=c(0.1, 1))


#' ### gene count density plots
ggplotly(plot_density_genes(gene_matrix_log2, main="gene count densities [raw]", xlab="log2(expression+1)"))
ggplotly(plot_density_genes(gene_matrix_vst, main="gene count densities [sf+vst]", xlab="sf+vst expression"))

ggplotly(plot_density_genes_cumulative(gene_matrix_log2, main="gene count densities [raw]", xlab="log2(expression+1)"))
ggplotly(plot_density_genes_cumulative(gene_matrix_vst, main="gene count densities [sf+vst]", xlab="sf+vst expression"))




