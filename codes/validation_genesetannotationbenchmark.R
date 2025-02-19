library(readxl)
library(msigdbr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(dplyr)


process_excel <- function(file_path) {
  sheets <- excel_sheets(file_path)
  gene_sets <- list()
  for (sheet in sheets) {
    df <- read_excel(file_path, sheet = sheet)
    for (column in names(df)) {
      gene_set_name <- paste(sheet, column, sep = "_")
      genes <- df[[column]][!is.na(df[[column]])]
      gene_sets[[gene_set_name]] <- genes
    }
  }
  return(gene_sets)
}
file_path <- "/Volumes/she4/knowledgegraph/data/validation/tirosh_MPs.xlsx"
result <- process_excel(file_path)
result_T = result[87:108]
result_B = result[1:12]
result_M = result[74:86]
result_tcr = result[109:112]
result_nk = result[113:116]

msigdb_genesets <- msigdbr(species = "Homo sapiens",category="C2") 
msigdb_genesetsH <- msigdbr(species = "Homo sapiens",category="H") 
msigdb_genesets = rbind(msigdb_genesets,msigdb_genesetsH)
msigdb_genesets = msigdb_genesets[which(msigdb_genesets$gs_subcat%in%c("","CP:REACTOME","CP:KEGG","CP:WIKIPATHWAYS")),]
length(unique(msigdb_genesets$gs_name))

universe_entrez <- mapIds(org.Hs.eg.db, keys = as.character(msigdb_genesets$entrez_gene), column = "SYMBOL", keytype = "ENTREZID")
msigdb_genesets$symbol = unname(universe_entrez)

res_T = list()
for(i in 1:length(result_T)){
  res = enricher(gene = unlist(result_T[i]),
                 universe = NULL,
                 TERM2GENE=msigdb_genesets[,c(3,4)],
                 pvalueCutoff = 0.05)
  res_T[[i]] = res[which(res@result$p.adjust<0.05),]
  print(i)
}
names(res_T) = names(result_T)
saveRDS(res_T,"/Volumes/she4/knowledgegraph/result/tirosh_T_standardannotation_small.rds")

res_B = list()
for(i in 1:length(result_B)){
  res = enricher(gene = unlist(result_B[i]),
                 universe = NULL,
                 TERM2GENE=msigdb_genesets[,c(3,4)],
                 pvalueCutoff = 0.05)
  res_B[[i]] = res[which(res@result$p.adjust<0.05),]
  print(i)
}
names(res_B) = names(result_B)
saveRDS(res_B,"/Volumes/she4/knowledgegraph/result/tirosh_B_standardannotation_small.rds")

res_M = list()
for(i in 1:length(result_M)){
  res = enricher(gene = unlist(result_M[i]),
                 universe = NULL,
                 TERM2GENE=msigdb_genesets[,c(3,4)],
                 pvalueCutoff = 0.05)
  res_M[[i]] = res[which(res@result$p.adjust<0.05),]
  print(i)
}
names(res_M) = names(result_M)
saveRDS(res_M,"/Volumes/she4/knowledgegraph/result/tirosh_M_standardannotation_small.rds")

res_tcr = list()
for(i in 1:length(result_tcr)){
  res = enricher(gene = unlist(result_tcr[i]),
                 universe = NULL,
                 TERM2GENE=msigdb_genesets[,c(3,4)],
                 pvalueCutoff = 0.05)
  res_tcr[[i]] = res[which(res@result$p.adjust<0.05),]
  print(i)
}
names(res_tcr) = names(result_tcr)
saveRDS(res_tcr,"/Volumes/she4/knowledgegraph/result/xianli_tcr_standardannotation_small.rds")

res_nk = list()
for(i in 1:length(result_nk)){
  res = enricher(gene = unlist(result_nk[i]),
                 universe = NULL,
                 TERM2GENE=msigdb_genesets[,c(3,4)],
                 pvalueCutoff = 0.05)
  res_nk[[i]] = res[which(res@result$p.adjust<0.05),]
  print(i)
}
names(res_nk) = names(result_nk)
saveRDS(res_nk,"/Volumes/she4/knowledgegraph/result/nk_standardannotation_small.rds")

library(ggplot2)
library(reshape2)
count = read.csv("/Volumes/she4/knowledgegraph/result/genesetannotation.csv")
count = count[,c(1,2,4,5)]
count_long = melt(count,id.vars=c("Cell.type","Gene.set"))
names(count_long)[3:4] = c("Method","Count")

pdf("/Volumes/she4/knowledgegraph/result/genesetannotation_compared.pdf",height=11,width=6)
ggplot(count_long, aes(x = Gene.set, y = Count, fill = Method)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = c("ORA_small" = "#FF9999", "KG" = "#66C2A5")) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    panel.grid.major = element_line(color = "gray90"),
    panel.grid.minor = element_line(color = "gray95"),
    legend.position = "bottom",
    legend.box = "vertical",
    legend.margin = margin(t = 20),
    strip.text.y = element_text(angle = 0)  # Make facet labels horizontal
  ) +
  facet_grid(Cell.type ~ ., scales = "free", space = "free") +
  labs(
    x = "Gene Sets of different immune cell types",
    y = "# of enriched terms",
    title = "Enrichment results",
    fill = "Method"
  ) +
  coord_flip()  
dev.off()

count = read.csv("/Volumes/she4/knowledgegraph/result/genesetannotation.csv")
count = count[,c(1,2,6,7)]
count_long = melt(count,id.vars=c("Cell.type","Gene.set"))
names(count_long)[3:4] = c("Method","Length")
pdf("/Volumes/she4/knowledgegraph/result/genesetannotation_phrase_length_compared.pdf",height=11,width=6)
ggplot(count_long, aes(x = Gene.set, y = Length, fill = Method)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = c("ORA_small_phrase_length" = "#FF9999", "KG_phrase_length" = "#66C2A5")) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    panel.grid.major = element_line(color = "gray90"),
    panel.grid.minor = element_line(color = "gray95"),
    legend.position = "bottom",
    legend.box = "vertical",
    legend.margin = margin(t = 20),
    strip.text.y = element_text(angle = 0)  # Make facet labels horizontal
  ) +
  facet_grid(Cell.type ~ ., scales = "free", space = "free") +
  labs(
    x = "Gene Sets of different immune cell types",
    y = "Average annotation string length",
    title = "Enrichment results",
    fill = "Method"
  ) +
  coord_flip()
dev.off()

### Create Table 4
library(dplyr)
library(tidyr)
csv_files <- list.files("/Volumes/she4/knowledgegraph/result/",pattern = "^pagerank_.*\\.csv$",full.names = T)
processed_data <- list()
for (file in csv_files) {
  data <- read.csv(file)
  sorted_data <- data %>%
    arrange(desc(Score))
  pathway_col <- sorted_data$Pathway
  processed_data[[file]] <- pathway_col
}
max_length <- max(sapply(processed_data, length))
processed_data <- lapply(processed_data, function(x) {
  length(x) <- max_length
  return(x)
})
result_df <- as.data.frame(processed_data)
colnames(result_df) <- gsub("\\.csv$", "", basename(csv_files))
write.csv(result_df, "~/Desktop/combined_pathways.csv", row.names = FALSE)

### Create Table 5 - ORA annotation
annoB = readRDS("/Volumes/she4/knowledgegraph/result/tirosh_B_standardannotation_small.rds")
annoT = readRDS("/Volumes/she4/knowledgegraph/result/tirosh_T_standardannotation_small.rds")
annoM = readRDS("/Volumes/she4/knowledgegraph/result/tirosh_M_standardannotation_small.rds")
annoNK = readRDS("/Volumes/she4/knowledgegraph/result/nk_standardannotation_small.rds")
annoTCR = readRDS("/Volumes/she4/knowledgegraph/result/xianli_tcr_standardannotation_small.rds")
anno = c(annoB,annoT,annoM,annoNK,annoTCR)
list = list()
for(i in names(anno)){
  list[[i]] = anno[[i]]$Description
}
max_length <- max(sapply(list, length))
list <- lapply(list, function(x) {
  length(x) <- max_length
  return(x)
})
result_df <- as.data.frame(list)
colnames(result_df) <- names(anno)
write.csv(result_df, "~/desktop/ORA_genesetannotation_small.csv", row.names = FALSE)

### phrase length
char_counts <- lapply(res_tcr, function(df) {
  nchar(df$Description)
})
mean_lengths <- round(sapply(char_counts, mean),0)
mean_lengths

files <- list.files("/Volumes/she4/knowledgegraph/result/", pattern="pagerank", full.names=TRUE)
results <- sapply(files, function(f) {
  df <- read.csv(f)
  round(mean(nchar(df$Pathway)),0)
})

names(results) <- basename(files)
print(results)
