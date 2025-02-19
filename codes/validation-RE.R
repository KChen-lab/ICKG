library(caret)
# library(KEGGgraph)
# library(KEGGREST)
# library(dplyr)
# library(reactome.db)
library(org.Hs.eg.db)
# library(graphite)

#################################### (-----------------------------------------) ####################################
#################################### Relationship between genes and gene ####################################
#################################### KEGG ####################################

############ (1) Load the kegg database ############

pathway_ids <- keggList("pathway",organism='hsa')
pathway_ids = names(pathway_ids)

for(i in 1:length(pathway_ids)){
  print(i)
  pathway_id <- pathway_ids[i]

  kgml_url <- paste0("https://rest.kegg.jp/get/", pathway_id, "/kgml")
  kgml_file <- paste0("/Volumes/she4/knowledgegraph/data/validation/kegg-relationship/",pathway_id, ".kgml")
  download.file(kgml_url, kgml_file)
  kgml_graph <- parseKGML(kgml_file)
  edges <- edges(kgml_graph)
  interactions <- data.frame(
    source = sapply(edges, function(edge) edge@entry1ID),
    target = sapply(edges, function(edge) edge@entry2ID),
    interaction_type = sapply(edges, function(edge) {
      if (length(edge@subtype) > 0) {
        # Extract subtype names, combine multiple subtypes with commas if they exist
        paste(sapply(edge@subtype, function(subtype) subtype@name), collapse = ", ")
      } else {
        # Assign NA if there is no subtype (no relationship defined)
        NA
      }
    })
  )

  # Filter out rows where interaction_type is NA (meaning no relationship)
  interactions <- interactions[!is.na(interactions$interaction_type), ]

  # Step 2: Extract the node information from the parsed KGML
  entries <- nodes(kgml_graph)

  # Step 3: Create a lookup table for node IDs and the associated genes
  node_gene_lookup <- sapply(entries, function(entry) {
    entry_id <- entry@entryID  # KEGG node ID
    entry_name <- entry@name  # KEGG gene names (comma-separated)
    gene_symbols <- gsub("hsa:", "", entry_name)  # Remove the 'hsa:' prefix to keep Entrez IDs
    return(gene_symbols)
  })

  # Convert the list into a named vector for easy lookup
  node_gene_lookup <- setNames(node_gene_lookup, sapply(entries, function(entry) entry@entryID))

  # Step 4: Map the source and target nodes to their respective genes
  map_node_to_genes <- function(node_id) {
    if (node_id %in% names(node_gene_lookup)) {
      return(node_gene_lookup[[node_id]])  # Return gene IDs for the node
    } else {
      return(NA)  # Return NA if no match found
    }
  }

  # Apply the mapping to the source and target columns
  interactions$source_gene <- lapply(interactions$source, map_node_to_genes)
  interactions$target_gene <- lapply(interactions$target, map_node_to_genes)
  # interactions$source_gene <- as.character(interactions$source_gene)
  # interactions$target_gene <- as.character(interactions$target_gene)

  # Step 5: Split the source and target genes into individual rows
  expanded_interactions <- do.call(rbind, lapply(1:nrow(interactions), function(i) {
    # Split the source and target genes into individual vectors
    source_genes <- unlist(interactions$source_gene[i])
    target_genes <- unlist(interactions$target_gene[i])

    # Create all possible source-target pairs
    expand.grid(source_gene = source_genes, target_gene = target_genes,
                interaction_type = interactions$interaction_type[i])
  }))

  # Step 6: Convert Entrez IDs to gene symbols using org.Hs.eg.db
  expanded_interactions$source_gene = as.character(expanded_interactions$source_gene)
  expanded_interactions$target_gene = as.character(expanded_interactions$target_gene)
  # Convert Entrez IDs to gene symbols for the source and target genes
  expanded_interactions$source_gene_symbol <- sapply(expanded_interactions$source_gene, function(id) {
    if (id %in% valid_ids) {
      mapIds(org.Hs.eg.db, keys = id, column = "SYMBOL", keytype = "ENTREZID", multiVals = "first")
    } else {
      NA  # Return NA if not a valid ENTREZID
    }
  })
  expanded_interactions$target_gene_symbol <- sapply(expanded_interactions$target_gene, function(id) {
    if (id %in% valid_ids) {
      mapIds(org.Hs.eg.db, keys = id, column = "SYMBOL", keytype = "ENTREZID", multiVals = "first")
    } else {
      NA  # Return NA if not a valid ENTREZID
    }
  })

  # Optional: Save the result to a CSV file
  write.csv(expanded_interactions,
            paste0("/Volumes/she4/knowledgegraph/data/validation/kegg-relationship/",pathway_id,"_expanded_interactions_with_gene_symbols.csv"),
            row.names = FALSE)

}

csv_files <- list.files("/Volumes/she4/knowledgegraph/data/validation/kegg-relationship/",pattern = "*symbols.csv",full.names = T)
combined_data <- data.frame()
for (file in csv_files) {
  data <- read.csv(file)
  if(nrow(data)==0){next}
  selected_data <- data[,c("source_gene_symbol", "target_gene_symbol", "interaction_type")]
  combined_data <- rbind(combined_data, selected_data)
}
# Remove duplicate rows and rows with any NA values
final_data <- combined_data %>%
  distinct() %>%
  na.omit()
# Check the first few rows of the final combined data
# final_data$interaction_group <- ifelse(
#   grepl("activation|phosphorylation|binding|association|expression|ubiquitination|state change|compound",
#         final_data$interaction_type, ignore.case = TRUE), "activate",
#   ifelse(grepl("inhibition|repression|dephosphorylation|dissociation|indirect effect|missing interaction",
#                final_data$interaction_type, ignore.case = TRUE), "inhibit",
#          "other")
# )
final_data$interaction_group <- ifelse(final_data$interaction_type%in%c("activation","activation, indirect",
                                                                        "activation, indirect effect"),"activate",
                                       ifelse(final_data$interaction_type%in%c("inhibition","inhibition, indirect", "inhibition, indirect effect", "inhibition, repression", "repression, indirect effect","repression"),"inhibit","other"))
write.csv(final_data, "/Volumes/she4/knowledgegraph/data/validation/kegg-relationship/combined_kegg_relationship.csv",
          row.names = FALSE)
final_data_cut = final_data[which(final_data$interaction_group!="other"),]
write.csv(final_data_cut, "/Volumes/she4/knowledgegraph/data/validation/kegg-relationship/combined_kegg_relationship_cut.csv",
          row.names = FALSE)

############ (2) Load the constructed KG ############
final_data = read.csv("/rsrch4/home/bcb/she4/combined_kegg_relationship.csv")
KG_nodes = read.csv("/rsrch4/home/bcb/she4/knowledgegraph/result/neo4j_normalized/M_nodes_pmid_gene_normalized.csv") # T
KG_nodes = KG_nodes[which(KG_nodes$node_type=="gene"),]
KG = read.csv("/rsrch4/home/bcb/she4/knowledgegraph/result/neo4j_normalized/M_relationships_pmid_gene_normalized.csv")
KG = KG[which(KG$X.START_ID%in%KG_nodes$name&KG$X.END_ID%in%KG_nodes$name),]

# cd_to_gene_name <- function(cd_genes) {
#   symbols <- mapIds(org.Hs.eg.db,
#                     keys = cd_genes,
#                     column = "SYMBOL",
#                     keytype = "ALIAS",
#                     multiVals = "first")
#   return(symbols)
# }
# KG$X.START_ID = toupper(KG$X.START_ID)
# KG$X.END_ID = toupper(KG$X.END_ID)
# KG$X.START_ID <- ifelse(is.na(unname(cd_to_gene_name(KG$X.START_ID))),KG$X.START_ID,unname(cd_to_gene_name(KG$X.START_ID)))
# KG$X.END_ID <- ifelse(is.na(unname(cd_to_gene_name(KG$X.END_ID))),KG$X.END_ID,unname(cd_to_gene_name(KG$X.END_ID)))

############ (3) Find the overlap and calculate accuracy (first need to find all synonyms in the groundtruth database) ############

get_gene_synonyms <- function(gene_symbol) {
  tryCatch({
    synonyms <- AnnotationDbi::select(org.Hs.eg.db, 
                                      keys = gene_symbol,
                                      columns = c("SYMBOL", "ALIAS"),
                                      keytype = "SYMBOL")
    unique(c(synonyms$SYMBOL, synonyms$ALIAS)[!is.na(c(synonyms$SYMBOL, synonyms$ALIAS))])
  }, error = function(e) {
    if (grepl("None of the keys entered are valid keys for 'SYMBOL'", e$message)) {
      warning(paste("Invalid SYMBOL key:", gene_symbol, ". Returning empty vector."))
      return(gene_symbol)  # Return an empty vector
    } else {
      stop(e)  # Re-throw the error if it's not the one we're expecting
    }
  })
}

final_data_source_gene_synonym_list = readRDS("/rsrch4/home/bcb/she4/knowledgegraph/data/validation/kegg-relationship/source_gene_synonyms.rds")
final_data_target_gene_synonym_list = readRDS("/rsrch4/home/bcb/she4/knowledgegraph/data/validation/kegg-relationship/target_gene_synonyms.rds")

ordered_final_data = NULL
order_KG = NULL
for(i in 1:nrow(KG)){
  start = NULL
  for(r in 1:length(final_data_source_gene_synonym_list)){
    start = c(start,sum(grepl(final_data_source_gene_synonym_list[[r]],pattern = KG$X.START_ID[i])))
  }
  start_index = which(start!=0) ## these two rows in the ground truth have matching source genes
  
  end = NULL
  for(r in 1:length(final_data_target_gene_synonym_list)){
    end = c(end,sum(grepl(final_data_target_gene_synonym_list[[r]],pattern = KG$X.END_ID[i])))
  }
  end_index = which(end!=0) ## these two rows in the ground truth have matching target genes
  
  row = intersect(start_index,end_index)
  if(length(row)==0){next}
  if(length(row)>1){row = row[1]}
  order_KG = rbind(order_KG,KG[i,])
  ordered_final_data = rbind(ordered_final_data,final_data[row,])
  print(i)
}
print(dim(ordered_final_data))
print(dim(order_KG))

xtab_full <- table(factor(ordered_final_data$interaction_group, levels=c("activate", "inhibit")),
                   factor(order_KG$relationship, levels=c("activate", "inhibit")))
confusionMatrix(xtab_full)
tab = confusionMatrix(xtab_full)
print(tab)
print(fisher.test(tab$table,alternative = "two.sided"))

