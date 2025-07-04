library(data.table)
gtf_path = "D:/Study/KI/Projects/###/Analysis/useful_files/gencode.v37lift37.annotation.gtf.gz"
#
#
gtf_file = fread(gtf_path, sep = "\t")
unique_link = "---ThisIsAUniqueLink---"
unique_ID_Name = sapply(unique(sapply(as.matrix(as.data.frame(gtf_file[, 9])), function(x){
  tmp = unlist(strsplit(x, 'gene_id \"', fixed = T))[2]
  GeneID = unlist(strsplit(tmp, '\"', fixed = T))[1]
  tmp = unlist(strsplit(x, 'gene_name \"', fixed = T))[2]
  GeneName = unlist(strsplit(tmp, '\"', fixed = T))[1]
  return(paste0(GeneID, unique_link, GeneName))
})), function(x){
  str_fragments = unlist(strsplit(x, unique_link, fixed = T))
  return(c(str_fragments[1], str_fragments[2]))
})
mapping_ID_Name = unique_ID_Name[2, ]
names(mapping_ID_Name) = unique_ID_Name[1, ]
saveRDS(mapping_ID_Name, "D:/Study/KI/Projects/###/Analysis/useful_files/GRCh37_mapping_ensemblID_genename.rds")
mask_unique_Name = !duplicated(unique_ID_Name[2, ])
mapping_Name_ID = unique_ID_Name[1, mask_unique_Name]
names(mapping_Name_ID) = unique_ID_Name[2, mask_unique_Name]
saveRDS(mapping_Name_ID, "D:/Study/KI/Projects/###/Analysis/useful_files/GRCh37_mapping_genename_ensemblID.rds")
multipleID_genes = unique(unique_ID_Name[2, !mask_unique_Name])
saveRDS(multipleID_genes, "D:/Study/KI/Projects/###/Analysis/useful_files/GRCh37_multipleID_genes.rds")
