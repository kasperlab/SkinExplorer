GRCh38v103_multipleID_genes = readRDS("D:/Study/KI/Projects/###/Analysis/useful_files/GRCh38v103_multipleID_genes.rds")
alias_gene_path = "D:/Study/KI/Projects/###/Analysis/useful_files/hgnc_complete_set.txt"
alias_gene_0 = data.table::fread(alias_gene_path, sep = "\t")
alias_gene_1 = as.matrix(as.data.frame(alias_gene_0[, c("symbol", "alias_symbol", "prev_symbol", "locus_group", "ensembl_gene_id")]))
alias_gene_1_main = alias_gene_1[alias_gene_1[, "locus_group"] %in% c("protein-coding gene", "non-coding RNA", "pseudogene"), ]
alias_gene_1_other = alias_gene_1[alias_gene_1[, "locus_group"] == "other", ]
alias_gene_list = list()
main_alias = c(unlist(strsplit(alias_gene_1_main[, "alias_symbol"], "|", fixed = T)), unlist(strsplit(alias_gene_1_main[, "prev_symbol"], "|", fixed = T)))
duplicated_main_alias = main_alias[duplicated(main_alias)]
for(ii in seq(nrow(alias_gene_1_main))){
  alias = setdiff(c(unlist(strsplit(alias_gene_1_main[ii, "alias_symbol"], "|", fixed = T)), unlist(strsplit(alias_gene_1_main[ii, "prev_symbol"], "|", fixed = T))), c(alias_gene_1_main[, "symbol"], duplicated_main_alias))
  if(length(alias) > 0){
    alias_gene_list[[alias_gene_1_main[ii, "symbol"]]] = alias
  }
  if(ii %% 1000 == 0){
    print(ii)
  }
}
all_genename = unique(c(names(alias_gene_list), unlist(alias_gene_list)))
for(ii in seq(nrow(alias_gene_1_other))){
  alias = setdiff(c(unlist(strsplit(alias_gene_1_other[ii, "alias_symbol"], "|", fixed = T)), unlist(strsplit(alias_gene_1_other[ii, "prev_symbol"], "|", fixed = T))), all_genename)
  if(length(alias) > 0){
    alias_gene_list[[alias_gene_1_other[ii, "symbol"]]] = alias
  }
  if(ii %% 1000 == 0){
    print(ii)
  }
}
#
from_c = c()
to_c = c()
for(ii in names(alias_gene_list)){
  this_alias = c(ii, alias_gene_list[[ii]])
  from_c = c(from_c, this_alias)
  to_c = c(to_c, rep(ii, length(this_alias)))
}
from_c_exclude = c(from_c[duplicated(from_c)], GRCh38v103_multipleID_genes)
from_c_mask = !from_c %in% from_c_exclude
from_c_filtered = from_c[from_c_mask]
from_c_filtered_seurat = gsub("_", "-", from_c_filtered, fixed = T)
mask_seurat_unique = !from_c_filtered_seurat %in% from_c_filtered
from_c_filtered_both = c(from_c_filtered, from_c_filtered_seurat[mask_seurat_unique])
to_c_filtered = to_c[from_c_mask]
to_c_filtered_both = c(to_c_filtered, to_c_filtered[mask_seurat_unique])
to_c_filtered_seurat = gsub("_", "-", to_c_filtered, fixed = T)
to_c_filtered_seurat_both = c(to_c_filtered_seurat, to_c_filtered_seurat[mask_seurat_unique])
names(to_c_filtered_both) = from_c_filtered_both
names(to_c_filtered_seurat_both) = from_c_filtered_both
saveRDS(to_c_filtered_both, "D:/Study/KI/Projects/###/Analysis/useful_files/hgnc_mapping_alias_genename.rds")
saveRDS(to_c_filtered_seurat_both, "D:/Study/KI/Projects/###/Analysis/useful_files/hgnc_mapping_alias_genename_seurat.rds")
#
# Ensembl Gene ID
#
ensembl_gene_id_for_alias_gene_list = list()
for(ii in seq(nrow(alias_gene_1))){
  gene_names = alias_gene_1[ii, "symbol"]
  if(length(gene_names) > 0){
    ensembl_gene_id_for_alias_gene_list[[alias_gene_1[ii, "ensembl_gene_id"]]] = gene_names
  }
  if(ii %% 1000 == 0){
    print(ii)
  }
}
from_c = c()
to_c = c()
for(ii in names(ensembl_gene_id_for_alias_gene_list)){
  these_gene_names = c(ii, ensembl_gene_id_for_alias_gene_list[[ii]])
  from_c = c(from_c, these_gene_names)
  to_c = c(to_c, rep(ii, length(these_gene_names)))
}
from_c_exclude = from_c[duplicated(from_c)]
from_c_mask = !from_c %in% from_c_exclude & from_c != ""
to_c_filtered = to_c[from_c_mask]
from_c_filtered = from_c[from_c_mask]
names(to_c_filtered) = from_c_filtered
# Seurat case
from_c_filtered_seurat = gsub("_", "-", from_c_filtered, fixed = T)
to_c_filtered_seurat = gsub("_", "-", to_c_filtered, fixed = T)
names(to_c_filtered_seurat) = from_c_filtered_seurat
#
to_c_filtered = c(to_c_filtered, to_c_filtered_seurat[!names(to_c_filtered_seurat) %in% names(to_c_filtered)])
saveRDS(to_c_filtered, "D:/Study/KI/Projects/###/Analysis/useful_files/hgnc_mapping_genename_ensemblID.rds")
#
# Gene ID -> Gene Name
#
ensembl_gene_id_c = as.character(alias_gene_1[, "ensembl_gene_id"])
symbol_c = as.character(alias_gene_1[, "symbol"])
ensembl_gene_id_valid = ensembl_gene_id_c != ""
from_c = ensembl_gene_id_c[ensembl_gene_id_valid]
to_c = symbol_c[ensembl_gene_id_valid]
names(to_c) = from_c
saveRDS(to_c, "D:/Study/KI/Projects/###/Analysis/useful_files/hgnc_mapping_ensemblID_genename.rds")
to_c_seurat = gsub("_", "-", to_c, fixed = T)
saveRDS(to_c_seurat, "D:/Study/KI/Projects/###/Analysis/useful_files/hgnc_mapping_ensemblID_genename_seurat.rds")




