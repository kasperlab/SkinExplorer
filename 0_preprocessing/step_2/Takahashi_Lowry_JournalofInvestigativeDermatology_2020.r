r_utilities = "D:/Study/KI/Projects/###/Analysis/code_library/utilities.r"
py_utilities = "D:/Study/KI/Projects/###/Analysis/code_library/utilities.py"
source(r_utilities)
library(reticulate)
py_run_file(py_utilities, convert = F)
######
X_GeneName_path = "D:/Study/KI/Projects/###/Analysis/useful_files/X_GeneName_Homo_sapiens.GRCh38.103.rds"
Y_GeneName_path = "D:/Study/KI/Projects/###/Analysis/useful_files/Y_GeneName_Homo_sapiens.GRCh38.103.rds"
#
hgnc_mapping_alias_genename = readRDS("D:/Study/KI/Projects/###/Analysis/useful_files/hgnc_mapping_alias_genename_seurat.rds")
hgnc_mapping_genename_ensemblID = readRDS("D:/Study/KI/Projects/###/Analysis/useful_files/hgnc_mapping_genename_ensemblID.rds")
GRCh38v103_mapping_genename_ensemblID = readRDS("D:/Study/KI/Projects/###/Analysis/useful_files/GRCh38v103_mapping_genename_ensemblID.rds")
GRCh37_mapping_genename_ensemblID = readRDS("D:/Study/KI/Projects/###/Analysis/useful_files/GRCh37_mapping_genename_ensemblID.rds")
IntegratedAnnotation_2_metadata = readRDS("D:/Study/KI/Projects/###/MetaStudies/Analysis/Integration/0_preprocessing/IntegratedAnnotation_2_metadata.rds")
IndividualAnalysis_metadata_list = readRDS("D:/Study/KI/Projects/###/MetaStudies/Analysis/Integration/0_preprocessing/IndividualAnalysis_metadata_list.rds")
#
# IO - Setting
#
py$output_dir = "D:/Study/KI/Projects/###/MetaStudies/Data/Skin/Human/data_preparation_HCA_v1_0/Takahashi_Lowry_JournalofInvestigativeDermatology_2020"
py_run_string("os.makedirs(output_dir, exist_ok=True)")
h5ad_ALL_file = "GSE129611_ALL.h5ad"
metadata_ALL_file = "GSE129611_ALL_metadata.tsv"
######
input_folder = "D:/Study/KI/Projects/###/MetaStudies/Data/Skin/Human/Takahashi_Lowry_JournalofInvestigativeDermatology_2020/GSE129611_RAW"
lib_id_c = c()
lib_id_repo_c = c()
reference_genome_c = c()
py_run_string("GEX_dict = {}")
for(this_file in list.files(input_folder, pattern = ".gz")){
  this_path = paste0(input_folder, "/", this_file)
  str_fragments = unlist(strsplit(this_file, "_", fixed = T))
  this_lib_repo = str_fragments[1]
  this_lib = unlist(strsplit(paste(str_fragments[-1], collapse = "_"), "-?(Rie|data)", perl = T))[1]
  tmp_df = as.data.frame(fread(this_path))
  py$tmp_df = tmp_df[, -1]
  py$tmp_var_df = data.frame(row.names = gsub("hg19_", "", as.character(tmp_df[, 1]), fixed = T))
  py_run_string(paste0("GEX_dict['", this_lib, "'] = ad.AnnData(X=scipy.sparse.csr_matrix(tmp_df.values.transpose()).astype(np.float32), obs=pd.DataFrame(index='", this_lib, "_' + tmp_df.columns.astype(str).values), var=tmp_var_df)"))
  lib_id_c = c(lib_id_c, rep(this_lib, eval(parse(text = paste0("py$GEX_dict$'", this_lib, "'$n_obs")))))
  lib_id_repo_c = c(lib_id_repo_c, rep(this_lib_repo, eval(parse(text = paste0("py$GEX_dict$'", this_lib, "'$n_obs")))))
  if(all(grepl("hg19_", as.character(tmp_df[, 1]), fixed = T))){
    reference_genome_c = c(reference_genome_c, rep("GRCh38", eval(parse(text = paste0("py$GEX_dict$'", this_lib, "'$n_obs")))))
  }else{
    if(any(grepl("hg19_", as.character(tmp_df[, 1]), fixed = T))){
      stop("Some (but not all) genes have hg19 tag!")
    }
    reference_genome_c = c(reference_genome_c, rep("GRCh37", eval(parse(text = paste0("py$GEX_dict$'", this_lib, "'$n_obs")))))
  }
  print(this_lib)
  rm(tmp_df)
  py_run_string("del tmp_df")
  py_run_string("del tmp_var_df")
  py_run_string("gc.collect()")
  gc()
}
#
py_run_string("adata_raw = ad.concat(GEX_dict, axis=0, join='outer')")
py_run_string("del GEX_dict")
GeneName = featureMapMerge("adata_raw", "adata_genename", hgnc_mapping_alias_genename)
py_run_string("del adata_raw")
GeneName_to_GeneID = mapping_GeneName_GeneID(GeneName, hgnc_mapping_genename_ensemblID)
GeneName_to_GeneID = mapping_GeneName_GeneID(GeneName_to_GeneID, GRCh38v103_mapping_genename_ensemblID)
GeneName_to_GeneID = mapping_GeneName_GeneID(GeneName_to_GeneID, GRCh37_mapping_genename_ensemblID)
names(GeneName_to_GeneID) = GeneName
GeneID = featureMapMerge("adata_genename", "adata", GeneName_to_GeneID)
py_run_string("del adata_genename")
py$var_data_frame = data.frame(row.names = GeneID, feature_is_filtered=F, feature_biotype="gene", feature_reference="NCBITaxon:9606", feature_name=names(GeneID))
py_run_string(paste0("var_data_frame.loc[var_data_frame.index.str.startswith('ERCC-'), 'feature_biotype'] = 'spike-in'"))
py_run_string(paste0("adata.var = var_data_frame.copy()"))
py_run_string("del var_data_frame")
py_run_string("gc.collect()")
#
# Metadata
#
HCA_metadata_list = load_HCA_metadata_xlsx("D:/Study/KI/Projects/###/MetaStudies/Data/Skin/Human/HCA_metadata/Takahashi_2020_HCA_tier 1_metadata.xlsx")
HCA_metadata_list_original = HCA_metadata_list # Keep the original for checking later
rownames(HCA_metadata_list$DN) = HCA_metadata_list$DN[, "donor_id"]
HCA_metadata_list$DN[is.na(HCA_metadata_list$DN)] = "unknown"
rownames(HCA_metadata_list$SP) = HCA_metadata_list$SP[, "library_id"]
HCA_metadata_list$SP[is.na(HCA_metadata_list$SP)] = "unknown"
#
donor_id_c = HCA_metadata_list$SP[lib_id_c, "donor_id"]
# Sex Identification
X_GeneName = readRDS(X_GeneName_path)
Y_GeneName = readRDS(Y_GeneName_path)
adata_X = py$adata$X
dimnames(adata_X) = list(py$adata$obs_names$values, py$adata$var[["feature_name"]])
ALL_expression = rowSums(adata_X)
Y_expression = rowSums(adata_X[, which(colnames(adata_X) %in% Y_GeneName)])
X_expression = rowSums(adata_X[, which(colnames(adata_X) %in% X_GeneName)])
Y_ratio = Y_expression / ALL_expression
X_ratio = X_expression / ALL_expression
Y_ratio_aggregate = aggregate(Y_ratio, by = list(donor_id_c), FUN = function(x){
  return(mean(x, na.rm = T))
})
X_ratio_aggregate = aggregate(X_ratio, by = list(donor_id_c), FUN = function(x){
  return(mean(x, na.rm = T))
})
Y_X_ratio = Y_ratio_aggregate$x / X_ratio_aggregate$x
names(Y_X_ratio) = Y_ratio_aggregate$Group.1
pdf(paste0(py$output_dir, '/DonorID_Sex_barplot.pdf'), width = 10, height = 5)
par(mar=c(11, 4, 4, 4))
barplot(Y_X_ratio, las=2)
dev.off()
#
### Study level
py$adata$obs["study_id"] = "Takahashi_Lowry_JournalofInvestigativeDermatology_2020"
### Experiment level
py$adata$obs["experiment_id"] = as.character(py$adata$obs[["study_id"]])
py$adata$obs["institute"] = HCA_metadata_list$SP[lib_id_c, "institute"]
py$adata$obs["author_batch_notes"] = HCA_metadata_list$SP[lib_id_c, "author_batch_notes"]
### Donor background level
py$adata$obs["organism"] = "Homo sapiens"
py$adata$obs["organism_ontology_term_id"] = HCA_metadata_list$DN[donor_id_c, "organism_ontology_term_id"]
py$adata$obs["donor_id"] = donor_id_c
py$adata$obs["donor_id_running"] = paste0(py$adata$obs[["experiment_id"]], "___", lib_id_c)
py$adata$obs["sex_ontology_term"] = HCA_metadata_list$DN[donor_id_c, "sex_ontology_term"]
py$adata$obs["sex_ontology_term_id"] = HCA_metadata_list$DN[donor_id_c, "sex_ontology_term_id"]
py$adata$obs["ethnicity_1"] = "unknown"
py$adata$obs["ethnicity_2"] = "unknown"
py$adata$obs["ethnicity_details"] = "unknown"
py$adata$obs["self_reported_ethnicity_ontology_term_id"] = "unknown"
py$adata$obs["genotype"] = "unknown"
### Donor status level
py$adata$obs["age_years"] = "unknown"
py$adata$obs["age_range"] = HCA_metadata_list$SP[lib_id_c, "age_range"]
py$adata$obs["development_stage_ontology_term"] = "unknown"
py$adata$obs["development_stage_ontology_term_id"] = HCA_metadata_list$SP[lib_id_c, "development_stage_ontology_term_id"]
py$adata$obs["disease"] = "normal"
py$adata$obs["disease_ontology_term_id"] = HCA_metadata_list$SP[lib_id_c, "disease_ontology_term_id"]
py$adata$obs["disease_status"] = "not applicable"
py$adata$obs["treatment_status"] = "unknown"
#
py$adata$obs["smoking_status"] = "unknown"
py$adata$obs["smoking_history"] = "unknown"
py$adata$obs["bmi"] = "unknown"
py$adata$obs["systemic_disorder"] = "unknown"
py$adata$obs["autoimmune_component"] = "unknown"
py$adata$obs["menopause_status"] = "unknown"
py$adata$obs["puperty_status"] = "unknown"
py$adata$obs["sun_exposure"] = "unknown"
py$adata$obs["skin_tone"] = "unknown"
py$adata$obs["skin_care"] = "unknown"
py$adata$obs["manner_of_death"] = HCA_metadata_list$DN[donor_id_c, "manner_of_death"]
### Sample level
py$adata$obs["sample_id"] = HCA_metadata_list$SP[lib_id_c, "sample_id"]
py$adata$obs["sample_id_running"] = as.character(py$adata$obs[["donor_id_running"]])
py$adata$obs["sample_source"] = HCA_metadata_list$SP[lib_id_c, "sample_source"]
py$adata$obs["sample_collection_method"] = HCA_metadata_list$SP[lib_id_c, "sample_collection_method"]
py$adata$obs["sampled_site_condition"] = HCA_metadata_list$SP[lib_id_c, "sampled_site_condition"]
py$adata$obs["anatomical_region_level_1"] = "Head"
py$adata$obs["anatomical_region_level_2"] = "Scalp"
py$adata$obs["anatomical_region_level_3"] = "unknown"
metadata_AR = as.character(py$adata$obs[["anatomical_region_level_3"]])
metadata_AR[metadata_AR %in% c("unknown", "not applicable")] = NA
metadata_AR[is.na(metadata_AR)] = as.character(py$adata$obs[["anatomical_region_level_2"]])[is.na(metadata_AR)]
metadata_AR[metadata_AR %in% c("unknown", "not applicable")] = NA
metadata_AR[is.na(metadata_AR)] = as.character(py$adata$obs[["anatomical_region_level_1"]])[is.na(metadata_AR)]
metadata_AR[metadata_AR %in% c("unknown", "not applicable")] = NA
metadata_AR[is.na(metadata_AR)] = "unknown"
py$adata$obs["anatomical_region"] = metadata_AR
py$adata$obs["tissue_type"] = HCA_metadata_list$SP[lib_id_c, "tissue_type"]
py$adata$obs["tissue_free_text"] = HCA_metadata_list$SP[lib_id_c, "tissue_free_text"]
py$adata$obs["tissue_ontology_term"] = HCA_metadata_list$SP[lib_id_c, "tissue_ontology_term"]
py$adata$obs["tissue_ontology_term_id"] = HCA_metadata_list$SP[lib_id_c, "tissue_ontology_term_id"]
py$adata$obs["skin_tissue"] = "unknown"
py$adata$obs["sample_collection_site"] = HCA_metadata_list$SP[lib_id_c, "sample_collection_site"]
py$adata$obs["sample_cultured"] = "No"
### Sequencing level
py$adata$obs["sequencing_type"] = "GEX"
py$adata$obs["sample_collection_year"] = HCA_metadata_list$SP[lib_id_c, "sample_collection_year"]
py$adata$obs["sample_collection_relative_time_point"] = HCA_metadata_list$SP[lib_id_c, "sample_collection_relative_time_point"]
py$adata$obs["library_id"] = lib_id_c
py$adata$obs["library_id_repository"] = HCA_metadata_list$SP[lib_id_c, "library_id_repository"]
py$adata$obs["library_preparation_batch"] = HCA_metadata_list$SP[lib_id_c, "library_preparation_batch"]
py$adata$obs["library_sequencing_run"] = HCA_metadata_list$SP[lib_id_c, "library_sequencing_run"]
py$adata$obs["sample_preservation_method"] = HCA_metadata_list$SP[lib_id_c, "sample_preservation_method"]
py$adata$obs["protocol_url"] = "unknown"
py$adata$obs["dissociation_protocol"] = HCA_metadata_list$SP[lib_id_c, "dissociation_protocol"]
py$adata$obs["cell_enrichment"] = HCA_metadata_list$SP[lib_id_c, "cell_enrichment"]
py$adata$obs["cell_viability_percentage"] = HCA_metadata_list$SP[lib_id_c, "cell_viability_percentage"]
py$adata$obs["cell_number_loaded"] = HCA_metadata_list$SP[lib_id_c, "cell_number_loaded"]
py$adata$obs["suspension_type"] = HCA_metadata_list$SP[lib_id_c, "suspension_type"]
py$adata$obs["assay_ontology_term"] = dplyr::case_match(
  lib_id_c,
  "DS1" ~ "Drop-seq",
  "DS2" ~ "Drop-seq",
  "DS3" ~ "Drop-seq",
  "10X1" ~ "10x 3' v2",
  "10X2" ~ "10x 3' v2",
  .default = lib_id_c
)
py$adata$obs["assay_ontology_term_id"] = dplyr::case_match(
  lib_id_c,
  "DS1" ~ "EFO:0008722",
  "DS2" ~ "EFO:0008722",
  "DS3" ~ "EFO:0008722",
  "10X1" ~ "EFO:0009899",
  "10X2" ~ "EFO:0009899",
  .default = lib_id_c
)
py$adata$obs["sequenced_fragment"] = "3 prime tag"
py$adata$obs["sequencing_platform"] = "Illumina HiSeq 2500"
### Preprocessing level
py$adata$obs["reference_genome"] = reference_genome_c
py$adata$obs["alignment_software"] = dplyr::case_match(
  lib_id_c,
  "DS1" ~ "HiSat2",
  "DS2" ~ "HiSat2",
  "DS3" ~ "HiSat2",
  "10X1" ~ "cell ranger 2.2",
  "10X2" ~ "cell ranger 2.2",
  .default = lib_id_c
)
py$adata$obs["intron_inclusion"] = "unknown"
### Analysis Level
py$adata$obs["is_primary_data"] = FALSE
py$adata$obs["author_cell_type_1"] = "unknown"
py$adata$obs["author_cell_type_2"] = "unknown"
py$adata$obs["author_cell_type_3"] = "unknown"
py$adata$obs["author_cell_type_4"] = "unknown"
py$adata$obs["author_cell_type_5"] = "unknown"
py$adata$obs["author_cell_type_6"] = "unknown"
### Previous Barcode
this_dataset_individual = "Takahashi 2020"
this_dataset_integrated = "Takahashi_Lowry_JournalofInvestigativeDermatology_2020"
individual_barcode_used = rownames(IndividualAnalysis_metadata_list[[this_dataset_individual]])
names(individual_barcode_used) = sapply(paste0(IndividualAnalysis_metadata_list[[this_dataset_individual]][, "SampleID"], "_", rownames(IndividualAnalysis_metadata_list[[this_dataset_individual]])), function(x){
  str_fragments = unlist(strsplit(x, ".", fixed = T))
  if(length(str_fragments) > 1){
    result = paste(str_fragments[-length(str_fragments)], collapse = ".")
  }else{
    result = x
  }
  return(result)
})
individual_barcode_all = individual_barcode_used[paste0(py$adata$obs[["study_id"]], "___", py$adata$obs_names$values)]
individual_barcode_all[is.na(individual_barcode_all)] = "unknown"
py$adata$obs["individual_barcode"] = individual_barcode_all
if(sum(individual_barcode_all != "unknown") != length(individual_barcode_used)){
  stop("Individual barcodes didn't match well!")
}else{
  print("Pass: Individual")
}
###
integrated_barcode_used = rownames(IntegratedAnnotation_2_metadata)[sapply(rownames(IntegratedAnnotation_2_metadata), function(x){
  str_fragments = unlist(strsplit(x, "---", fixed = T))
  return(str_fragments[length(str_fragments)])
}) == this_dataset_integrated]
names(integrated_barcode_used) = sapply(paste0(IntegratedAnnotation_2_metadata[integrated_barcode_used, "SampleID"], "_", sapply(integrated_barcode_used, function(x){
  str_fragments = unlist(strsplit(x, "---", fixed = T))
  return(str_fragments[-length(str_fragments)])
})), function(x){
  str_fragments = unlist(strsplit(x, ".", fixed = T))
  if(length(str_fragments) > 1){
    result = paste(str_fragments[-length(str_fragments)], collapse = ".")
  }else{
    result = x
  }
  return(result)
})
integrated_barcode_all = integrated_barcode_used[paste0(py$adata$obs[["study_id"]], "___", py$adata$obs_names$values)]
integrated_barcode_all[is.na(integrated_barcode_all)] = "unknown"
py$adata$obs["integrated_barcode"] = integrated_barcode_all
if(sum(integrated_barcode_all != "unknown") != length(integrated_barcode_used)){
  stop("Integrated barcodes didn't match well!")
}else{
  print("Pass: Integrated")
}
#
# Embeddings
#
#
# Put unique identifier for this experiment
#
py_run_string("adata.obs_names = adata.obs['experiment_id'].astype(str) + '___' + adata.obs_names.astype(str)")
#
# Check before writing
#
print(py$adata)
print(py$adata$obs)
print(py$adata$var)
#
# Write
#
py_run_string("sc.pp.filter_cells(adata, min_genes=1)")
py_run_string(paste0("adata.write(output_dir + '/", h5ad_ALL_file, "')"))
py_run_string(paste0("adata.obs.to_csv(output_dir + '/", metadata_ALL_file, "', sep='\t', index=True)"))
py_run_string("print(output_dir)")
#
metadata_matrix = as.matrix(py$adata$obs)
write.table(x = HCA_metadata_usage(HCA_metadata=HCA_metadata_list_original$DT, ADATA_metadata=metadata_matrix,
                                   metadata_rowname_column="assay_ontology_term_id"),
            file = paste0(py$output_dir, "/", "HCA_metadata_usage_DT.txt"), sep = "\t", col.names=T, row.names=F, quote=T)
write.table(x = HCA_metadata_usage(HCA_metadata=HCA_metadata_list_original$DN, ADATA_metadata=metadata_matrix,
                                   metadata_rowname_column="donor_id"),
            file = paste0(py$output_dir, "/", "HCA_metadata_usage_DN.txt"), sep = "\t", col.names=T, row.names=F, quote=T)
write.table(x = HCA_metadata_usage(HCA_metadata=HCA_metadata_list_original$SP, ADATA_metadata=metadata_matrix,
                                   metadata_rowname_column="sample_id"),
            file = paste0(py$output_dir, "/", "HCA_metadata_usage_SP.txt"), sep = "\t", col.names=T, row.names=F, quote=T)

