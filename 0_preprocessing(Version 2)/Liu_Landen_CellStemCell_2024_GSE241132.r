r_utilities = "D:/Study/KI/Projects/###/Analysis/code_library/utilities.r"
py_utilities = "D:/Study/KI/Projects/###/Analysis/code_library/utilities.py"
source(r_utilities)
library(reticulate)
py_run_file(py_utilities, convert = F)
######
X_GeneID_path = "D:/Study/KI/Projects/###/Analysis/useful_files/X_GeneID_Homo_sapiens.GRCh38.103.rds"
Y_GeneID_path = "D:/Study/KI/Projects/###/Analysis/useful_files/Y_GeneID_Homo_sapiens.GRCh38.103.rds"
#
hgnc_mapping_ensemblID_genename = readRDS("D:/Study/KI/Projects/###/Analysis/useful_files/hgnc_mapping_ensemblID_genename_seurat.rds")
IntegratedAnnotation_2_metadata = readRDS("D:/Study/KI/Projects/###/MetaStudies/Analysis/Integration/0_preprocessing/IntegratedAnnotation_2_metadata.rds")
IndividualAnalysis_metadata_list = readRDS("D:/Study/KI/Projects/###/MetaStudies/Analysis/Integration/0_preprocessing/IndividualAnalysis_metadata_list.rds")
#
# IO - Setting
#
py$output_dir = "D:/Study/KI/Projects/###/MetaStudies/Data/Skin/Human/data_preparation_HCA_v1_0/Liu_Landen_CellStemCell_2024"
py_run_string("os.makedirs(output_dir, exist_ok=True)")
h5ad_ALL_file = "GSE241132_ALL.h5ad"
metadata_ALL_file = "GSE241132_ALL_metadata.tsv"
h5ad_HC_file = "GSE241132_HC.h5ad"
metadata_HC_file = "GSE241132_HC_metadata.tsv"
h5ad_Ning_file = "GSE241132_Ning.h5ad"
metadata_Ning_file = "GSE241132_Ning_metadata.tsv"
######
input_folder = "D:/Study/KI/Projects/###/MetaStudies/Data/Skin/Human/Liu_Landen_CellStemCell_2024/GSE241132_RAW"
sample_folders = list.dirs(input_folder, recursive = F)
py_run_string("GEX_dict = {}")
for(ii in sample_folders){
  ii_splits = unlist(strsplit(ii, "[_/]", perl = T))
  this_sample = ii_splits[length(ii_splits)]
  this_path = paste0(ii, "/", this_sample)
  ReadMtx_py("tmp_adata", data_dir = this_path, prefix = paste0(this_sample, "_"), feature_column_name = c("gene_ids", "feature_name", "feature_types"))
  py_run_string(paste0("tmp_adata.obs_names = '", this_sample, "_' + tmp_adata.obs_names.astype(str).values"))
  py$tmp_adata$obs["library_id_repository"] = rep(ii_splits[length(ii_splits) - 1], py$tmp_adata$n_obs)
  py_run_string(paste0("GEX_dict['", this_sample, "'] = tmp_adata"))
  print(this_sample)
}
#
py_run_string("adata = ad.concat(GEX_dict, axis=0, join='outer')")
py_run_string("SameVar = AssertSameDictAdataVar(GEX_dict)")
if(py$SameVar){
  py_run_string("adata.var = tmp_adata.var.copy()")
}
py_run_string("del tmp_adata")
py_run_string("del GEX_dict")
GeneID = py$adata$var_names$values
py$var_data_frame = data.frame(row.names = GeneID, feature_is_filtered=F, feature_biotype="gene", feature_reference="NCBITaxon:9606", feature_name=mapping_GeneName_GeneID(GeneID, hgnc_mapping_ensemblID_genename))
py_run_string(paste0("var_data_frame.loc[var_data_frame.index.str.startswith('ERCC-'), 'feature_biotype'] = 'spike-in'"))
py_run_string(paste0("adata.var = var_data_frame.copy()"))
py_run_string("del var_data_frame")
py_run_string("gc.collect()")
#
GSE241132_metadata_path = "D:/Study/KI/Projects/###/MetaStudies/Data/Skin/Human/Liu_Landen_CellStemCell_2024/GSE241132_cell_metadata.txt.gz"
metadata_raw = as.matrix(read.table(gzfile(GSE241132_metadata_path), header = T, row.names = 1, sep = "\t"))
author_cell_type_1 = metadata_raw[, "newMainCellTypes"][py$adata$obs_names$values]
author_cell_type_1[is.na(author_cell_type_1)] = "unknown"
author_cell_type_2 = metadata_raw[, "oldCellTypes_InATACproject"][py$adata$obs_names$values]
author_cell_type_2[is.na(author_cell_type_2)] = "unknown"
author_cell_type_3 = metadata_raw[, "newCellTypes"][py$adata$obs_names$values]
author_cell_type_3[is.na(author_cell_type_3)] = "unknown"
#
# Metadata
#
sample_id_running_c = sapply(py$adata$obs_names$values, function(x){
  return(unlist(strsplit(x, "_"))[1])
})
day_c = sapply(sample_id_running_c, function(x){
  return(unlist(strsplit(x, "D"))[2])
})
lib_id_c = as.character(py$adata$obs[["library_id_repository"]])
#
HCA_metadata_list = load_HCA_metadata_xlsx("D:/Study/KI/Projects/###/MetaStudies/Data/Skin/Human/HCA_metadata/Liu_2024_HCA_tier 1_metadata.xlsx")
HCA_metadata_list_original = HCA_metadata_list # Keep the original for checking later
HCA_metadata_list_Hao = load_HCA_metadata_xlsx("D:/Study/KI/Projects/###/MetaStudies/Data/Skin/Human/HCA_metadata/Liu_2024_HCA_tier 1_metadata_Hao.xlsx")
rownames(HCA_metadata_list$DN) = HCA_metadata_list$DN[, "donor_id"]
HCA_metadata_list$DN[is.na(HCA_metadata_list$DN)] = "unknown"
HCA_metadata_list$SP = HCA_metadata_list_Hao$SP
rownames(HCA_metadata_list$SP) = HCA_metadata_list$SP[, "library_id_repository"]
HCA_metadata_list$SP[is.na(HCA_metadata_list$SP)] = "unknown"
#
donor_id_c = HCA_metadata_list$SP[lib_id_c, "donor_id"]
# Sex Identification
X_GeneID = readRDS(X_GeneID_path)
Y_GeneID = readRDS(Y_GeneID_path)
adata_X = py$adata$X
dimnames(adata_X) = list(py$adata$obs_names$values, py$adata$var_names$values)
ALL_expression = rowSums(adata_X)
Y_expression = rowSums(adata_X[, which(colnames(adata_X) %in% Y_GeneID)])
X_expression = rowSums(adata_X[, which(colnames(adata_X) %in% X_GeneID)])
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
pdf(paste0(py$output_dir, '/GSE241132_DonorID_Sex_barplot.pdf'), width = 10, height = 5)
par(mar=c(11, 4, 4, 4))
barplot(Y_X_ratio, las=2)
dev.off()
#
### Study level
py$adata$obs["study_id"] = "Liu_Landen_CellStemCell_2024"
### Experiment level
py_run_string("adata.obs['experiment_id'] = adata.obs['study_id'].astype(str) + '_' + 'GSE241132'")
py$adata$obs["institute"] = "Karolinska Institutet"
py$adata$obs["author_batch_notes"] = HCA_metadata_list$SP[lib_id_c, "author_batch_notes"]
### Donor background level
py$adata$obs["organism"] = "Homo sapiens"
py$adata$obs["organism_ontology_term_id"] = HCA_metadata_list$DN[donor_id_c, "organism_ontology_term_id"]
py$adata$obs["donor_id"] = donor_id_c
py$adata$obs["donor_id_running"] = paste0(py$adata$obs[["experiment_id"]], "___", donor_id_c)
py$adata$obs["sex_ontology_term"] = HCA_metadata_list$DN[donor_id_c, "sex_ontology_term"]
py$adata$obs["sex_ontology_term_id"] = HCA_metadata_list$DN[donor_id_c, "sex_ontology_term_id"]
py$adata$obs["ethnicity_1"] = "European Ancestry"
py$adata$obs["ethnicity_2"] = "Caucasian"
py$adata$obs["ethnicity_details"] = "Caucasian"
py$adata$obs["self_reported_ethnicity_ontology_term_id"] = "unknown"
py$adata$obs["genotype"] = "unknown"
### Donor status level
py$adata$obs["age_years"] = dplyr::case_match(
  donor_id_c,
  "PWH26" ~ "29",
  "PWH27" ~ "22",
  "PWH28" ~ "24",
  .default = donor_id_c
)
age_years_character = py$adata$obs[["age_years"]]
age_years_numeric = as.numeric(age_years_character)
development_stage_ontology_term_id = rep(NA, length(age_years_character))
development_stage_ontology_term_id[age_years_numeric >= 0 & age_years_numeric <= 14] = "HsapDv:0000264"
development_stage_ontology_term_id[age_years_numeric >= 15 & age_years_numeric <= 19] = "HsapDv:0000268"
development_stage_ontology_term_id[age_years_numeric >= 20 & age_years_numeric <= 29] = "HsapDv:0000237"
development_stage_ontology_term_id[age_years_numeric >= 30 & age_years_numeric <= 39] = "HsapDv:0000238"
development_stage_ontology_term_id[age_years_numeric >= 40 & age_years_numeric <= 49] = "HsapDv:0000239"
development_stage_ontology_term_id[age_years_numeric >= 50 & age_years_numeric <= 59] = "HsapDv:0000240"
development_stage_ontology_term_id[age_years_numeric >= 60 & age_years_numeric <= 69] = "HsapDv:0000241"
development_stage_ontology_term_id[age_years_numeric >= 70 & age_years_numeric <= 79] = "HsapDv:0000242"
development_stage_ontology_term_id[age_years_numeric >= 80 & age_years_numeric <= 89] = "HsapDv:0000243"
development_stage_ontology_term_id[age_years_numeric >= 90 & age_years_numeric <= 99] = "HsapDv:0000244"
development_stage_ontology_term_id[age_years_numeric >= 100 & age_years_numeric <= 109] = "HsapDv:0000245"
development_stage_ontology_term_id[is.na(development_stage_ontology_term_id)] = "unknown"
development_stage = dplyr::case_match(
  development_stage_ontology_term_id,
  "HsapDv:0000264" ~ "pediatric stage",
  "HsapDv:0000268" ~ "15-19 year-old",
  "HsapDv:0000237" ~ "third decade stage",
  "HsapDv:0000238" ~ "fourth decade stage",
  "HsapDv:0000239" ~ "fifth decade stage",
  "HsapDv:0000240" ~ "sixth decade stage",
  "HsapDv:0000241" ~ "seventh decade stage",
  "HsapDv:0000242" ~ "eighth decade stage",
  "HsapDv:0000243" ~ "ninth decade stage",
  "HsapDv:0000244" ~ "tenth decade stage",
  "HsapDv:0000245" ~ "eleventh decade stage",
  .default = development_stage_ontology_term_id
)
age_range_c = dplyr::case_match(
  development_stage_ontology_term_id,
  "HsapDv:0000264" ~ "0-14",
  "HsapDv:0000268" ~ "15-19",
  "HsapDv:0000237" ~ "20-29",
  "HsapDv:0000238" ~ "30-39",
  "HsapDv:0000239" ~ "40-49",
  "HsapDv:0000240" ~ "50-59",
  "HsapDv:0000241" ~ "60-69",
  "HsapDv:0000242" ~ "70-79",
  "HsapDv:0000243" ~ "80-89",
  "HsapDv:0000244" ~ "90-99",
  "HsapDv:0000245" ~ "100-109",
  .default = development_stage_ontology_term_id
)
py$adata$obs["age_range"] = age_range_c
py$adata$obs["development_stage_ontology_term"] = development_stage
py$adata$obs["development_stage_ontology_term_id"] = HCA_metadata_list$SP[lib_id_c, "development_stage_ontology_term_id"]
py$adata$obs["disease"] = "normal"
py$adata$obs["disease_ontology_term_id"] = HCA_metadata_list$SP[lib_id_c, "disease_ontology_term_id"]
py$adata$obs["disease_status"] = "not applicable"
py$adata$obs["treatment_status"] = dplyr::case_match(
  day_c,
  "0" ~ "No",
  "1" ~ "Wound 1D",
  "7" ~ "Wound 7D",
  "30" ~ "Wound 30D",
  .default = day_c
)
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
py$adata$obs["sample_id_running"] = paste0(py$adata$obs[["experiment_id"]], "___", sample_id_running_c)
py$adata$obs["sample_source"] = HCA_metadata_list$SP[lib_id_c, "sample_source"]
py$adata$obs["sample_collection_method"] = HCA_metadata_list$SP[lib_id_c, "sample_collection_method"]
py$adata$obs["sampled_site_condition"] = dplyr::case_match(
  day_c,
  "0" ~ "healthy",
  "1" ~ "wound",
  "7" ~ "wound",
  "30" ~ "wound",
  .default = day_c
)
py$adata$obs["anatomical_region_level_1"] = "Torso"
py$adata$obs["anatomical_region_level_2"] = "Back"
py$adata$obs["anatomical_region_level_3"] = "Lower back"
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
py$adata$obs["skin_tissue"] = "Full Skin"
py$adata$obs["sample_collection_site"] = HCA_metadata_list$SP[lib_id_c, "sample_collection_site"]
py$adata$obs["sample_cultured"] = "No"
### Sequencing level
py$adata$obs["sequencing_type"] = "GEX"
py$adata$obs["sample_collection_year"] = HCA_metadata_list$SP[lib_id_c, "sample_collection_year"]
py$adata$obs["sample_collection_relative_time_point"] = HCA_metadata_list$SP[lib_id_c, "sample_collection_relative_time_point"]
py$adata$obs["library_id"] = HCA_metadata_list$SP[lib_id_c, "library_id"]
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
py$adata$obs["assay_ontology_term"] = "10x 3' v3"
py$adata$obs["assay_ontology_term_id"] = "EFO:0009922"
py$adata$obs["sequenced_fragment"] = "3 prime tag"
py$adata$obs["sequencing_platform"] = "Illumina NovaSeq 6000" # S4
### Preprocessing level
py$adata$obs["reference_genome"] = "GRCh38"
py$adata$obs["alignment_software"] = "cell ranger 5.0.1"
py$adata$obs["intron_inclusion"] = "unknown"
### Analysis Level
py$adata$obs["is_primary_data"] = FALSE
py$adata$obs["author_cell_type_1"] = author_cell_type_1
py$adata$obs["author_cell_type_2"] = author_cell_type_2
py$adata$obs["author_cell_type_3"] = author_cell_type_3
py$adata$obs["author_cell_type_4"] = "unknown"
py$adata$obs["author_cell_type_5"] = "unknown"
py$adata$obs["author_cell_type_6"] = "unknown"
### Previous Barcode
this_dataset_individual = "Liu 2024"
this_dataset_integrated = "Ning"
individual_barcode_used = rownames(IndividualAnalysis_metadata_list[[this_dataset_individual]])
names(individual_barcode_used) = individual_barcode_used
individual_barcode_all = individual_barcode_used[py$adata$obs_names$values]
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
names(integrated_barcode_used) = sapply(integrated_barcode_used, function(x){
  str_fragments = unlist(strsplit(x, "---", fixed = T))
  return(str_fragments[-length(str_fragments)])
})
integrated_barcode_all = integrated_barcode_used[py$adata$obs_names$values]
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
# Filtering
#
py$adata_ALL = py$adata[py$adata$obs_names$values %in% paste0("Liu_Landen_CellStemCell_2024_GSE241132___", rownames(metadata_raw))]
py$adata_Ning = py$adata[as.character(py$adata$obs[["individual_barcode"]]) != "unknown", ]
#
# Write
#
py_run_string("sc.pp.filter_cells(adata_ALL, min_genes=1)")
py_run_string(paste0("adata_ALL.write(output_dir + '/", h5ad_ALL_file, "')"))
py_run_string(paste0("adata_ALL.obs.to_csv(output_dir + '/", metadata_ALL_file, "', sep='\t', index=True)"))
py_run_string("adata_HC = adata_ALL[adata_ALL.obs['sampled_site_condition'] == 'healthy']")
py_run_string(paste0("adata_HC.write(output_dir + '/", h5ad_HC_file, "')"))
py_run_string(paste0("adata_HC.obs.to_csv(output_dir + '/", metadata_HC_file, "', sep='\t', index=True)"))
#
py_run_string("sc.pp.filter_cells(adata_Ning, min_genes=1)")
py_run_string(paste0("adata_Ning.write(output_dir + '/", h5ad_Ning_file, "')"))
py_run_string(paste0("adata_Ning.obs.to_csv(output_dir + '/", metadata_Ning_file, "', sep='\t', index=True)"))
py_run_string("print(output_dir)")
#
metadata_matrix = as.matrix(py$adata_ALL$obs)
write.table(x = HCA_metadata_usage(HCA_metadata=HCA_metadata_list_original$DT, ADATA_metadata=metadata_matrix,
                                   metadata_rowname_column="assay_ontology_term_id"),
            file = paste0(py$output_dir, "/GSE241132_", "HCA_metadata_usage_DT.txt"), sep = "\t", col.names=T, row.names=F, quote=T)
write.table(x = HCA_metadata_usage(HCA_metadata=HCA_metadata_list_original$DN, ADATA_metadata=metadata_matrix,
                                   metadata_rowname_column="donor_id"),
            file = paste0(py$output_dir, "/GSE241132_", "HCA_metadata_usage_DN.txt"), sep = "\t", col.names=T, row.names=F, quote=T)
write.table(x = HCA_metadata_usage(HCA_metadata=HCA_metadata_list_original$SP, ADATA_metadata=metadata_matrix,
                                   metadata_rowname_column="sample_id"),
            file = paste0(py$output_dir, "/GSE241132_", "HCA_metadata_usage_SP.txt"), sep = "\t", col.names=T, row.names=F, quote=T)

