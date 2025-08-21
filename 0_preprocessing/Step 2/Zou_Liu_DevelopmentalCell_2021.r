r_utilities = "/home/haoy/projects/Analysis/code_library/utilities.r"
py_utilities = "/home/haoy/projects/Analysis/code_library/utilities.py"
source(r_utilities)
library(reticulate)
py_run_file(py_utilities, convert = F)
######
X_GeneID_path = "/home/haoy/projects/Analysis/useful_files/X_GeneID_Homo_sapiens.GRCh38.103.rds"
Y_GeneID_path = "/home/haoy/projects/Analysis/useful_files/Y_GeneID_Homo_sapiens.GRCh38.103.rds"
#
hgnc_mapping_ensemblID_genename = readRDS("/home/haoy/projects/Analysis/useful_files/hgnc_mapping_ensemblID_genename_seurat.rds")
IntegratedAnnotation_2_metadata = readRDS("/home/haoy/projects/MetaStudiesAnalysis/Integration/Output/Integration/v67_gauss/IntegratedAnnotation_2_metadata.rds")
IndividualAnalysis_metadata_list = readRDS("/home/haoy/projects/MetaStudiesAnalysis/Integration/Output/Integration/v67_gauss/IndividualAnalysis_metadata_list.rds")
#
# IO - Setting
#
py$output_dir = "/home/haoy/projects/MetaStudiesAnalysis/Integration/Data/Skin/Human/data_preparation_HCA_v1_0/Zou_Liu_DevelopmentalCell_2021"
py_run_string("os.makedirs(output_dir, exist_ok=True)")
h5ad_ALL_file = "HRA000395_ALL.h5ad"
metadata_ALL_file = "HRA000395_ALL_metadata.tsv"
######
input_folder = "/home/haoy/projects/MetaStudiesAnalysis/Integration/Data/Zou_Liu_DevelopmentalCell_2021"
sample_id_c = c()
py_run_string("GEX_dict = {}")
for(this_sample in grep("HRS", list.dirs(path = input_folder, full.names = F, recursive = F), value = T)){
  ReadMtx_py("tmp_adata", data_dir = paste0(input_folder, "/", this_sample, "/outs/filtered_feature_bc_matrix"), prefix = "", feature_column_name = c("gene_ids", "feature_name", "feature_types"))
  py_run_string(paste0("tmp_adata.obs_names = '", this_sample, "_' + tmp_adata.obs_names.astype(str).values"))
  py_run_string(paste0("GEX_dict['", this_sample, "'] = tmp_adata"))
  sample_id_c = c(sample_id_c, rep(this_sample, py$tmp_adata$n_obs))
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
# Metadata
#
HCA_metadata_list = load_HCA_metadata_xlsx("/home/haoy/projects/MetaStudiesAnalysis/Integration/Data/HCA_metadata/Zou_2021_HCA_tier 1_metadata.xlsx")
HCA_metadata_list_original = HCA_metadata_list # Keep the original for checking later
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
Y_ratio_aggregate = aggregate(Y_ratio, by = list(sample_id_c), FUN = function(x){
  return(mean(x, na.rm = T))
})
X_ratio_aggregate = aggregate(X_ratio, by = list(sample_id_c), FUN = function(x){
  return(mean(x, na.rm = T))
})
Y_X_ratio = Y_ratio_aggregate$x / X_ratio_aggregate$x
names(Y_X_ratio) = Y_ratio_aggregate$Group.1
pdf(paste0(py$output_dir, '/SampleID_Sex_barplot.pdf'), width = 10, height = 5)
par(mar=c(11, 4, 4, 4))
barplot(Y_X_ratio, las=2)
dev.off()
#
### Study level
py$adata$obs["study_id"] = "Zou_Liu_DevelopmentalCell_2021"
### Experiment level
py$adata$obs["experiment_id"] = as.character(py$adata$obs[["study_id"]])
py$adata$obs["institute"] = "unknown"
py$adata$obs["author_batch_notes"] = "unknown"
### Donor background level
py$adata$obs["organism"] = "Homo sapiens"
py$adata$obs["organism_ontology_term_id"] = "NCBITaxon:9606"
py$adata$obs["donor_id"] = dplyr::case_match(
  sample_id_c,
  "HRS118996" ~ "HRI077736",
  "HRS118997" ~ "HRI077737",
  "HRS118998" ~ "HRI077738",
  "HRS118999" ~ "HRI077739",
  "HRS119000" ~ "HRI077740",
  "HRS119001" ~ "HRI077741",
  "HRS119002" ~ "HRI077742",
  "HRS119003" ~ "HRI077743",
  "HRS119004" ~ "HRI077744",
  .default = sample_id_c
)
py$adata$obs["donor_id_running"] = paste0(py$adata$obs[["experiment_id"]], "___", dplyr::case_match(
  sample_id_c,
  "HRS118996" ~ "Skin-Y-18",
  "HRS118997" ~ "Skin-Y-22",
  "HRS118998" ~ "Skin-Y-23",
  "HRS118999" ~ "Skin-M-44",
  "HRS119000" ~ "Skin-M-47",
  "HRS119001" ~ "Skin-M-48",
  "HRS119002" ~ "Skin-O-70",
  "HRS119003" ~ "Skin-O-73",
  "HRS119004" ~ "Skin-O-76",
  .default = sample_id_c
))
py$adata$obs["sex_ontology_term"] = "female"
py$adata$obs["sex_ontology_term_id"] = "PATO:0000383"
py$adata$obs["ethnicity_1"] = "Asian ancestry"
py$adata$obs["ethnicity_2"] = "East Asian ancestry"
py$adata$obs["ethnicity_details"] = "Han Chinese"
py$adata$obs["self_reported_ethnicity_ontology_term_id"] = "HANCESTRO:0027"
py$adata$obs["genotype"] = "unknown"
### Donor status level
py$adata$obs["age_years"] = dplyr::case_match(
  sample_id_c,
  "HRS118996" ~ "18",
  "HRS118997" ~ "22",
  "HRS118998" ~ "23",
  "HRS118999" ~ "44",
  "HRS119000" ~ "47",
  "HRS119001" ~ "48",
  "HRS119002" ~ "70",
  "HRS119003" ~ "73",
  "HRS119004" ~ "76",
  .default = sample_id_c
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
py$adata$obs["development_stage_ontology_term_id"] = development_stage_ontology_term_id
py$adata$obs["disease"] = "normal"
py$adata$obs["disease_ontology_term_id"] = "PATO:0000461"
py$adata$obs["disease_status"] = "not applicable"
py$adata$obs["treatment_status"] = "No"
#
py$adata$obs["smoking_status"] = "No"
py$adata$obs["smoking_history"] = "No"
py$adata$obs["bmi"] = "unknown"
py$adata$obs["systemic_disorder"] = "unknown"
py$adata$obs["autoimmune_component"] = "unknown"
py$adata$obs["menopause_status"] = "unknown"
py$adata$obs["puperty_status"] = "unknown"
py$adata$obs["sun_exposure"] = "unknown"
py$adata$obs["skin_tone"] = "unknown"
py$adata$obs["skin_care"] = "unknown"
py$adata$obs["manner_of_death"] = 'not applicable'
### Sample level
py$adata$obs["sample_id"] = sample_id_c
py$adata$obs["sample_id_running"] = as.character(py$adata$obs[["donor_id_running"]])
py$adata$obs["sample_source"] = "surgical donor"
py$adata$obs["sample_collection_method"] = "surgical resection"
py$adata$obs["sampled_site_condition"] = "healthy"
py$adata$obs["anatomical_region_level_1"] = "Head"
py$adata$obs["anatomical_region_level_2"] = "Face"
py$adata$obs["anatomical_region_level_3"] = "Eyelids"
metadata_AR = as.character(py$adata$obs[["anatomical_region_level_3"]])
metadata_AR[metadata_AR %in% c("unknown", "not applicable")] = NA
metadata_AR[is.na(metadata_AR)] = as.character(py$adata$obs[["anatomical_region_level_2"]])[is.na(metadata_AR)]
metadata_AR[metadata_AR %in% c("unknown", "not applicable")] = NA
metadata_AR[is.na(metadata_AR)] = as.character(py$adata$obs[["anatomical_region_level_1"]])[is.na(metadata_AR)]
metadata_AR[metadata_AR %in% c("unknown", "not applicable")] = NA
metadata_AR[is.na(metadata_AR)] = "unknown"
py$adata$obs["anatomical_region"] = metadata_AR
py$adata$obs["tissue_type"] = "tissue"
py$adata$obs["tissue_free_text"] = "unknown"
py$adata$obs["tissue_ontology_term"] = "unknown"
py$adata$obs["tissue_ontology_term_id"] = "unknown"
py$adata$obs["skin_tissue"] = "unknown"
py$adata$obs["sample_collection_site"] = "unknown"
py$adata$obs["sample_cultured"] = "No"
### Sequencing level
py$adata$obs["sequencing_type"] = "GEX"
py$adata$obs["sample_collection_year"] = "unknown"
py$adata$obs["sample_collection_relative_time_point"] = "unknown"
py$adata$obs["library_id"] = "unknown"
py$adata$obs["library_id_repository"] = "unknown"
py$adata$obs["library_preparation_batch"] = "unknown"
py$adata$obs["library_sequencing_run"] = "unknown"
py$adata$obs["sample_preservation_method"] = "unknown"
py$adata$obs["protocol_url"] = "unknown"
py$adata$obs["dissociation_protocol"] = "unknown"
py$adata$obs["cell_enrichment"] = "unknown"
py$adata$obs["cell_viability_percentage"] = "unknown"
py$adata$obs["cell_number_loaded"] = "unknown"
py$adata$obs["suspension_type"] = "cell"
py$adata$obs["assay_ontology_term"] = "10x 3' v2"
py$adata$obs["assay_ontology_term_id"] = "EFO:0009899"
py$adata$obs["sequenced_fragment"] = "3 prime tag"
py$adata$obs["sequencing_platform"] = "Illumina NovaSeq 6000"
### Preprocessing level
py$adata$obs["reference_genome"] = "GRCh38" # Regenerated by ourselves
py$adata$obs["alignment_software"] = "cell ranger 7.1.0" # Regenerated by ourselves
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
this_dataset_individual = "Zou 2020"
this_dataset_integrated = "Zou_Liu_DevelopmentalCell_2020"
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
individual_barcode_all = individual_barcode_used[stringi::stri_replace_all_fixed(str = paste0(py$adata$obs[["study_id"]], "___", py$adata$obs_names$values),
                                                                                 pattern = c("Zou_Liu_DevelopmentalCell_2021", "HRS118996", "HRS118997", "HRS118998", "HRS118999", "HRS119000", "HRS119001", "HRS119002", "HRS119003", "HRS119004"),
                                                                                 replacement = c("Zou_Liu_DevelopmentalCell_2020", "Skin-Y-18", "Skin-Y-22", "Skin-Y-23", "Skin-M-44", "Skin-M-47", "Skin-M-48", "Skin-O-70", "Skin-O-73", "Skin-O-76"),
                                                                                 vectorize_all = FALSE)]
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
integrated_barcode_all = integrated_barcode_used[stringi::stri_replace_all_fixed(str = paste0(py$adata$obs[["study_id"]], "___", py$adata$obs_names$values),
                                                                                 pattern = c("Zou_Liu_DevelopmentalCell_2021", "HRS118996", "HRS118997", "HRS118998", "HRS118999", "HRS119000", "HRS119001", "HRS119002", "HRS119003", "HRS119004"),
                                                                                 replacement = c("Zou_Liu_DevelopmentalCell_2020", "Skin-Y-18", "Skin-Y-22", "Skin-Y-23", "Skin-M-44", "Skin-M-47", "Skin-M-48", "Skin-O-70", "Skin-O-73", "Skin-O-76"),
                                                                                 vectorize_all = FALSE)]
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

