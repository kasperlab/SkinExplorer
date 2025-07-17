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
py$output_dir = "D:/Study/KI/Projects/###/MetaStudies/Data/Skin/Human/data_preparation_HCA_v1_0/Rindler_Brunner_MolecularCancer_2021"
py_run_string("os.makedirs(output_dir, exist_ok=True)")
h5ad_ALL_file = "GSE173205_ALL.h5ad"
metadata_ALL_file = "GSE173205_ALL_metadata.tsv"
h5ad_HC_file = "GSE173205_HC.h5ad"
metadata_HC_file = "GSE173205_HC_metadata.tsv"
######
input_folder = "D:/Study/KI/Projects/###/MetaStudies/Data/Skin/Human/Rindler_Brunner_MolecularCancer_2021/GSE173205_RAW"
prefix_unique = sapply(list.files(input_folder, pattern = "_matrix.mtx.gz"), function(x){
  str_fragments = unlist(strsplit(x, "_", fixed = T))
  return(paste(str_fragments[-length(str_fragments)], collapse = "_"))
})
lib_unique = sapply(prefix_unique, function(x){
  return(unlist(strsplit(x, "_", fixed = T))[1])
})
#
HCA_metadata_list = load_HCA_metadata_xlsx("D:/Study/KI/Projects/###/MetaStudies/Data/Skin/Human/HCA_metadata/Rindler_2021_HCA_tier 1_metadata.xlsx")
HCA_metadata_list_original = HCA_metadata_list # Keep the original for checking later
HCA_metadata_list_Hao = load_HCA_metadata_xlsx("D:/Study/KI/Projects/###/MetaStudies/Data/Skin/Human/HCA_metadata/Rindler_2021_HCA_tier 1_metadata_Hao.xlsx")
HCA_metadata_list$SP = HCA_metadata_list_Hao$SP
HCA_metadata_list$SP[, "cell_number_loaded"] = HCA_metadata_list_original$SP[, "cell_number_loaded"]
rownames(HCA_metadata_list$SP) = HCA_metadata_list$SP[, "library_id_repository"]
HCA_metadata_list$SP[is.na(HCA_metadata_list$SP)] = "unknown"
rownames(HCA_metadata_list$DN) = HCA_metadata_list$DN[, "donor_id"]
HCA_metadata_list$DN[is.na(HCA_metadata_list$DN)] = "unknown"
#
lib_sample_mapping = HCA_metadata_list$SP[lib_unique, "sample_id"]
names(lib_sample_mapping) = lib_unique
#
lib_id_c = c()
sample_id_running_c = c()
py_run_string("GEX_dict = {}")
for(this_prefix in prefix_unique){
  this_lib = unlist(strsplit(this_prefix, "_", fixed = T))[1]
  this_sample = lib_sample_mapping[this_lib]
  ReadMtx_py("tmp_adata", data_dir = input_folder, prefix = paste0(this_prefix, "_"), feature_column_name = c("gene_ids", "feature_name", "feature_types"))
  py_run_string(paste0("tmp_adata.obs_names = '", this_sample, "_' + tmp_adata.obs_names.astype(str).values"))
  py_run_string(paste0("GEX_dict['", this_sample, "'] = tmp_adata"))
  lib_id_c = c(lib_id_c, rep(this_lib, py$tmp_adata$n_obs))
  sample_id_running_c = c(sample_id_running_c, rep(this_sample, py$tmp_adata$n_obs))
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
pdf(paste0(py$output_dir, '/DonorID_Sex_barplot.pdf'), width = 10, height = 5)
par(mar=c(11, 4, 4, 4))
barplot(Y_X_ratio, las=2)
dev.off()
#
### Study level
py$adata$obs["study_id"] = "Rindler_Brunner_MolecularCancer_2021"
### Experiment level
py$adata$obs["experiment_id"] = as.character(py$adata$obs[["study_id"]])
py$adata$obs["institute"] = HCA_metadata_list$SP[lib_id_c, "institute"]
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
  sample_id_running_c,
  "MF309_patch" ~ "76",
  "MF309_tumor" ~ "76",
  "MF309_follow_up" ~ "76",
  "MF311_patch" ~ "74",
  "MF311_plaque" ~ "74",
  "MF311_follow_up" ~ "75",
  "MF312_patch" ~ "55",
  "MF312_plaque" ~ "55",
  "MF312_follow_up" ~ "56",
  "MF318_patch" ~ "55",
  "MF318_plaque" ~ "55",
  "MF309_nonlesional" ~ "76",
  "MF311_nonlesional" ~ "74",
  "P65_lesional" ~ "53",
  "P65_nonlesional" ~ "53",
  "P73_lesional" ~ "82",
  "P73_nonlesional" ~ "82",
  "P84_lesional" ~ "75",
  "P84_nonlesional" ~ "75",
  "P90_lesional" ~ "75",
  "P90_nonlesional" ~ "75",
  "P107_early_stage" ~ "39",
  "P138_early_stage" ~ "47",
  "P112_HC" ~ "51",
  "P115_HC" ~ "48",
  "P116_HC" ~ "57",
  "P121_HC" ~ "44",
  .default = sample_id_running_c
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
py$adata$obs["disease"] = dplyr::case_match(
  sample_id_running_c,
  "MF309_patch" ~ "mycosis fungoides",
  "MF309_tumor" ~ "mycosis fungoides",
  "MF309_follow_up" ~ "mycosis fungoides",
  "MF311_patch" ~ "mycosis fungoides",
  "MF311_plaque" ~ "mycosis fungoides",
  "MF311_follow_up" ~ "mycosis fungoides",
  "MF312_patch" ~ "mycosis fungoides",
  "MF312_plaque" ~ "mycosis fungoides",
  "MF312_follow_up" ~ "mycosis fungoides",
  "MF318_patch" ~ "mycosis fungoides",
  "MF318_plaque" ~ "mycosis fungoides",
  "MF309_nonlesional" ~ "mycosis fungoides",
  "MF311_nonlesional" ~ "mycosis fungoides",
  "P65_lesional" ~ "mycosis fungoides",
  "P65_nonlesional" ~ "mycosis fungoides",
  "P73_lesional" ~ "folliculotropic mycosis fungoides",
  "P73_nonlesional" ~ "folliculotropic mycosis fungoides",
  "P84_lesional" ~ "mycosis fungoides",
  "P84_nonlesional" ~ "mycosis fungoides",
  "P90_lesional" ~ "mycosis fungoides",
  "P90_nonlesional" ~ "mycosis fungoides",
  "P107_early_stage" ~ "mycosis fungoides",
  "P138_early_stage" ~ "mycosis fungoides",
  "P112_HC" ~ "normal",
  "P115_HC" ~ "normal",
  "P116_HC" ~ "normal",
  "P121_HC" ~ "normal",
  .default = sample_id_running_c
)
py$adata$obs["disease_ontology_term_id"] = HCA_metadata_list$SP[lib_id_c, "disease_ontology_term_id"]
py$adata$obs["disease_status"] = dplyr::case_match(
  sample_id_running_c,
  "MF309_patch" ~ "IVA1(T3N0M0B2)",
  "MF309_tumor" ~ "IVA1(T3N0M0B2)",
  "MF309_follow_up" ~ "IVA1(T3N0M0B2)",
  "MF311_patch" ~ "IVA1(T2N1M0B2)",
  "MF311_plaque" ~ "IVA1(T2N1M0B2)",
  "MF311_follow_up" ~ "IVA1(T4N1M0B2)",
  "MF312_patch" ~ "IIB(T2N0M0B0) - Recurrent",
  "MF312_plaque" ~ "IIB(T2N0M0B0) - Recurrent",
  "MF312_follow_up" ~ "IIB(T1N0M0B0) - Recurrent",
  "MF318_patch" ~ "IIB(T3N0M0B0)",
  "MF318_plaque" ~ "IIB(T3N0M0B0)",
  "MF309_nonlesional" ~ "IVA1(T3N0M0B2)",
  "MF311_nonlesional" ~ "IVA1(T2N1M0B2)",
  "P65_lesional" ~ "IB(T2N0M0B0)",
  "P65_nonlesional" ~ "IB(T2N0M0B0)",
  "P73_lesional" ~ "IIB(T3N0M0B0)",
  "P73_nonlesional" ~ "IIB(T3N0M0B0)",
  "P84_lesional" ~ "IIB(T3N0M0B1)",
  "P84_nonlesional" ~ "IIB(T3N0M0B1)",
  "P90_lesional" ~ "IB(T2N0M0B0)",
  "P90_nonlesional" ~ "IB(T2N0M0B0)",
  "P107_early_stage" ~ "IA(T1N0M0B0)",
  "P138_early_stage" ~ "IA(T1N0M0B0)",
  "P112_HC" ~ "not applicable",
  "P115_HC" ~ "not applicable",
  "P116_HC" ~ "not applicable",
  "P121_HC" ~ "not applicable",
  .default = sample_id_running_c
)
py$adata$obs["treatment_status"] = dplyr::case_match(
  sample_id_running_c,
  "MF309_patch" ~ "Previous:GCS(Topical),IFN alpha;Ongoing:ECP",
  "MF309_tumor" ~ "Previous:GCS(Topical),IFN alpha;Ongoing:ECP",
  "MF309_follow_up" ~ "Previous:Brentuximab;Ongoing:ECP",
  "MF311_patch" ~ "Previous:Alemtuzumab,Bexarotene,GCS(Topical/Systemic),IFN alpha,NB-UVB,PUVA,Radiotherapy;Ongoing:ECP)",
  "MF311_plaque" ~ "Previous:Alemtuzumab,Bexarotene,GCS(Topical/Systemic),IFN alpha,NB-UVB,PUVA,Radiotherapy;Ongoing:ECP)",
  "MF311_follow_up" ~ "Previous:Alemtuzumab,Bexarotene,GCS(Topical/Systemic),IFN alpha,NB-UVB,PUVA,Radiotherapy;Ongoing:ECP)",
  "MF312_patch" ~ "Previous:GCS(Topical),NB-UVB,Radiotherapy;Ongoing:Acitretin,ECP,IFN alpha",
  "MF312_plaque" ~ "Previous:GCS(Topical),NB-UVB,Radiotherapy;Ongoing:Acitretin,ECP,IFN alpha",
  "MF312_follow_up" ~ "Previous:GCS(Topical),NB-UVB,Radiotherapy;Ongoing:Chlormethine,ECP",
  "MF318_patch" ~ "Previous:CNI,GCS(Topical),NB-UVB;Ongoing:No",
  "MF318_plaque" ~ "Previous:CNI,GCS(Topical),NB-UVB;Ongoing:No",
  "MF309_nonlesional" ~ "Previous:Brentuximab;Ongoing:ECP",
  "MF311_nonlesional" ~ "Previous:Alemtuzumab,Bexarotene,GCS(Topical/Systemic),IFN alpha,NB-UVB,PUVA,Radiotherapy;Ongoing:ECP)",
  "P65_lesional" ~ "Previous:GCS(Topical);Ongoing:No",
  "P65_nonlesional" ~ "Previous:GCS(Topical);Ongoing:No",
  "P73_lesional" ~ "Previous:GCS(Topical),PUVA,Re-PUVA,Radiotherapy;Ongoing:No",
  "P73_nonlesional" ~ "Previous:GCS(Topical),PUVA,Re-PUVA,Radiotherapy;Ongoing:No",
  "P84_lesional" ~ "Previous:GCS(Topical),NB-UVB,;Ongoing:No",
  "P84_nonlesional" ~ "Previous:GCS(Topical),NB-UVB,;Ongoing:PUVA",
  "P90_lesional" ~ "Previous:Acitretin,Bexarotene,GCS(Topical),NB-UVB,PUVA;Ongoing:No",
  "P90_nonlesional" ~ "Previous:Acitretin,Bexarotene,GCS(Topical),NB-UVB,PUVA;Ongoing:No",
  "P107_early_stage" ~ "Previous:GCS(Topical);Ongoing:No",
  "P138_early_stage" ~ "Previous:GCS(Topical);Ongoing:No",
  "P112_HC" ~ "No",
  "P115_HC" ~ "No",
  "P116_HC" ~ "No",
  "P121_HC" ~ "No",
  .default = sample_id_running_c
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
py$adata$obs["sampled_site_condition"] = HCA_metadata_list$SP[lib_id_c, "sampled_site_condition"]
py$adata$obs["anatomical_region_level_1"] = dplyr::case_match(
  sample_id_running_c,
  "MF309_patch" ~ "Torso",
  "MF309_tumor" ~ "Torso",
  "MF309_follow_up" ~ "Torso",
  "MF311_patch" ~ "Torso",
  "MF311_plaque" ~ "Torso",
  "MF311_follow_up" ~ "Extremities",
  "MF312_patch" ~ "Extremities",
  "MF312_plaque" ~ "Extremities",
  "MF312_follow_up" ~ "Extremities",
  "MF318_patch" ~ "Extremities",
  "MF318_plaque" ~ "Extremities",
  "MF309_nonlesional" ~ "Torso",
  "MF311_nonlesional" ~ "Extremities",
  "P65_lesional" ~ "Extremities",
  "P65_nonlesional" ~ "Extremities",
  "P73_lesional" ~ "Extremities",
  "P73_nonlesional" ~ "Extremities",
  "P84_lesional" ~ "Extremities",
  "P84_nonlesional" ~ "Extremities",
  "P90_lesional" ~ "Extremities",
  "P90_nonlesional" ~ "Extremities",
  "P107_early_stage" ~ "Torso",
  "P138_early_stage" ~ "Extremities",
  "P112_HC" ~ "unknown",
  "P115_HC" ~ "unknown",
  "P116_HC" ~ "unknown",
  "P121_HC" ~ "unknown",
  .default = sample_id_running_c
)
py$adata$obs["anatomical_region_level_2"] = dplyr::case_match(
  sample_id_running_c,
  "MF309_patch" ~ "unknown",
  "MF309_tumor" ~ "unknown",
  "MF309_follow_up" ~ "unknown",
  "MF311_patch" ~ "unknown",
  "MF311_plaque" ~ "unknown ",
  "MF311_follow_up" ~ "Arm",
  "MF312_patch" ~ "Arm",
  "MF312_plaque" ~ "Arm",
  "MF312_follow_up" ~ "Arm",
  "MF318_patch" ~ "Leg",
  "MF318_plaque" ~ "Leg",
  "MF309_nonlesional" ~ "unknown",
  "MF311_nonlesional" ~ "Leg",
  "P65_lesional" ~ "Arm",
  "P65_nonlesional" ~ "Arm",
  "P73_lesional" ~ "Leg",
  "P73_nonlesional" ~ "Leg",
  "P84_lesional" ~ "Arm",
  "P84_nonlesional" ~ "Arm",
  "P90_lesional" ~ "Arm",
  "P90_nonlesional" ~ "Arm",
  "P107_early_stage" ~ "Back",
  "P138_early_stage" ~ "Leg",
  "P112_HC" ~ "unknown",
  "P115_HC" ~ "unknown",
  "P116_HC" ~ "unknown",
  "P121_HC" ~ "unknown",
  .default = sample_id_running_c
)
py$adata$obs["anatomical_region_level_3"] = dplyr::case_match(
  sample_id_running_c,
  "MF309_patch" ~ "unknown",
  "MF309_tumor" ~ "unknown",
  "MF309_follow_up" ~ "unknown",
  "MF311_patch" ~ "unknown",
  "MF311_plaque" ~ "unknown ",
  "MF311_follow_up" ~ "unknown",
  "MF312_patch" ~ "unknown",
  "MF312_plaque" ~ "unknown",
  "MF312_follow_up" ~ "unknown",
  "MF318_patch" ~ "Thigh",
  "MF318_plaque" ~ "Thigh",
  "MF309_nonlesional" ~ "unknown",
  "MF311_nonlesional" ~ "unknown",
  "P65_lesional" ~ "unknown",
  "P65_nonlesional" ~ "unknown",
  "P73_lesional" ~ "unknown",
  "P73_nonlesional" ~ "unknown",
  "P84_lesional" ~ "unknown",
  "P84_nonlesional" ~ "unknown",
  "P90_lesional" ~ "unknown",
  "P90_nonlesional" ~ "unknown",
  "P107_early_stage" ~ "Lower back",
  "P138_early_stage" ~ "Thigh",
  "P112_HC" ~ "unknown",
  "P115_HC" ~ "unknown",
  "P116_HC" ~ "unknown",
  "P121_HC" ~ "unknown",
  .default = sample_id_running_c
)
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
py$adata$obs["library_id"] = HCA_metadata_list$SP[lib_id_c, "library_id"]
py$adata$obs["library_id_repository"] = lib_id_c
py$adata$obs["library_preparation_batch"] = HCA_metadata_list$SP[lib_id_c, "library_preparation_batch"]
py$adata$obs["library_sequencing_run"] = HCA_metadata_list$SP[lib_id_c, "library_sequencing_run"]
py$adata$obs["sample_preservation_method"] = HCA_metadata_list$SP[lib_id_c, "sample_preservation_method"]
py$adata$obs["protocol_url"] = "unknown"
py$adata$obs["dissociation_protocol"] = HCA_metadata_list$SP[lib_id_c, "dissociation_protocol"]
py$adata$obs["cell_enrichment"] = HCA_metadata_list$SP[lib_id_c, "cell_enrichment"]
py$adata$obs["cell_viability_percentage"] = HCA_metadata_list$SP[lib_id_c, "cell_viability_percentage"]
py$adata$obs["cell_number_loaded"] = HCA_metadata_list$SP[lib_id_c, "cell_number_loaded"]
py$adata$obs["suspension_type"] = HCA_metadata_list$SP[lib_id_c, "suspension_type"]
py$adata$obs["assay_ontology_term"] = "10x 5' v1"
py$adata$obs["assay_ontology_term_id"] = "EFO:0011025"
py$adata$obs["sequenced_fragment"] = "5 prime tag"
py$adata$obs["sequencing_platform"] = "Illumina NovaSeq 6000"
### Preprocessing level
py$adata$obs["reference_genome"] = "GRCh38"
py$adata$obs["alignment_software"] = "cell ranger 3.0.2"
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
this_dataset_individual = "Rindler 2021"
this_dataset_integrated = "Rindler_Brunner_MolecularCancer_2021"
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
individual_barcode_all = individual_barcode_used[paste0(py$adata$obs[["study_id"]], "___", stringi::stri_replace_all_fixed(str = py$adata$obs_names$values,
                                                                                                                           pattern = c("P112_HC", "P115_HC", "P116_HC", "P121_HC"),
                                                                                                                           replacement = c("P112", "P115", "P116", "P121"),
                                                                                                                           vectorize_all = FALSE))]
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
integrated_barcode_all = integrated_barcode_used[paste0(py$adata$obs[["study_id"]], "___", stringi::stri_replace_all_fixed(str = py$adata$obs_names$values,
                                                                                                                           pattern = c("P112_HC", "P115_HC", "P116_HC", "P121_HC"),
                                                                                                                           replacement = c("P112", "P115", "P116", "P121"),
                                                                                                                           vectorize_all = FALSE))]
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
py_run_string("adata_HC = adata[adata.obs['sampled_site_condition'] == 'healthy']")
py_run_string(paste0("adata_HC.write(output_dir + '/", h5ad_HC_file, "')"))
py_run_string(paste0("adata_HC.obs.to_csv(output_dir + '/", metadata_HC_file, "', sep='\t', index=True)"))
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

