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
IntegratedAnnotation_2_metadata = readRDS("D:/Study/KI/Projects/###/MetaStudies/Analysis/Integration/0_preprocessing/IntegratedAnnotation_2_metadata.rds")
IndividualAnalysis_metadata_list = readRDS("D:/Study/KI/Projects/###/MetaStudies/Analysis/Integration/0_preprocessing/IndividualAnalysis_metadata_list.rds")
#
# IO - Setting
#
py$output_dir = "D:/Study/KI/Projects/###/MetaStudies/Data/Skin/Human/data_preparation_HCA_v1_0/Belote_Torres_NatureCellBiology_2021"
py_run_string("os.makedirs(output_dir, exist_ok=True)")
h5ad_ALL_file = "GSE151091_ALL.h5ad"
metadata_ALL_file = "GSE151091_ALL_metadata.tsv"
h5ad_POSTPARTUM_file = "GSE151091_POSTPARTUM.h5ad"
metadata_POSTPARTUM_file = "GSE151091_POSTPARTUM_metadata.tsv"
h5ad_Belote_file = "GSE151091_Belote.h5ad"
metadata_Belote_file = "GSE151091_Belote_metadata.tsv"
######
csv_path = "D:/Study/KI/Projects/###/MetaStudies/Data/Skin/Human/Belote_Torres_NatureCellBiology_2021/GSE151091_raw_matrix.csv.gz"
######
this_df = as.data.frame(fread(csv_path))
rownames(this_df) = this_df[, 1]
this_obj = CreateSeuratObject(as.matrix(this_df[, -1]), min.cells = 1)
rm(this_df)
raw_counts = t(this_obj[["RNA"]]$counts)
rm(this_obj)
py$adata_raw = py$ad$AnnData(X=raw_counts, obs=data.frame(row.names = rownames(raw_counts)), var=data.frame(row.names = colnames(raw_counts)))
rm(raw_counts)
GeneName = featureMapMerge("adata_raw", "adata_genename", hgnc_mapping_alias_genename)
py_run_string("del adata_raw")
GeneName_to_GeneID = mapping_GeneName_GeneID(GeneName, hgnc_mapping_genename_ensemblID)
GeneName_to_GeneID = mapping_GeneName_GeneID(GeneName_to_GeneID, GRCh38v103_mapping_genename_ensemblID)
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
GSE151091_metadata_path = "D:/Study/KI/Projects/###/MetaStudies/Data/Skin/Human/Belote_Torres_NatureCellBiology_2021/GSE151091_Metadata.csv.gz"
metadata_raw = as.matrix(read.table(GSE151091_metadata_path, header = T, row.names = 1, sep = ","))
GSE151091_HighResLouvain_path = "D:/Study/KI/Projects/###/MetaStudies/Data/Skin/Human/Belote_Torres_NatureCellBiology_2021/GSE151091_HighResLouvain.csv.gz"
HighResLouvain = as.matrix(read.table(GSE151091_HighResLouvain_path, header = T, row.names = 2, sep = ","))
GSE151091_LowResLouvain_path = "D:/Study/KI/Projects/###/MetaStudies/Data/Skin/Human/Belote_Torres_NatureCellBiology_2021/GSE151091_LowResLouvain.csv.gz"
LowResLouvain = as.matrix(read.table(GSE151091_LowResLouvain_path, header = T, row.names = 2, sep = ","))
#
donor_id_running_c = sapply(py$adata$obs_names$values, function(x){
  str_fragments = unlist(strsplit(x, "_"))
  return(paste(str_fragments[seq(length(str_fragments) - 2)], collapse = "_"))
})
#
anatomical_region_level_1_c = dplyr::case_match(
  donor_id_running_c,
  "9.5WK02" ~ "Extremities", # leg,arm
  "10WK03" ~ "Extremities", # arm,palm,leg
  "12WKM01" ~ "Extremities", # sole,palm,leg
  "12WK05" ~ "Extremities", # sole,leg,palm,arm
  "16WKM04" ~ "Extremities", # palm
  "18WKM06" ~ "Extremities", # sole,leg,palm
  "FS030_LM" ~ "Torso", # foreskin
  "FS043_LM" ~ "Torso", # foreskin
  "A1021M" ~ "Extremities", # ankle or upper foot
  "A1038LM" ~ "Extremities", # arch,calf,heel,shin
  "A1012M" ~ "Extremities", # arm
  "A1022M" ~ "Extremities", # arm
  "A1015LM" ~ "Extremities", # arm
  "A1025L" ~ "Extremities", # left thigh
  "A1016LM" ~ "Extremities", # anterior forearm
  "A1020LM" ~ "Extremities", # arm
  "A1033M" ~ "Extremities", # thigh
  "A1011L" ~ "Extremities", # arm
  "A1026L" ~ "Extremities", # right thigh
  "A1014L" ~ "Extremities", # thigh
  "A1046M" ~ "Extremities", # shin,heel,arch,calf
  "A1017LM" ~ "Extremities", # medial calf
  .default = donor_id_running_c
)
anatomical_region_level_2_c = dplyr::case_match(
  donor_id_running_c,
  "9.5WK02" ~ "unknown", # leg,arm
  "10WK03" ~ "unknown", # arm,palm,leg
  "12WKM01" ~ "unknown", # sole,palm,leg
  "12WK05" ~ "unknown", # sole,leg,palm,arm
  "16WKM04" ~ "Hand", # palm
  "18WKM06" ~ "unknown", # sole,leg,palm
  "FS030_LM" ~ "Genitalia", # foreskin
  "FS043_LM" ~ "Genitalia", # foreskin
  "A1021M" ~ "Foot", # ankle or upper foot
  "A1038LM" ~ "unknown", # arch,calf,heel,shin
  "A1012M" ~ "Arm", # arm
  "A1022M" ~ "Arm", # arm
  "A1015LM" ~ "Arm", # arm
  "A1025L" ~ "Leg", # left thigh
  "A1016LM" ~ "Arm", # anterior forearm
  "A1020LM" ~ "Arm", # arm
  "A1033M" ~ "Leg", # thigh
  "A1011L" ~ "Arm", # arm
  "A1026L" ~ "Leg", # right thigh
  "A1014L" ~ "Leg", # thigh
  "A1046M" ~ "unknown", # shin,heel,arch,calf
  "A1017LM" ~ "Leg", # medial calf
  .default = donor_id_running_c
)
anatomical_region_level_3_c = dplyr::case_match(
  donor_id_running_c,
  "9.5WK02" ~ "unknown", # leg,arm
  "10WK03" ~ "unknown", # arm,palm,leg
  "12WKM01" ~ "unknown", # sole,palm,leg
  "12WK05" ~ "unknown", # sole,leg,palm,arm
  "16WKM04" ~ "Palm", # palm
  "18WKM06" ~ "unknown", # sole,leg,palm
  "FS030_LM" ~ "Foreskin", # foreskin
  "FS043_LM" ~ "Foreskin", # foreskin
  "A1021M" ~ "unknown", # ankle or upper foot
  "A1038LM" ~ "unknown", # arch,calf,heel,shin
  "A1012M" ~ "unknown", # arm
  "A1022M" ~ "unknown", # arm
  "A1015LM" ~ "unknown", # arm
  "A1025L" ~ "Thigh", # left thigh
  "A1016LM" ~ "Forearm", # anterior forearm
  "A1020LM" ~ "unknown", # arm
  "A1033M" ~ "Thigh", # thigh
  "A1011L" ~ "unknown", # arm
  "A1026L" ~ "Thigh", # right thigh
  "A1014L" ~ "Thigh", # thigh
  "A1046M" ~ "unknown", # shin,heel,arch,calf
  "A1017LM" ~ "Calf", # medial calf
  .default = donor_id_running_c
)
#
metadata_anatomical_region_level_1_c = dplyr::case_match(
  metadata_raw[, "anatomical_location"],
  "ankle or upper foot" ~ "Extremities",
  "anterior forearm" ~ "Extremities",
  "arch" ~ "Extremities",
  "arm" ~ "Extremities",
  "calf" ~ "Extremities",
  "foreskin" ~ "Torso",
  "heel" ~ "Extremities",
  "left thigh" ~ "Extremities",
  "leg" ~ "Extremities",
  "medial calf" ~ "Extremities",
  "palm" ~ "Extremities",
  "right thigh" ~ "Extremities",
  "shin" ~ "Extremities",
  "sole" ~ "Extremities",
  "thigh" ~ "Extremities",
  .default = metadata_raw[, "anatomical_location"]
)
metadata_anatomical_region_level_2_c = dplyr::case_match(
  metadata_raw[, "anatomical_location"],
  "ankle or upper foot" ~ "Foot",
  "anterior forearm" ~ "Arm",
  "arch" ~ "Foot",
  "arm" ~ "Arm",
  "calf" ~ "Leg",
  "foreskin" ~ "Foreskin",
  "heel" ~ "Foot",
  "left thigh" ~ "Leg",
  "leg" ~ "Leg",
  "medial calf" ~ "Leg",
  "palm" ~ "Hand",
  "right thigh" ~ "Leg",
  "shin" ~ "Leg",
  "sole" ~ "Foot",
  "thigh" ~ "Leg",
  .default = metadata_raw[, "anatomical_location"]
)
metadata_anatomical_region_level_3_c = dplyr::case_match(
  metadata_raw[, "anatomical_location"],
  "ankle or upper foot" ~ "unknown",
  "anterior forearm" ~ "Forearm",
  "arch" ~ "Plantar",
  "arm" ~ "unknown",
  "calf" ~ "Calf",
  "foreskin" ~ "not applicable",
  "heel" ~ "Heel",
  "left thigh" ~ "Thigh",
  "leg" ~ "unknown",
  "medial calf" ~ "Calf",
  "palm" ~ "Palm",
  "right thigh" ~ "Thigh",
  "shin" ~ "Shin",
  "sole" ~ "Plantar",
  "thigh" ~ "Thigh",
  .default = metadata_raw[, "anatomical_location"]
)
metadata_cell_enrichment_c = stringi::stri_replace_all_fixed(str = metadata_raw[, "gate_label"],
                                                             pattern = c("a6high", "a6low", "a6mid", "ckitpos", "null"),
                                                             replacement = c("CD11c−KIT−ITGA6_High", "CD11c−KIT−ITGA6_Low", "CD11c−KIT−ITGA6_Middle", "CD11c−KIT+", "No"),
                                                             vectorize_all = FALSE)
#
names(anatomical_region_level_1_c) = py$adata$obs_names$values
names(metadata_anatomical_region_level_1_c) = rownames(metadata_raw)
anatomical_region_level_1_c[rownames(metadata_raw)] = metadata_anatomical_region_level_1_c
names(anatomical_region_level_2_c) = py$adata$obs_names$values
names(metadata_anatomical_region_level_2_c) = rownames(metadata_raw)
anatomical_region_level_2_c[rownames(metadata_raw)] = metadata_anatomical_region_level_2_c
names(anatomical_region_level_3_c) = py$adata$obs_names$values
names(metadata_anatomical_region_level_3_c) = rownames(metadata_raw)
anatomical_region_level_3_c[rownames(metadata_raw)] = metadata_anatomical_region_level_3_c
#
names(metadata_cell_enrichment_c) = rownames(metadata_raw)
cell_enrichment_c = metadata_cell_enrichment_c[py$adata$obs_names$values]
cell_enrichment_c[is.na(cell_enrichment_c)] = "unknown"
#
author_cell_type_1 = metadata_raw[, "class_3"][py$adata$obs_names$values]
author_cell_type_1[is.na(author_cell_type_1)] = "unknown"
#
HCA_metadata_list = load_HCA_metadata_xlsx("D:/Study/KI/Projects/###/MetaStudies/Data/Skin/Human/HCA_metadata/Belote_2021_HCA_tier 1_metadata.xlsx")
HCA_metadata_list_original = HCA_metadata_list # Keep the original for checking later
HCA_metadata_list_DN = HCA_metadata_list$DN[, -1]
rownames(HCA_metadata_list_DN) = HCA_metadata_list$DN[, 1]
HCA_metadata_list_LIB = matrix(nrow = 0, ncol = ncol(HCA_metadata_list$DN) + ncol(HCA_metadata_list$SP), dimnames = list(c(), c(colnames(HCA_metadata_list$DN), colnames(HCA_metadata_list$SP))))
for(ii in seq(nrow(HCA_metadata_list$SP))){
  this_SP_MD = HCA_metadata_list$SP[ii, ]
  this_LIB = unlist(strsplit(this_SP_MD$library_id, "||", fixed = T))
  this_LIB_SP_MD = this_SP_MD[rep(1, length(this_LIB)), ]
  this_LIB_DN_MD = HCA_metadata_list_DN[this_LIB_SP_MD[, "donor_id"], ]
  this_LIB_MD = cbind(this_LIB_DN_MD, this_LIB_SP_MD)
  rownames(this_LIB_MD) = this_LIB
  HCA_metadata_list_LIB = rbind(HCA_metadata_list_LIB, this_LIB_MD)
}
HCA_metadata_list_LIB_reordered = HCA_metadata_list_LIB[py$adata$obs_names$values, ]
HCA_metadata_list_LIB_reordered[is.na(HCA_metadata_list_LIB_reordered)] = "unknown"
rownames(HCA_metadata_list_LIB_reordered) = py$adata$obs_names$values
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
Y_ratio_aggregate = aggregate(Y_ratio, by = list(donor_id_running_c), FUN = function(x){
  return(mean(x, na.rm = T))
})
X_ratio_aggregate = aggregate(X_ratio, by = list(donor_id_running_c), FUN = function(x){
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
py$adata$obs["study_id"] = "Belote_Torres_NatureCellBiology_2021"
### Experiment level
py$adata$obs["experiment_id"] = as.character(py$adata$obs[["study_id"]])
py$adata$obs["institute"] = "unknown"
py$adata$obs["author_batch_notes"] = "unknown"
### Donor background level
py$adata$obs["organism"] = "Homo sapiens"
py$adata$obs["organism_ontology_term_id"] = "NCBITaxon:9606"
py$adata$obs["donor_id"] = HCA_metadata_list_LIB_reordered[, "donor_id"]
py$adata$obs["donor_id_running"] = paste0(py$adata$obs[["experiment_id"]], "___", donor_id_running_c)
py$adata$obs["sex_ontology_term"] = dplyr::case_match(
  donor_id_running_c,
  "9.5WK02" ~ "unknown",
  "10WK03" ~ "unknown",
  "12WKM01" ~ "male",
  "12WK05" ~ "unknown",
  "16WKM04" ~ "unknown",
  "18WKM06" ~ "unknown",
  "FS030_LM" ~ "male",
  "FS043_LM" ~ "male",
  "A1021M" ~ "female",
  "A1038LM" ~ "male",
  "A1012M" ~ "female",
  "A1022M" ~ "male",
  "A1015LM" ~ "female",
  "A1025L" ~ "male",
  "A1016LM" ~ "male",
  "A1020LM" ~ "female",
  "A1033M" ~ "male",
  "A1011L" ~ "male",
  "A1026L" ~ "male",
  "A1014L" ~ "female",
  "A1046M" ~ "female",
  "A1017LM" ~ "male",
  .default = donor_id_running_c
)
py$adata$obs["sex_ontology_term_id"] = dplyr::case_match(
  donor_id_running_c,
  "9.5WK02" ~ "unknown",
  "10WK03" ~ "unknown",
  "12WKM01" ~ "PATO:0000384",
  "12WK05" ~ "unknown",
  "16WKM04" ~ "unknown",
  "18WKM06" ~ "unknown",
  "FS030_LM" ~ "PATO:0000384",
  "FS043_LM" ~ "PATO:0000384",
  "A1021M" ~ "PATO:0000383",
  "A1038LM" ~ "PATO:0000384",
  "A1012M" ~ "PATO:0000383",
  "A1022M" ~ "PATO:0000384",
  "A1015LM" ~ "PATO:0000383",
  "A1025L" ~ "PATO:0000384",
  "A1016LM" ~ "PATO:0000384",
  "A1020LM" ~ "PATO:0000383",
  "A1033M" ~ "PATO:0000384",
  "A1011L" ~ "PATO:0000384",
  "A1026L" ~ "PATO:0000384",
  "A1014L" ~ "PATO:0000383",
  "A1046M" ~ "PATO:0000383",
  "A1017LM" ~ "PATO:0000384",
  .default = donor_id_running_c
)
py$adata$obs["ethnicity_1"] = dplyr::case_match(
  donor_id_running_c,
  "9.5WK02" ~ "unknown",
  "10WK03" ~ "unknown",
  "12WKM01" ~ "unknown",
  "12WK05" ~ "unknown",
  "16WKM04" ~ "unknown",
  "18WKM06" ~ "unknown",
  "FS030_LM" ~ "unknown",
  "FS043_LM" ~ "unknown",
  "A1021M" ~ "Asian ancestry",
  "A1038LM" ~ "European ancestry",
  "A1012M" ~ "Asian ancestry",
  "A1022M" ~ "Latin or Admixed American ancestry",
  "A1015LM" ~ "Asian ancestry",
  "A1025L" ~ "European ancestry",
  "A1016LM" ~ "European ancestry",
  "A1020LM" ~ "Asian ancestry",
  "A1033M" ~ "Latin or Admixed American ancestry",
  "A1011L" ~ "European ancestry",
  "A1026L" ~ "European ancestry",
  "A1014L" ~ "European ancestry",
  "A1046M" ~ "Latin or Admixed American ancestry",
  "A1017LM" ~ "Latin or Admixed American ancestry",
  .default = donor_id_running_c
)
py$adata$obs["ethnicity_2"] = "unknown"
py$adata$obs["ethnicity_details"] = dplyr::case_match(
  donor_id_running_c,
  "9.5WK02" ~ "unknown",
  "10WK03" ~ "unknown",
  "12WKM01" ~ "unknown",
  "12WK05" ~ "unknown",
  "16WKM04" ~ "unknown",
  "18WKM06" ~ "unknown",
  "FS030_LM" ~ "unknown",
  "FS043_LM" ~ "unknown",
  "A1021M" ~ "Asian",
  "A1038LM" ~ "White",
  "A1012M" ~ "Asian",
  "A1022M" ~ "Hispanic/Latinx",
  "A1015LM" ~ "Asian",
  "A1025L" ~ "White",
  "A1016LM" ~ "White",
  "A1020LM" ~ "Asian",
  "A1033M" ~ "Hispanic/Latinx",
  "A1011L" ~ "White",
  "A1026L" ~ "White",
  "A1014L" ~ "White",
  "A1046M" ~ "Hispanic/Latinx",
  "A1017LM" ~ "Hispanic/Latinx",
  .default = donor_id_running_c
)
py$adata$obs["self_reported_ethnicity_ontology_term_id"] = "unknown"
py$adata$obs["genotype"] = "unknown"
### Donor status level
py$adata$obs["age_years"] = dplyr::case_match(
  donor_id_running_c,
  "9.5WK02" ~ "FW_9.5",
  "10WK03" ~ "FW_10",
  "12WKM01" ~ "FW_12",
  "12WK05" ~ "FW_12",
  "16WKM04" ~ "FW_16",
  "18WKM06" ~ "FW_18",
  "FS030_LM" ~ "0",
  "FS043_LM" ~ "0",
  "A1021M" ~ "24",
  "A1038LM" ~ "35",
  "A1012M" ~ "37",
  "A1022M" ~ "42",
  "A1015LM" ~ "52",
  "A1025L" ~ "56",
  "A1016LM" ~ "58",
  "A1020LM" ~ "60",
  "A1033M" ~ "61",
  "A1011L" ~ "65",
  "A1026L" ~ "66",
  "A1014L" ~ "68",
  "A1046M" ~ "77",
  "A1017LM" ~ "81",
  .default = donor_id_running_c
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
development_stage_ontology_term_id_HCA = HCA_metadata_list_LIB_reordered[, "development_stage_ontology_term_id"]
mask_unknown = development_stage_ontology_term_id_HCA == "unknown"
development_stage_ontology_term_id_HCA[mask_unknown] = development_stage_ontology_term_id[mask_unknown]
py$adata$obs["development_stage_ontology_term_id"] = development_stage_ontology_term_id_HCA
py$adata$obs["disease"] = "normal"
py$adata$obs["disease_ontology_term_id"] = "PATO:0000461"
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
py$adata$obs["skin_tone"] = dplyr::case_match(
  donor_id_running_c,
  "9.5WK02" ~ "unknown",
  "10WK03" ~ "unknown",
  "12WKM01" ~ "unknown",
  "12WK05" ~ "unknown",
  "16WKM04" ~ "unknown",
  "18WKM06" ~ "unknown",
  "FS030_LM" ~ "4",
  "FS043_LM" ~ "4",
  "A1021M" ~ "5",
  "A1038LM" ~ "4",
  "A1012M" ~ "5",
  "A1022M" ~ "5",
  "A1015LM" ~ "4",
  "A1025L" ~ "3",
  "A1016LM" ~ "4",
  "A1020LM" ~ "4",
  "A1033M" ~ "5",
  "A1011L" ~ "3",
  "A1026L" ~ "3",
  "A1014L" ~ "3",
  "A1046M" ~ "5",
  "A1017LM" ~ "4",
  .default = donor_id_running_c
)
# according to https://www.researchgate.net/figure/Skin-type-and-tanning-ability-based-on-the-Fitzpatrick-skin-pigmentation-scale_tbl1_285056396
# L = 3, LM = 4, M = 5
py$adata$obs["skin_care"] = "unknown"
py$adata$obs["manner_of_death"] = dplyr::case_when(
  grepl("FW", age_years_character, fixed = T) ~ "unknown",
  .default = "not applicable"
)
### Sample level
py$adata$obs["sample_id"] = HCA_metadata_list_LIB_reordered[, "sample_id"]
py$adata$obs["sample_id_running"] = as.character(py$adata$obs[["donor_id_running"]])
py$adata$obs["sample_source"] = dplyr::case_when(
  grepl("FW", age_years_character, fixed = T) ~ "postmortem donor",
  .default = "surgical donor"
)
py$adata$obs["sample_collection_method"] = "surgical resection"
py$adata$obs["sampled_site_condition"] = "healthy"
py$adata$obs["anatomical_region_level_1"] = anatomical_region_level_1_c
py$adata$obs["anatomical_region_level_2"] = anatomical_region_level_2_c
py$adata$obs["anatomical_region_level_3"] = anatomical_region_level_3_c
metadata_AR = as.character(py$adata$obs[["anatomical_region_level_3"]])
metadata_AR[metadata_AR %in% c("unknown", "not applicable")] = NA
metadata_AR[is.na(metadata_AR)] = as.character(py$adata$obs[["anatomical_region_level_2"]])[is.na(metadata_AR)]
metadata_AR[metadata_AR %in% c("unknown", "not applicable")] = NA
metadata_AR[is.na(metadata_AR)] = as.character(py$adata$obs[["anatomical_region_level_1"]])[is.na(metadata_AR)]
metadata_AR[metadata_AR %in% c("unknown", "not applicable")] = NA
metadata_AR[is.na(metadata_AR)] = "unknown"
py$adata$obs["anatomical_region"] = metadata_AR
py$adata$obs["tissue_type"] = "tissue"
py$adata$obs["tissue_free_text"] = HCA_metadata_list_LIB_reordered[, "tissue_free_text"]
py$adata$obs["tissue_ontology_term"] = HCA_metadata_list_LIB_reordered[, "tissue_ontology_term"]
py$adata$obs["tissue_ontology_term_id"] = HCA_metadata_list_LIB_reordered[, "tissue_ontology_term_id"]
py$adata$obs["skin_tissue"] = "Epidermis"
py$adata$obs["sample_collection_site"] = "unknown"
py$adata$obs["sample_cultured"] = "No"
### Sequencing level
py$adata$obs["sequencing_type"] = "GEX"
py$adata$obs["sample_collection_year"] = "unknown"
py$adata$obs["sample_collection_relative_time_point"] = "unknown"
py$adata$obs["library_id"] = py$adata$obs_names$values
py$adata$obs["library_id_repository"] = "unknown"
py$adata$obs["library_preparation_batch"] = HCA_metadata_list_LIB_reordered[, "library_preparation_batch"]
py$adata$obs["library_sequencing_run"] = "unknown"
py$adata$obs["sample_preservation_method"] = "frozen at -80C" # Plates were stored at −80°C until used for sorting.
py$adata$obs["protocol_url"] = "unknown"
py$adata$obs["dissociation_protocol"] = "unknown"
py$adata$obs["cell_enrichment"] = cell_enrichment_c
py$adata$obs["cell_viability_percentage"] = "unknown"
py$adata$obs["cell_number_loaded"] = "unknown"
py$adata$obs["suspension_type"] = "cell"
py$adata$obs["assay_ontology_term"] = "Smart-seq2"
py$adata$obs["assay_ontology_term_id"] = "EFO:0008931"
py$adata$obs["sequenced_fragment"] = "full length"
py$adata$obs["sequencing_platform"] = "Illumina NovaSeq 6000"
### Preprocessing level
py$adata$obs["reference_genome"] = "GRCh38"
py$adata$obs["alignment_software"] = "cell ranger 3.0.2"
py$adata$obs["intron_inclusion"] = "unknown"
### Analysis Level
py$adata$obs["is_primary_data"] = FALSE
py$adata$obs["author_cell_type_1"] = author_cell_type_1
py$adata$obs["author_cell_type_2"] = "unknown"
py$adata$obs["author_cell_type_3"] = "unknown"
py$adata$obs["author_cell_type_4"] = "unknown"
py$adata$obs["author_cell_type_5"] = "unknown"
py$adata$obs["author_cell_type_6"] = "unknown"
### Previous Barcode
this_dataset_individual = "Belote 2021"
this_dataset_integrated = "Belote_Torres_NatureCellBiology_2021"
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
  return(paste(str_fragments[-length(str_fragments)], collapse = "---"))
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
py$adata_ALL = py$adata[py$adata$obs_names$values %in% paste0("Belote_Torres_NatureCellBiology_2021___", rownames(metadata_raw))]
py$adata_Belote = py$adata[as.character(py$adata$obs[["individual_barcode"]]) != "unknown", ]
#
# Write
#
py_run_string("sc.pp.filter_cells(adata_ALL, min_genes=1)")
py_run_string(paste0("adata_ALL.write(output_dir + '/", h5ad_ALL_file, "')"))
py_run_string(paste0("adata_ALL.obs.to_csv(output_dir + '/", metadata_ALL_file, "', sep='\t', index=True)"))
py_run_string("adata_POSTPARTUM = adata_ALL[adata_ALL.obs['age_years'].str.contains('FW', regex=False)]")
py_run_string(paste0("adata_POSTPARTUM.write(output_dir + '/", h5ad_POSTPARTUM_file, "')"))
py_run_string(paste0("adata_POSTPARTUM.obs.to_csv(output_dir + '/", metadata_POSTPARTUM_file, "', sep='\t', index=True)"))
#
py_run_string("sc.pp.filter_cells(adata_Belote, min_genes=1)")
py_run_string(paste0("adata_Belote.write(output_dir + '/", h5ad_Belote_file, "')"))
py_run_string(paste0("adata_Belote.obs.to_csv(output_dir + '/", metadata_Belote_file, "', sep='\t', index=True)"))
py_run_string("print(output_dir)")
#
metadata_matrix = as.matrix(py$adata_ALL$obs)
write.table(x = HCA_metadata_usage(HCA_metadata=HCA_metadata_list_original$DT, ADATA_metadata=metadata_matrix,
                                   metadata_rowname_column="assay_ontology_term_id"),
            file = paste0(py$output_dir, "/", "HCA_metadata_usage_DT.txt"), sep = "\t", col.names=T, row.names=F, quote=T)
write.table(x = HCA_metadata_usage(HCA_metadata=HCA_metadata_list_original$DN, ADATA_metadata=metadata_matrix,
                                   metadata_rowname_column="donor_id"),
            file = paste0(py$output_dir, "/", "HCA_metadata_usage_DN.txt"), sep = "\t", col.names=T, row.names=F, quote=T)
write.table(x = HCA_metadata_usage(HCA_metadata=HCA_metadata_list_original$SP, ADATA_metadata=metadata_matrix,
                                   metadata_rowname_column="sample_id"),
            file = paste0(py$output_dir, "/", "HCA_metadata_usage_SP.txt"), sep = "\t", col.names=T, row.names=F, quote=T)

