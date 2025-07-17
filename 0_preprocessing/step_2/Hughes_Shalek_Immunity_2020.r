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
py$output_dir = "D:/Study/KI/Projects/###/MetaStudies/Data/Skin/Human/data_preparation_HCA_v1_0/Hughes_Shalek_Immunity_2020"
py_run_string("os.makedirs(output_dir, exist_ok=True)")
h5ad_ALL_file = "GSE150672_ALL.h5ad"
metadata_ALL_file = "GSE150672_ALL_metadata.tsv"
h5ad_HC_file = "GSE150672_HC.h5ad"
metadata_HC_file = "GSE150672_HC_metadata.tsv"
######
csv_path = "D:/Study/KI/Projects/###/MetaStudies/Data/Skin/Human/Hughes_Shalek_Immunity_2020/GSE150672_Skin_Expression_counts.csv.gz"
######
this_df = as.data.frame(fread(csv_path))
rownames(this_df) = this_df[, 1]
this_obj = CreateSeuratObject(as.matrix(this_df[, -1]), min.cells = 1)
rm(this_df)
######
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
lib_id_c = sapply(py$adata$obs_names$values, function(x){
  return(unlist(strsplit(x, "_"))[1])
})
#
GSE150672_metadata_path = "D:/Study/KI/Projects/###/MetaStudies/Data/Skin/Human/Hughes_Shalek_Immunity_2020/1-s2.0-S107476132030409X-mmc3.txt"
metadata_raw = as.matrix(read.table(GSE150672_metadata_path, header = T, row.names = 1, sep = "\t"))[py$adata$obs_names$values, ]
metadata = metadata_raw[py$adata$obs_names$values, ]
author_cell_type_1 = metadata[, "CellType"]
author_cell_type_2 = metadata[, "Specific_CellType"]
author_cell_type_3 = metadata[, "SingleR_Labels"]
umap_all = t(apply(metadata[, c("UMAP1", "UMAP2")], 1, as.numeric))
umap_sample = t(apply(metadata[, c("Sample_UMAP1", "Sample_UMAP2")], 1, as.numeric))
tsne_sample = t(apply(metadata[, c("Sample_TSNE1", "Sample_TSNE2")], 1, as.numeric))
rm(metadata)
#
HCA_metadata_list = load_HCA_metadata_xlsx("D:/Study/KI/Projects/###/MetaStudies/Data/Skin/Human/HCA_metadata/Hughes_2020_HCA_tier 1_metadata.xlsx")
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
HCA_metadata_list_LIB[is.na(HCA_metadata_list_LIB)] = "unknown"
#
donor_id_c = HCA_metadata_list_LIB[lib_id_c, "donor_id"]
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
py$adata$obs["study_id"] = "Hughes_Shalek_Immunity_2020"
### Experiment level
py$adata$obs["experiment_id"] = as.character(py$adata$obs[["study_id"]])
py$adata$obs["institute"] = HCA_metadata_list_LIB[lib_id_c, "institute"]
py$adata$obs["author_batch_notes"] = HCA_metadata_list_LIB[lib_id_c, "author_batch_notes"]
### Donor background level
py$adata$obs["organism"] = "Homo sapiens"
py$adata$obs["organism_ontology_term_id"] = HCA_metadata_list_LIB[lib_id_c, "organism_ontology_term_id"]
py$adata$obs["donor_id"] = donor_id_c
py$adata$obs["donor_id_running"] = paste0(py$adata$obs[["experiment_id"]], "___", donor_id_c)
py$adata$obs["sex_ontology_term"] = HCA_metadata_list_LIB[lib_id_c, "sex_ontology_term"]
py$adata$obs["sex_ontology_term_id"] = HCA_metadata_list_LIB[lib_id_c, "sex_ontology_term_id"]
py$adata$obs["ethnicity_1"] = dplyr::case_match(
  donor_id_c,
  "Acne1" ~ "Latin or Admixed American ancestry",
  "Acne2" ~ "Latin or Admixed American ancestry",
  "Acne3" ~ "European ancestry",
  "Acne4" ~ "European ancestry",
  "Alopecia1" ~ "Asian ancestry/European ancestry",
  "GA1" ~ "European ancestry",
  "GA2" ~ "European ancestry",
  "Leprosy1" ~ "Latin or Admixed American ancestry",
  "Leprosy2" ~ "unknown",
  "Leprosy3" ~ "unknown",
  "Leprosy4" ~ "unknown",
  "Normal1" ~ "European ancestry",
  "Normal2" ~ "unknown",
  "Normal3" ~ "Asian ancestry",
  "Psoriasis1" ~ "European ancestry",
  "Psoriasis2" ~ "European ancestry",
  "Psoriasis3" ~ "European ancestry",
  "Psoriasis4" ~ "European ancestry",
  "Psoriasis5" ~ "European ancestry",
  .default = donor_id_c
)
py$adata$obs["ethnicity_2"] = "unknown"
py$adata$obs["ethnicity_details"] = dplyr::case_match(
  donor_id_c,
  "Acne1" ~ "Hispanic - White",
  "Acne2" ~ "Hispanic - White",
  "Acne3" ~ "Non-hispanic - White",
  "Acne4" ~ "Non-hispanic - White",
  "Alopecia1" ~ "Non-Hispanic - Asian/White",
  "GA1" ~ "Non-hispanic - White",
  "GA2" ~ "Non-hispanic - White",
  "Leprosy1" ~ "Hispanic - White",
  "Leprosy2" ~ "unknown",
  "Leprosy3" ~ "unknown",
  "Leprosy4" ~ "unknown",
  "Normal1" ~ "Non-hispanic - White",
  "Normal2" ~ "unknown",
  "Normal3" ~ "Non-hispanic - Asian",
  "Psoriasis1" ~ "Non-hispanic - White",
  "Psoriasis2" ~ "Non-hispanic - White",
  "Psoriasis3" ~ "Non-hispanic - White",
  "Psoriasis4" ~ "Non-hispanic - White",
  "Psoriasis5" ~ "Non-hispanic - White",
  .default = donor_id_c
)
py$adata$obs["self_reported_ethnicity_ontology_term_id"] = dplyr::case_match(
  donor_id_c,
  "Acne1" ~ "HANCESTRO:0014",
  "Acne2" ~ "HANCESTRO:0014",
  "Acne3" ~ "HANCESTRO:0005",
  "Acne4" ~ "HANCESTRO:0005",
  "Alopecia1" ~ "HANCESTRO:0008/HANCESTRO:0005",
  "GA1" ~ "HANCESTRO:0005",
  "GA2" ~ "HANCESTRO:0005",
  "Leprosy1" ~ "HANCESTRO:0014",
  "Leprosy2" ~ "unknown",
  "Leprosy3" ~ "unknown",
  "Leprosy4" ~ "unknown",
  "Normal1" ~ "HANCESTRO:0005",
  "Normal2" ~ "unknown",
  "Normal3" ~ "HANCESTRO:0008",
  "Psoriasis1" ~ "HANCESTRO:0005",
  "Psoriasis2" ~ "HANCESTRO:0005",
  "Psoriasis3" ~ "HANCESTRO:0005",
  "Psoriasis4" ~ "HANCESTRO:0005",
  "Psoriasis5" ~ "HANCESTRO:0005",
  .default = donor_id_c
)
py$adata$obs["genotype"] = "unknown"
### Donor status level
py$adata$obs["age_years"] = dplyr::case_match(
  donor_id_c,
  "Acne1" ~ "46",
  "Acne2" ~ "29",
  "Acne3" ~ "29",
  "Acne4" ~ "26",
  "Alopecia1" ~ "27",
  "GA1" ~ "61",
  "GA2" ~ "58",
  "Leprosy1" ~ "30",
  "Leprosy2" ~ "59",
  "Leprosy3" ~ "38",
  "Leprosy4" ~ "71",
  "Normal1" ~ "63",
  "Normal2" ~ "66-70",
  "Normal3" ~ "56-60",
  "Psoriasis1" ~ "56",
  "Psoriasis2" ~ "63",
  "Psoriasis3" ~ "61",
  "Psoriasis4" ~ "25",
  "Psoriasis5" ~ "35",
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
#
development_stage_ontology_term_id[age_years_character == "56-60"] = "HsapDv:0000240"
development_stage_ontology_term_id[age_years_character == "66-70"] = "HsapDv:0000241"
#
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
py$adata$obs["development_stage_ontology_term_id"] = HCA_metadata_list_LIB[lib_id_c, "development_stage_ontology_term_id"]
py$adata$obs["disease"] = dplyr::case_match(
  donor_id_c,
  "Acne1" ~ "acne",
  "Acne2" ~ "acne",
  "Acne3" ~ "acne",
  "Acne4" ~ "acne",
  "Alopecia1" ~ "alopecia",
  "GA1" ~ "granuloma annulare",
  "GA2" ~ "granuloma annulare",
  "Leprosy1" ~ "leprosy",
  "Leprosy2" ~ "leprosy",
  "Leprosy3" ~ "leprosy",
  "Leprosy4" ~ "leprosy",
  "Normal1" ~ "normal",
  "Normal2" ~ "normal",
  "Normal3" ~ "normal",
  "Psoriasis1" ~ "psoriasis",
  "Psoriasis2" ~ "psoriasis",
  "Psoriasis3" ~ "psoriasis",
  "Psoriasis4" ~ "psoriasis",
  "Psoriasis5" ~ "psoriasis",
  .default = donor_id_c
)
py$adata$obs["disease_ontology_term_id"] = HCA_metadata_list_LIB[lib_id_c, "disease_ontology_term_id"]
py$adata$obs["disease_status"] = dplyr::case_match(
  donor_id_c,
  "Acne1" ~ "unknown",
  "Acne2" ~ "unknown",
  "Acne3" ~ "unknown",
  "Acne4" ~ "unknown",
  "Alopecia1" ~ "unknown",
  "GA1" ~ "unknown",
  "GA2" ~ "unknown",
  "Leprosy1" ~ "unknown",
  "Leprosy2" ~ "unknown",
  "Leprosy3" ~ "unknown",
  "Leprosy4" ~ "unknown",
  "Normal1" ~ "not applicable",
  "Normal2" ~ "not applicable",
  "Normal3" ~ "not applicable",
  "Psoriasis1" ~ "unknown",
  "Psoriasis2" ~ "unknown",
  "Psoriasis3" ~ "unknown",
  "Psoriasis4" ~ "unknown",
  "Psoriasis5" ~ "unknown",
  .default = donor_id_c
)
py$adata$obs["treatment_status"] = "No"
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
py$adata$obs["manner_of_death"] = HCA_metadata_list_LIB[lib_id_c, "manner_of_death"]
### Sample level
py$adata$obs["sample_id"] = HCA_metadata_list_LIB[lib_id_c, "sample_id"]
py$adata$obs["sample_id_running"] = as.character(py$adata$obs[["donor_id_running"]])
py$adata$obs["sample_source"] = HCA_metadata_list_LIB[lib_id_c, "sample_source"]
py$adata$obs["sample_collection_method"] = HCA_metadata_list_LIB[lib_id_c, "sample_collection_method"]
py$adata$obs["sampled_site_condition"] = HCA_metadata_list_LIB[lib_id_c, "sampled_site_condition"]
py$adata$obs["anatomical_region_level_1"] = "unknown"
py$adata$obs["anatomical_region_level_2"] = "unknown"
py$adata$obs["anatomical_region_level_3"] = "unknown"
metadata_AR = as.character(py$adata$obs[["anatomical_region_level_3"]])
metadata_AR[metadata_AR %in% c("unknown", "not applicable")] = NA
metadata_AR[is.na(metadata_AR)] = as.character(py$adata$obs[["anatomical_region_level_2"]])[is.na(metadata_AR)]
metadata_AR[metadata_AR %in% c("unknown", "not applicable")] = NA
metadata_AR[is.na(metadata_AR)] = as.character(py$adata$obs[["anatomical_region_level_1"]])[is.na(metadata_AR)]
metadata_AR[metadata_AR %in% c("unknown", "not applicable")] = NA
metadata_AR[is.na(metadata_AR)] = "unknown"
py$adata$obs["anatomical_region"] = metadata_AR
py$adata$obs["tissue_type"] = HCA_metadata_list_LIB[lib_id_c, "tissue_type"]
py$adata$obs["tissue_free_text"] = HCA_metadata_list_LIB[lib_id_c, "tissue_free_text"]
py$adata$obs["tissue_ontology_term"] = HCA_metadata_list_LIB[lib_id_c, "tissue_ontology_term"]
py$adata$obs["tissue_ontology_term_id"] = HCA_metadata_list_LIB[lib_id_c, "tissue_ontology_term_id"]
py$adata$obs["skin_tissue"] = "unknown"
py$adata$obs["sample_collection_site"] = HCA_metadata_list_LIB[lib_id_c, "sample_collection_site"]
py$adata$obs["sample_cultured"] = "No"
### Sequencing level
py$adata$obs["sequencing_type"] = "GEX"
py$adata$obs["sample_collection_year"] = HCA_metadata_list_LIB[lib_id_c, "sample_collection_year"]
py$adata$obs["sample_collection_relative_time_point"] = HCA_metadata_list_LIB[lib_id_c, "sample_collection_relative_time_point"]
py$adata$obs["library_id"] = HCA_metadata_list_LIB[lib_id_c, "library_id"]
py$adata$obs["library_id_repository"] = dplyr::case_match(
  donor_id_c,
  "Acne1" ~ "GSM4693360",
  "Acne2" ~ "GSM4693361",
  "Acne3" ~ "GSM4693362",
  "Acne4" ~ "GSM4693363",
  "Alopecia1" ~ "GSM4693364",
  "GA1" ~ "GSM4693365",
  "GA2" ~ "GSM4693366",
  "Leprosy1" ~ "GSM4693367",
  "Leprosy2" ~ "GSM4693368",
  "Leprosy3" ~ "GSM4693369",
  "Leprosy4" ~ "GSM4693370",
  "Normal1" ~ "GSM4693371",
  "Normal2" ~ "GSM4693372",
  "Normal3" ~ "GSM4693373",
  "Psoriasis1" ~ "GSM4693374",
  "Psoriasis2" ~ "GSM4693375",
  "Psoriasis3" ~ "GSM4693376||GSM4693377",
  "Psoriasis4" ~ "GSM4693378||GSM4693379",
  "Psoriasis5" ~ "GSM4693380||GSM4693381",
  .default = donor_id_c
)
py$adata$obs["library_preparation_batch"] = HCA_metadata_list_LIB[lib_id_c, "library_preparation_batch"]
py$adata$obs["library_sequencing_run"] = HCA_metadata_list_LIB[lib_id_c, "library_sequencing_run"]
py$adata$obs["sample_preservation_method"] = HCA_metadata_list_LIB[lib_id_c, "sample_preservation_method"]
py$adata$obs["protocol_url"] = "unknown"
py$adata$obs["dissociation_protocol"] = HCA_metadata_list_LIB[lib_id_c, "dissociation_protocol"]
py$adata$obs["cell_enrichment"] = HCA_metadata_list_LIB[lib_id_c, "cell_enrichment"]
py$adata$obs["cell_viability_percentage"] = HCA_metadata_list_LIB[lib_id_c, "cell_viability_percentage"]
py$adata$obs["cell_number_loaded"] = "unknown"
py$adata$obs["suspension_type"] = "cell"
py$adata$obs["assay_ontology_term"] = "Seq-Well S3"
py$adata$obs["assay_ontology_term_id"] = "EFO:0030019"
py$adata$obs["sequenced_fragment"] = "3 prime tag"
py$adata$obs["sequencing_platform"] = "Illumina NextSeq 500"
### Preprocessing level
py$adata$obs["reference_genome"] = "GRCh37"
py$adata$obs["alignment_software"] = "unknown"
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
this_dataset_individual = "Hughes 2020"
this_dataset_integrated = "Hughes_Shalek_Immunity_2020"
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
py$adata$obsm["Original_UMAP"] = umap_all
py$adata$obsm["Original_UMAP_Sample"] = umap_sample
py$adata$obsm["Original_TSNE_Sample"] = tsne_sample
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

