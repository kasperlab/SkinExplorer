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
py$output_dir = "D:/Study/KI/Projects/###/MetaStudies/Data/Skin/Human/data_preparation_HCA_v1_0/Reynolds_Greenleaf_NatureGenetics_2023"
py_run_string("os.makedirs(output_dir, exist_ok=True)")
h5ad_ALL_file = "GSE212447_ALL.h5ad"
metadata_ALL_file = "GSE212447_ALL_metadata.tsv"
h5ad_HC_file = "GSE212447_HC.h5ad"
metadata_HC_file = "GSE212447_HC_metadata.tsv"
######
input_folder = "D:/Study/KI/Projects/###/MetaStudies/Data/Skin/Human/Reynolds_Greenleaf_NatureGenetics_2023/GSE212447_RAW"
sample_prefixes = c("GSM6532919_AA2", "GSM6532920_AA4", "GSM6532921_AA7", "GSM6532922_AA8", "GSM6532923_C_PB1", "GSM6532924_C_PB2", "GSM6532925_C_PB3", "GSM6532926_C_SD1", "GSM6532927_C_SD2", "GSM6532928_C_SD3")
sample_id_c = c()
library_id_repository_c = c()
py_run_string("GEX_dict = {}")
for(ii in sample_prefixes){
  ii_splits = unlist(strsplit(ii, "_"))
  this_sample = paste(ii_splits[-1], collapse = "_")
  ReadMtx_py("tmp_adata", data_dir = input_folder, prefix = paste0(ii, "_"), feature_column_name = c("gene_ids", "feature_name", "feature_types"))
  py_run_string(paste0("tmp_adata.obs_names = '", this_sample, "_' + tmp_adata.obs_names.astype(str).values"))
  py_run_string(paste0("GEX_dict['", this_sample, "'] = tmp_adata"))
  sample_id_c = c(sample_id_c, rep(this_sample, py$tmp_adata$n_obs))
  library_id_repository_c = c(library_id_repository_c, rep(ii_splits[1], py$tmp_adata$n_obs))
  print(ii)
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
py$adata$obs["study_id"] = "Reynolds_Greenleaf_NatureGenetics_2023"
### Experiment level
py$adata$obs["experiment_id"] = as.character(py$adata$obs[["study_id"]])
py$adata$obs["institute"] = "unknown"
py$adata$obs["author_batch_notes"] = "unknown"
### Donor background level
py$adata$obs["organism"] = "Homo sapiens"
py$adata$obs["organism_ontology_term_id"] = "NCBITaxon:9606"
py$adata$obs["donor_id"] = sample_id_c
py$adata$obs["donor_id_running"] = paste0(py$adata$obs[["experiment_id"]], "___", sample_id_c)
py$adata$obs["sex_ontology_term"] = dplyr::case_match(
  sample_id_c,
  "AA2" ~ "male",
  "AA4" ~ "female",
  "AA7" ~ "female",
  "AA8" ~ "female",
  "C_PB1" ~ "male",
  "C_PB2" ~ "female",
  "C_PB3" ~ "female",
  "C_SD1" ~ "female",
  "C_SD2" ~ "female",
  "C_SD3" ~ "female",
  .default = sample_id_c
) # S Table 1
py$adata$obs["sex_ontology_term_id"] = dplyr::case_match(
  sample_id_c,
  "AA2" ~ "PATO:0000384",
  "AA4" ~ "PATO:0000383",
  "AA7" ~ "PATO:0000383",
  "AA8" ~ "PATO:0000383",
  "C_PB1" ~ "PATO:0000384",
  "C_PB2" ~ "PATO:0000383",
  "C_PB3" ~ "PATO:0000383",
  "C_SD1" ~ "PATO:0000383",
  "C_SD2" ~ "PATO:0000383",
  "C_SD3" ~ "PATO:0000383",
  .default = sample_id_c
)
py$adata$obs["ethnicity_1"] = "unknown"
py$adata$obs["ethnicity_2"] = "unknown"
py$adata$obs["ethnicity_details"] = "unknown"
py$adata$obs["self_reported_ethnicity_ontology_term_id"] = "unknown"
py$adata$obs["genotype"] = "unknown"
### Donor status level
py$adata$obs["age_years"] = dplyr::case_match(
  sample_id_c,
  "AA2" ~ "25",
  "AA4" ~ "40",
  "AA7" ~ "35",
  "AA8" ~ "55",
  "C_PB1" ~ "30",
  "C_PB2" ~ "20",
  "C_PB3" ~ "25",
  "C_SD1" ~ "45",
  "C_SD2" ~ "45",
  "C_SD3" ~ "65",
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
py$adata$obs["disease"] = dplyr::case_when(
  grepl("AA", sample_id_c, fixed = T) ~ "alopecia areata",
  grepl("C_PB", sample_id_c, fixed = T) ~ "normal",
  grepl("C_SD", sample_id_c, fixed = T) ~ "normal",
  .default = sample_id_c
)
py$adata$obs["disease_ontology_term_id"] = dplyr::case_when(
  grepl("AA", sample_id_c, fixed = T) ~ "MONDO:0005340",
  grepl("C_PB", sample_id_c, fixed = T) ~ "PATO:0000461",
  grepl("C_SD", sample_id_c, fixed = T) ~ "PATO:0000461",
  .default = sample_id_c
)
py$adata$obs["disease_status"] = dplyr::case_when(
  grepl("C_PB", sample_id_c, fixed = T) ~ "not applicable",
  grepl("C_SD", sample_id_c, fixed = T) ~ "not applicable",
  grepl("AA", sample_id_c, fixed = T) ~ "unknown",
  .default = sample_id_c
)
py$adata$obs["treatment_status"] = dplyr::case_when(
  grepl("C_PB", sample_id_c, fixed = T) ~ "No",
  grepl("C_SD", sample_id_c, fixed = T) ~ "No",
  grepl("AA", sample_id_c, fixed = T) ~ "unknown",
  .default = sample_id_c
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
py$adata$obs["manner_of_death"] = 'not applicable'
### Sample level
py$adata$obs["sample_id"] = sample_id_c
py$adata$obs["sample_id_running"] = as.character(py$adata$obs[["donor_id_running"]])
py$adata$obs["sample_source"] = "surgical donor"
py$adata$obs["sample_collection_method"] = dplyr::case_when(
  grepl("C_PB", sample_id_c, fixed = T) ~ "biopsy",
  grepl("C_SD", sample_id_c, fixed = T) ~ "surgical resection",
  grepl("AA", sample_id_c, fixed = T) ~ "biopsy",
  .default = sample_id_c
)
py$adata$obs["sampled_site_condition"] = dplyr::case_when(
  grepl("C_PB", sample_id_c, fixed = T) ~ "healthy",
  grepl("C_SD", sample_id_c, fixed = T) ~ "healthy",
  grepl("AA", sample_id_c, fixed = T) ~ "diseased",
  .default = sample_id_c
)
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
py$adata$obs["tissue_type"] = "tissue"
py$adata$obs["tissue_free_text"] = "unknown"
py$adata$obs["tissue_ontology_term"] = "skin of scalp"
py$adata$obs["tissue_ontology_term_id"] = "UBERON:8300000"
py$adata$obs["skin_tissue"] = "unknown"
py$adata$obs["sample_collection_site"] = "unknown"
py$adata$obs["sample_cultured"] = "No"
### Sequencing level
py$adata$obs["sequencing_type"] = "GEX"
py$adata$obs["sample_collection_year"] = "unknown"
py$adata$obs["sample_collection_relative_time_point"] = "unknown"
py$adata$obs["library_id"] = "unknown"
py$adata$obs["library_id_repository"] = library_id_repository_c
py$adata$obs["library_preparation_batch"] = "unknown"
py$adata$obs["library_sequencing_run"] = "unknown"
py$adata$obs["sample_preservation_method"] = dplyr::case_match(
  sample_id_c,
  "AA2" ~ "fresh",
  "AA4" ~ "fresh",
  "AA7" ~ "frozen at -80C",
  "AA8" ~ "frozen at -80C",
  "C_PB1" ~ "fresh",
  "C_PB2" ~ "fresh",
  "C_PB3" ~ "fresh",
  "C_SD1" ~ "fresh",
  "C_SD2" ~ "fresh",
  "C_SD3" ~ "fresh",
  .default = sample_id_c
)
py$adata$obs["protocol_url"] = "unknown"
py$adata$obs["dissociation_protocol"] = "unknown"
py$adata$obs["cell_enrichment"] = "unknown"
py$adata$obs["cell_viability_percentage"] = "unknown"
py$adata$obs["cell_number_loaded"] = "unknown"
py$adata$obs["suspension_type"] = "cell"
py$adata$obs["assay_ontology_term"] = "10x 3' v3"
py$adata$obs["assay_ontology_term_id"] = "EFO:0009922"
py$adata$obs["sequenced_fragment"] = "3 prime tag"
py$adata$obs["sequencing_platform"] = "Illumina NextSeq 550"
### Preprocessing level
py$adata$obs["reference_genome"] = "GRCh38"
py$adata$obs["alignment_software"] = "cell ranger 3.1.0"
py$adata$obs["intron_inclusion"] = "unknown"
### Analysis Level
py$adata$obs["is_primary_data"] = FALSE
py$adata$obs["author_cell_type_1"] = "unknown"
py$adata$obs["author_cell_type_2"] = "unknown"
py$adata$obs["author_cell_type_3"] = "unknown"
py$adata$obs["author_cell_type_4"] = "unknown"
py$adata$obs["author_cell_type_5"] = "unknown"
py$adata$obs["author_cell_type_6"] = "unknown"
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

