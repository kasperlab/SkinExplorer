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
py$output_dir = "D:/Study/KI/Projects/###/MetaStudies/Data/Skin/Human/data_preparation_HCA_v1_0/Cheng_Cho_CellReports_2018"
py_run_string("os.makedirs(output_dir, exist_ok=True)")
h5ad_ALL_file = "Cheng_Cho_CellReports_2018_ALL.h5ad"
metadata_ALL_file = "Cheng_Cho_CellReports_2018_ALL_metadata.tsv"
h5ad_HC_file = "Cheng_Cho_CellReports_2018_HC.h5ad"
metadata_HC_file = "Cheng_Cho_CellReports_2018_HC_metadata.tsv"
######
input_folder = "D:/Study/KI/Projects/###/MetaStudies/Data/Internal/Hao/Cheng_2018_raw_data"
sample_unique = c("UCSFSKIN1_ABD4epi", "UCSFSKIN2_BRST41epi", "UCSFSKIN3_BRST53epi", "UCSFSKIN4_FORE12epi", "UCSFSKIN5_FORE8epi", "UCSFSKIN6_FORE9epi", "UCSFSKIN7_PSO48epi", "UCSFSKIN8_PSO49epi", "UCSFSKIN9_PSO14epi", "UCSFSKIN10_SCALP11epi", "UCSFSKIN11_SCALP26epi", "UCSFSKIN12_SCALP32epi")
sample_id_running_c = c()
py_run_string("GEX_dict = {}")
for(this_sample in sample_unique){
  tmp_df = as.data.frame(fread(paste0(input_folder, "/", this_sample, ".txt")))
  py$tmp_df = tmp_df[, -1]
  py$tmp_var_df = data.frame(row.names = as.character(tmp_df[, 1]))
  py_run_string(paste0("GEX_dict['", this_sample, "'] = ad.AnnData(X=scipy.sparse.csr_matrix(tmp_df.values.transpose()).astype(np.float32), obs=pd.DataFrame(index='", this_sample, "_' + tmp_df.columns.astype(str).values), var=tmp_var_df)"))
  sample_id_running_c = c(sample_id_running_c, rep(this_sample, eval(parse(text = paste0("py$GEX_dict$'", this_sample, "'$n_obs")))))
  print(this_sample)
  rm(tmp_df)
  py_run_string("del tmp_df")
  py_run_string("del tmp_var_df")
  py_run_string("gc.collect()")
  gc()
}
py_run_string("adata = ad.concat(GEX_dict, axis=0, join='outer')")
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
Y_ratio_aggregate = aggregate(Y_ratio, by = list(sample_id_running_c), FUN = function(x){
  return(mean(x, na.rm = T))
})
X_ratio_aggregate = aggregate(X_ratio, by = list(sample_id_running_c), FUN = function(x){
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
py$adata$obs["study_id"] = "Cheng_Cho_CellReports_2018"
### Experiment level
py$adata$obs["experiment_id"] = as.character(py$adata$obs[["study_id"]])
py$adata$obs["institute"] = "unknown"
py$adata$obs["author_batch_notes"] = "unknown"
### Donor background level
py$adata$obs["organism"] = "Homo sapiens"
py$adata$obs["organism_ontology_term_id"] = "NCBITaxon:9606"
py$adata$obs["donor_id"] = "unknown"
py$adata$obs["donor_id_running"] = paste0(py$adata$obs[["experiment_id"]], "___", sample_id_running_c)
py$adata$obs["sex_ontology_term"] = dplyr::case_match(
  sample_id_running_c,
  "UCSFSKIN1_ABD4epi" ~ "female",
  "UCSFSKIN2_BRST41epi" ~ "female",
  "UCSFSKIN3_BRST53epi" ~ "female",
  "UCSFSKIN4_FORE12epi" ~ "male",
  "UCSFSKIN5_FORE8epi" ~ "male",
  "UCSFSKIN6_FORE9epi" ~ "male",
  "UCSFSKIN7_PSO48epi" ~ "male",
  "UCSFSKIN8_PSO49epi" ~ "female",
  "UCSFSKIN9_PSO14epi" ~ "female",
  "UCSFSKIN10_SCALP11epi" ~ "male",
  "UCSFSKIN11_SCALP26epi" ~ "female",
  "UCSFSKIN12_SCALP32epi" ~ "female",
  .default = sample_id_running_c
) # S Table 1
py$adata$obs["sex_ontology_term_id"] = dplyr::case_match(
  sample_id_running_c,
  "UCSFSKIN1_ABD4epi" ~ "PATO:0000383",
  "UCSFSKIN2_BRST41epi" ~ "PATO:0000383",
  "UCSFSKIN3_BRST53epi" ~ "PATO:0000383",
  "UCSFSKIN4_FORE12epi" ~ "PATO:0000384",
  "UCSFSKIN5_FORE8epi" ~ "PATO:0000384",
  "UCSFSKIN6_FORE9epi" ~ "PATO:0000384",
  "UCSFSKIN7_PSO48epi" ~ "PATO:0000384",
  "UCSFSKIN8_PSO49epi" ~ "PATO:0000383",
  "UCSFSKIN9_PSO14epi" ~ "PATO:0000383",
  "UCSFSKIN10_SCALP11epi" ~ "PATO:0000384",
  "UCSFSKIN11_SCALP26epi" ~ "PATO:0000383",
  "UCSFSKIN12_SCALP32epi" ~ "PATO:0000383",
  .default = sample_id_running_c
)
py$adata$obs["ethnicity_1"] = "unknown"
py$adata$obs["ethnicity_2"] = "unknown"
py$adata$obs["ethnicity_details"] = "unknown"
py$adata$obs["self_reported_ethnicity_ontology_term_id"] = "unknown"
py$adata$obs["genotype"] = "unknown"
### Donor status level
py$adata$obs["age_years"] = "unknown"
py$adata$obs["age_range"] = "unknown"
py$adata$obs["development_stage_ontology_term"] = "unknown"
py$adata$obs["development_stage_ontology_term_id"] = "unknown"
py$adata$obs["disease"] = dplyr::case_match(
  sample_id_running_c,
  "UCSFSKIN1_ABD4epi" ~ "normal",
  "UCSFSKIN2_BRST41epi" ~ "normal",
  "UCSFSKIN3_BRST53epi" ~ "normal",
  "UCSFSKIN4_FORE12epi" ~ "normal",
  "UCSFSKIN5_FORE8epi" ~ "normal",
  "UCSFSKIN6_FORE9epi" ~ "normal",
  "UCSFSKIN7_PSO48epi" ~ "psoriasis",
  "UCSFSKIN8_PSO49epi" ~ "psoriasis",
  "UCSFSKIN9_PSO14epi" ~ "psoriasis",
  "UCSFSKIN10_SCALP11epi" ~ "normal",
  "UCSFSKIN11_SCALP26epi" ~ "normal",
  "UCSFSKIN12_SCALP32epi" ~ "normal",
  .default = sample_id_running_c
)
py$adata$obs["disease_ontology_term_id"] = dplyr::case_match(
  sample_id_running_c,
  "UCSFSKIN1_ABD4epi" ~ "PATO:0000461",
  "UCSFSKIN2_BRST41epi" ~ "PATO:0000461",
  "UCSFSKIN3_BRST53epi" ~ "PATO:0000461",
  "UCSFSKIN4_FORE12epi" ~ "PATO:0000461",
  "UCSFSKIN5_FORE8epi" ~ "PATO:0000461",
  "UCSFSKIN6_FORE9epi" ~ "PATO:0000461",
  "UCSFSKIN7_PSO48epi" ~ "MONDO:0005083",
  "UCSFSKIN8_PSO49epi" ~ "MONDO:0005083",
  "UCSFSKIN9_PSO14epi" ~ "MONDO:0005083",
  "UCSFSKIN10_SCALP11epi" ~ "PATO:0000461",
  "UCSFSKIN11_SCALP26epi" ~ "PATO:0000461",
  "UCSFSKIN12_SCALP32epi" ~ "PATO:0000461",
  .default = sample_id_running_c
)
py$adata$obs["disease_status"] = dplyr::case_match(
  sample_id_running_c,
  "UCSFSKIN1_ABD4epi" ~ "not applicable",
  "UCSFSKIN2_BRST41epi" ~ "not applicable",
  "UCSFSKIN3_BRST53epi" ~ "not applicable",
  "UCSFSKIN4_FORE12epi" ~ "not applicable",
  "UCSFSKIN5_FORE8epi" ~ "not applicable",
  "UCSFSKIN6_FORE9epi" ~ "not applicable",
  "UCSFSKIN7_PSO48epi" ~ "unknown",
  "UCSFSKIN8_PSO49epi" ~ "unknown",
  "UCSFSKIN9_PSO14epi" ~ "unknown",
  "UCSFSKIN10_SCALP11epi" ~ "not applicable",
  "UCSFSKIN11_SCALP26epi" ~ "not applicable",
  "UCSFSKIN12_SCALP32epi" ~ "not applicable",
  .default = sample_id_running_c
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
py$adata$obs["manner_of_death"] = 'not applicable'
### Sample level
py$adata$obs["sample_id"] = "unknown"
py$adata$obs["sample_id_running"] = as.character(py$adata$obs[["donor_id_running"]])
py$adata$obs["sample_source"] = "surgical donor"
py$adata$obs["sample_collection_method"] = dplyr::case_match(
  sample_id_running_c,
  "UCSFSKIN1_ABD4epi" ~ "surgical resection",
  "UCSFSKIN2_BRST41epi" ~ "surgical resection",
  "UCSFSKIN3_BRST53epi" ~ "surgical resection",
  "UCSFSKIN4_FORE12epi" ~ "surgical resection",
  "UCSFSKIN5_FORE8epi" ~ "surgical resection",
  "UCSFSKIN6_FORE9epi" ~ "surgical resection",
  "UCSFSKIN7_PSO48epi" ~ "unknown",
  "UCSFSKIN8_PSO49epi" ~ "unknown",
  "UCSFSKIN9_PSO14epi" ~ "unknown",
  "UCSFSKIN10_SCALP11epi" ~ "surgical resection",
  "UCSFSKIN11_SCALP26epi" ~ "surgical resection",
  "UCSFSKIN12_SCALP32epi" ~ "surgical resection",
  .default = sample_id_running_c
)
py$adata$obs["sampled_site_condition"] = dplyr::case_match(
  sample_id_running_c,
  "UCSFSKIN1_ABD4epi" ~ "healthy",
  "UCSFSKIN2_BRST41epi" ~ "healthy",
  "UCSFSKIN3_BRST53epi" ~ "healthy",
  "UCSFSKIN4_FORE12epi" ~ "healthy",
  "UCSFSKIN5_FORE8epi" ~ "healthy",
  "UCSFSKIN6_FORE9epi" ~ "healthy",
  "UCSFSKIN7_PSO48epi" ~ "diseased",
  "UCSFSKIN8_PSO49epi" ~ "diseased",
  "UCSFSKIN9_PSO14epi" ~ "diseased",
  "UCSFSKIN10_SCALP11epi" ~ "healthy",
  "UCSFSKIN11_SCALP26epi" ~ "healthy",
  "UCSFSKIN12_SCALP32epi" ~ "healthy",
  .default = sample_id_running_c
)
py$adata$obs["anatomical_region_level_1"] = dplyr::case_match(
  sample_id_running_c,
  "UCSFSKIN1_ABD4epi" ~ "Torso",
  "UCSFSKIN2_BRST41epi" ~ "Torso",
  "UCSFSKIN3_BRST53epi" ~ "Torso",
  "UCSFSKIN4_FORE12epi" ~ "Torso",
  "UCSFSKIN5_FORE8epi" ~ "Torso",
  "UCSFSKIN6_FORE9epi" ~ "Torso",
  "UCSFSKIN7_PSO48epi" ~ "Torso",
  "UCSFSKIN8_PSO49epi" ~ "Torso",
  "UCSFSKIN9_PSO14epi" ~ "Torso",
  "UCSFSKIN10_SCALP11epi" ~ "Head",
  "UCSFSKIN11_SCALP26epi" ~ "Head",
  "UCSFSKIN12_SCALP32epi" ~ "Head",
  .default = sample_id_running_c
)
py$adata$obs["anatomical_region_level_2"] = dplyr::case_match(
  sample_id_running_c,
  "UCSFSKIN1_ABD4epi" ~ "Abdominal",
  "UCSFSKIN2_BRST41epi" ~ "Breast",
  "UCSFSKIN3_BRST53epi" ~ "Breast",
  "UCSFSKIN4_FORE12epi" ~ "Genitalia",
  "UCSFSKIN5_FORE8epi" ~ "Genitalia",
  "UCSFSKIN6_FORE9epi" ~ "Genitalia",
  "UCSFSKIN7_PSO48epi" ~ "unknown",
  "UCSFSKIN8_PSO49epi" ~ "unknown",
  "UCSFSKIN9_PSO14epi" ~ "unknown",
  "UCSFSKIN10_SCALP11epi" ~ "Scalp",
  "UCSFSKIN11_SCALP26epi" ~ "Scalp",
  "UCSFSKIN12_SCALP32epi" ~ "Scalp",
  .default = sample_id_running_c
)
py$adata$obs["anatomical_region_level_3"] = dplyr::case_match(
  sample_id_running_c,
  "UCSFSKIN1_ABD4epi" ~ "unknown",
  "UCSFSKIN2_BRST41epi" ~ "not applicable",
  "UCSFSKIN3_BRST53epi" ~ "not applicable",
  "UCSFSKIN4_FORE12epi" ~ "Foreskin",
  "UCSFSKIN5_FORE8epi" ~ "Foreskin",
  "UCSFSKIN6_FORE9epi" ~ "Foreskin",
  "UCSFSKIN7_PSO48epi" ~ "unknown",
  "UCSFSKIN8_PSO49epi" ~ "unknown",
  "UCSFSKIN9_PSO14epi" ~ "unknown",
  "UCSFSKIN10_SCALP11epi" ~ "unknown",
  "UCSFSKIN11_SCALP26epi" ~ "unknown",
  "UCSFSKIN12_SCALP32epi" ~ "unknown",
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
py$adata$obs["tissue_type"] = "tissue"
py$adata$obs["tissue_free_text"] = "unknown"
py$adata$obs["tissue_ontology_term"] = "unknown"
py$adata$obs["tissue_ontology_term_id"] = "unknown"
py$adata$obs["skin_tissue"] = "Epidermis"
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
py$adata$obs["sequencing_platform"] = dplyr::case_match(
  sample_id_running_c,
  "UCSFSKIN1_ABD4epi" ~ "Illumina HiSeq 2500",
  "UCSFSKIN2_BRST41epi" ~ "Illumina HiSeq 2500, Illumina NovaSeq 6000",
  "UCSFSKIN3_BRST53epi" ~ "Illumina NovaSeq 6000",
  "UCSFSKIN4_FORE12epi" ~ "Illumina HiSeq 2500",
  "UCSFSKIN5_FORE8epi" ~ "Illumina HiSeq 2500",
  "UCSFSKIN6_FORE9epi" ~ "Illumina HiSeq 2500",
  "UCSFSKIN7_PSO48epi" ~ "Illumina HiSeq 2500, Illumina NovaSeq 6000",
  "UCSFSKIN8_PSO49epi" ~ "Illumina HiSeq 2500, Illumina NovaSeq 6000",
  "UCSFSKIN9_PSO14epi" ~ "Illumina HiSeq 2500/4000",
  "UCSFSKIN10_SCALP11epi" ~ "Illumina HiSeq 4000",
  "UCSFSKIN11_SCALP26epi" ~ "Illumina HiSeq 2500",
  "UCSFSKIN12_SCALP32epi" ~ "Illumina HiSeq 2500",
  .default = sample_id_running_c
)
### Preprocessing level
py$adata$obs["reference_genome"] = "GRCh38"
py$adata$obs["alignment_software"] = "cell ranger 2.0.2"
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
this_dataset_individual = "Cheng 2018"
this_dataset_integrated = "Cheng_Cho_CellReports_2018"
individual_barcode_used = rownames(IndividualAnalysis_metadata_list[[this_dataset_individual]])
names(individual_barcode_used) = sapply(rownames(IndividualAnalysis_metadata_list[[this_dataset_individual]]), function(x){
  str_fragments = unlist(strsplit(x, "---", fixed = T))
  return(paste0(str_fragments[2], "_", str_fragments[1]))
})
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
  return(paste0(str_fragments[2], "_", str_fragments[1]))
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
# Write
#
py_run_string("sc.pp.filter_cells(adata, min_genes=1)")
py_run_string(paste0("adata.write(output_dir + '/", h5ad_ALL_file, "')"))
py_run_string(paste0("adata.obs.to_csv(output_dir + '/", metadata_ALL_file, "', sep='\t', index=True)"))
py_run_string("adata_HC = adata[adata.obs['sampled_site_condition'] == 'healthy']")
py_run_string(paste0("adata_HC.write(output_dir + '/", h5ad_HC_file, "')"))
py_run_string(paste0("adata_HC.obs.to_csv(output_dir + '/", metadata_HC_file, "', sep='\t', index=True)"))
py_run_string("print(output_dir)")

