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
py$output_dir = "D:/Study/KI/Projects/###/MetaStudies/Data/Skin/Human/data_preparation_HCA_v1_0/He_Guttman_JournalofAllergyandClinicalImmunology_2020"
py_run_string("os.makedirs(output_dir, exist_ok=True)")
h5ad_ALL_file = "GSE147424_ALL.h5ad"
metadata_ALL_file = "GSE147424_ALL_metadata.tsv"
h5ad_HC_file = "GSE147424_HC.h5ad"
metadata_HC_file = "GSE147424_HC_metadata.tsv"
######
input_folder = "D:/Study/KI/Projects/###/MetaStudies/Data/Skin/Human/He_Guttman_JournalofAllergyandClinicalImmunology_2020"
sample_unique = c("S1_LS", "S2_LS", "S3_NL", "S4_H", "S5_LS", "S6_H", "S7_LS", "S8_H", "S9_H", "S10_H", "S11_NL", "S12_H", "S13_H", "S14_NL", "S15_NL", "S16_NL", "S17_H")
names(sample_unique) = c("GSM4430459", "GSM4430460", "GSM4430461", "GSM4430462", "GSM4430463", "GSM4430464", "GSM4430465", "GSM4430466", "GSM4430467", "GSM4430468", "GSM4430469", "GSM4430470", "GSM4430471", "GSM4430472", "GSM4430473", "GSM4430474", "GSM4430475")
sample_txts = list.files(input_folder, pattern = ".clean.data.txt.gz")
sample_id_running_c = c()
py_run_string("GEX_dict = {}")
for(this_txt in sample_txts){
  this_sample = sample_unique[unlist(strsplit(this_txt, "_"))[1]]
  tmp_df = as.data.frame(fread(paste0(input_folder, "/", this_txt)))
  tmp_df_without_1 = tmp_df[, -1]
  colnames(tmp_df_without_1) = sapply(colnames(tmp_df_without_1), function(x){
    return(paste(unlist(strsplit(x, "_"))[-1], collapse = "_"))
  })
  raw_counts_mat = apply(expm1(tmp_df_without_1) / 10000, 2, function(x){
    all_values = sort(unique(x))
    interval = min(all_values[-1] - all_values[-length(all_values)])
    return(x / interval)
  })
  py$tmp_df = as.data.frame(raw_counts_mat)
  py$tmp_var_df = data.frame(row.names = as.character(tmp_df[, 1]))
  py_run_string(paste0("GEX_dict['", this_sample, "'] = ad.AnnData(X=scipy.sparse.csr_matrix(tmp_df.values.transpose()).astype(np.float32), obs=pd.DataFrame(index='", this_sample, "_' + tmp_df.columns.astype(str).values), var=tmp_var_df)"))
  sample_id_running_c = c(sample_id_running_c, rep(this_sample, eval(parse(text = paste0("py$GEX_dict$'", this_sample, "'$n_obs")))))
  print(this_sample)
  rm(tmp_df)
  rm(tmp_df_without_1)
  rm(raw_counts_mat)
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
py$adata$obs["study_id"] = "He_Guttman_JournalofAllergyandClinicalImmunology_2020"
### Experiment level
py$adata$obs["experiment_id"] = as.character(py$adata$obs[["study_id"]])
py$adata$obs["institute"] = "unknown"
py$adata$obs["author_batch_notes"] = "unknown"
### Donor background level
py$adata$obs["organism"] = "Homo sapiens"
py$adata$obs["organism_ontology_term_id"] = "unknown"
py$adata$obs["donor_id"] = "unknown"
py$adata$obs["donor_id_running"] = paste0(py$adata$obs[["experiment_id"]], "___", sample_id_running_c)
py$adata$obs["sex_ontology_term"] = dplyr::case_match(
  sample_id_running_c,
  "S1_LS" ~ "male",
  "S2_LS" ~ "female",
  "S3_NL" ~ "male",
  "S4_H" ~ "male",
  "S5_LS" ~ "female",
  "S6_H" ~ "female",
  "S7_LS" ~ "male",
  "S8_H" ~ "male",
  "S9_H" ~ "male",
  "S10_H" ~ "female",
  "S11_NL" ~ "male",
  "S12_H" ~ "female",
  "S13_H" ~ "male",
  "S14_NL" ~ "female",
  "S15_NL" ~ "female",
  "S16_NL" ~ "male",
  "S17_H" ~ "female",
  .default = sample_id_running_c
) # By checking X/Y expression
py$adata$obs["sex_ontology_term_id"] = dplyr::case_match(
  sample_id_running_c,
  "S1_LS" ~ "PATO:0000384",
  "S2_LS" ~ "PATO:0000383",
  "S3_NL" ~ "PATO:0000384",
  "S4_H" ~ "PATO:0000384",
  "S5_LS" ~ "PATO:0000383",
  "S6_H" ~ "PATO:0000383",
  "S7_LS" ~ "PATO:0000384",
  "S8_H" ~ "PATO:0000384",
  "S9_H" ~ "PATO:0000384",
  "S10_H" ~ "PATO:0000383",
  "S11_NL" ~ "PATO:0000384",
  "S12_H" ~ "PATO:0000383",
  "S13_H" ~ "PATO:0000384",
  "S14_NL" ~ "PATO:0000383",
  "S15_NL" ~ "PATO:0000383",
  "S16_NL" ~ "PATO:0000384",
  "S17_H" ~ "PATO:0000383",
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
  "S1_LS" ~ "atopic eczema",
  "S2_LS" ~ "atopic eczema",
  "S3_NL" ~ "atopic eczema",
  "S4_H" ~ "normal",
  "S5_LS" ~ "atopic eczema",
  "S6_H" ~ "normal",
  "S7_LS" ~ "atopic eczema",
  "S8_H" ~ "normal",
  "S9_H" ~ "normal",
  "S10_H" ~ "normal",
  "S11_NL" ~ "atopic eczema",
  "S12_H" ~ "normal",
  "S13_H" ~ "normal",
  "S14_NL" ~ "atopic eczema",
  "S15_NL" ~ "atopic eczema",
  "S16_NL" ~ "atopic eczema",
  "S17_H" ~ "normal",
  .default = sample_id_running_c
)
py$adata$obs["disease_ontology_term_id"] = dplyr::case_match(
  sample_id_running_c,
  "S1_LS" ~ "MONDO:0004980",
  "S2_LS" ~ "MONDO:0004980",
  "S3_NL" ~ "MONDO:0004980",
  "S4_H" ~ "PATO:0000461",
  "S5_LS" ~ "MONDO:0004980",
  "S6_H" ~ "PATO:0000461",
  "S7_LS" ~ "MONDO:0004980",
  "S8_H" ~ "PATO:0000461",
  "S9_H" ~ "PATO:0000461",
  "S10_H" ~ "PATO:0000461",
  "S11_NL" ~ "MONDO:0004980",
  "S12_H" ~ "PATO:0000461",
  "S13_H" ~ "PATO:0000461",
  "S14_NL" ~ "MONDO:0004980",
  "S15_NL" ~ "MONDO:0004980",
  "S16_NL" ~ "MONDO:0004980",
  "S17_H" ~ "PATO:0000461",
  .default = sample_id_running_c
)
py$adata$obs["disease_status"] = dplyr::case_match(
  sample_id_running_c,
  "S1_LS" ~ "unknown",
  "S2_LS" ~ "unknown",
  "S3_NL" ~ "unknown",
  "S4_H" ~ "not applicable",
  "S5_LS" ~ "unknown",
  "S6_H" ~ "not applicable",
  "S7_LS" ~ "unknown",
  "S8_H" ~ "not applicable",
  "S9_H" ~ "not applicable",
  "S10_H" ~ "not applicable",
  "S11_NL" ~ "unknown",
  "S12_H" ~ "not applicable",
  "S13_H" ~ "not applicable",
  "S14_NL" ~ "unknown",
  "S15_NL" ~ "unknown",
  "S16_NL" ~ "unknown",
  "S17_H" ~ "not applicable",
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
py$adata$obs["sample_collection_method"] = "biopsy"
py$adata$obs["sampled_site_condition"] = dplyr::case_match(
  sample_id_running_c,
  "S1_LS" ~ "diseased",
  "S2_LS" ~ "diseased",
  "S3_NL" ~ "healthy",
  "S4_H" ~ "healthy",
  "S5_LS" ~ "diseased",
  "S6_H" ~ "healthy",
  "S7_LS" ~ "diseased",
  "S8_H" ~ "healthy",
  "S9_H" ~ "healthy",
  "S10_H" ~ "healthy",
  "S11_NL" ~ "healthy",
  "S12_H" ~ "healthy",
  "S13_H" ~ "healthy",
  "S14_NL" ~ "healthy",
  "S15_NL" ~ "healthy",
  "S16_NL" ~ "healthy",
  "S17_H" ~ "healthy",
  .default = sample_id_running_c
)
py$adata$obs["anatomical_region_level_1"] = "Extremities"
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
py$adata$obs["sample_preservation_method"] = dplyr::case_match(
  sample_id_running_c,
  "S1_LS" ~ "frozen",
  "S2_LS" ~ "frozen",
  "S3_NL" ~ "frozen",
  "S4_H" ~ "frozen",
  "S5_LS" ~ "frozen",
  "S6_H" ~ "frozen",
  "S7_LS" ~ "frozen",
  "S8_H" ~ "fresh",
  "S9_H" ~ "frozen",
  "S10_H" ~ "frozen",
  "S11_NL" ~ "frozen",
  "S12_H" ~ "frozen",
  "S13_H" ~ "frozen",
  "S14_NL" ~ "frozen",
  "S15_NL" ~ "frozen",
  "S16_NL" ~ "frozen",
  "S17_H" ~ "frozen",
  .default = sample_id_running_c
)
py$adata$obs["protocol_url"] = "unknown"
py$adata$obs["dissociation_protocol"] = "unknown"
py$adata$obs["cell_enrichment"] = "unknown"
py$adata$obs["cell_viability_percentage"] = "unknown"
py$adata$obs["cell_number_loaded"] = "unknown"
py$adata$obs["suspension_type"] = "cell"
py$adata$obs["assay_ontology_term"] = "10x 3' v2"
py$adata$obs["assay_ontology_term_id"] = "EFO:0009899"
py$adata$obs["sequenced_fragment"] = "3 prime tag"
py$adata$obs["sequencing_platform"] = "Illumina HiSeq 2500"
### Preprocessing level
py$adata$obs["reference_genome"] = "GRCh38"
py$adata$obs["alignment_software"] = "STAR aligner"
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
this_dataset_individual = "He 2020"
this_dataset_integrated = "He_Guttman_JournalofAllergyandClinicalImmunology_2020"
individual_barcode_used = rownames(IndividualAnalysis_metadata_list[[this_dataset_individual]])
names(individual_barcode_used) = individual_barcode_used
individual_barcode_all = individual_barcode_used[sapply(py$adata$obs_names$values, function(x){
  str_fragments = unlist(strsplit(x, "_", fixed = T))
  return(paste0(str_fragments[1], "_", str_fragments[3]))
})]
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
integrated_barcode_all = integrated_barcode_used[sapply(py$adata$obs_names$values, function(x){
  str_fragments = unlist(strsplit(x, "_", fixed = T))
  return(paste0(str_fragments[1], "_", str_fragments[3]))
})]
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

