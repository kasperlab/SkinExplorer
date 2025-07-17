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
py$output_dir = "D:/Study/KI/Projects/###/MetaStudies/Data/Skin/Human/data_preparation_HCA_v1_0/Gaydosik_Fuschiotti_ClinicalCancerResearch_2019"
py_run_string("os.makedirs(output_dir, exist_ok=True)")
h5ad_ALL_file = "GSE128531_ALL.h5ad"
metadata_ALL_file = "GSE128531_ALL_metadata.tsv"
h5ad_HC_file = "GSE128531_HC.h5ad"
metadata_HC_file = "GSE128531_HC_metadata.tsv"
######
input_folder = "D:/Study/KI/Projects/###/MetaStudies/Data/Skin/Human/Gaydosik_Fuschiotti_ClinicalCancerResearch_2019"
sample_unique = c("CTCL-2", "CTCL-5", "CTCL-6", "CTCL-8", "CTCL-12", "HC-1", "HC-2", "HC-3", "HC-4")
names(sample_unique) = c("GSM3679033", "GSM3679034", "GSM3679035", "GSM3679036", "GSM3679037", "GSM3679038", "GSM3679039", "GSM3679040", "GSM3679041")
sample_csvs = list.files(input_folder, pattern = "csv.gz")
sample_id_running_c = c()
py_run_string("GEX_dict = {}")
for(this_csv in sample_csvs){
  this_sample = sample_unique[unlist(strsplit(this_csv, "_"))[1]]
  tmp_df = as.data.frame(fread(paste0(input_folder, "/", this_csv)))
  tmp_df_without_1 = tmp_df[, -1]
  colnames(tmp_df_without_1) = sapply(colnames(tmp_df_without_1), function(x){
    return(paste(unlist(strsplit(x, "_"))[-1], collapse = "_"))
  })
  py$tmp_df = tmp_df_without_1
  py$tmp_var_df = data.frame(row.names = as.character(tmp_df[, 1]))
  py_run_string(paste0("GEX_dict['", this_sample, "'] = ad.AnnData(X=scipy.sparse.csr_matrix(tmp_df.values.transpose()).astype(np.float32), obs=pd.DataFrame(index='", this_sample, "_' + tmp_df.columns.astype(str).values), var=tmp_var_df)"))
  sample_id_running_c = c(sample_id_running_c, rep(this_sample, eval(parse(text = paste0("py$GEX_dict$'", this_sample, "'$n_obs")))))
  print(this_sample)
  rm(tmp_df)
  rm(tmp_df_without_1)
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
HCA_metadata_list = load_HCA_metadata_xlsx("D:/Study/KI/Projects/###/MetaStudies/Data/Skin/Human/HCA_metadata/Gaydosik_2019_HCA_tier 1_metadata.xlsx")
HCA_metadata_list_original = HCA_metadata_list # Keep the original for checking later
rownames(HCA_metadata_list$DN) = c("CTCL-2", "CTCL-5", "CTCL-6", "CTCL-8", "CTCL-12", "HC-1", "HC-2", "HC-3", "HC-4")
HCA_metadata_list$DN[is.na(HCA_metadata_list$DN)] = "unknown"
rownames(HCA_metadata_list$SP) = c("CTCL-2", "CTCL-5", "CTCL-6", "CTCL-8", "CTCL-12", "HC-1", "HC-2", "HC-3", "HC-4")
HCA_metadata_list$SP[is.na(HCA_metadata_list$SP)] = "unknown"
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
py$adata$obs["study_id"] = "Gaydosik_Fuschiotti_ClinicalCancerResearch_2019"
### Experiment level
py$adata$obs["experiment_id"] = as.character(py$adata$obs[["study_id"]])
py$adata$obs["institute"] = HCA_metadata_list$SP[sample_id_running_c, "institute"]
py$adata$obs["author_batch_notes"] = HCA_metadata_list$SP[sample_id_running_c, "author_batch_notes"]
### Donor background level
py$adata$obs["organism"] = "Homo sapiens"
py$adata$obs["organism_ontology_term_id"] = HCA_metadata_list$DN[sample_id_running_c, "organism_ontology_term_id"]
py$adata$obs["donor_id"] = HCA_metadata_list$DN[sample_id_running_c, "donor_id"]
py$adata$obs["donor_id_running"] = paste0(py$adata$obs[["experiment_id"]], "___", sample_id_running_c)
py$adata$obs["sex_ontology_term"] = HCA_metadata_list$DN[sample_id_running_c, "sex_ontology_term"]
py$adata$obs["sex_ontology_term_id"] = HCA_metadata_list$DN[sample_id_running_c, "sex_ontology_term_id"]
py$adata$obs["ethnicity_1"] = "unknown"
py$adata$obs["ethnicity_2"] = "unknown"
py$adata$obs["ethnicity_details"] = "unknown"
py$adata$obs["self_reported_ethnicity_ontology_term_id"] = "unknown"
py$adata$obs["genotype"] = "unknown"
### Donor status level
py$adata$obs["age_years"] = dplyr::case_match(
  sample_id_running_c,
  "CTCL-2" ~ "68",
  "CTCL-5" ~ "76",
  "CTCL-6" ~ "83",
  "CTCL-8" ~ "46",
  "CTCL-12" ~ "77",
  "HC-1" ~ "64",
  "HC-2" ~ "48",
  "HC-3" ~ "54",
  "HC-4" ~ "61",
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
py$adata$obs["development_stage_ontology_term_id"] = HCA_metadata_list$SP[sample_id_running_c, "development_stage_ontology_term_id"]
py$adata$obs["disease"] = dplyr::case_when(
  grepl("HC", sample_id_running_c, fixed = T) ~ "normal",
  grepl("CTCL", sample_id_running_c, fixed = T) ~ "primary cutaneous T-cell lymphoma",
  .default = sample_id_running_c
)
py$adata$obs["disease_ontology_term_id"] = HCA_metadata_list$SP[sample_id_running_c, "disease_ontology_term_id"]
py$adata$obs["disease_status"] = dplyr::case_match(
  sample_id_running_c,
  "CTCL-2" ~ "IVA(T4NXM0B2)",
  "CTCL-5" ~ "IIB(T3N0M0B0)",
  "CTCL-6" ~ "IIB(T3NXM0B0)",
  "CTCL-8" ~ "IIB(T3N0M0B0)",
  "CTCL-12" ~ "IVA(T4NXB0M0)",
  "HC-1" ~ "not applicable",
  "HC-2" ~ "not applicable",
  "HC-3" ~ "not applicable",
  "HC-4" ~ "not applicable",
  .default = sample_id_running_c
)
py$adata$obs["treatment_status"] = dplyr::case_match(
  sample_id_running_c,
  "CTCL-2" ~ "Previous:NB-UVB;Ongoing:Bexarotene,GCS(Topical),Romidepsin",
  "CTCL-5" ~ "Previous:Bexarotene,GCS(Topical),LEB,Methotrexate,Nitrogen Mustard,Pralatrexate,PUVA,Radiotherapy;Ongoing:No",
  "CTCL-6" ~ "Previous:Bexarotene,LEB,Nitrogen Mustard,Pralatrexate;Ongoing:No",
  "CTCL-8" ~ "Previous:Methotrexate,NB-UVB,PUVA,Vorinostat;Ongoing:Bexarotene,GCS(Topical),IFN alpha,Nitrogen Mustard",
  "CTCL-12" ~ "Previous:Bexarotene,Triamcinolone,Valchlor;Ongoing:Romidepsin",
  "HC-1" ~ "unknown",
  "HC-2" ~ "unknown",
  "HC-3" ~ "unknown",
  "HC-4" ~ "unknown",
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
py$adata$obs["manner_of_death"] = dplyr::case_match(
  sample_id_running_c,
  "CTCL-2" ~ "not applicable",
  "CTCL-5" ~ "4",
  "CTCL-6" ~ "4",
  "CTCL-8" ~ "not applicable",
  "CTCL-12" ~ "not applicable",
  "HC-1" ~ "not applicable",
  "HC-2" ~ "not applicable",
  "HC-3" ~ "not applicable",
  "HC-4" ~ "not applicable",
  .default = sample_id_running_c
)
### Sample level
py$adata$obs["sample_id"] = HCA_metadata_list$SP[sample_id_running_c, "sample_id"]
py$adata$obs["sample_id_running"] = as.character(py$adata$obs[["donor_id_running"]])
py$adata$obs["sample_source"] = HCA_metadata_list$SP[sample_id_running_c, "sample_source"]
py$adata$obs["sample_collection_method"] = HCA_metadata_list$SP[sample_id_running_c, "sample_collection_method"]
py$adata$obs["sampled_site_condition"] = HCA_metadata_list$SP[sample_id_running_c, "sampled_site_condition"]
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
py$adata$obs["tissue_type"] = HCA_metadata_list$SP[sample_id_running_c, "tissue_type"]
py$adata$obs["tissue_free_text"] = HCA_metadata_list$SP[sample_id_running_c, "tissue_free_text"]
py$adata$obs["tissue_ontology_term"] = HCA_metadata_list$SP[sample_id_running_c, "tissue_ontology_term"]
py$adata$obs["tissue_ontology_term_id"] = HCA_metadata_list$SP[sample_id_running_c, "tissue_ontology_term_id"]
py$adata$obs["skin_tissue"] = "unknown"
py$adata$obs["sample_collection_site"] = HCA_metadata_list$SP[sample_id_running_c, "sample_collection_site"]
py$adata$obs["sample_cultured"] = "No"
### Sequencing level
py$adata$obs["sequencing_type"] = "GEX"
py$adata$obs["sample_collection_year"] = dplyr::case_match(
  sample_id_running_c,
  "CTCL-2" ~ "2017",
  "CTCL-5" ~ "2017",
  "CTCL-6" ~ "2017",
  "CTCL-8" ~ "2017",
  "CTCL-12" ~ "2018",
  "HC-1" ~ "unknown",
  "HC-2" ~ "unknown",
  "HC-3" ~ "unknown",
  "HC-4" ~ "unknown",
  .default = sample_id_running_c
)
py$adata$obs["sample_collection_relative_time_point"] = HCA_metadata_list$SP[sample_id_running_c, "sample_collection_relative_time_point"]
py$adata$obs["library_id"] = HCA_metadata_list$SP[sample_id_running_c, "library_id"]
py$adata$obs["library_id_repository"] = HCA_metadata_list$SP[sample_id_running_c, "library_id_repository"]
py$adata$obs["library_preparation_batch"] = HCA_metadata_list$SP[sample_id_running_c, "library_preparation_batch"]
py$adata$obs["library_sequencing_run"] = HCA_metadata_list$SP[sample_id_running_c, "library_sequencing_run"]
py$adata$obs["sample_preservation_method"] = HCA_metadata_list$SP[sample_id_running_c, "sample_preservation_method"]
py$adata$obs["protocol_url"] = "unknown"
py$adata$obs["dissociation_protocol"] = HCA_metadata_list$SP[sample_id_running_c, "dissociation_protocol"]
py$adata$obs["cell_enrichment"] = HCA_metadata_list$SP[sample_id_running_c, "cell_enrichment"]
py$adata$obs["cell_viability_percentage"] = HCA_metadata_list$SP[sample_id_running_c, "cell_viability_percentage"]
py$adata$obs["cell_number_loaded"] = HCA_metadata_list$SP[sample_id_running_c, "cell_number_loaded"]
py$adata$obs["suspension_type"] = HCA_metadata_list$SP[sample_id_running_c, "suspension_type"]
py$adata$obs["assay_ontology_term"] = "10x 3' v2"
py$adata$obs["assay_ontology_term_id"] = "EFO:0009899"
py$adata$obs["sequenced_fragment"] = "3 prime tag"
py$adata$obs["sequencing_platform"] = "Illumina NextSeq 500"
### Preprocessing level
py$adata$obs["reference_genome"] = "GRCh38"
py$adata$obs["alignment_software"] = "cell ranger"
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
this_dataset_individual = "Gaydosik 2019"
this_dataset_integrated = "Gaydosik_Fuschiotti_ClinicalCancerResearch_2019"
individual_barcode_used = rownames(IndividualAnalysis_metadata_list[[this_dataset_individual]])
names(individual_barcode_used) = individual_barcode_used
individual_barcode_all = individual_barcode_used[stringi::stri_replace_all_fixed(str = py$adata$obs_names$values,
                                                                                 pattern = c("HC-1", "HC-2", "HC-3", "HC-4"),
                                                                                 replacement = c("SC50nor", "SC68nor", "SC124nor", "SC125nor"),
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
names(integrated_barcode_used) = sapply(integrated_barcode_used, function(x){
  str_fragments = unlist(strsplit(x, "---", fixed = T))
  return(paste(str_fragments[-length(str_fragments)], collapse = "---"))
})
integrated_barcode_all = integrated_barcode_used[stringi::stri_replace_all_fixed(str = py$adata$obs_names$values,
                                                                                 pattern = c("HC-1", "HC-2", "HC-3", "HC-4"),
                                                                                 replacement = c("SC50nor", "SC68nor", "SC124nor", "SC125nor"),
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

