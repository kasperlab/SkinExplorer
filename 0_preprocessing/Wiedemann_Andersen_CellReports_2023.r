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
py$output_dir = "D:/Study/KI/Projects/###/MetaStudies/Data/Skin/Human/data_preparation_HCA_v1_0/Wiedemann_Andersen_CellReports_2023"
py_run_string("os.makedirs(output_dir, exist_ok=True)")
h5ad_ALL_file = "GSE202352_ALL.h5ad"
metadata_ALL_file = "GSE202352_ALL_metadata.tsv"
######
input_folder = "D:/Study/KI/Projects/###/MetaStudies/Data/Skin/Human/Wiedemann_Andersen_CellReports_2023/GSE202352_RAW"
sample_folders = list.dirs(input_folder, recursive = F, full.names = F)
lib_id_repo_c = c()
donor_id_running_c = c()
bodysite_c = c()
skin_tissue_c = c()
sample_id_running_c = c()
py_run_string("GEX_dict = {}")
for(this_folder in sample_folders){
  str_fragments = unlist(strsplit(gsub(".filtered_feature_bc_matrix", "", this_folder, fixed = T), "_", fixed = T))
  this_lib_repo = str_fragments[1]
  this_donor_id_running = str_fragments[4]
  this_bodysite = str_fragments[5]
  if(length(str_fragments) == 5){
    this_skin_tissue = "Epidermis + Dermis"
  }else{
    if(length(str_fragments) != 6){
      stop("Error: str_fragments")
    }else{
      if(str_fragments[6] == "dermis"){
        this_skin_tissue = "Dermis"
      }else if(str_fragments[6] == "epi"){
        this_skin_tissue = "Epidermis"
      }else{
        stop("Error: str_fragments[6]")
      }
    }
  }
  this_sample_id_running = stringi::stri_replace_all_fixed(str = paste0(this_donor_id_running, "_", this_bodysite, "_", this_skin_tissue),
                                                           pattern = c("Epidermis", " + ", "Dermis"),
                                                           replacement = c("E", "", "D"),
                                                           vectorize_all = FALSE)
  ReadMtx_py("tmp_adata", data_dir = paste0(input_folder, "/", this_folder), prefix = "", feature_column_name = c("gene_ids", "feature_name", "feature_types"))
  py_run_string(paste0("tmp_adata.obs_names = '", this_lib_repo, "_' + tmp_adata.obs_names.astype(str).values"))
  py_run_string(paste0("GEX_dict['", this_lib_repo, "'] = tmp_adata"))
  lib_id_repo_c = c(lib_id_repo_c, rep(this_lib_repo, py$tmp_adata$n_obs))
  donor_id_running_c = c(donor_id_running_c, rep(this_donor_id_running, py$tmp_adata$n_obs))
  bodysite_c = c(bodysite_c, rep(this_bodysite, py$tmp_adata$n_obs))
  skin_tissue_c = c(skin_tissue_c, rep(this_skin_tissue, py$tmp_adata$n_obs))
  sample_id_running_c = c(sample_id_running_c, rep(this_sample_id_running, py$tmp_adata$n_obs))
  print(this_sample_id_running)
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
py$adata$obs["study_id"] = "Wiedemann_Andersen_CellReports_2023"
### Experiment level
py$adata$obs["experiment_id"] = as.character(py$adata$obs[["study_id"]])
py$adata$obs["institute"] = "unknown"
py$adata$obs["author_batch_notes"] = "unknown"
### Donor background level
py$adata$obs["organism"] = "Homo sapiens"
py$adata$obs["organism_ontology_term_id"] = "NCBITaxon:9606"
py$adata$obs["donor_id"] = "unknown"
py$adata$obs["donor_id_running"] = paste0(py$adata$obs[["experiment_id"]], "___", donor_id_running_c)
py$adata$obs["sex_ontology_term"] = dplyr::case_match(
  donor_id_running_c,
  "Subject1" ~ "male",
  "Subject2" ~ "male",
  "Subject3" ~ "male",
  "Subject4" ~ "female",
  .default = sample_id_running_c
) # by inferring, 3 male, 1 female, between the ages of 39 and 65
py$adata$obs["sex_ontology_term_id"] = dplyr::case_match(
  donor_id_running_c,
  "Subject1" ~ "PATO:0000384",
  "Subject2" ~ "PATO:0000384",
  "Subject3" ~ "PATO:0000384",
  "Subject4" ~ "PATO:0000383",
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
py$adata$obs["disease"] = "normal"
py$adata$obs["disease_ontology_term_id"] = "PATO:0000461"
py$adata$obs["disease_status"] = "not applicable"
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
py$adata$obs["sample_id_running"] = as.character(paste0(py$adata$obs[["experiment_id"]], "___", sample_id_running_c))
py$adata$obs["sample_source"] = "surgical donor"
py$adata$obs["sample_collection_method"] = "biopsy"
py$adata$obs["sampled_site_condition"] = "healthy"
py$adata$obs["anatomical_region_level_1"] = dplyr::case_match(
  bodysite_c,
  "hip" ~ "Torso",
  "palm" ~ "Extremities",
  "sole" ~ "Extremities",
  .default = bodysite_c
)
py$adata$obs["anatomical_region_level_2"] = dplyr::case_match(
  bodysite_c,
  "hip" ~ "Hip",
  "palm" ~ "Hand",
  "sole" ~ "Foot",
  .default = bodysite_c
)
py$adata$obs["anatomical_region_level_3"] = dplyr::case_match(
  bodysite_c,
  "hip" ~ "not applicable",
  "palm" ~ "Palm",
  "sole" ~ "Plantar",
  .default = bodysite_c
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
py$adata$obs["skin_tissue"] = skin_tissue_c
py$adata$obs["sample_collection_site"] = "unknown"
py$adata$obs["sample_cultured"] = "No"
### Sequencing level
py$adata$obs["sequencing_type"] = "GEX"
py$adata$obs["sample_collection_year"] = "unknown"
py$adata$obs["sample_collection_relative_time_point"] = "unknown"
py$adata$obs["library_id"] = "unknown"
py$adata$obs["library_id_repository"] = lib_id_repo_c
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
py$adata$obs["reference_genome"] = "GRCh37"
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
this_dataset_individual = "Wiedemann 2023"
this_dataset_integrated = "Wiedemann_Andersen_CellReports_2023"
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
py_run_string("print(output_dir)")

