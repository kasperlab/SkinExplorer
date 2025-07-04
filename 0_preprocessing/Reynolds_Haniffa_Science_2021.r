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
#
# IO - Setting
#
py$output_dir = "D:/Study/KI/Projects/###/MetaStudies/Data/Skin/Human/data_preparation_HCA_v1_0/Reynolds_Haniffa_Science_2021"
py_run_string("os.makedirs(output_dir, exist_ok=True)")
h5ad_10X_file = "Reynolds_Haniffa_Science_2021_10X.h5ad"
metadata_10X_file = "Reynolds_Haniffa_Science_2021_10X_metadata.tsv"
h5ad_10X_POSTPARTUM_file = "Reynolds_Haniffa_Science_2021_10X_POSTPARTUM.h5ad"
metadata_10X_POSTPARTUM_file = "Reynolds_Haniffa_Science_2021_10X_POSTPARTUM_metadata.tsv"
h5ad_10X_POSTPARTUM_NORMAL_file = "Reynolds_Haniffa_Science_2021_10X_POSTPARTUM_NORMAL.h5ad"
metadata_10X_POSTPARTUM_NORMAL_file = "Reynolds_Haniffa_Science_2021_10X_POSTPARTUM_NORMAL_metadata.tsv"
######
py$postpartum_h5ad = "D:/Study/KI/Projects/###/MetaStudies/Data/Skin/Human/Reynolds_Haniffa_Science_2021/zenodo.4569496/submission_210120.h5ad"
py$fetal_h5ad = "D:/Study/KI/Projects/###/MetaStudies/Data/Skin/Human/Reynolds_Haniffa_Science_2021/zenodo.4569496/fetal_submission.h5ad"
#
# postpartum
#
py_run_string("GEX_dict = {}")
py_run_string("GEX_dict['postpartum'] = sc.read_h5ad(postpartum_h5ad)")
postpartum_var = as.matrix(py$GEX_dict$postpartum$var)
colnames_postpartum_var = colnames(postpartum_var)
for(ii in c("gene_ids", "feature_types")){
  var_ids = grep(ii, colnames_postpartum_var, value = T)
  var_ids_1 = var_ids[1]
  if(sum(sapply(var_ids[-1], function(x){
    return(sum(postpartum_var[, var_ids_1] != postpartum_var[, x]))
  })) == 0){
    for(jj in var_ids[-1]){
      py_run_string(paste0("del GEX_dict['postpartum'].var['", jj, "']"))
    }
  }
}
py_run_string("GEX_dict['postpartum'].var['feature_name'] = GEX_dict['postpartum'].var_names.values.copy()")
py_run_string("GEX_dict['postpartum'].var['gene_ids'] = GEX_dict['postpartum'].var['gene_ids-SKN8090524'].copy()")
py_run_string("del GEX_dict['postpartum'].var['gene_ids-SKN8090524']")
py_run_string("GEX_dict['postpartum'].var['feature_types'] = GEX_dict['postpartum'].var['feature_types-SKN8090524'].copy()")
py_run_string("del GEX_dict['postpartum'].var['feature_types-SKN8090524']")
py_run_string("GEX_dict['postpartum'].var_names = GEX_dict['postpartum'].var['gene_ids'].copy()")
py_run_string(paste0("print(GEX_dict['postpartum'].var)"))
#
# fetal
#
py_run_string("GEX_dict['fetal'] = sc.read_h5ad(fetal_h5ad)")
fetal_var = as.matrix(py$GEX_dict$fetal$var)
colnames_fetal_var = colnames(fetal_var)
for(ii in c("gene_ids", "feature_types")){
  var_ids = grep(ii, colnames_fetal_var, value = T)
  var_ids_1 = var_ids[1]
  if(sum(sapply(var_ids[-1], function(x){
    return(sum(fetal_var[, var_ids_1] != fetal_var[, x]))
  })) == 0){
    for(jj in var_ids[-1]){
      py_run_string(paste0("del GEX_dict['fetal'].var['", jj, "']"))
    }
  }
}
py_run_string("GEX_dict['fetal'].var['feature_name'] = GEX_dict['fetal'].var_names.values.copy()")
py_run_string("GEX_dict['fetal'].var['gene_ids'] = GEX_dict['fetal'].var['gene_ids-4834STDY7002879'].copy()")
py_run_string("del GEX_dict['fetal'].var['gene_ids-4834STDY7002879']")
py_run_string("GEX_dict['fetal'].var['feature_types'] = GEX_dict['fetal'].var['feature_types-4834STDY7002879'].copy()")
py_run_string("del GEX_dict['fetal'].var['feature_types-4834STDY7002879']")
py_run_string("GEX_dict['fetal'].var_names = GEX_dict['fetal'].var['gene_ids'].copy()")
py_run_string(paste0("print(GEX_dict['fetal'].var)"))
#
# Merge
#
py_run_string("adata = ad.concat(GEX_dict, axis=0, join='outer')")
py_run_string("SameVar = AssertSameDictAdataVar(GEX_dict)")
if(py$SameVar){
  py_run_string("adata.var = GEX_dict['postpartum'].var.copy()")
}
py_run_string("del GEX_dict")
GeneID = py$adata$var_names$values
py$var_data_frame = data.frame(row.names = GeneID, feature_is_filtered=F, feature_biotype="gene", feature_reference="NCBITaxon:9606", feature_name=mapping_GeneName_GeneID(GeneID, hgnc_mapping_ensemblID_genename))
py_run_string(paste0("var_data_frame.loc[var_data_frame.index.str.startswith('ERCC-'), 'feature_biotype'] = 'spike-in'"))
py_run_string(paste0("adata.var = var_data_frame.copy()"))
py_run_string("del var_data_frame")
py_run_string("gc.collect()")
py_run_string("adata.obs_names = adata.obs_names.to_series().str.split('-').apply(lambda x: x[len(x) - 1] + '_' + '-'.join(x[:(len(x) - 1)])).copy()")
#
# Remake the adata file to make sure it's clean.
#
py_run_string("adata = ad.AnnData(X=adata.X.copy(),
                                  obs=adata.obs.copy(),
                                  var=adata.var.copy())")
#
# Metadata
#
sdrf_E_MTAB_8142_df = as.data.frame(fread("D:/Study/KI/Projects/###/MetaStudies/Data/Skin/Human/Reynolds_Haniffa_Science_2021/E-MTAB-8142/E-MTAB-8142.sdrf.txt", header = T))
sdrf_E_MTAB_8142_df_single_matrix = as.matrix(sdrf_E_MTAB_8142_df[sdrf_E_MTAB_8142_df[, "Comment[read_type]"] == "single", ])
rownames(sdrf_E_MTAB_8142_df_single_matrix) = sdrf_E_MTAB_8142_df_single_matrix[, "Source Name"]
metadata_negriscrep_df = as.data.frame(fread("D:/Study/KI/Projects/###/MetaStudies/Analysis/Integration/KC_benchmarking/metadata_negriscrep.csv", header = T))
metadata_negriscrep_mat = as.matrix(metadata_negriscrep_df[, -1])
rownames(metadata_negriscrep_mat) = as.character(metadata_negriscrep_df[, 1])
umap_coordinates_VN_1_df = as.data.frame(fread("D:/Study/KI/Projects/###/MetaStudies/Analysis/Integration/KC_benchmarking/umap_coordinates_VN_1.csv", header = T))
umap_coordinates_VN_1 = as.matrix(umap_coordinates_VN_1_df[, -1])
rownames(umap_coordinates_VN_1) = as.character(umap_coordinates_VN_1_df[, 1])
#
cell_ids_old = sapply(py$adata$obs_names$values, function(x){
  str_fragments = unlist(strsplit(x, "_"))
  return(paste0(str_fragments[2], "-", str_fragments[1]))
})
sample_id_running_c = sapply(py$adata$obs_names$values, function(x){
  return(unlist(strsplit(x, "_", fixed = T))[1])
})
mask_samples_have_metadata = sample_id_running_c %in% rownames(sdrf_E_MTAB_8142_df_single_matrix)
#
donor_id_running_c = as.character(py$adata$obs[["donor_id"]])
donor_id_running_c[mask_samples_have_metadata] = toupper(sdrf_E_MTAB_8142_df_single_matrix[sample_id_running_c[mask_samples_have_metadata], "Characteristics[individual]"])
sex_c = stringi::stri_replace_all_fixed(str = as.character(py$adata$obs[["Sex"]]),
                                        pattern = c("Female", "Male"),
                                        replacement = c("female", "male"),
                                        vectorize_all = FALSE)
print(sum(is.na(sex_c)))
age_c = dplyr::case_when(
  grepl("F", donor_id_running_c, fixed = T) ~ "up to 8 FW",
  .default = as.character(py$adata$obs[["Age"]])
)
age_c[is.na(age_c)] = "unknown"
disease_c = dplyr::case_when(
  grepl("E", donor_id_running_c, fixed = T) ~ "atopic eczema",
  grepl("F", donor_id_running_c, fixed = T) ~ "normal",
  grepl("P", donor_id_running_c, fixed = T) ~ "psoriasis",
  grepl("S", donor_id_running_c, fixed = T) ~ "normal",
  .default = donor_id_running_c
)
sample_source_c = dplyr::case_when(
  grepl("E", donor_id_running_c, fixed = T) ~ "surgical donor",
  grepl("F", donor_id_running_c, fixed = T) ~ "postmortem donor",
  grepl("P", donor_id_running_c, fixed = T) ~ "surgical donor",
  grepl("S", donor_id_running_c, fixed = T) ~ "surgical donor",
  .default = donor_id_running_c
)
sample_collection_method_c = dplyr::case_when(
  grepl("E", donor_id_running_c, fixed = T) ~ "biopsy",
  grepl("F", donor_id_running_c, fixed = T) ~ "biopsy",
  grepl("P", donor_id_running_c, fixed = T) ~ "biopsy",
  grepl("S", donor_id_running_c, fixed = T) ~ "surgical resection",
  .default = donor_id_running_c
)
site_c = as.character(py$adata$obs[["Site"]])
sampled_site_condition_c = dplyr::case_when(
  grepl("E", donor_id_running_c, fixed = T) & grepl("lesion", site_c, fixed = T) ~ "diseased",
  grepl("E", donor_id_running_c, fixed = T) & grepl("non_lesion", site_c, fixed = T) ~ "healthy",
  grepl("F", donor_id_running_c, fixed = T) ~ "healthy",
  grepl("P", donor_id_running_c, fixed = T) & grepl("lesion", site_c, fixed = T) ~ "diseased",
  grepl("P", donor_id_running_c, fixed = T) & grepl("non_lesion", site_c, fixed = T) ~ "healthy",
  grepl("S", donor_id_running_c, fixed = T) ~ "healthy",
  .default = donor_id_running_c
)
tissue_free_text_c = as.character(py$adata$obs[["Location"]])
tissue_free_text_c[is.na(tissue_free_text_c)] = "unknown"
skin_tissue_c = as.character(py$adata$obs[["Tissue"]])
skin_tissue_c[is.na(skin_tissue_c)] = "unknown"
cell_enrichment_c = as.character(py$adata$obs[["Enrichment"]])
cell_enrichment_c[is.na(cell_enrichment_c)] = "unknown"
author_cell_type_1 = as.character(py$adata$obs[["full_clustering"]])
author_cell_type_1[is.na(author_cell_type_1)] = as.character(py$adata$obs[["anno_final"]])[is.na(author_cell_type_1)]
author_cell_type_1[is.na(author_cell_type_1)] = "unknown"
#
mask_cells_have_metadata = cell_ids_old %in% rownames(metadata_negriscrep_mat)
author_cell_type_2 = rep("unknown", py$adata$n_obs)
author_cell_type_2[mask_cells_have_metadata] = metadata_negriscrep_mat[cell_ids_old[mask_cells_have_metadata], "clusters"]
#
umap_VN_1 = matrix(nrow = py$adata$n_obs, ncol = 2)
umap_VN_1[mask_cells_have_metadata, ] = umap_coordinates_VN_1[cell_ids_old[mask_cells_have_metadata], ]
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
py$adata$obs["study_id"] = "Reynolds_Haniffa_Science_2021"
### Experiment level
py$adata$obs["experiment_id"] = as.character(py$adata$obs[["study_id"]])
py$adata$obs["institute"] = "unknown"
py$adata$obs["author_batch_notes"] = "unknown"
### Donor background level
py$adata$obs["organism"] = "Homo sapiens"
py$adata$obs["organism_ontology_term_id"] = "unknown"
py$adata$obs["donor_id"] = "unknown"
py$adata$obs["donor_id_running"] = paste0(py$adata$obs[["experiment_id"]], "___", donor_id_running_c)
py$adata$obs["sex_ontology_term"] = sex_c
py$adata$obs["sex_ontology_term_id"] = stringi::stri_replace_all_fixed(str = sex_c,
                                                                       pattern = c("female", "male"),
                                                                       replacement = c("PATO:0000383", "PATO:0000384"),
                                                                       vectorize_all = FALSE)
py$adata$obs["ethnicity_1"] = "unknown"
py$adata$obs["ethnicity_2"] = "unknown"
py$adata$obs["ethnicity_details"] = "unknown"
py$adata$obs["self_reported_ethnicity_ontology_term_id"] = "unknown"
py$adata$obs["genotype"] = "unknown"
### Donor status level
py$adata$obs["age_years"] = age_c
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
py$adata$obs["disease"] = disease_c
py$adata$obs["disease_ontology_term_id"] = dplyr::case_match(
  disease_c,
  "atopic eczema" ~ "MONDO:0004980",
  "psoriasis" ~ "MONDO:0005083",
  "normal" ~ "PATO:0000461",
  .default = disease_c
)
py$adata$obs["disease_status"] = dplyr::case_match(
  disease_c,
  "atopic eczema" ~ "unknown",
  "psoriasis" ~ "unknown",
  "normal" ~ "not applicable",
  .default = disease_c
)
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
py$adata$obs["skin_tone"] = "unknown"
py$adata$obs["skin_care"] = "unknown"
py$adata$obs["manner_of_death"] = "unknown"
### Sample level
py$adata$obs["sample_id"] = "unknown"
py$adata$obs["sample_id_running"] = paste0(py$adata$obs[["experiment_id"]], "___", sample_id_running_c)
py$adata$obs["sample_source"] = sample_source_c
py$adata$obs["sample_collection_method"] = sample_collection_method_c
py$adata$obs["sampled_site_condition"] = sampled_site_condition_c
py$adata$obs["anatomical_region_level_1"] = dplyr::case_match(
  tissue_free_text_c,
  "Breast" ~ "Torso",
  "Lower_back" ~ "Torso",
  "Trunk_Limb" ~ "unknown",
  .default = tissue_free_text_c
)
py$adata$obs["anatomical_region_level_2"] = dplyr::case_match(
  tissue_free_text_c,
  "Breast" ~ "Breast",
  "Lower_back" ~ "Back",
  "Trunk_Limb" ~ "unknown",
  .default = tissue_free_text_c
)
py$adata$obs["anatomical_region_level_3"] = dplyr::case_match(
  tissue_free_text_c,
  "Breast" ~ "not applicable",
  "Lower_back" ~ "Lower back",
  "Trunk_Limb" ~ "unknown",
  .default = tissue_free_text_c
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
py$adata$obs["library_id_repository"] = "unknown"
py$adata$obs["library_preparation_batch"] = "unknown"
py$adata$obs["library_sequencing_run"] = "unknown"
py$adata$obs["sample_preservation_method"] = "unknown"
py$adata$obs["protocol_url"] = "unknown"
py$adata$obs["dissociation_protocol"] = "unknown"
py$adata$obs["cell_enrichment"] = cell_enrichment_c
py$adata$obs["cell_viability_percentage"] = "unknown"
py$adata$obs["cell_number_loaded"] = "unknown"
py$adata$obs["suspension_type"] = "cell"
py$adata$obs["assay_ontology_term"] = "10x 3' v2"
py$adata$obs["assay_ontology_term_id"] = "EFO:0009899"
py$adata$obs["sequenced_fragment"] = "3 prime tag"
py$adata$obs["sequencing_platform"] = "Illumina HiSeq 4000"
### Preprocessing level
py$adata$obs["reference_genome"] = "GRCh38"
py$adata$obs["alignment_software"] = "cell ranger 3.1.0"
py$adata$obs["intron_inclusion"] = "unknown"
### Analysis Level
py$adata$obs["is_primary_data"] = FALSE
py$adata$obs["author_cell_type_1"] = author_cell_type_1
py$adata$obs["author_cell_type_2"] = author_cell_type_2
py$adata$obs["author_cell_type_3"] = "unknown"
py$adata$obs["author_cell_type_4"] = "unknown"
py$adata$obs["author_cell_type_5"] = "unknown"
py$adata$obs["author_cell_type_6"] = "unknown"
#
# Embeddings
#
py$adata$obsm["VN_UMAP"] = umap_VN_1
py_run_string("sc.settings.verbosity = 3")
py_run_string("sc.logging.print_header()")
py_run_string("sc.settings.set_figure_params(dpi=300, dpi_save=300, facecolor='white', fontsize=10)")
py_run_string("sc.settings.autoshow = False")
py_run_string("figure_dir = output_dir + '/figures'")
py_run_string("sc.settings.figdir = figure_dir")
dir.create(py$figure_dir, showWarnings = F, recursive = T)
#
adata_str = "adata"
visualization_embedding = "VN_UMAP"
#
py_run_string(paste0("tmp_adata = ad.AnnData(obs=", adata_str, ".obs.copy(), obsm={'", visualization_embedding, "': ", adata_str, ".obsm['", visualization_embedding, "'].copy()})"))
py_run_string(paste0("sc.pl.embedding(tmp_adata, '", visualization_embedding, "', color=['author_cell_type_2'], size=8, frameon=False, ncols=1, save='_author_cell_type_2.pdf')"))
py_run_string("del tmp_adata")
py_run_string("gc.collect()")
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
py_run_string(paste0("adata.write(output_dir + '/", h5ad_10X_file, "')"))
py_run_string(paste0("adata.obs.to_csv(output_dir + '/", metadata_10X_file, "', sep='\t', index=True)"))
py_run_string("adata_10X_POSTPARTUM = adata[adata.obs['age_years'] != 'up to 8 FW']")
py_run_string(paste0("adata_10X_POSTPARTUM.write(output_dir + '/", h5ad_10X_POSTPARTUM_file, "')"))
py_run_string(paste0("adata_10X_POSTPARTUM.obs.to_csv(output_dir + '/", metadata_10X_POSTPARTUM_file, "', sep='\t', index=True)"))
py_run_string("adata_10X_POSTPARTUM_NORMAL = adata_10X_POSTPARTUM[adata_10X_POSTPARTUM.obs['disease'] == 'normal']")
py_run_string(paste0("adata_10X_POSTPARTUM_NORMAL.write(output_dir + '/", h5ad_10X_POSTPARTUM_NORMAL_file, "')"))
py_run_string(paste0("adata_10X_POSTPARTUM_NORMAL.obs.to_csv(output_dir + '/", metadata_10X_POSTPARTUM_NORMAL_file, "', sep='\t', index=True)"))
py_run_string("print(output_dir)")
