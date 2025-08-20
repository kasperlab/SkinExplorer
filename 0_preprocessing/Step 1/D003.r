r_utilities = "D:/Study/KI/Projects/###/Analysis/code_library/utilities.r"
py_utilities = "D:/Study/KI/Projects/###/Analysis/code_library/utilities.py"
source(r_utilities)
library(reticulate)
py_run_file(py_utilities, convert = F)
######
X_gene_path = "D:/Study/KI/Projects/###/Analysis/useful_files/X_Gene_Homo_sapiens.GRCh38.103.rds"
Y_gene_path = "D:/Study/KI/Projects/###/Analysis/useful_files/Y_Gene_Homo_sapiens.GRCh38.103.rds"
alias_gene_c_seurat_path = "D:/Study/KI/Projects/###/Analysis/useful_files/hgnc_complete_set_c_seurat.rds"
alias_gene_c_seurat = readRDS(alias_gene_c_seurat_path)
######
py$output_dir = "D:/Study/KI/Projects/###/MetaStudies/Analysis/Reynolds_Haniffa_Science_2021/Output/convert_output_v20241021"
dir.create(py$output_dir, showWarnings = F, recursive = T)
###
py$healthy_path = "D:/Study/KI/Projects/###/MetaStudies/Analysis/Reynolds_Haniffa_Science_2021/Input/healthy.h5ad"
py$all_path = "D:/Study/KI/Projects/###/MetaStudies/Analysis/Reynolds_Haniffa_Science_2021/Input/submission_210120.h5ad"
###
#
py_run_string("adata_healthy = sc.read(healthy_path)")
py_run_string("adata_all = sc.read(all_path)")
py_run_string("index_healthy = pd.Series(range(adata_all.X.shape[0]), index=adata_all.obs.index).reindex(adata_healthy.obs.index).values
")
py_run_string("newkeys = np.setdiff1d(adata_all.obs.columns.values, adata_healthy.obs.columns.values)")
py_run_string("metadata_healthy_merged = adata_healthy.obs.merge(adata_all.obs.iloc[index_healthy, :].loc[:, newkeys], how='left', on='index')")
py_run_string("metadata_healthy_merged_sorted = metadata_healthy_merged.reindex(adata_healthy.obs_names)")
py_run_string("adata = ad.AnnData(adata_all.X[index_healthy, :])")
py_run_string("adata.var.index = adata_all.var.index.values")
py_run_string("adata.obs.index = metadata_healthy_merged_sorted.index.values")
py_run_string("del healthy_path")
py_run_string("del all_path")
py_run_string("del adata_all")
py_run_string("del adata_healthy")
py_run_string("del index_healthy")
py_run_string("del newkeys")
py_run_string("del metadata_healthy_merged")
######
# metadata
cell_number = py$adata$n_obs
experiment_id = rep("Reynolds_Haniffa_Science_2021", cell_number)
donor_id = as.character(py$metadata_healthy_merged_sorted[, "donor_id"])
subject_type = rep("Alive - Healthy", cell_number)
sample_id = as.character(py$metadata_healthy_merged_sorted[, "sample_id"])
sample_type = rep("Surgical resection", cell_number)
sampled_site_condition = as.character(py$metadata_healthy_merged_sorted[, "Site"])
tissue = as.character(py$metadata_healthy_merged_sorted[, "Tissue"])
biological_unit = rep("Cells", cell_number)
library_platform = rep("10x_3'_v2", cell_number) # https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-8142/sdrf?full=true
strand_sequence = rep("3'", cell_number) # https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-8142/sdrf?full=true
sequencing_platform = rep("Illumina HiSeq 4000", cell_number) # P-MTAB-87569
reference_genome = rep("GRCh38", cell_number)
sample_status = rep("NA", cell_number)
sample_cultured = rep("No", cell_number)
anatomical_region_level1 = rep("NA", cell_number)
anatomical_region_level2 = rep("NA", cell_number)
anatomical_region_level3 = rep("NA", cell_number)
ethnicity_1 = rep("NA", cell_number)
ethnicity_2 = rep("NA", cell_number)
ethnicity_detail = rep("NA", cell_number)
age = as.character(py$metadata_healthy_merged_sorted[, "Age"])
sex = as.character(py$metadata_healthy_merged_sorted[, "Sex"])
original_annotation = as.character(py$metadata_healthy_merged_sorted[, "full_clustering"])
cell_enrichment = as.character(py$metadata_healthy_merged_sorted[, "Enrichment"])
######
age[age == "NaN"] = "NA"
sex[tolower(sex) == "female"] = "F"
sex[tolower(sex) == "male"] = "M"
bodysite = tolower(as.character(py$metadata_healthy_merged_sorted[, "Location"]))
this_bodysite = "breast"
anatomical_region_level1[bodysite == this_bodysite] = "Torso"
anatomical_region_level2[bodysite == this_bodysite] = "Breast"
anatomical_region_level3[bodysite == this_bodysite] = "NA" # String for scanpy
this_sampled_site_condition = "non_lesion"
sampled_site_condition[sampled_site_condition == this_sampled_site_condition] = "Healthy"
stage = as.character(py$metadata_healthy_merged_sorted[, "stage"])
###
py$adata$obs[["ExperimentID"]] = experiment_id
py$adata$obs[["DonorID"]] = paste0(experiment_id, "___", donor_id)
py$adata$obs[["Age"]] = age
py$adata$obs[["Sex"]] = sex
py$adata$obs[["Ethnicity1"]] = ethnicity_1
py$adata$obs[["Ethnicity2"]] = ethnicity_2
py$adata$obs[["EthnicityDetail"]] = ethnicity_detail
py$adata$obs[["SubjectType"]] = subject_type
py$adata$obs[["SampleID"]] = sample_id
py$adata$obs[["SampleType"]] = sample_type
py$adata$obs[["SampledSiteCondition"]] = sampled_site_condition
py$adata$obs[["Tissue"]] = tissue
py$adata$obs[["BiologicalUnit"]] = biological_unit
py$adata$obs[["LibraryPlatform"]] = library_platform
py$adata$obs[["StrandSequence"]] = strand_sequence
py$adata$obs[["SequencingPlatform"]] = sequencing_platform
py$adata$obs[["ReferenceGenome"]] = reference_genome
py$adata$obs[["SampleStatus"]] = sample_status
py$adata$obs[["SampleCultured"]] = sample_cultured
py$adata$obs[["AnatomicalRegionLevel1"]] = anatomical_region_level1
py$adata$obs[["AnatomicalRegionLevel2"]] = anatomical_region_level2
py$adata$obs[["AnatomicalRegionLevel3"]] = anatomical_region_level3
py$adata$obs[["OriginalAnnotation"]] = original_annotation
py$adata$obs[["CellEnrichment"]] = cell_enrichment
######
write.table(py$metadata_healthy_merged_sorted, paste0(py$output_dir, "/metadata_healthy_merged_sorted.tsv"), sep = "\t", row.names = T, col.names = T, quote = F)
#
py$adata = unify_genename(py$adata, alias_gene_c_seurat)
#
metadata_negriscrep_dt = fread("D:/Study/KI/Projects/###/MetaStudies/Analysis/Integration/KC_benchmarking/metadata_negriscrep.csv", header = T)
metadata_negriscrep_df = as.data.frame(metadata_negriscrep_dt)
#
py$negriscrep_cells = as.character(metadata_negriscrep_df[, "cellnames"])
py_run_string("adata_negriscrep = adata[negriscrep_cells]")
py$adata_negriscrep$obs[["Annotation_negriscrep"]] = as.character(metadata_negriscrep_df[, "clusters"])
#
umap_coordinates_VN_1_dt = fread("D:/Study/KI/Projects/###/MetaStudies/Analysis/Integration/KC_benchmarking/umap_coordinates_VN_1.csv", header = T)
umap_coordinates_VN_1_df = as.data.frame(metadata_negriscrep_dt)
#
py$negriscrep_cells = as.character(metadata_negriscrep_df[, "cellnames"])
py_run_string("adata_negriscrep = adata[negriscrep_cells]")
py$adata_negriscrep$obs[["Annotation_negriscrep"]] = as.character(metadata_negriscrep_df[, "clusters"])
#
# Write
#
py_run_string("adata.write(output_dir + '/healthy.h5ad')")
py_run_string("adata.obs.to_csv(output_dir + '/healthy_metadata.txt', sep='\t', index=True)")
py_run_string("adata_negriscrep.write(output_dir + '/negriscrep.h5ad')")
py_run_string("adata_negriscrep.obs.to_csv(output_dir + '/negriscrep_metadata.txt', sep='\t', index=True)")
py_run_string("print(output_dir)")


