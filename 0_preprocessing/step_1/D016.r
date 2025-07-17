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
py$h5_path_1 = "D:/Study/KI/Projects/###/MetaStudies/Data/Skin/Human/Boothby_Rosenblum_Nature_2021/GSE183031_RAW/GSM5549903_HC01_filtered_feature_bc_matrix.h5"
py$h5_path_2 = "D:/Study/KI/Projects/###/MetaStudies/Data/Skin/Human/Boothby_Rosenblum_Nature_2021/GSE183031_RAW/GSM5549904_HC02_filtered_feature_bc_matrix.h5"
py$h5_path_3 = "D:/Study/KI/Projects/###/MetaStudies/Data/Skin/Human/Boothby_Rosenblum_Nature_2021/GSE183031_RAW/GSM5549905_HC03_filtered_feature_bc_matrix.h5"
###
output_dir = "D:/Study/KI/Projects/###/MetaStudies/Data/Skin/Human/Boothby_Rosenblum_Nature_2021/convert_output"
dir.create(output_dir, showWarnings = F, recursive = T)
loom_path = paste0(output_dir, "/GSE183031_HC.loom")
metadata_path = paste0(output_dir, "/GSE183031_HC_metadata.tsv")
###
# HC
py_run_string("adata_1 = sc.read_10x_h5(h5_path_1)")
py$adata_1 = unify_genename(py$adata_1, alias_gene_c_seurat, verbose = T)
py_run_string("adata_2 = sc.read_10x_h5(h5_path_2)")
py$adata_2 = unify_genename(py$adata_2, alias_gene_c_seurat, verbose = T)
py_run_string("adata_3 = sc.read_10x_h5(h5_path_3)")
py$adata_3 = unify_genename(py$adata_3, alias_gene_c_seurat, verbose = T)
######
py_run_string("adata = ad.concat({
              'HC01': adata_1,
              'HC02': adata_2,
              'HC03': adata_3
              }, axis=0, join='outer', index_unique='---')")
py_run_string("adata.var_names_make_unique()")
py_run_string("sc.pp.filter_cells(adata, min_genes=1)")
py_run_string("sc.pp.filter_genes(adata, min_cells=1)")
this_obj = py$adata$obs
######
# metadata
cell_number = py$adata$n_obs
experiment_id = rep("Boothby_Rosenblum_Nature_2021", cell_number)
donor_id = c(rep("HC01", py$adata_1$n_obs),
             rep("HC02", py$adata_2$n_obs),
             rep("HC03", py$adata_3$n_obs))
subject_type = rep("Alive - Disease", cell_number)
sample_id = paste0(experiment_id, "___", donor_id)
sample_type = rep("Surgical resection", cell_number)
sampled_site_condition = rep("Disease adjacent", cell_number)
tissue = rep(NA, cell_number)
biological_unit = rep("Cells", cell_number)
library_platform = rep("10x_3'_v3", cell_number)
strand_sequence = rep("3'", cell_number)
sequencing_platform = rep("Illumina NovaSeq 6000", cell_number)
reference_genome = rep("GRCh38", cell_number)
sample_status = rep(NA, cell_number)
sample_cultured = rep("No", cell_number)
anatomical_region_level1 = rep(NA, cell_number)
anatomical_region_level2 = rep(NA, cell_number)
anatomical_region_level3 = rep(NA, cell_number)
ethnicity_1 = rep(NA, cell_number)
ethnicity_2 = rep(NA, cell_number)
ethnicity_detail = rep(NA, cell_number)
age = rep(NA, cell_number)
sex = rep(NA, cell_number)
original_annotation = rep(NA, cell_number)
######
this_donor = "HC01"
age[donor_id == this_donor] = 66
sex[donor_id == this_donor] = "M"
anatomical_region_level1[donor_id == this_donor] = "Head"
anatomical_region_level2[donor_id == this_donor] = "Face"
anatomical_region_level3[donor_id == this_donor] = "Cheek"
this_donor = "HC02"
age[donor_id == this_donor] = 67
sex[donor_id == this_donor] = "M"
anatomical_region_level1[donor_id == this_donor] = "Head"
anatomical_region_level2[donor_id == this_donor] = "Face"
anatomical_region_level3[donor_id == this_donor] = "Forehead"
this_donor = "HC03"
age[donor_id == this_donor] = 64
sex[donor_id == this_donor] = "M"
anatomical_region_level1[donor_id == this_donor] = "Extremities"
anatomical_region_level2[donor_id == this_donor] = "Leg"
anatomical_region_level3[donor_id == this_donor] = NA
######
this_obj[["ExperimentID"]] = experiment_id
this_obj[["DonorID"]] = paste0(experiment_id, "___", donor_id)
this_obj[["Age"]] = age
this_obj[["Sex"]] = sex
this_obj[["Ethnicity1"]] = ethnicity_1
this_obj[["Ethnicity2"]] = ethnicity_2
this_obj[["EthnicityDetail"]] = ethnicity_detail
this_obj[["SubjectType"]] = subject_type
this_obj[["SampleID"]] = sample_id
this_obj[["SampleType"]] = sample_type
this_obj[["SampledSiteCondition"]] = sampled_site_condition
this_obj[["Tissue"]] = tissue
this_obj[["BiologicalUnit"]] = biological_unit
this_obj[["LibraryPlatform"]] = library_platform
this_obj[["StrandSequence"]] = strand_sequence
this_obj[["SequencingPlatform"]] = sequencing_platform
this_obj[["ReferenceGenome"]] = reference_genome
this_obj[["SampleStatus"]] = sample_status
this_obj[["SampleCultured"]] = sample_cultured
this_obj[["AnatomicalRegionLevel1"]] = anatomical_region_level1
this_obj[["AnatomicalRegionLevel2"]] = anatomical_region_level2
this_obj[["AnatomicalRegionLevel3"]] = anatomical_region_level3
this_obj[["OriginalAnnotation"]] = original_annotation
######
py$output_path = loom_path
py_run_string("loompy.create(output_path, adata.X.T, {'Gene': adata.var_names.values}, {'CellID': adata.obs_names.values})")
write.table(as.matrix(this_obj), metadata_path, sep = "\t", row.names = T, col.names = T, quote = F)



