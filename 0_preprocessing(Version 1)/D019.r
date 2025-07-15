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
py$h5_path_1 = "D:/Study/KI/Projects/###/MetaStudies/Data/Skin/Human/Kim_Nagao_NatureMedicine_2020/GSE132802_RAW/GSM3892577_SKIN_HV1_F1_filtered_gene_bc_matrices_h5.h5"
py$h5_path_2 = "D:/Study/KI/Projects/###/MetaStudies/Data/Skin/Human/Kim_Nagao_NatureMedicine_2020/GSE132802_RAW/GSM3892578_SKIN_HV1_F2_filtered_gene_bc_matrices_h5.h5"
py$h5_path_3 = "D:/Study/KI/Projects/###/MetaStudies/Data/Skin/Human/Kim_Nagao_NatureMedicine_2020/GSE132802_RAW/GSM3892579_SKIN_HV2_filtered_gene_bc_matrices_h5.h5"
py$h5_path_4 = "D:/Study/KI/Projects/###/MetaStudies/Data/Skin/Human/Kim_Nagao_NatureMedicine_2020/GSE132802_RAW/GSM3892580_SKIN_HV3_filtered_gene_bc_matrices_h5.h5"
py$h5_path_5 = "D:/Study/KI/Projects/###/MetaStudies/Data/Skin/Human/Kim_Nagao_NatureMedicine_2020/GSE132802_RAW/GSM3892581_SKIN_HV4_filtered_gene_bc_matrices_h5.h5"
py$h5_path_6 = "D:/Study/KI/Projects/###/MetaStudies/Data/Skin/Human/Kim_Nagao_NatureMedicine_2020/GSE132802_RAW/GSM3892582_SKIN_HV5_filtered_gene_bc_matrices_h5.h5"
###
output_dir = "D:/Study/KI/Projects/###/MetaStudies/Data/Skin/Human/Kim_Nagao_NatureMedicine_2020/convert_output"
dir.create(output_dir, showWarnings = F, recursive = T)
loom_path = paste0(output_dir, "/GSE132802_HC.loom")
metadata_path = paste0(output_dir, "/GSE132802_HC_metadata.tsv")
###
# HC
py_run_string("adata_1 = sc.read_10x_h5(h5_path_1)")
py$adata_1 = unify_genename(py$adata_1, alias_gene_c_seurat, verbose = T)
py_run_string("adata_2 = sc.read_10x_h5(h5_path_2)")
py$adata_2 = unify_genename(py$adata_2, alias_gene_c_seurat, verbose = T)
py_run_string("adata_3 = sc.read_10x_h5(h5_path_3)")
py$adata_3 = unify_genename(py$adata_3, alias_gene_c_seurat, verbose = T)
py_run_string("adata_4 = sc.read_10x_h5(h5_path_4)")
py$adata_4 = unify_genename(py$adata_4, alias_gene_c_seurat, verbose = T)
py_run_string("adata_5 = sc.read_10x_h5(h5_path_5)")
py$adata_5 = unify_genename(py$adata_5, alias_gene_c_seurat, verbose = T)
py_run_string("adata_6 = sc.read_10x_h5(h5_path_6)")
py$adata_6 = unify_genename(py$adata_6, alias_gene_c_seurat, verbose = T)
######
py_run_string("adata = ad.concat({
              'SKIN_HV1_duplicate1': adata_1,
              'SKIN_HV1_duplicate2': adata_2,
              'SKIN_HV2': adata_3,
              'SKIN_HV3': adata_4,
              'SKIN_HV4': adata_5,
              'SKIN_HV5': adata_6
              }, axis=0, join='outer', index_unique='---')")
py_run_string("adata.var_names_make_unique()")
py_run_string("sc.pp.filter_cells(adata, min_genes=1)")
py_run_string("sc.pp.filter_genes(adata, min_cells=1)")
this_obj = py$adata$obs
######
# metadata
cell_number = py$adata$n_obs
experiment_id = rep("Kim_Nagao_NatureMedicine_2020", cell_number)
donor_id = sapply(py$adata$obs_names$values, function(x){
  unlist(strsplit(x, "---"))[2]
})
subject_type = rep("Alive - Healthy", cell_number)
sample_id = paste0(experiment_id, "___", donor_id)
sample_type = rep("Biopsy", cell_number)
sampled_site_condition = rep("Healthy", cell_number)
tissue = rep(NA, cell_number)
biological_unit = rep("Cells", cell_number)
library_platform = rep(NA, cell_number)
strand_sequence = rep(NA, cell_number)
sequencing_platform = rep("Illumina HiSeq 3000", cell_number)
reference_genome = rep("GRCh38", cell_number)
sample_status = rep(NA, cell_number)
sample_cultured = rep("No", cell_number)
anatomical_region_level1 = rep("Torso", cell_number)
anatomical_region_level2 = rep("Flank", cell_number)
anatomical_region_level3 = rep(NA, cell_number)
ethnicity_1 = rep(NA, cell_number)
ethnicity_2 = rep(NA, cell_number)
ethnicity_detail = rep(NA, cell_number)
age = rep(NA, cell_number)
sex = rep(NA, cell_number)
original_annotation = rep(NA, cell_number)
######
sex[donor_id %in% c("SKIN_HV1_duplicate1", "SKIN_HV1_duplicate2")] = "M"
library_platform[donor_id %in% c("SKIN_HV1_duplicate1", "SKIN_HV1_duplicate2")] = "10x_3'_v2"
strand_sequence[donor_id %in% c("SKIN_HV1_duplicate1", "SKIN_HV1_duplicate2")] = "3'"
donor_id[donor_id %in% c("SKIN_HV1_duplicate1", "SKIN_HV1_duplicate2")] = "SKIN_HV1"
this_donor = "SKIN_HV2"
sex[donor_id == this_donor] = "F"
library_platform[donor_id == this_donor] = "10x_3'_v2"
strand_sequence[donor_id == this_donor] = "3'"
this_donor = "SKIN_HV3"
sex[donor_id == this_donor] = "F"
library_platform[donor_id == this_donor] = "10x_3'_v2"
strand_sequence[donor_id == this_donor] = "3'"
this_donor = "SKIN_HV4"
sex[donor_id == this_donor] = "M"
library_platform[donor_id == this_donor] = "10x_5'"
strand_sequence[donor_id == this_donor] = "5'"
this_donor = "SKIN_HV5"
sex[donor_id == this_donor] = "M"
library_platform[donor_id == this_donor] = "10x_5'"
strand_sequence[donor_id == this_donor] = "5'"
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



