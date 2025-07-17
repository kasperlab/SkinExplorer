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
output_dir = "D:/Study/KI/Projects/###/MetaStudies/Data/Skin/Human/Tabib_Lafyatis_NatureCommunications_2021/convert_output"
dir.create(output_dir, showWarnings = F, recursive = T)
loom_path = paste0(output_dir, "/GSE138669_HC.loom")
metadata_path = paste0(output_dir, "/GSE138669_HC_metadata.tsv")
######
folder = "D:/Study/KI/Projects/###/MetaStudies/Data/Skin/Human/Tabib_Lafyatis_NatureCommunications_2021/GSE138669_RAW"
file_list = c("GSM4115868_SC1", "GSM4115870_SC4", "GSM4115872_SC18", "GSM4115874_SC32", "GSM4115875_SC33", "GSM4115876_SC34", "GSM4115878_SC50", "GSM4115880_SC68", "GSM4115885_SC124", "GSM4115886_SC125")
######
# Cell annotaion from Tabib_Lafyatis_JournalofInvestigativeDermatology_2017
metadata_csv_path = "D:/Study/KI/Projects/###/MetaStudies/Data/Skin/Human/Tabib_Lafyatis_JournalofInvestigativeDermatology_2017/Skin_6Control_Metadata/Skin_6Control_Metadata.csv"
metadata_dt = fread(metadata_csv_path)
metadata_df = as.data.frame(metadata_dt)
rownames(metadata_df) = sapply(metadata_df[, 1], function(x){
  str_fragments = unlist(strsplit(x, "control_"))
  return(paste0(gsub(".", "-", str_fragments[2], fixed = T), "---", str_fragments[1]))
})
metadata_df = metadata_df[, -1]
######
# HC
min_expressed_gene = 200
py_run_string("concat_dict = {}")
for(ii in file_list){
  str_fragments = unlist(strsplit(ii, "[_]"))
  this_donor = str_fragments[2]
  this_path = paste0(folder, "/", ii, "raw_feature_bc_matrix.h5")
  py_run_string(paste0("tmp_adata = sc.read_10x_h5('", this_path, "')"))
  py_run_string(paste0("tmp_adata = tmp_adata[(tmp_adata.X > 0).sum(1) >= ", min_expressed_gene, "]"))
  py_run_string("sc.pp.filter_cells(tmp_adata, min_genes=1)")
  py_run_string("sc.pp.filter_genes(tmp_adata, min_cells=1)")
  py$tmp_adata = unify_genename(py$tmp_adata, alias_gene_c_seurat, verbose = T)
  py_run_string(paste0("concat_dict['", this_donor, "'] = tmp_adata.copy()"))
  py_run_string("del tmp_adata")
}
py_run_string("adata = ad.concat(concat_dict, axis=0, join='outer', index_unique='---')")
py_run_string("del concat_dict")
py_run_string("adata.var_names_make_unique()")
py_run_string("sc.pp.filter_cells(adata, min_genes=1)")
py_run_string("sc.pp.filter_genes(adata, min_cells=1)")
this_obj = py$adata$obs
######
# metadata
cell_number = py$adata$n_obs
cellID = py$adata$obs_names$values
experiment_id = rep("Tabib_Lafyatis_NatureCommunications_2021", cell_number)
donor_id = sapply(cellID, function(x){
  unlist(strsplit(x, "---"))[2]
})
subject_type = rep("Alive - Healthy", cell_number)
sample_id = paste0(experiment_id, "___", donor_id)
sample_type = rep("Biopsy", cell_number)
sampled_site_condition = rep("Healthy", cell_number)
tissue = rep(NA, cell_number)
biological_unit = rep("Cells", cell_number)
library_platform = rep(NA, cell_number)
strand_sequence = rep("3'", cell_number)
sequencing_platform = rep("Illumina NextSeq 500", cell_number)
reference_genome = rep("GRCh38", cell_number)
sample_status = rep(NA, cell_number)
sample_cultured = rep("No", cell_number)
anatomical_region_level1 = rep("Extremities", cell_number)
anatomical_region_level2 = rep("Arm", cell_number)
anatomical_region_level3 = rep("Forearm", cell_number)
ethnicity_1 = rep(NA, cell_number)
ethnicity_2 = rep(NA, cell_number)
ethnicity_detail = rep(NA, cell_number)
age = rep(NA, cell_number)
sex = rep(NA, cell_number)
original_annotation = rep(NA, cell_number)
######
this_donor = "SC1"
age[donor_id == this_donor] = 63
sex[donor_id == this_donor] = "M"
ethnicity_1[donor_id == this_donor] = "European Ancestry"
ethnicity_detail[donor_id == this_donor] = "White"
library_platform[donor_id == this_donor] = "10x_3'_v1"
this_donor = "SC4"
age[donor_id == this_donor] = 54
sex[donor_id == this_donor] = "M"
ethnicity_1[donor_id == this_donor] = "European Ancestry"
ethnicity_detail[donor_id == this_donor] = "White"
library_platform[donor_id == this_donor] = "10x_3'_v1"
this_donor = "SC18"
age[donor_id == this_donor] = 66
sex[donor_id == this_donor] = "F"
ethnicity_1[donor_id == this_donor] = "European Ancestry"
ethnicity_detail[donor_id == this_donor] = "White"
library_platform[donor_id == this_donor] = "10x_3'_v1"
this_donor = "SC32"
age[donor_id == this_donor] = 23
sex[donor_id == this_donor] = "F"
ethnicity_1[donor_id == this_donor] = "Asian"
library_platform[donor_id == this_donor] = "10x_3'_v1"
this_donor = "SC33"
age[donor_id == this_donor] = 62
sex[donor_id == this_donor] = "F"
ethnicity_1[donor_id == this_donor] = "European Ancestry"
ethnicity_detail[donor_id == this_donor] = "White"
library_platform[donor_id == this_donor] = "10x_3'_v1"
this_donor = "SC34"
age[donor_id == this_donor] = 24
sex[donor_id == this_donor] = "M"
ethnicity_1[donor_id == this_donor] = "European Ancestry"
ethnicity_detail[donor_id == this_donor] = "White"
library_platform[donor_id == this_donor] = "10x_3'_v1"
this_donor = "SC50"
age[donor_id == this_donor] = 64
sex[donor_id == this_donor] = "M"
ethnicity_1[donor_id == this_donor] = "European Ancestry"
ethnicity_detail[donor_id == this_donor] = "White"
library_platform[donor_id == this_donor] = "10x_3'_v2"
this_donor = "SC68"
age[donor_id == this_donor] = 48
sex[donor_id == this_donor] = "F"
ethnicity_1[donor_id == this_donor] = "European Ancestry"
ethnicity_detail[donor_id == this_donor] = "White"
library_platform[donor_id == this_donor] = "10x_3'_v2"
this_donor = "SC124"
age[donor_id == this_donor] = 54
sex[donor_id == this_donor] = "M"
ethnicity_1[donor_id == this_donor] = "European Ancestry"
ethnicity_detail[donor_id == this_donor] = "White"
library_platform[donor_id == this_donor] = "10x_3'_v2"
this_donor = "SC125"
age[donor_id == this_donor] = 61
sex[donor_id == this_donor] = "M"
ethnicity_1[donor_id == this_donor] = "African Ancestry"
ethnicity_detail[donor_id == this_donor] = "Black"
library_platform[donor_id == this_donor] = "10x_3'_v2"
######
names(original_annotation) = cellID
original_annotation = as.character(metadata_df[cellID, "res.0.6"])
original_annotation[original_annotation == "0"] = "Fibroblast_0"
original_annotation[original_annotation == "1"] = "Keratinocyte_1"
original_annotation[original_annotation == "2"] = "EndothelialCell_2"
original_annotation[original_annotation == "3"] = "Fibroblast_3"
original_annotation[original_annotation == "4"] = "Fibroblast_4"
original_annotation[original_annotation == "5"] = "Keratinocyte_5"
original_annotation[original_annotation == "6"] = "Pericyte_6"
original_annotation[original_annotation == "7"] = "Keratinocyte_7"
original_annotation[original_annotation == "8"] = "Macrophage/DC_8"
original_annotation[original_annotation == "9"] = "Lymphocyte_9"
original_annotation[original_annotation == "10"] = "Pericyte_10"
original_annotation[original_annotation == "11"] = "Keratinocyte_11"
original_annotation[original_annotation == "12"] = "SecretoryEpith_12"
original_annotation[original_annotation == "13"] = "SmoothMuscle_13"
original_annotation[original_annotation == "14"] = "Keratinocyte_14"
original_annotation[original_annotation == "15"] = "Melanocyte_15"
original_annotation[original_annotation == "16"] = "NeuralCell_16"
original_annotation[original_annotation == "17"] = "CornifiedEnv_17"
original_annotation[original_annotation == "18"] = "BCell_18"
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



