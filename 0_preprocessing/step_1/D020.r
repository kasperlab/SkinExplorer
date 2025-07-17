r_utilities = "D:/Study/KI/Projects/###/Analysis/code_library/utilities.r"
source(r_utilities)
######
X_gene_path = "D:/Study/KI/Projects/###/Analysis/useful_files/X_Gene_Homo_sapiens.GRCh38.103.rds"
Y_gene_path = "D:/Study/KI/Projects/###/Analysis/useful_files/Y_Gene_Homo_sapiens.GRCh38.103.rds"
alias_gene_c_seurat_path = "D:/Study/KI/Projects/###/Analysis/useful_files/hgnc_complete_set_c_seurat.rds"
alias_gene_c_seurat = readRDS(alias_gene_c_seurat_path)
######
# P112
mtx_path_1 = "D:/Study/KI/Projects/###/MetaStudies/Data/Skin/Human/Rindler_Brunner_MolecularCancer_2021/GSE173205_RAW/GSM5534590_P112_HC_matrix.mtx.gz"
barcodes_path_1 = "D:/Study/KI/Projects/###/MetaStudies/Data/Skin/Human/Rindler_Brunner_MolecularCancer_2021/GSE173205_RAW/GSM5534590_P112_HC_barcodes.tsv.gz"
genes_path_1 = "D:/Study/KI/Projects/###/MetaStudies/Data/Skin/Human/Rindler_Brunner_MolecularCancer_2021/GSE173205_RAW/GSM5534590_P112_HC_features.tsv.gz"
# P115
mtx_path_2 = "D:/Study/KI/Projects/###/MetaStudies/Data/Skin/Human/Rindler_Brunner_MolecularCancer_2021/GSE173205_RAW/GSM5534591_P115_HC_matrix.mtx.gz"
barcodes_path_2 = "D:/Study/KI/Projects/###/MetaStudies/Data/Skin/Human/Rindler_Brunner_MolecularCancer_2021/GSE173205_RAW/GSM5534591_P115_HC_barcodes.tsv.gz"
genes_path_2 = "D:/Study/KI/Projects/###/MetaStudies/Data/Skin/Human/Rindler_Brunner_MolecularCancer_2021/GSE173205_RAW/GSM5534591_P115_HC_features.tsv.gz"
# P116
mtx_path_3 = "D:/Study/KI/Projects/###/MetaStudies/Data/Skin/Human/Rindler_Brunner_MolecularCancer_2021/GSE173205_RAW/GSM5534592_P116_HC_matrix.mtx.gz"
barcodes_path_3 = "D:/Study/KI/Projects/###/MetaStudies/Data/Skin/Human/Rindler_Brunner_MolecularCancer_2021/GSE173205_RAW/GSM5534592_P116_HC_barcodes.tsv.gz"
genes_path_3 = "D:/Study/KI/Projects/###/MetaStudies/Data/Skin/Human/Rindler_Brunner_MolecularCancer_2021/GSE173205_RAW/GSM5534592_P116_HC_features.tsv.gz"
# P121
mtx_path_4 = "D:/Study/KI/Projects/###/MetaStudies/Data/Skin/Human/Rindler_Brunner_MolecularCancer_2021/GSE173205_RAW/GSM5534593_P121_HC_matrix.mtx.gz"
barcodes_path_4 = "D:/Study/KI/Projects/###/MetaStudies/Data/Skin/Human/Rindler_Brunner_MolecularCancer_2021/GSE173205_RAW/GSM5534593_P121_HC_barcodes.tsv.gz"
genes_path_4 = "D:/Study/KI/Projects/###/MetaStudies/Data/Skin/Human/Rindler_Brunner_MolecularCancer_2021/GSE173205_RAW/GSM5534593_P121_HC_features.tsv.gz"
### Output
output_dir = "D:/Study/KI/Projects/###/MetaStudies/Data/Skin/Human/Rindler_Brunner_MolecularCancer_2021/convert_output"
dir.create(output_dir, showWarnings = F, recursive = T)
loom_path = paste0(output_dir, "/GSE173205_HC.loom")
metadata_path = paste0(output_dir, "/GSE173205_HC_metadata.tsv")
###
expression_matrix_1 = ReadMtx(mtx = mtx_path_1, cells = barcodes_path_1, features = genes_path_1)
expression_matrix_2 = ReadMtx(mtx = mtx_path_2, cells = barcodes_path_2, features = genes_path_2)
expression_matrix_3 = ReadMtx(mtx = mtx_path_3, cells = barcodes_path_3, features = genes_path_3)
expression_matrix_4 = ReadMtx(mtx = mtx_path_4, cells = barcodes_path_4, features = genes_path_4)
gene_bc_list = list(expression_matrix_1, expression_matrix_2, expression_matrix_3, expression_matrix_4)
cell_id_list = list(colnames(expression_matrix_1), colnames(expression_matrix_2), colnames(expression_matrix_3), colnames(expression_matrix_4))
this_matrix = merge_matrices(gene_bc_list, is_sparse = T)
rm(gene_bc_list)
this_obj = CreateSeuratObject(this_matrix, min.cells = 1)
this_obj = unify_genename_seurat(this_obj, alias_gene_c_seurat, verbose = T)
rm(this_matrix)
######
# metadata
#View(this_obj_0@meta.data)
cell_number = ncol(this_obj)
experiment_id = rep("Rindler_Brunner_MolecularCancer_2021", cell_number)
donor_id = c(rep("P112", length(cell_id_list[[1]])), 
             rep("P115", length(cell_id_list[[2]])), 
             rep("P116", length(cell_id_list[[3]])), 
             rep("P121", length(cell_id_list[[4]])))
subject_type = rep("Alive - Healthy", cell_number)
sample_id = paste0(experiment_id, "___", donor_id)
sample_type = rep("Biopsy", cell_number)
sampled_site_condition = rep("Healthy", cell_number)
tissue = rep(NA, cell_number)
biological_unit = rep("Cells", cell_number)
library_platform = rep("10x_5'_v1.1", cell_number)
strand_sequence = rep("5'", cell_number)
sequencing_platform = rep("Illumina NovaSeq 6000 SP", cell_number)
reference_genome = rep("GRCh38", cell_number)
sample_status = rep("Fresh", cell_number) # For each described sample, one single 6 mm punch biopsy was taken and processed immediately.
sample_cultured = rep("No", cell_number)
anatomical_region_level1 = rep(NA, cell_number)
anatomical_region_level2 = rep(NA, cell_number)
anatomical_region_level3 = rep(NA, cell_number)
ethnicity_1 = rep("European Ancestry", cell_number)
ethnicity_2 = rep("Caucasian", cell_number)
ethnicity_detail = rep(NA, cell_number)
age = rep(NA, cell_number)
sex = rep(NA, cell_number)
original_annotation = rep(NA, cell_number)
######
this_donor = "P112"
age[donor_id == this_donor] = "51"
sex[donor_id == this_donor] = "F"
this_donor = "P115"
age[donor_id == this_donor] = "48"
sex[donor_id == this_donor] = "M"
this_donor = "P116"
age[donor_id == this_donor] = "57"
sex[donor_id == this_donor] = "F"
this_donor = "P121"
age[donor_id == this_donor] = "44"
sex[donor_id == this_donor] = "F"
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
if(file.exists(loom_path)){
  unlink(loom_path)
}
save_h5(loom_path, t(as_matrix(this_obj@assays$RNA@counts)))
write.table(this_obj@meta.data, metadata_path, sep = "\t", row.names = T, col.names = T, quote = F)
