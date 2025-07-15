r_utilities = "D:/Study/KI/Projects/###/Analysis/code_library/utilities.r"
source(r_utilities)
######
X_gene_path = "D:/Study/KI/Projects/###/Analysis/useful_files/X_Gene_Homo_sapiens.GRCh38.103.rds"
Y_gene_path = "D:/Study/KI/Projects/###/Analysis/useful_files/Y_Gene_Homo_sapiens.GRCh38.103.rds"
alias_gene_c_seurat_path = "D:/Study/KI/Projects/###/Analysis/useful_files/hgnc_complete_set_c_seurat.rds"
alias_gene_c_seurat = readRDS(alias_gene_c_seurat_path)
######
# 1
mtx_path_1 = "D:/Study/KI/Projects/###/MetaStudies/Data/Skin/Human/Wang_Atwood_NatureCommunications_2020/GSM4431901_Donor1_ForeskinEpidermis_10XGenomics_matrix.mtx.gz" #Donor 1
barcodes_path_1 = "D:/Study/KI/Projects/###/MetaStudies/Data/Skin/Human/Wang_Atwood_NatureCommunications_2020/GSM4431901_Donor1_ForeskinEpidermis_10XGenomics_barcodes.tsv.gz" #Donor 1
genes_path_1 = "D:/Study/KI/Projects/###/MetaStudies/Data/Skin/Human/Wang_Atwood_NatureCommunications_2020/GSM4431901_Donor1_ForeskinEpidermis_10XGenomics_genes.tsv.gz" #Donor 1
# 2
mtx_path_2 = "D:/Study/KI/Projects/###/MetaStudies/Data/Skin/Human/Wang_Atwood_NatureCommunications_2020/GSM4431902_Donor2_ForeskinEpidermis_10XGenomics_matrix.mtx.gz" #Donor 2
barcodes_path_2 = "D:/Study/KI/Projects/###/MetaStudies/Data/Skin/Human/Wang_Atwood_NatureCommunications_2020/GSM4431902_Donor2_ForeskinEpidermis_10XGenomics_barcodes.tsv.gz" #Donor 2
genes_path_2 = "D:/Study/KI/Projects/###/MetaStudies/Data/Skin/Human/Wang_Atwood_NatureCommunications_2020/GSM4431902_Donor2_ForeskinEpidermis_10XGenomics_genes.tsv.gz" #Donor 2
# 3
mtx_path_3 = "D:/Study/KI/Projects/###/MetaStudies/Data/Skin/Human/Wang_Atwood_NatureCommunications_2020/GSM4431903_Donor3_ForeskinEpidermis_10XGenomics_matrix.mtx.gz" #Donor 3
barcodes_path_3 = "D:/Study/KI/Projects/###/MetaStudies/Data/Skin/Human/Wang_Atwood_NatureCommunications_2020/GSM4431903_Donor3_ForeskinEpidermis_10XGenomics_barcodes.tsv.gz" #Donor 3
genes_path_3 = "D:/Study/KI/Projects/###/MetaStudies/Data/Skin/Human/Wang_Atwood_NatureCommunications_2020/GSM4431903_Donor3_ForeskinEpidermis_10XGenomics_genes.tsv.gz" #Donor 3
# 4
mtx_path_4 = "D:/Study/KI/Projects/###/MetaStudies/Data/Skin/Human/Wang_Atwood_NatureCommunications_2020/GSM4431904_Donor4_ForeskinEpidermis_10XGenomics_matrix.mtx.gz" #Donor 4
barcodes_path_4 = "D:/Study/KI/Projects/###/MetaStudies/Data/Skin/Human/Wang_Atwood_NatureCommunications_2020/GSM4431904_Donor4_ForeskinEpidermis_10XGenomics_barcodes.tsv.gz" #Donor 4
genes_path_4 = "D:/Study/KI/Projects/###/MetaStudies/Data/Skin/Human/Wang_Atwood_NatureCommunications_2020/GSM4431904_Donor4_ForeskinEpidermis_10XGenomics_genes.tsv.gz" #Donor 4
# 5
mtx_path_5 = "D:/Study/KI/Projects/###/MetaStudies/Data/Skin/Human/Wang_Atwood_NatureCommunications_2020/GSM4431905_Donor5_ForeskinEpidermis_10XGenomics_matrix.mtx.gz" #Donor 5
barcodes_path_5 = "D:/Study/KI/Projects/###/MetaStudies/Data/Skin/Human/Wang_Atwood_NatureCommunications_2020/GSM4431905_Donor5_ForeskinEpidermis_10XGenomics_barcodes.tsv.gz" #Donor 5
genes_path_5 = "D:/Study/KI/Projects/###/MetaStudies/Data/Skin/Human/Wang_Atwood_NatureCommunications_2020/GSM4431905_Donor5_ForeskinEpidermis_10XGenomics_genes.tsv.gz" #Donor 5
######
### Output
output_dir = "D:/Study/KI/Projects/###/MetaStudies/Data/Skin/Human/Wang_Atwood_NatureCommunications_2020/convert_output"
dir.create(output_dir, showWarnings = F, recursive = T)
loom_path = paste0(output_dir, "/GSE147482.loom")
metadata_path = paste0(output_dir, "/GSE147482_metadata.tsv")
###
expression_matrix_1 = ReadMtx(mtx = mtx_path_1, cells = barcodes_path_1, features = genes_path_1)
expression_matrix_2 = ReadMtx(mtx = mtx_path_2, cells = barcodes_path_2, features = genes_path_2)
expression_matrix_3 = ReadMtx(mtx = mtx_path_3, cells = barcodes_path_3, features = genes_path_3)
expression_matrix_4 = ReadMtx(mtx = mtx_path_4, cells = barcodes_path_4, features = genes_path_4)
expression_matrix_5 = ReadMtx(mtx = mtx_path_5, cells = barcodes_path_5, features = genes_path_5)
gene_bc_list = list(expression_matrix_1, expression_matrix_2, expression_matrix_3, expression_matrix_4, expression_matrix_5)
cell_id_list = list(colnames(expression_matrix_1), colnames(expression_matrix_2), colnames(expression_matrix_3), colnames(expression_matrix_4), colnames(expression_matrix_5))
this_matrix = merge_matrices(gene_bc_list, is_sparse = T)
rm(gene_bc_list)
this_obj = CreateSeuratObject(this_matrix, min.cells = 1)
this_obj = unify_genename_seurat(this_obj, alias_gene_c_seurat, verbose = T)
rm(this_matrix)
######
# metadata
cell_number = ncol(this_obj)
experiment_id = rep("Wang_Atwood_NatureCommunications_2020", cell_number)
donor_id = c(rep("Donor1", length(cell_id_list[[1]])), 
             rep("Donor2", length(cell_id_list[[2]])), 
             rep("Donor3", length(cell_id_list[[3]])), 
             rep("Donor4", length(cell_id_list[[4]])), 
             rep("Donor5", length(cell_id_list[[5]])))
subject_type = rep("Alive - Healthy", cell_number)
sample_id = paste0(experiment_id, "___", donor_id)
sample_type = rep("Surgical resection", cell_number)
sampled_site_condition = rep("Healthy", cell_number)
tissue = rep("Epidermis", cell_number)
biological_unit = rep("Cells", cell_number)
library_platform = rep("10x_3'_v2", cell_number)
strand_sequence = rep("3'", cell_number)
sequencing_platform = rep("Illumina NextSeq 500", cell_number)
reference_genome = rep("GRCh38", cell_number)
sample_status = rep(NA, cell_number)
sample_cultured = rep("No", cell_number)
anatomical_region_level1 = rep("Torso", cell_number)
anatomical_region_level2 = rep("Foreskin", cell_number)
anatomical_region_level3 = rep(NA, cell_number)
ethnicity_1 = rep(NA, cell_number)
ethnicity_2 = rep(NA, cell_number)
ethnicity_detail = rep(NA, cell_number)
age = rep("Neonatal", cell_number)
sex = rep("M", cell_number)
original_annotation = rep(NA, cell_number)
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
