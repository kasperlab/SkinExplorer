r_utilities = "D:/Study/KI/Projects/###/Analysis/code_library/utilities.r"
source(r_utilities)
######
X_gene_path = "D:/Study/KI/Projects/###/Analysis/useful_files/X_Gene_Homo_sapiens.GRCh38.103.rds"
Y_gene_path = "D:/Study/KI/Projects/###/Analysis/useful_files/Y_Gene_Homo_sapiens.GRCh38.103.rds"
alias_gene_c_seurat_path = "D:/Study/KI/Projects/###/Analysis/useful_files/hgnc_complete_set_c_seurat.rds"
alias_gene_c_seurat = readRDS(alias_gene_c_seurat_path)
######
output_dir = "D:/Study/KI/Projects/###/MetaStudies/Data/Skin/Human/Dunlap_Rao_JCIInsight_2022/convert_output"
dir.create(output_dir, showWarnings = F, recursive = T)
loom_path = paste0(output_dir, "/GSE186476_HC.loom")
metadata_path = paste0(output_dir, "/GSE186476_HC_metadata.tsv")
######
folder = "D:/Study/KI/Projects/###/MetaStudies/Data/Skin/Human/Dunlap_Rao_JCIInsight_2022/GSE186476_RAW"
file_list = paste0("Healthy_D", sprintf("%.2d", seq(14)), "_counts")
######
gene_bc_list = list()
cell_id_list = list()
donor_id = c()
for(ii in file_list){
  str_fragments = unlist(strsplit(ii, "[_]"))
  this_donor = str_fragments[2]
  this_path = paste0(folder, "/", ii)
  gene_bc_list[[this_donor]] = ReadMtx_1(data_dir = this_path)
  cell_id_list[[this_donor]] = paste0(this_donor, "_", colnames(gene_bc_list[[this_donor]]))
  donor_id = c(donor_id, rep(this_donor, length(cell_id_list[[this_donor]])))
}
this_matrix = merge_matrices(gene_bc_list, is_sparse = T)
cell_id = unlist(cell_id_list)
colnames(this_matrix) = cell_id
rm(cell_id_list)
rm(cell_id)
rm(gene_bc_list)
this_obj = CreateSeuratObject(this_matrix, min.cells = 1)
this_obj = unify_genename_seurat(this_obj, alias_gene_c_seurat, verbose = T)
rm(this_matrix)
######
# metadata
#View(this_obj_0@meta.data)
cell_number = ncol(this_obj)
experiment_id = rep("Dunlap_Rao_JCIInsight_2022", cell_number)
subject_type = rep("Alive - Healthy", cell_number)
sample_type = rep("Biopsy", cell_number)
sampled_site_condition = rep("Healthy", cell_number)
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
X_Genes = readRDS(X_gene_path)
Y_Genes = readRDS(Y_gene_path)
Y_X_ratio = colSums(this_obj[["RNA"]]$counts[which(rownames(this_obj[["RNA"]]$counts) %in% Y_Genes), ]) / colSums(this_obj[["RNA"]]$counts[which(rownames(this_obj[["RNA"]]$counts) %in% X_Genes), ])
this_obj[["Y_X_ratio"]] = Y_X_ratio
this_obj[["DonorID"]] = donor_id
VlnPlot(this_obj, features = "Y_X_ratio", group.by = "DonorID")
ggsave(paste0(output_dir, "/Y_X_ratio.pdf"), width = 15, height = 8)
######
sex[donor_id %in% c("D01", "D02", "D07", "D09", "D10", "D13")] = "F" # Inference by ChrY-ChrX ratio
sex[donor_id %in% c("D03", "D04", "D05", "D06", "D08", "D11", "D12", "D14")] = "M" # Inference by ChrY-ChrX ratio
######
this_obj[["ExperimentID"]] = experiment_id
this_obj[["DonorID"]] = paste0(experiment_id, "___", donor_id)
this_obj[["Age"]] = age
this_obj[["Sex"]] = sex
this_obj[["Ethnicity1"]] = ethnicity_1
this_obj[["Ethnicity2"]] = ethnicity_2
this_obj[["EthnicityDetail"]] = ethnicity_detail
this_obj[["SubjectType"]] = subject_type
this_obj[["SampleID"]] = paste0(experiment_id, "___", donor_id)
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
save_h5(loom_path, t(as_matrix(this_obj[["RNA"]]$counts)))
write.table(this_obj@meta.data, metadata_path, sep = "\t", row.names = T, col.names = T, quote = F)
