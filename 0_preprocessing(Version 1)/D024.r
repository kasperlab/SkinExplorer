r_utilities = "D:/Study/KI/Projects/###/Analysis/code_library/utilities.r"
source(r_utilities)
######
X_gene_path = "D:/Study/KI/Projects/###/Analysis/useful_files/X_Gene_Homo_sapiens.GRCh38.103.rds"
Y_gene_path = "D:/Study/KI/Projects/###/Analysis/useful_files/Y_Gene_Homo_sapiens.GRCh38.103.rds"
alias_gene_c_seurat_path = "D:/Study/KI/Projects/###/Analysis/useful_files/hgnc_complete_set_c_seurat.rds"
alias_gene_c_seurat = readRDS(alias_gene_c_seurat_path)
######
output_dir = "D:/Study/KI/Projects/###/MetaStudies/Data/Skin/Human/Xu_Chen_Nature_2022/convert_output"
dir.create(output_dir, showWarnings = F, recursive = T)
loom_path = paste0(output_dir, "/OMIX691_HC.loom")
metadata_path = paste0(output_dir, "/OMIX691_HC_metadata.tsv")
######
folder = "D:/Study/KI/Projects/###/MetaStudies/Data/Skin/Human/Xu_Chen_Nature_2022"
file_list = c("OMIX691-20-01", "OMIX691-20-02", "OMIX691-20-03", "OMIX691-20-04", "OMIX691-20-05")
######
gene_bc_list = list()
cell_id_list = list()
donor_id = c()
for(ii in file_list){
  str_fragments = unlist(strsplit(ii, "[-]"))
  this_donor = paste0("H", str_fragments[3])
  this_path = paste0(folder, "/", ii)
  gene_bc_list[[this_donor]] = ReadMtx(mtx = paste0(this_path, "/matrix.mtx"), cells = paste0(this_path, "/barcodes.tsv"), features = paste0(this_path, "/genes.tsv"))
  cell_id_list[[this_donor]] = paste0(this_donor, "_", colnames(gene_bc_list[[this_donor]]))
  donor_id = c(donor_id, rep(this_donor, length(cell_id_list[[this_donor]])))
}
this_matrix = merge_matrices(gene_bc_list, is_sparse = T)
rm(gene_bc_list)
this_obj = CreateSeuratObject(this_matrix, min.cells = 1)
this_obj = unify_genename_seurat(this_obj, alias_gene_c_seurat, verbose = T)
rm(this_matrix)
######
# metadata
cell_number = ncol(this_obj)
experiment_id = rep("Xu_Chen_Nature_2022", cell_number)
subject_type = rep("Alive - Healthy", cell_number)
sample_id = paste0(experiment_id, "___", donor_id)
sample_type = rep("Biopsy", cell_number)
sampled_site_condition = rep("Healthy", cell_number)
tissue = rep("Epidermis", cell_number)
biological_unit = rep("Cells", cell_number)
library_platform = rep("10x_3'_v2", cell_number)
strand_sequence = rep("3'", cell_number)
sequencing_platform = rep("Illumina NextSeq 500", cell_number)
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
this_donor = "H01"
age[donor_id == this_donor] = 78
sex[donor_id == this_donor] = "M"
anatomical_region_level1[donor_id == this_donor] = "Extremities"
anatomical_region_level2[donor_id == this_donor] = "Leg"
anatomical_region_level3[donor_id == this_donor] = "Knee"
this_donor = "H02"
age[donor_id == this_donor] = 41
sex[donor_id == this_donor] = "F"
anatomical_region_level1[donor_id == this_donor] = "Torso"
anatomical_region_level2[donor_id == this_donor] = "Hip"
anatomical_region_level3[donor_id == this_donor] = NA
this_donor = "H03"
age[donor_id == this_donor] = 18
sex[donor_id == this_donor] = "F"
anatomical_region_level1[donor_id == this_donor] = "Head"
anatomical_region_level2[donor_id == this_donor] = "Face"
anatomical_region_level3[donor_id == this_donor] = NA
this_donor = "H04"
age[donor_id == this_donor] = 33
sex[donor_id == this_donor] = "M"
anatomical_region_level1[donor_id == this_donor] = "Torso"
anatomical_region_level2[donor_id == this_donor] = "Chest"
anatomical_region_level3[donor_id == this_donor] = NA
this_donor = "H05"
age[donor_id == this_donor] = 41
sex[donor_id == this_donor] = "F"
anatomical_region_level1[donor_id == this_donor] = "Extremities"
anatomical_region_level2[donor_id == this_donor] = "Arm"
anatomical_region_level3[donor_id == this_donor] = "Shoulder"
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
