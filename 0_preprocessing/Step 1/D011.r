r_utilities = "D:/Study/KI/Projects/###/Analysis/code_library/utilities.r"
source(r_utilities)
######
X_gene_path = "D:/Study/KI/Projects/###/Analysis/useful_files/X_Gene_Homo_sapiens.GRCh38.103.rds"
Y_gene_path = "D:/Study/KI/Projects/###/Analysis/useful_files/Y_Gene_Homo_sapiens.GRCh38.103.rds"
alias_gene_c_seurat_path = "D:/Study/KI/Projects/###/Analysis/useful_files/hgnc_complete_set_c_seurat.rds"
alias_gene_c_seurat = readRDS(alias_gene_c_seurat_path)
###
txt_path_1 = "D:/Study/KI/Projects/###/MetaStudies/Data/Skin/Human/He_Guttman_JournalofAllergyandClinicalImmunology_2020/GSM4430462_MS.sample4.clean.data.txt.gz" #S4_H
txt_path_2 = "D:/Study/KI/Projects/###/MetaStudies/Data/Skin/Human/He_Guttman_JournalofAllergyandClinicalImmunology_2020/GSM4430464_MS.sample6.clean.data.txt.gz" #S6_H
txt_path_3 = "D:/Study/KI/Projects/###/MetaStudies/Data/Skin/Human/He_Guttman_JournalofAllergyandClinicalImmunology_2020/GSM4430466_MS.sample8.clean.data.txt.gz" #S8_H
txt_path_4 = "D:/Study/KI/Projects/###/MetaStudies/Data/Skin/Human/He_Guttman_JournalofAllergyandClinicalImmunology_2020/GSM4430467_MS.sample9.clean.data.txt.gz" #S9_H
txt_path_5 = "D:/Study/KI/Projects/###/MetaStudies/Data/Skin/Human/He_Guttman_JournalofAllergyandClinicalImmunology_2020/GSM4430468_MS.sample10.clean.data.txt.gz" #S10_H
txt_path_6 = "D:/Study/KI/Projects/###/MetaStudies/Data/Skin/Human/He_Guttman_JournalofAllergyandClinicalImmunology_2020/GSM4430470_MS.sample12.clean.data.txt.gz" #S12_H
txt_path_7 = "D:/Study/KI/Projects/###/MetaStudies/Data/Skin/Human/He_Guttman_JournalofAllergyandClinicalImmunology_2020/GSM4430471_MS.sample13.clean.data.txt.gz" #S13_H
txt_path_8 = "D:/Study/KI/Projects/###/MetaStudies/Data/Skin/Human/He_Guttman_JournalofAllergyandClinicalImmunology_2020/GSM4430475_MS.sample17.clean.data.txt.gz" #S17_H
###
output_dir = "D:/Study/KI/Projects/###/MetaStudies/Data/Skin/Human/He_Guttman_JournalofAllergyandClinicalImmunology_2020/convert_output"
dir.create(output_dir, showWarnings = F, recursive = T)
loom_path = paste0(output_dir, "/GSE147424_HC.loom")
metadata_path = paste0(output_dir, "/GSE147424_HC_metadata.tsv")
###
cell_id_list = list()
gene_bc_list = list()
# HC
for(ii in c(txt_path_1,
            txt_path_2,
            txt_path_3,
            txt_path_4,
            txt_path_5,
            txt_path_6,
            txt_path_7,
            txt_path_8)){
  this_dt = fread(ii)
  this_df = as.data.frame(this_dt)
  rownames(this_df) = this_df[, 1]
  gene_bc_list[[ii]] = as.matrix(this_df[, -1])
  cell_id_list[[ii]] = colnames(gene_bc_list[[ii]])
}
this_matrix_normalized = merge_matrices(gene_bc_list, is_sparse = T)
rm(gene_bc_list)
this_matrix_percent = expm1(this_matrix_normalized) / 10000
raw_counts_mat = apply(this_matrix_percent, 2, function(x){
  all_values = sort(unique(x))
  interval = min(all_values[-1] - all_values[-length(all_values)])
  return(x / interval)
})
raw_counts = as(raw_counts_mat, "CsparseMatrix")
this_obj = CreateSeuratObject(raw_counts, min.cells = 1)
this_obj = unify_genename_seurat(this_obj, alias_gene_c_seurat, verbose = T)
rm(this_matrix_normalized)
rm(this_matrix_percent)
rm(raw_counts_mat)
######
# metadata
#View(this_obj_0@meta.data)
cell_number = ncol(this_obj)
experiment_id = rep("He_Guttman_JournalofAllergyandClinicalImmunology_2020", cell_number)
donor_id = c(rep("S4_H", length(cell_id_list[[txt_path_1]])),
             rep("S6_H", length(cell_id_list[[txt_path_2]])),
             rep("S8_H", length(cell_id_list[[txt_path_3]])),
             rep("S9_H", length(cell_id_list[[txt_path_4]])),
             rep("S10_H", length(cell_id_list[[txt_path_5]])),
             rep("S12_H", length(cell_id_list[[txt_path_6]])),
             rep("S13_H", length(cell_id_list[[txt_path_7]])),
             rep("S17_H", length(cell_id_list[[txt_path_8]])))
subject_type = rep("Alive - Healthy", cell_number)
sample_id = paste0(experiment_id, "___", donor_id)
sample_type = rep("Biopsy", cell_number)
sampled_site_condition = rep("Healthy", cell_number)
tissue = rep(NA, cell_number)
biological_unit = rep("Cells", cell_number)
library_platform = rep("10x_3'_v2", cell_number)
strand_sequence = rep("3'", cell_number)
sequencing_platform = rep("Illumina HiSeq 2500", cell_number)
reference_genome = rep("GRCh38", cell_number)
sample_status = rep(NA, cell_number)
sample_cultured = rep("No", cell_number)
anatomical_region_level1 = rep("Extremities", cell_number)
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
Y_X_ratio = colSums(raw_counts[which(rownames(raw_counts) %in% Y_Genes), ]) / colSums(raw_counts[which(rownames(raw_counts) %in% X_Genes), ])
this_obj[["Y_X_ratio"]] = Y_X_ratio
this_obj[["DonorID"]] = donor_id
VlnPlot(this_obj, features = "Y_X_ratio", group.by = "DonorID")
ggsave(paste0(output_dir, "/Y_X_ratio.pdf"), width = 15, height = 8)
rm(raw_counts)
######
this_donor = "S4_H"
sex[donor_id == this_donor] = "M"
sample_status[donor_id == this_donor] = "Frozen"
this_donor = "S6_H"
sex[donor_id == this_donor] = "F"
sample_status[donor_id == this_donor] = "Frozen"
this_donor = "S8_H"
sex[donor_id == this_donor] = "M"
sample_status[donor_id == this_donor] = "Fresh"
this_donor = "S9_H"
sex[donor_id == this_donor] = "M"
sample_status[donor_id == this_donor] = "Frozen"
this_donor = "S10_H"
sex[donor_id == this_donor] = "F"
sample_status[donor_id == this_donor] = "Frozen"
this_donor = "S12_H"
sex[donor_id == this_donor] = "F"
sample_status[donor_id == this_donor] = "Frozen"
this_donor = "S13_H"
sex[donor_id == this_donor] = "M"
sample_status[donor_id == this_donor] = "Frozen"
this_donor = "S17_H"
sex[donor_id == this_donor] = "F"
sample_status[donor_id == this_donor] = "Frozen"
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
