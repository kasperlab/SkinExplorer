r_utilities = "D:/Study/KI/Projects/###/Analysis/code_library/utilities.r"
source(r_utilities)
######
X_gene_path = "D:/Study/KI/Projects/###/Analysis/useful_files/X_Gene_Homo_sapiens.GRCh38.103.rds"
Y_gene_path = "D:/Study/KI/Projects/###/Analysis/useful_files/Y_Gene_Homo_sapiens.GRCh38.103.rds"
alias_gene_c_seurat_path = "D:/Study/KI/Projects/###/Analysis/useful_files/hgnc_complete_set_c_seurat.rds"
alias_gene_c_seurat = readRDS(alias_gene_c_seurat_path)
######
csv_path_1 = "D:/Study/KI/Projects/###/MetaStudies/Data/Skin/Human/Gaydosik_Fuschiotti_ClinicalCancerResearch_2019/GSM3679038_Labeled_SC50_011917_SK_NOR_GRCh38raw.csv.gz" #HC1
csv_path_2 = "D:/Study/KI/Projects/###/MetaStudies/Data/Skin/Human/Gaydosik_Fuschiotti_ClinicalCancerResearch_2019/GSM3679039_Labeled_SC68_051517_SK_NOR_GRCh38raw.csv.gz" #HC2
csv_path_3 = "D:/Study/KI/Projects/###/MetaStudies/Data/Skin/Human/Gaydosik_Fuschiotti_ClinicalCancerResearch_2019/GSM3679040_Labeled_SC124_080317_SK_NOR_GRCh38raw.csv.gz" #HC3
csv_path_4 = "D:/Study/KI/Projects/###/MetaStudies/Data/Skin/Human/Gaydosik_Fuschiotti_ClinicalCancerResearch_2019/GSM3679041_Labeled_SC125_080317_SK_NOR_GRCh38raw.csv.gz" #HC4
###
output_dir = "D:/Study/KI/Projects/###/MetaStudies/Data/Skin/Human/Gaydosik_Fuschiotti_ClinicalCancerResearch_2019/convert_output"
dir.create(output_dir, showWarnings = F, recursive = T)
loom_path = paste0(output_dir, "/GSE128531_HC.loom")
metadata_path = paste0(output_dir, "/GSE128531_HC_metadata.tsv")
###
cell_id_list = list()
gene_bc_list = list()
# HC
for(ii in c(csv_path_1,
            csv_path_2,
            csv_path_3,
            csv_path_4)){
  this_dt = fread(ii)
  this_df = as.data.frame(this_dt)
  rownames(this_df) = this_df[, 1]
  gene_bc_list[[ii]] = as.matrix(this_df[, -1])
  cell_id_list[[ii]] = colnames(gene_bc_list[[ii]])
}
this_matrix = merge_matrices(gene_bc_list, is_sparse = T)
rm(gene_bc_list)
this_obj = CreateSeuratObject(this_matrix, min.cells = 1)
this_obj = unify_genename_seurat(this_obj, alias_gene_c_seurat, verbose = T)
rm(this_matrix)
######
# metadata
#View(this_obj_0@meta.data)
cell_number = ncol(this_obj)
experiment_id = rep("Gaydosik_Fuschiotti_ClinicalCancerResearch_2019", cell_number)
donor_id = c(rep("HC-1", length(cell_id_list[[csv_path_1]])), 
             rep("HC-2", length(cell_id_list[[csv_path_2]])), 
             rep("HC-3", length(cell_id_list[[csv_path_3]])), 
             rep("HC-4", length(cell_id_list[[csv_path_4]])))
subject_type = rep(NA, cell_number)
sample_id = paste0(experiment_id, "___", donor_id)
sample_type = rep("Biopsy", cell_number)
sampled_site_condition = rep("Healthy", cell_number)
tissue = rep(NA, cell_number)
biological_unit = rep("Cells", cell_number)
library_platform = rep("10x_3'_v2", cell_number)
strand_sequence = rep("3'", cell_number)
sequencing_platform = rep("Illumina HiSeq 4000", cell_number)
reference_genome = rep("GRCh38", cell_number)
sample_status = rep("Fresh", cell_number) # immediately after resection, no longer than 1 h
sample_cultured = rep("No", cell_number)
anatomical_region_level1 = rep("Torso", cell_number)
anatomical_region_level2 = rep("Groins", cell_number)
anatomical_region_level3 = rep(NA, cell_number)
ethnicity_1 = rep(NA, cell_number)
ethnicity_2 = rep(NA, cell_number)
ethnicity_detail = rep(NA, cell_number)
age = rep(NA, cell_number)
sex = rep(NA, cell_number)
original_annotation = rep(NA, cell_number)
######
this_donor = "HC-1"
age[donor_id == this_donor] = 64
sex[donor_id == this_donor] = "M"
this_donor = "HC-2"
age[donor_id == this_donor] = 48
sex[donor_id == this_donor] = "F"
this_donor = "HC-3"
age[donor_id == this_donor] = 54
sex[donor_id == this_donor] = "M"
this_donor = "HC-4"
age[donor_id == this_donor] = 61
sex[donor_id == this_donor] = "M"
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
