r_utilities = "D:/Study/KI/Projects/###/Analysis/code_library/utilities.r"
source(r_utilities)
######
X_gene_path = "D:/Study/KI/Projects/###/Analysis/useful_files/X_Gene_Homo_sapiens.GRCh38.103.rds"
Y_gene_path = "D:/Study/KI/Projects/###/Analysis/useful_files/Y_Gene_Homo_sapiens.GRCh38.103.rds"
alias_gene_c_seurat_path = "D:/Study/KI/Projects/###/Analysis/useful_files/hgnc_complete_set_c_seurat.rds"
alias_gene_c_seurat = readRDS(alias_gene_c_seurat_path)
######
output_dir = "D:/Study/KI/Projects/###/MetaStudies/Data/Skin/Human/Rojahn_Brunner_JournalofAllergyandClinicalImmunology_2020/convert_output"
dir.create(output_dir, showWarnings = F, recursive = T)
loom_path = paste0(output_dir, "/GSE153760_HC.loom")
metadata_path = paste0(output_dir, "/GSE153760_HC_metadata.tsv")
######
folder = "D:/Study/KI/Projects/###/MetaStudies/Data/Skin/Human/Rojahn_Brunner_JournalofAllergyandClinicalImmunology_2020/GSE153760_RAW"
file_list = c("GSM4653863_HC1", "GSM4653864_HC2", "GSM4653865_HC3", "GSM4653866_HC4", "GSM4653867_HC5", "GSM4653868_HC6", "GSM4653869_HC7")
######
gene_bc_list = list()
cell_id_list = list()
donor_id = c()
for(ii in file_list){
  str_fragments = unlist(strsplit(ii, "[_]"))
  this_donor = str_fragments[2]
  this_path = paste0(folder, "/", ii)
  gene_bc_list[[this_donor]] = ReadMtx(mtx = paste0(this_path, "_matrix.mtx.gz"), cells = paste0(this_path, "_barcodes.tsv.gz"), features = paste0(this_path, "_features.tsv.gz"))
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
experiment_id = rep("Rojahn_Brunner_JournalofAllergyandClinicalImmunology_2020", cell_number)
subject_type = rep("Alive - Healthy", cell_number)
sample_id = paste0(experiment_id, "___", donor_id)
sample_type = rep(NA, cell_number)
sampled_site_condition = rep("Healthy", cell_number)
tissue = rep("Epidermis", cell_number)
biological_unit = rep("Cells", cell_number)
library_platform = rep(NA, cell_number)
strand_sequence = rep("3'", cell_number)
sequencing_platform = rep("Illumina HiSeq 3000", cell_number)
reference_genome = rep("GRCh38", cell_number)
sample_status = rep(NA, cell_number)
sample_cultured = rep("No", cell_number)
anatomical_region_level1 = rep("Extremities", cell_number)
anatomical_region_level2 = rep("Arm", cell_number)
anatomical_region_level3 = rep("Elbow", cell_number)
ethnicity_1 = rep("European Ancestry", cell_number)
ethnicity_2 = rep("Caucasian", cell_number)
ethnicity_detail = rep(NA, cell_number)
age = rep(NA, cell_number)
sex = rep(NA, cell_number)
original_annotation = rep(NA, cell_number)
######
this_donor = "HC1"
age[donor_id == this_donor] = 40
sex[donor_id == this_donor] = "M"
library_platform[donor_id == this_donor] = "10x_3'_v2"
sample_type[donor_id == this_donor] = "Suction blister"
this_donor = "HC2"
age[donor_id == this_donor] = 42
sex[donor_id == this_donor] = "F"
library_platform[donor_id == this_donor] = "10x_3'_v2"
sample_type[donor_id == this_donor] = "Suction blister"
this_donor = "HC3"
age[donor_id == this_donor] = 47
sex[donor_id == this_donor] = "F"
library_platform[donor_id == this_donor] = "10x_3'_v2"
sample_type[donor_id == this_donor] = "Suction blister"
this_donor = "HC4"
age[donor_id == this_donor] = 49
sex[donor_id == this_donor] = "F"
library_platform[donor_id == this_donor] = "10x_3'_v2"
sample_type[donor_id == this_donor] = "Suction blister"
this_donor = "HC5"
age[donor_id == this_donor] = 39
sex[donor_id == this_donor] = "M"
library_platform[donor_id == this_donor] = "10x_3'_v3"
sample_type[donor_id == this_donor] = "Suction blister"
this_donor = "HC6"
age[donor_id == this_donor] = 27
sex[donor_id == this_donor] = "F"
library_platform[donor_id == this_donor] = "10x_3'_v3"
sample_type[donor_id == this_donor] = "Biopsy"
this_donor = "HC7"
age[donor_id == this_donor] = 42
sex[donor_id == this_donor] = "F"
library_platform[donor_id == this_donor] = "10x_3'_v3"
sample_type[donor_id == this_donor] = "Biopsy"
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
