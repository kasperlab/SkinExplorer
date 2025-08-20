r_utilities = "D:/Study/KI/Projects/###/Analysis/code_library/utilities.r"
source(r_utilities)
######
X_gene_path = "D:/Study/KI/Projects/###/Analysis/useful_files/X_Gene_Homo_sapiens.GRCh38.103.rds"
Y_gene_path = "D:/Study/KI/Projects/###/Analysis/useful_files/Y_Gene_Homo_sapiens.GRCh38.103.rds"
alias_gene_c_seurat_path = "D:/Study/KI/Projects/###/Analysis/useful_files/hgnc_complete_set_c_seurat.rds"
alias_gene_c_seurat = readRDS(alias_gene_c_seurat_path)
######
output_dir = "D:/Study/KI/Projects/###/MetaStudies/Data/Skin/Human/Theocharidis_Bhasin_NatureCommunications_2022/convert_output"
dir.create(output_dir, showWarnings = F, recursive = T)
loom_path = paste0(output_dir, "/GSE165816_HC.loom")
metadata_path = paste0(output_dir, "/GSE165816_HC_metadata.tsv")
######
folder = "D:/Study/KI/Projects/###/MetaStudies/Data/Skin/Human/Theocharidis_Bhasin_NatureCommunications_2022/GSE165816_RAW"
file_list = c(
              #"GSM5050521_G1", # levonorgestrel, tretinoin, norethindrone acetate, clobetasol, minocycline, celecoxib
              #"GSM5050540_G16", "GSM5050542_G18", # Pre-diabetes, premature ventricular contractions, colonic polyps
              "GSM5050552_G28", "GSM5050553_G29",
              #"GSM5050555_G31", "GSM5050556_G32", "GSM5050560_G36", # Hypertension, angina
              #"GSM5050574_G50", # Hypertension, heart murmur, hypothyroidism
              "GSM5050534_G10", "GSM5050548_G24",
              #"GSM5050567_G43", # Hypertension, pre-diabetes, hyperlipidemia, history of iron deficiency, history of thrombotic thrombocytopenic purpura
              #"GSM5050538_G14", # Hypothyroidism
              #"GSM5050568_G44", # Asthma, hypertension
              "GSM5050564_G40")
######
gene_bc_list = list()
sample_id = c()
for(ii in file_list){
  str_fragments = unlist(strsplit(ii, "[_]"))
  this_sample = str_fragments[2]
  this_path = paste0(folder, "/", ii, "counts.csv.gz")
  this_dt = fread(this_path)
  this_df = as.data.frame(this_dt)
  rownames(this_df) = this_df[, 1]
  this_matrix = as.matrix(this_df[, -1])
  this_cell_id = paste0(this_sample, "_", colnames(this_matrix))
  colnames(this_matrix) = this_cell_id
  sample_id = c(sample_id, rep(this_sample, length(this_cell_id)))
  this_seurat_obj = CreateSeuratObject(this_matrix, min.cells = 1)
  this_seurat_obj = unify_genename_seurat(this_seurat_obj, alias_gene_c_seurat, verbose = T)
  gene_bc_list[[ii]] = as_matrix(this_seurat_obj[["RNA"]]$counts)
  rm(this_dt, this_df, this_matrix, this_seurat_obj)
  print(ii)
}
gene_bc_matrix = merge_matrices(gene_bc_list, is_sparse = T)
rm(gene_bc_list)
this_obj = CreateSeuratObject(gene_bc_matrix, min.cells = 1)
rm(gene_bc_matrix)
######
# metadata
cell_number = ncol(this_obj)
experiment_id = rep("Theocharidis_Bhasin_NatureCommunications_2022", cell_number)
donor_id = rep(NA, cell_number)
subject_type = rep("Alive - Healthy", cell_number)
sample_type = rep(NA, cell_number)
sampled_site_condition = rep("Healthy", cell_number)
tissue = rep(NA, cell_number)
biological_unit = rep("Cells", cell_number)
library_platform = rep("10x_3'", cell_number)
strand_sequence = rep("3'", cell_number)
sequencing_platform = rep("Illumina NovaSeq 6000 S4", cell_number)
reference_genome = rep("GRCh38", cell_number)
sample_status = rep(NA, cell_number)
sample_cultured = rep("No", cell_number)
anatomical_region_level1 = rep("Extremities", cell_number)
anatomical_region_level2 = rep(NA, cell_number)
anatomical_region_level3 = rep(NA, cell_number)
ethnicity_1 = rep("European Ancestry", cell_number)
ethnicity_2 = rep(NA, cell_number)
ethnicity_detail = rep("White", cell_number)
age = rep(NA, cell_number)
sex = rep(NA, cell_number)
original_annotation = rep(NA, cell_number)
######
this_sample = "G28"
age[sample_id == this_sample] = 74
sex[sample_id == this_sample] = "M"
donor_id[sample_id == this_sample] = "H03"
anatomical_region_level2[sample_id == this_sample] = "Foot"
anatomical_region_level3[sample_id == this_sample] = "Dorsal foot"
sample_type[sample_id == this_sample] = "Surgical resection"
this_sample = "G29"
age[sample_id == this_sample] = 74
sex[sample_id == this_sample] = "M"
donor_id[sample_id == this_sample] = "H03"
anatomical_region_level2[sample_id == this_sample] = "Arm"
anatomical_region_level3[sample_id == this_sample] = "Forearm"
sample_type[sample_id == this_sample] = "Biopsy"
this_sample = "G10"
age[sample_id == this_sample] = 59
sex[sample_id == this_sample] = "M"
donor_id[sample_id == this_sample] = "H06"
anatomical_region_level2[sample_id == this_sample] = "Foot"
anatomical_region_level3[sample_id == this_sample] = "Dorsal foot"
sample_type[sample_id == this_sample] = "Surgical resection"
this_sample = "G24"
age[sample_id == this_sample] = 59
sex[sample_id == this_sample] = "M"
donor_id[sample_id == this_sample] = "H06"
anatomical_region_level2[sample_id == this_sample] = "Foot"
anatomical_region_level3[sample_id == this_sample] = "Dorsal foot"
sample_type[sample_id == this_sample] = "Surgical resection"
this_sample = "G40"
age[sample_id == this_sample] = 34
sex[sample_id == this_sample] = "F"
donor_id[sample_id == this_sample] = "H10"
anatomical_region_level2[sample_id == this_sample] = "Foot"
anatomical_region_level3[sample_id == this_sample] = "Dorsal foot"
sample_type[sample_id == this_sample] = "Surgical resection"
######
this_obj[["ExperimentID"]] = experiment_id
this_obj[["DonorID"]] = paste0(experiment_id, "___", donor_id)
this_obj[["Age"]] = age
this_obj[["Sex"]] = sex
this_obj[["Ethnicity1"]] = ethnicity_1
this_obj[["Ethnicity2"]] = ethnicity_2
this_obj[["EthnicityDetail"]] = ethnicity_detail
this_obj[["SubjectType"]] = subject_type
this_obj[["SampleID"]] = paste0(experiment_id, "___", sample_id)
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
