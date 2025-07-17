r_utilities = "D:/Study/KI/Projects/###/Analysis/code_library/utilities.r"
source(r_utilities)
######
X_gene_path = "D:/Study/KI/Projects/###/Analysis/useful_files/X_Gene_Homo_sapiens.GRCh38.103.rds"
Y_gene_path = "D:/Study/KI/Projects/###/Analysis/useful_files/Y_Gene_Homo_sapiens.GRCh38.103.rds"
alias_gene_c_seurat_path = "D:/Study/KI/Projects/###/Analysis/useful_files/hgnc_complete_set_c_seurat.rds"
alias_gene_c_seurat = readRDS(alias_gene_c_seurat_path)
######
csv_path = "D:/Study/KI/Projects/###/MetaStudies/Data/Skin/Human/Belote_Torres_NatureCellBiology_2021/GSE151091_raw_matrix.csv.gz"
######
output_dir = "D:/Study/KI/Projects/###/MetaStudies/Data/Skin/Human/Belote_Torres_NatureCellBiology_2021/convert_output"
dir.create(output_dir, showWarnings = F, recursive = T)
loom_path = paste0(output_dir, "/GSE151091_raw_matrix.loom")
metadata_path = paste0(output_dir, "/GSE151091_raw_matrix_metadata.tsv")
this_dt = fread(csv_path)
this_df = as.data.frame(this_dt)
rownames(this_df) = this_df[, 1]
this_matrix = as.matrix(this_df[, -1])
this_obj = CreateSeuratObject(this_matrix, min.cells = 1)
this_obj = unify_genename_seurat(this_obj, alias_gene_c_seurat, verbose = T)
######
# metadata
#View(this_obj_0@meta.data)
cell_number = ncol(this_obj)
experiment_id = rep("Belote_Torres_NatureCellBiology_2021", cell_number)
donor_id = sapply(colnames(this_obj), function(x){
  str_fragments = unlist(strsplit(x, "_"))
  return(paste(str_fragments[seq(length(str_fragments) - 2)], collapse = "_"))
})
subject_type = rep("Alive - Healthy", cell_number)
sample_id = paste0(experiment_id, "___", donor_id)
sample_type = rep("Surgical resection", cell_number)
sampled_site_condition = rep("Healthy", cell_number)
tissue = rep("Epidermis", cell_number)
biological_unit = rep("Cells", cell_number)
library_platform = rep("Smart-Seq2", cell_number)
strand_sequence = rep("Full", cell_number)
sequencing_platform = rep("Illumina NovaSeq 6000", cell_number) # GSE151091
reference_genome = rep("GRCh38", cell_number)
sample_status = rep("Fresh", cell_number) # Tissue dissociation was started on the same day as sample acquisition.
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
X_gene_path = "D:/Study/KI/Projects/###/Analysis/HumanSkinAtlas/X_Gene_Homo_sapiens.GRCh38.103.rds"
Y_gene_path = "D:/Study/KI/Projects/###/Analysis/HumanSkinAtlas/Y_Gene_Homo_sapiens.GRCh38.103.rds"
X_Genes = readRDS(X_gene_path)
Y_Genes = readRDS(Y_gene_path)
Y_X_ratio = colSums(this_matrix[which(rownames(this_matrix) %in% Y_Genes), ]) / colSums(this_matrix[which(rownames(this_matrix) %in% X_Genes), ])
this_obj[["Y_X_ratio"]] = Y_X_ratio
this_obj[["DonorID"]] = donor_id
VlnPlot(this_obj, features = "Y_X_ratio", group.by = "DonorID")
ggsave(paste0(output_dir, "/Y_X_ratio.pdf"), width = 15, height = 8)
rm(this_matrix)
######
this_donor = "9.5WK02"
anatomical_region_level1[donor_id == this_donor] = "Extremities"
anatomical_region_level2[donor_id == this_donor] = NA
anatomical_region_level3[donor_id == this_donor] = NA ### leg,arm
ethnicity_1[donor_id == this_donor] = "Unknown"
age[donor_id == this_donor] = "FW_9.5"
sex[donor_id == this_donor] = "F" # Inference by ChrY-ChrX ratio
this_donor = "10WK03"
anatomical_region_level1[donor_id == this_donor] = "Extremities"
anatomical_region_level2[donor_id == this_donor] = NA
anatomical_region_level3[donor_id == this_donor] = NA ### leg,arm,palm
ethnicity_1[donor_id == this_donor] = "Unknown"
age[donor_id == this_donor] = "FW_10"
sex[donor_id == this_donor] = "M" # Inference by ChrY-ChrX ratio
this_donor = "12WKM01"
anatomical_region_level1[donor_id == this_donor] = "Extremities"
anatomical_region_level2[donor_id == this_donor] = NA
anatomical_region_level3[donor_id == this_donor] = NA ### leg,palm,sole
ethnicity_1[donor_id == this_donor] = "Unknown"
age[donor_id == this_donor] = "FW_12"
sex[donor_id == this_donor] = "M"
donor_id[donor_id == this_donor] = "12WK01"
this_donor = "12WK05"
anatomical_region_level1[donor_id == this_donor] = "Extremities"
anatomical_region_level2[donor_id == this_donor] = NA
anatomical_region_level3[donor_id == this_donor] = NA ### sole,leg,palm,arm
ethnicity_1[donor_id == this_donor] = "Unknown"
age[donor_id == this_donor] = "FW_12"
sex[donor_id == this_donor] = "M" # Inference by ChrY-ChrX ratio
this_donor = "16WKM04"
anatomical_region_level1[donor_id == this_donor] = "Extremities"
anatomical_region_level2[donor_id == this_donor] = "Hand"
anatomical_region_level3[donor_id == this_donor] = "Palm" ### palm
ethnicity_1[donor_id == this_donor] = "Unknown"
age[donor_id == this_donor] = "FW_16"
sex[donor_id == this_donor] = "M" # Inference by ChrY-ChrX ratio
donor_id[donor_id == this_donor] = "16WK04"
this_donor = "18WKM06"
anatomical_region_level1[donor_id == this_donor] = "Extremities"
anatomical_region_level2[donor_id == this_donor] = NA
anatomical_region_level3[donor_id == this_donor] = NA ### leg,palm,sole
ethnicity_1[donor_id == this_donor] = "Unknown"
age[donor_id == this_donor] = "FW_18"
sex[donor_id == this_donor] = "M" # Inference by ChrY-ChrX ratio
donor_id[donor_id == this_donor] = "18WK06"
this_donor = "FS030_LM"
anatomical_region_level1[donor_id == this_donor] = "Torso"
anatomical_region_level2[donor_id == this_donor] = "Foreskin"
anatomical_region_level3[donor_id == this_donor] = NA
ethnicity_1[donor_id == this_donor] = "Unknown"
age[donor_id == this_donor] = 0
sex[donor_id == this_donor] = "M"
this_donor = "FS043_LM"
anatomical_region_level1[donor_id == this_donor] = "Torso"
anatomical_region_level2[donor_id == this_donor] = "Foreskin"
anatomical_region_level3[donor_id == this_donor] = NA
ethnicity_1[donor_id == this_donor] = "Unknown"
age[donor_id == this_donor] = 0
sex[donor_id == this_donor] = "M"
this_donor = "A1021M"
anatomical_region_level1[donor_id == this_donor] = "Extremities"
anatomical_region_level2[donor_id == this_donor] = "Foot"
anatomical_region_level3[donor_id == this_donor] = NA # ankle or upper foot
ethnicity_1[donor_id == this_donor] = "Asian"
age[donor_id == this_donor] = 24
sex[donor_id == this_donor] = "F"
donor_id[donor_id == this_donor] = "A1021"
this_donor = "A1038LM"
anatomical_region_level1[donor_id == this_donor] = "Extremities"
anatomical_region_level2[donor_id == this_donor] = NA
anatomical_region_level3[donor_id == this_donor] = NA # arch,calf,heel,shin
ethnicity_1[donor_id == this_donor] = "European Ancestry"
ethnicity_detail[donor_id == this_donor] = "White"
age[donor_id == this_donor] = 35
sex[donor_id == this_donor] = "M"
donor_id[donor_id == this_donor] = "A1038"
this_donor = "A1012M"
anatomical_region_level1[donor_id == this_donor] = "Extremities"
anatomical_region_level2[donor_id == this_donor] = "Arm"
anatomical_region_level3[donor_id == this_donor] = NA
ethnicity_1[donor_id == this_donor] = "Asian"
age[donor_id == this_donor] = 37
sex[donor_id == this_donor] = "F"
donor_id[donor_id == this_donor] = "A1012"
this_donor = "A1022M"
anatomical_region_level1[donor_id == this_donor] = "Extremities"
anatomical_region_level2[donor_id == this_donor] = "Arm"
anatomical_region_level3[donor_id == this_donor] = NA
ethnicity_1[donor_id == this_donor] = "Hispanic/Latinx"
age[donor_id == this_donor] = 42
sex[donor_id == this_donor] = "M"
donor_id[donor_id == this_donor] = "A1022"
this_donor = "A1015LM"
anatomical_region_level1[donor_id == this_donor] = "Extremities"
anatomical_region_level2[donor_id == this_donor] = "Arm"
anatomical_region_level3[donor_id == this_donor] = NA
ethnicity_1[donor_id == this_donor] = "Asian"
age[donor_id == this_donor] = 52
sex[donor_id == this_donor] = "F"
donor_id[donor_id == this_donor] = "A1015"
this_donor = "A1025L"
anatomical_region_level1[donor_id == this_donor] = "Extremities"
anatomical_region_level2[donor_id == this_donor] = "Leg"
anatomical_region_level3[donor_id == this_donor] = "Thigh" # left thigh
ethnicity_1[donor_id == this_donor] = "European Ancestry"
ethnicity_detail[donor_id == this_donor] = "White"
age[donor_id == this_donor] = 56
sex[donor_id == this_donor] = "M"
donor_id[donor_id == this_donor] = "A1025"
this_donor = "A1016LM"
anatomical_region_level1[donor_id == this_donor] = "Extremities"
anatomical_region_level2[donor_id == this_donor] = "Arm"
anatomical_region_level3[donor_id == this_donor] = "Forearm" # anterior forearm
ethnicity_1[donor_id == this_donor] = "European Ancestry"
ethnicity_detail[donor_id == this_donor] = "White"
age[donor_id == this_donor] = 58
sex[donor_id == this_donor] = "M"
donor_id[donor_id == this_donor] = "A1016"
this_donor = "A1020LM"
anatomical_region_level1[donor_id == this_donor] = "Extremities"
anatomical_region_level2[donor_id == this_donor] = "Arm"
anatomical_region_level3[donor_id == this_donor] = NA
ethnicity_1[donor_id == this_donor] = "Asian"
age[donor_id == this_donor] = 60
sex[donor_id == this_donor] = "F"
donor_id[donor_id == this_donor] = "A1020"
this_donor = "A1033M"
anatomical_region_level1[donor_id == this_donor] = "Extremities"
anatomical_region_level2[donor_id == this_donor] = "Leg"
anatomical_region_level3[donor_id == this_donor] = "Thigh"
ethnicity_1[donor_id == this_donor] = "Hispanic/Latinx"
age[donor_id == this_donor] = 61
sex[donor_id == this_donor] = "M"
donor_id[donor_id == this_donor] = "A1033"
this_donor = "A1011L"
anatomical_region_level1[donor_id == this_donor] = "Extremities"
anatomical_region_level2[donor_id == this_donor] = "Arm"
anatomical_region_level3[donor_id == this_donor] = NA
ethnicity_1[donor_id == this_donor] = "European Ancestry"
ethnicity_detail[donor_id == this_donor] = "White"
age[donor_id == this_donor] = 65
sex[donor_id == this_donor] = "M"
donor_id[donor_id == this_donor] = "A1011"
this_donor = "A1026L"
anatomical_region_level1[donor_id == this_donor] = "Extremities"
anatomical_region_level2[donor_id == this_donor] = "Leg"
anatomical_region_level3[donor_id == this_donor] = "Thigh" # right thigh
ethnicity_1[donor_id == this_donor] = "European Ancestry"
ethnicity_detail[donor_id == this_donor] = "White"
age[donor_id == this_donor] = 66
sex[donor_id == this_donor] = "M"
donor_id[donor_id == this_donor] = "A1026"
this_donor = "A1014L"
anatomical_region_level1[donor_id == this_donor] = "Extremities"
anatomical_region_level2[donor_id == this_donor] = "Leg"
anatomical_region_level3[donor_id == this_donor] = "Thigh"
ethnicity_1[donor_id == this_donor] = "European Ancestry"
ethnicity_detail[donor_id == this_donor] = "White"
age[donor_id == this_donor] = 68
sex[donor_id == this_donor] = "F"
donor_id[donor_id == this_donor] = "A1014"
this_donor = "A1046M"
anatomical_region_level1[donor_id == this_donor] = "Extremities"
anatomical_region_level2[donor_id == this_donor] = NA
anatomical_region_level3[donor_id == this_donor] = NA # shin,heel,arch,calf
ethnicity_1[donor_id == this_donor] = "Hispanic/Latinx"
age[donor_id == this_donor] = 77
sex[donor_id == this_donor] = "F"
donor_id[donor_id == this_donor] = "A1046"
this_donor = "A1017LM"
anatomical_region_level1[donor_id == this_donor] = "Extremities"
anatomical_region_level2[donor_id == this_donor] = "Leg"
anatomical_region_level3[donor_id == this_donor] = "Calf"
ethnicity_1[donor_id == this_donor] = "Hispanic/Latinx"
age[donor_id == this_donor] = 81
original_annotation[donor_id == this_donor] = "Melanocyte"
sex[donor_id == this_donor] = "M"
donor_id[donor_id == this_donor] = "A1017"
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
filter_cell =
  setdiff(seq(cell_number), grep("FW", this_obj$Age))
filtered_obj = this_obj[, filter_cell]
filter_gene = rowSums(filtered_obj) > 0
filtered_obj = filtered_obj[filter_gene, ]
######
if(file.exists(loom_path)){
  unlink(loom_path)
}
save_h5(loom_path, t(as_matrix(filtered_obj@assays$RNA@counts)))
write.table(filtered_obj@meta.data, metadata_path, sep = "\t", row.names = T, col.names = T, quote = F)
