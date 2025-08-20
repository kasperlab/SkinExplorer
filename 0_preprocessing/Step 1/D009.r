r_utilities = "D:/Study/KI/Projects/###/Analysis/code_library/utilities.r"
source(r_utilities)
######
X_gene_path = "D:/Study/KI/Projects/###/Analysis/useful_files/X_Gene_Homo_sapiens.GRCh38.103.rds"
Y_gene_path = "D:/Study/KI/Projects/###/Analysis/useful_files/Y_Gene_Homo_sapiens.GRCh38.103.rds"
alias_gene_c_seurat_path = "D:/Study/KI/Projects/###/Analysis/useful_files/hgnc_complete_set_c_seurat.rds"
alias_gene_c_seurat = readRDS(alias_gene_c_seurat_path)
######
# GSE144236_CAL27_path = "D:/Study/KI/Projects/###/MetaStudies/Data/Skin/Human/Ji_Khavari_Cell_2020/GSE144236_CAL27_counts.txt.gz"
# GSE144236_CAL27_vitro_path = "D:/Study/KI/Projects/###/MetaStudies/Data/Skin/Human/Ji_Khavari_Cell_2020/GSE144236_CAL27_vitro_counts.txt.gz"
GSE144236_cSCC_path = "D:/Study/KI/Projects/###/MetaStudies/Data/Skin/Human/Ji_Khavari_Cell_2020/GSE144236_cSCC_counts.txt.gz"
GSE144236_cSCC_rds_path = "D:/Study/KI/Projects/###/MetaStudies/Data/Skin/Human/Ji_Khavari_Cell_2020/GSE144236_cSCC_counts.rds"
GSE144236_patient_metadata_path = "D:/Study/KI/Projects/###/MetaStudies/Data/Skin/Human/Ji_Khavari_Cell_2020/GSE144236_patient_metadata_new.txt.gz"
# GSE144236_SCC13_path = "D:/Study/KI/Projects/###/MetaStudies/Data/Skin/Human/Ji_Khavari_Cell_2020/GSE144236_SCC13_counts.txt.gz"
# GSE144236_XG_TME_path = "D:/Study/KI/Projects/###/MetaStudies/Data/Skin/Human/Ji_Khavari_Cell_2020/GSE144236_XG_TME_counts.txt.gz"
# Output Dir
output_dir = "D:/Study/KI/Projects/###/MetaStudies/Data/Skin/Human/Ji_Khavari_Cell_2020/convert_output"
dir.create(output_dir, showWarnings = F, recursive = T)
loom_path = paste0(output_dir, "/GSE151091_raw_matrix.loom")
metadata_path = paste0(output_dir, "/GSE151091_raw_matrix_metadata.tsv")
######
# GSE144236_CAL27 = as.matrix(read.table(gzfile(GSE144236_CAL27_path), header = T, row.names = 1, sep = "\t"))
# GSE144236_CAL27_vitro = as.matrix(read.table(gzfile(GSE144236_CAL27_vitro_path), header = T, row.names = 1, sep = "\t"))
if(!file.exists(GSE144236_cSCC_rds_path)){
  GSE144236_cSCC = as.matrix(read.table(gzfile(GSE144236_cSCC_path), header = T, row.names = 1, sep = "\t"))
  saveRDS(GSE144236_cSCC, GSE144236_cSCC_rds_path)
}else{
  GSE144236_cSCC = readRDS(GSE144236_cSCC_rds_path)
}
raw_count = GSE144236_cSCC[-c(1, 2), ]
GSE144236_patient_metadata = as.matrix(read.table(gzfile(GSE144236_patient_metadata_path), header = T, row.names = 1, sep = "\t"))
# GSE144236_SCC13 = as.matrix(read.table(gzfile(GSE144236_SCC13_path), header = T, row.names = 1, sep = "\t"))
# GSE144236_XG_TME = as.matrix(read.table(gzfile(GSE144236_XG_TME_path), header = T, row.names = 1, sep = "\t"))
######
this_obj = CreateSeuratObject(raw_count, min.cells = 1)
this_obj = unify_genename_seurat(this_obj, alias_gene_c_seurat, verbose = T)
metadata = GSE144236_patient_metadata[colnames(this_obj), ]
######
# metadata
cell_number = ncol(this_obj)
experiment_id = rep("Ji_Khavari_Cell_2020", cell_number)
donor_id = metadata[, "patient"]
subject_type = rep("Alive - Disease", cell_number)
sample_id = paste0(experiment_id, "___", donor_id)
sample_type = rep("Surgical resection", cell_number)
sampled_site_condition = metadata[, "tum.norm"]
tissue = rep(NA, cell_number)
biological_unit = rep("Cells", cell_number)
library_platform = rep("10x_3'_v2", cell_number)
strand_sequence = rep("3'", cell_number)
sequencing_platform = rep("Illumina HiSeq 4000", cell_number)
reference_genome = rep("GRCh37", cell_number)
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
original_annotation = metadata[, "level3_celltype"]
######
sampled_site_condition[sampled_site_condition == "Normal"] = "Disease adjacent"
this_donor = "P1"
anatomical_region_level1[donor_id == this_donor] = "Extremities"
anatomical_region_level2[donor_id == this_donor] = "Hand"
anatomical_region_level3[donor_id == this_donor] = "Dorsal hand" ### right dorsal hand
age[donor_id == this_donor] = 83
sex[donor_id == this_donor] = "M"
this_donor = "P2"
anatomical_region_level1[donor_id == this_donor] = "Extremities"
anatomical_region_level2[donor_id == this_donor] = "Arm"
anatomical_region_level3[donor_id == this_donor] = "Forearm" ### left forearm
age[donor_id == this_donor] = 69
sex[donor_id == this_donor] = "M"
this_donor = "P3"
anatomical_region_level1[donor_id == this_donor] = "Head"
anatomical_region_level2[donor_id == this_donor] = "Scalp"
anatomical_region_level3[donor_id == this_donor] = "Vertex" ### vertex scalp
age[donor_id == this_donor] = 71
sex[donor_id == this_donor] = "F"
this_donor = "P4"
anatomical_region_level1[donor_id == this_donor] = "Extremities"
anatomical_region_level2[donor_id == this_donor] = "Hand"
anatomical_region_level3[donor_id == this_donor] = "Dorsal hand" ### left dorsal hand
age[donor_id == this_donor] = 75
sex[donor_id == this_donor] = "F"
this_donor = "P5"
anatomical_region_level1[donor_id == this_donor] = "Head"
anatomical_region_level2[donor_id == this_donor] = "Scalp"
anatomical_region_level3[donor_id == this_donor] = "Vertex" ### left vertex scalp
age[donor_id == this_donor] = 98
sex[donor_id == this_donor] = "M"
this_donor = "P6"
anatomical_region_level1[donor_id == this_donor] = "Torso"
anatomical_region_level2[donor_id == this_donor] = "Chest"
anatomical_region_level3[donor_id == this_donor] = NA ### left chest
age[donor_id == this_donor] = 96
sex[donor_id == this_donor] = "F"
this_donor = "P7"
anatomical_region_level1[donor_id == this_donor] = "Extremities"
anatomical_region_level2[donor_id == this_donor] = "Arm"
anatomical_region_level3[donor_id == this_donor] = "Shoulder" ### right shoulder
age[donor_id == this_donor] = 75
sex[donor_id == this_donor] = "M"
this_donor = "P8"
anatomical_region_level1[donor_id == this_donor] = "Head"
anatomical_region_level2[donor_id == this_donor] = "Face"
anatomical_region_level3[donor_id == this_donor] = "Forehead" ### left forehead
age[donor_id == this_donor] = 63
sex[donor_id == this_donor] = "M"
this_donor = "P9"
anatomical_region_level1[donor_id == this_donor] = "Extremities"
anatomical_region_level2[donor_id == this_donor] = "Arm"
anatomical_region_level3[donor_id == this_donor] = "Forearm" ### right forearm
age[donor_id == this_donor] = 87
sex[donor_id == this_donor] = "M"
this_donor = "P10"
anatomical_region_level1[donor_id == this_donor] = "Head"
anatomical_region_level2[donor_id == this_donor] = "Face"
anatomical_region_level3[donor_id == this_donor] = "Tragus" ### right tragus
age[donor_id == this_donor] = 71
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
split_list = SplitObject(this_obj, split.by = "SampledSiteCondition")
######
if(file.exists(loom_path)){
  unlink(loom_path)
}
save_h5(loom_path, t(as_matrix(split_list[["Disease adjacent"]]@assays$RNA@counts)))
write.table(split_list[["Disease adjacent"]]@meta.data, metadata_path, sep = "\t", row.names = T, col.names = T, quote = F)



