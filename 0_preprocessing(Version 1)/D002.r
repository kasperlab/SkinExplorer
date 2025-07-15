r_utilities = "D:/Study/KI/Projects/###/Analysis/code_library/utilities.r"
source(r_utilities)
######
X_gene_path = "D:/Study/KI/Projects/###/Analysis/useful_files/X_Gene_Homo_sapiens.GRCh38.103.rds"
Y_gene_path = "D:/Study/KI/Projects/###/Analysis/useful_files/Y_Gene_Homo_sapiens.GRCh38.103.rds"
alias_gene_c_seurat_path = "D:/Study/KI/Projects/###/Analysis/useful_files/hgnc_complete_set_c_seurat.rds"
alias_gene_c_seurat = readRDS(alias_gene_c_seurat_path)
######
seurat_obj = readRDS("D:/Study/KI/Projects/###/MetaStudies/Data/From Ning/NingGroup_Seurat_afterQC.rds")
###
output_dir = "/home/haoy/projects/MetaStudiesAnalysis/Integration/Data/Ning/convert_output"
dir.create(output_dir, showWarnings = F, recursive = T)
loom_path = paste0(output_dir, "/Ning.loom")
metadata_path = paste0(output_dir, "/Ning_metadata.tsv")
healthy_loom_path = paste0(output_dir, "/Ning_healthy.loom")
healthy_metadata_path = paste0(output_dir, "/Ning_healthy_metadata.tsv")
######
this_obj = CreateSeuratObject(seurat_obj@assays$RNA@counts, min.cells = 1)
this_obj = unify_genename_seurat(this_obj, alias_gene_c_seurat, verbose = T)
######
# metadata
cell_number = ncol(this_obj)
experiment_id = rep("Ning", cell_number)
donor_id = seurat_obj@meta.data$Patient
subject_type = rep("Alive - Healthy", cell_number)
sample_id = paste0(experiment_id, "___", seurat_obj@meta.data$orig.ident)
sample_type = rep("Biopsy", cell_number)
sampled_site_condition = seurat_obj@meta.data$Condition
tissue = rep(NA, cell_number) # Full Skin
biological_unit = rep("Cells", cell_number)
library_platform = rep("10x_3'_v3", cell_number)
strand_sequence = rep("3'", cell_number)
sequencing_platform = rep("Illumina NovaSeq 6000 S4", cell_number)
reference_genome = rep("GRCh38", cell_number)
sample_status = rep("Fresh", cell_number)
sample_cultured = rep("No", cell_number)
anatomical_region_level1 = rep("Torso", cell_number)
anatomical_region_level2 = rep("Back", cell_number)
anatomical_region_level3 = rep("Lower back", cell_number)
ethnicity_1 = rep("European Ancestry", cell_number)
ethnicity_2 = rep("Caucasian", cell_number)
ethnicity_detail = rep(NA, cell_number)
age = seurat_obj@meta.data$Age
sex = as.character(seurat_obj@meta.data$Gender)
original_annotation = rep(NA, cell_number)
######
sex[sex == "Male"] = "M"
sex[sex == "Female"] = "F"
sampled_site_condition[sampled_site_condition == "Skin"] = "Healthy"
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
save_h5(loom_path, t(as.matrix(this_obj@assays$RNA@counts)))
write.table(this_obj@meta.data, metadata_path, sep = "\t", row.names = T, col.names = T, quote = F)
######
split_list = SplitObject(this_obj, split.by = "SampledSiteCondition")
######
if(file.exists(healthy_loom_path)){
  unlink(healthy_loom_path)
}
save_h5(healthy_loom_path, t(as_matrix(split_list[["Healthy"]]@assays$RNA@counts)))
write.table(split_list[["Healthy"]]@meta.data, healthy_metadata_path, sep = "\t", row.names = T, col.names = T, quote = F)



