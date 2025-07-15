r_utilities = "D:/Study/KI/Projects/###/Analysis/code_library/utilities.r"
source(r_utilities)
######
X_gene_path = "D:/Study/KI/Projects/###/Analysis/useful_files/X_Gene_Homo_sapiens.GRCh38.103.rds"
Y_gene_path = "D:/Study/KI/Projects/###/Analysis/useful_files/Y_Gene_Homo_sapiens.GRCh38.103.rds"
alias_gene_c_seurat_path = "D:/Study/KI/Projects/###/Analysis/useful_files/hgnc_complete_set_c_seurat.rds"
alias_gene_c_seurat = readRDS(alias_gene_c_seurat_path)
######
GSE150672_path = "D:/Study/KI/Projects/###/MetaStudies/Data/Skin/Human/Hughes_Shalek_Immunity_2020/GSE150672_Skin_Expression_counts.csv.gz"
GSE150672_rds_path = "D:/Study/KI/Projects/###/MetaStudies/Data/Skin/Human/Hughes_Shalek_Immunity_2020/GSE150672_Skin_Expression_counts.rds"
GSE150672_metadata_path = "D:/Study/KI/Projects/###/MetaStudies/Data/Skin/Human/Hughes_Shalek_Immunity_2020/1-s2.0-S107476132030409X-mmc3.txt"
# Output Dir
output_dir = "D:/Study/KI/Projects/###/MetaStudies/Data/Skin/Human/Hughes_Shalek_Immunity_2020/convert_output"
dir.create(output_dir, showWarnings = F, recursive = T)
loom_path = paste0(output_dir, "/GSE150672_HC.loom")
metadata_path = paste0(output_dir, "/GSE150672_HC_metadata.tsv")
######
if(!file.exists(GSE150672_rds_path)){
  GSE150672 = as.matrix(read.table(gzfile(GSE150672_path), header = T, row.names = 1, sep = ","))
  saveRDS(GSE150672, GSE150672_rds_path)
}else{
  GSE150672 = readRDS(GSE150672_rds_path)
}
######
metadata = as.matrix(read.table(GSE150672_metadata_path, header = T, row.names = 1, sep = "\t"))
######
all_donor_id = sapply(colnames(GSE150672), function(x){
  unlist(strsplit(x, "_"))[1]
})
GSE150672_HC = GSE150672[, all_donor_id %in% c("Normal", "NW", "NS810")]
metadata_HC = metadata[colnames(GSE150672_HC), ]
######
this_obj = CreateSeuratObject(GSE150672_HC, min.cells = 1)
this_obj = unify_genename_seurat(this_obj, alias_gene_c_seurat, verbose = T)
######
# metadata
cell_number = ncol(this_obj)
experiment_id = rep("Hughes_Shalek_Immunity_2020", cell_number)
donor_id = metadata_HC[, "SampleCondition"]
subject_type = rep("Alive - Healthy", cell_number)
sample_id = paste0(experiment_id, "___", donor_id)
sample_type = rep("Biopsy", cell_number)
sampled_site_condition = rep("Healthy", cell_number)
tissue = rep(NA, cell_number)
biological_unit = rep("Cells", cell_number)
library_platform = rep("Seq-well S3", cell_number)
strand_sequence = rep(NA, cell_number)
sequencing_platform = rep("Illumina NextSeq 500", cell_number)
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
original_annotation = metadata_HC[, "Specific_CellType"]
######
this_donor = "Normal1"
sex[donor_id == this_donor] = "F"
age[donor_id == this_donor] = "63"
ethnicity_1[donor_id == this_donor] = "European Ancestry"
ethnicity_detail[donor_id == this_donor] = "White"
this_donor = "Normal2"
sex[donor_id == this_donor] = "F"
age[donor_id == this_donor] = "66-70"
ethnicity_1[donor_id == this_donor] = NA
ethnicity_detail[donor_id == this_donor] = NA
this_donor = "Normal3"
sex[donor_id == this_donor] = "M"
age[donor_id == this_donor] = "56-60"
ethnicity_1[donor_id == this_donor] = "Asian"
ethnicity_detail[donor_id == this_donor] = NA
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
