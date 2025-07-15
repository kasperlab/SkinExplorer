r_utilities = "D:/Study/KI/Projects/###/Analysis/code_library/utilities.r"
py_utilities = "D:/Study/KI/Projects/###/Analysis/code_library/utilities.py"
source(r_utilities)
library(reticulate)
py_run_file(py_utilities, convert = F)
######
X_gene_path = "D:/Study/KI/Projects/###/Analysis/useful_files/X_Gene_Homo_sapiens.GRCh38.103.rds"
Y_gene_path = "D:/Study/KI/Projects/###/Analysis/useful_files/Y_Gene_Homo_sapiens.GRCh38.103.rds"
alias_gene_c_seurat_path = "D:/Study/KI/Projects/###/Analysis/useful_files/hgnc_complete_set_c_seurat.rds"
alias_gene_c_seurat = readRDS(alias_gene_c_seurat_path)
######
###
UCSFSKIN1_ABD4epi_path = "D:/Study/KI/Projects/###/MetaStudies/Data/Internal/Hao/Cheng_2018_raw_data/UCSFSKIN1_ABD4epi.txt"
UCSFSKIN2_BRST41epi_path = "D:/Study/KI/Projects/###/MetaStudies/Data/Internal/Hao/Cheng_2018_raw_data/UCSFSKIN2_BRST41epi.txt"
UCSFSKIN3_BRST53epi_path = "D:/Study/KI/Projects/###/MetaStudies/Data/Internal/Hao/Cheng_2018_raw_data/UCSFSKIN3_BRST53epi.txt"
UCSFSKIN4_FORE12epi_path = "D:/Study/KI/Projects/###/MetaStudies/Data/Internal/Hao/Cheng_2018_raw_data/UCSFSKIN4_FORE12epi.txt"
UCSFSKIN5_FORE8epi_path = "D:/Study/KI/Projects/###/MetaStudies/Data/Internal/Hao/Cheng_2018_raw_data/UCSFSKIN5_FORE8epi.txt"
UCSFSKIN6_FORE9epi_path = "D:/Study/KI/Projects/###/MetaStudies/Data/Internal/Hao/Cheng_2018_raw_data/UCSFSKIN6_FORE9epi.txt"
#UCSFSKIN7_PSO48epi_path = "D:/Study/KI/Projects/###/MetaStudies/Data/Internal/Hao/Cheng_2018_raw_data/UCSFSKIN7_PSO48epi.txt"
#UCSFSKIN8_PSO49epi_path = "D:/Study/KI/Projects/###/MetaStudies/Data/Internal/Hao/Cheng_2018_raw_data/UCSFSKIN8_PSO49epi.txt"
#UCSFSKIN9_PSO14epi_path = "D:/Study/KI/Projects/###/MetaStudies/Data/Internal/Hao/Cheng_2018_raw_data/UCSFSKIN9_PSO14epi.txt"
UCSFSKIN10_SCALP11epi_path = "D:/Study/KI/Projects/###/MetaStudies/Data/Internal/Hao/Cheng_2018_raw_data/UCSFSKIN10_SCALP11epi.txt"
UCSFSKIN11_SCALP26epi_path = "D:/Study/KI/Projects/###/MetaStudies/Data/Internal/Hao/Cheng_2018_raw_data/UCSFSKIN11_SCALP26epi.txt"
UCSFSKIN12_SCALP32epi_path = "D:/Study/KI/Projects/###/MetaStudies/Data/Internal/Hao/Cheng_2018_raw_data/UCSFSKIN12_SCALP32epi.txt"
# Output Dir
output_dir = "D:/Study/KI/Projects/###/MetaStudies/Data/Skin/Human/Cheng_Cho_CellReports_2018/convert_output"
dir.create(output_dir, showWarnings = F, recursive = T)
# Use GTF
gtf_path = "D:/Study/KI/Projects/###/Analysis/useful_files/Homo_sapiens.GRCh38.103.gtf.gz"
rds_path = "D:/Study/KI/Projects/###/Analysis/useful_files/Homo_sapiens.GRCh38.103_ChrID_GeneID_TranscriptID_GeneName.rds"
txt_path = "D:/Study/KI/Projects/###/Analysis/useful_files/Homo_sapiens.GRCh38.103_ChrID_GeneID_TranscriptID_GeneName.txt"
GTF_info = extract_gtf_information(gtf_path, rds_path, txt_path)
X_Genes = unique(as.character(GTF_info[GTF_info[, "ChrID"] == "X", c("GeneID", "GeneName")]))
Y_Genes = unique(as.character(GTF_info[GTF_info[, "ChrID"] == "Y", c("GeneID", "GeneName")]))
######
donor = "UCSFSKIN1_ABD4epi"
tmp = as.data.frame(fread(UCSFSKIN1_ABD4epi_path))
loom_path = paste0(output_dir, "/", donor, ".loom")
metadata_path = paste0(output_dir, "/", donor, "_metadata.tsv")
# load raw data
this_matrix = as.matrix(tmp[, -1])
rownames(this_matrix) = as.character(tmp[, 1])
this_matrix = geneid_to_genename(this_matrix, GTF_info[, "GeneID"], GTF_info[, "GeneName"], verbose = T)
# Metadata
experiment_id = rep("Cheng_Cho_CellReports_2018", ncol(this_matrix))
donor_id = rep(donor, ncol(this_matrix))
metadata = data.frame(row.names = colnames(this_matrix),
                      ExperimentID = experiment_id,
                      DonorID = paste0(experiment_id, "___", donor_id),
                      Age = rep(NA, ncol(this_matrix)),
                      Sex = rep(NA, ncol(this_matrix)),
                      Ethnicity1 = rep(NA, ncol(this_matrix)),
                      Ethnicity2 = rep(NA, ncol(this_matrix)),
                      EthnicityDetail = rep(NA, ncol(this_matrix)),
                      SubjectType = rep("Alive - Healthy", ncol(this_matrix)),
                      SampleID = paste0(experiment_id, "___", donor_id),
                      SampleType = rep("Surgical resection", ncol(this_matrix)),
                      SampledSiteCondition = rep("Healthy", ncol(this_matrix)),
                      Tissue = rep("Epidermis", ncol(this_matrix)),
                      BiologicalUnit = rep("Cells", ncol(this_matrix)),
                      LibraryPlatform = rep("10x_3'_v2", ncol(this_matrix)),
                      StrandSequence = rep("3'", ncol(this_matrix)),
                      SequencingPlatform = rep("lllumina HiSeq 2500/4000, Illumina NovaSeq 6000 S2", ncol(this_matrix)),
                      ReferenceGenome = rep("GRCh38", ncol(this_matrix)),
                      SampleStatus = rep(NA, ncol(this_matrix)),
                      SampleCultured = rep("No", ncol(this_matrix)),
                      AnatomicalRegionLevel1 = rep("Torso", ncol(this_matrix)),
                      AnatomicalRegionLevel2 = rep("Abdominal", ncol(this_matrix)),
                      AnatomicalRegionLevel3 = rep(NA, ncol(this_matrix)),
                      OriginalAnnotation = rep(NA, ncol(this_matrix)),
                      Y_X_ratio = colSums(this_matrix[which(rownames(this_matrix) %in% Y_Genes), ]) / colSums(this_matrix[which(rownames(this_matrix) %in% X_Genes), ]))
if(file.exists(loom_path)){
  unlink(loom_path)
}
save_h5(loom_path, t(this_matrix))
write.table(metadata, metadata_path, sep = "\t", row.names = T, col.names = T, quote = F)
######
######
donor = "UCSFSKIN2_BRST41epi"
tmp = as.data.frame(fread(UCSFSKIN2_BRST41epi_path))
loom_path = paste0(output_dir, "/", donor, ".loom")
metadata_path = paste0(output_dir, "/", donor, "_metadata.tsv")
# load raw data
this_matrix = as.matrix(tmp[, -1])
rownames(this_matrix) = as.character(tmp[, 1])
this_matrix = geneid_to_genename(this_matrix, GTF_info[, "GeneID"], GTF_info[, "GeneName"], verbose = T)
# Metadata
donor_id = rep(donor, ncol(this_matrix))
experiment_id = rep("Cheng_Cho_CellReports_2018", ncol(this_matrix))
metadata = data.frame(row.names = colnames(this_matrix),
                      ExperimentID = experiment_id,
                      DonorID = paste0(experiment_id, "___", donor_id),
                      Age = rep(NA, ncol(this_matrix)),
                      Sex = rep(NA, ncol(this_matrix)),
                      Ethnicity1 = rep(NA, ncol(this_matrix)),
                      Ethnicity2 = rep(NA, ncol(this_matrix)),
                      EthnicityDetail = rep(NA, ncol(this_matrix)),
                      SubjectType = rep("Alive - Healthy", ncol(this_matrix)),
                      SampleID = paste0(experiment_id, "___", donor_id),
                      SampleType = rep("Surgical resection", ncol(this_matrix)),
                      SampledSiteCondition = rep("Healthy", ncol(this_matrix)),
                      Tissue = rep("Epidermis", ncol(this_matrix)),
                      BiologicalUnit = rep("Cells", ncol(this_matrix)),
                      LibraryPlatform = rep("10x_3'_v2", ncol(this_matrix)),
                      StrandSequence = rep("3'", ncol(this_matrix)),
                      SequencingPlatform = rep("lllumina HiSeq 2500/4000, Illumina NovaSeq 6000 S2", ncol(this_matrix)),
                      ReferenceGenome = rep("GRCh38", ncol(this_matrix)),
                      SampleStatus = rep(NA, ncol(this_matrix)),
                      SampleCultured = rep("No", ncol(this_matrix)),
                      AnatomicalRegionLevel1 = rep("Torso", ncol(this_matrix)),
                      AnatomicalRegionLevel2 = rep("Breast", ncol(this_matrix)),
                      AnatomicalRegionLevel3 = rep(NA, ncol(this_matrix)),
                      OriginalAnnotation = rep(NA, ncol(this_matrix)),
                      Y_X_ratio = colSums(this_matrix[which(rownames(this_matrix) %in% Y_Genes), ]) / colSums(this_matrix[which(rownames(this_matrix) %in% X_Genes), ]))
if(file.exists(loom_path)){
  unlink(loom_path)
}
save_h5(loom_path, t(this_matrix))
write.table(metadata, metadata_path, sep = "\t", row.names = T, col.names = T, quote = F)
######
######
donor = "UCSFSKIN3_BRST53epi"
tmp = as.data.frame(fread(UCSFSKIN3_BRST53epi_path))
loom_path = paste0(output_dir, "/", donor, ".loom")
metadata_path = paste0(output_dir, "/", donor, "_metadata.tsv")
# load raw data
this_matrix = as.matrix(tmp[, -1])
rownames(this_matrix) = as.character(tmp[, 1])
this_matrix = geneid_to_genename(this_matrix, GTF_info[, "GeneID"], GTF_info[, "GeneName"], verbose = T)
# Metadata
donor_id = rep(donor, ncol(this_matrix))
experiment_id = rep("Cheng_Cho_CellReports_2018", ncol(this_matrix))
metadata = data.frame(row.names = colnames(this_matrix),
                      ExperimentID = experiment_id,
                      DonorID = paste0(experiment_id, "___", donor_id),
                      Age = rep(NA, ncol(this_matrix)),
                      Sex = rep(NA, ncol(this_matrix)),
                      Ethnicity1 = rep(NA, ncol(this_matrix)),
                      Ethnicity2 = rep(NA, ncol(this_matrix)),
                      EthnicityDetail = rep(NA, ncol(this_matrix)),
                      SubjectType = rep("Alive - Healthy", ncol(this_matrix)),
                      SampleID = paste0(experiment_id, "___", donor_id),
                      SampleType = rep("Surgical resection", ncol(this_matrix)),
                      SampledSiteCondition = rep("Healthy", ncol(this_matrix)),
                      Tissue = rep("Epidermis", ncol(this_matrix)),
                      BiologicalUnit = rep("Cells", ncol(this_matrix)),
                      LibraryPlatform = rep("10x_3'_v2", ncol(this_matrix)),
                      StrandSequence = rep("3'", ncol(this_matrix)),
                      SequencingPlatform = rep("lllumina HiSeq 2500/4000, Illumina NovaSeq 6000 S2", ncol(this_matrix)),
                      ReferenceGenome = rep("GRCh38", ncol(this_matrix)),
                      SampleStatus = rep(NA, ncol(this_matrix)),
                      SampleCultured = rep("No", ncol(this_matrix)),
                      AnatomicalRegionLevel1 = rep("Torso", ncol(this_matrix)),
                      AnatomicalRegionLevel2 = rep("Breast", ncol(this_matrix)),
                      AnatomicalRegionLevel3 = rep(NA, ncol(this_matrix)),
                      OriginalAnnotation = rep(NA, ncol(this_matrix)),                      Y_X_ratio = colSums(this_matrix[which(rownames(this_matrix) %in% Y_Genes), ]) / colSums(this_matrix[which(rownames(this_matrix) %in% X_Genes), ]))
if(file.exists(loom_path)){
  unlink(loom_path)
}
save_h5(loom_path, t(this_matrix))
write.table(metadata, metadata_path, sep = "\t", row.names = T, col.names = T, quote = F)
######
######
donor = "UCSFSKIN4_FORE12epi"
tmp = as.data.frame(fread(UCSFSKIN4_FORE12epi_path))
loom_path = paste0(output_dir, "/", donor, ".loom")
metadata_path = paste0(output_dir, "/", donor, "_metadata.tsv")
# load raw data
this_matrix = as.matrix(tmp[, -1])
rownames(this_matrix) = as.character(tmp[, 1])
this_matrix = geneid_to_genename(this_matrix, GTF_info[, "GeneID"], GTF_info[, "GeneName"], verbose = T)
# Metadata
donor_id = rep(donor, ncol(this_matrix))
experiment_id = rep("Cheng_Cho_CellReports_2018", ncol(this_matrix))
metadata = data.frame(row.names = colnames(this_matrix),
                      ExperimentID = experiment_id,
                      DonorID = paste0(experiment_id, "___", donor_id),
                      Age = rep(NA, ncol(this_matrix)),
                      Sex = rep(NA, ncol(this_matrix)),
                      Ethnicity1 = rep(NA, ncol(this_matrix)),
                      Ethnicity2 = rep(NA, ncol(this_matrix)),
                      EthnicityDetail = rep(NA, ncol(this_matrix)),
                      SubjectType = rep("Alive - Healthy", ncol(this_matrix)),
                      SampleID = paste0(experiment_id, "___", donor_id),
                      SampleType = rep("Surgical resection", ncol(this_matrix)),
                      SampledSiteCondition = rep("Healthy", ncol(this_matrix)),
                      Tissue = rep("Epidermis", ncol(this_matrix)),
                      BiologicalUnit = rep("Cells", ncol(this_matrix)),
                      LibraryPlatform = rep("10x_3'_v2", ncol(this_matrix)),
                      StrandSequence = rep("3'", ncol(this_matrix)),
                      SequencingPlatform = rep("lllumina HiSeq 2500/4000, Illumina NovaSeq 6000 S2", ncol(this_matrix)),
                      ReferenceGenome = rep("GRCh38", ncol(this_matrix)),
                      SampleStatus = rep(NA, ncol(this_matrix)),
                      SampleCultured = rep("No", ncol(this_matrix)),
                      AnatomicalRegionLevel1 = rep("Torso", ncol(this_matrix)),
                      AnatomicalRegionLevel2 = rep("Foreskin", ncol(this_matrix)),
                      AnatomicalRegionLevel3 = rep(NA, ncol(this_matrix)),
                      OriginalAnnotation = rep(NA, ncol(this_matrix)),
                      Y_X_ratio = colSums(this_matrix[which(rownames(this_matrix) %in% Y_Genes), ]) / colSums(this_matrix[which(rownames(this_matrix) %in% X_Genes), ]))
if(file.exists(loom_path)){
  unlink(loom_path)
}
save_h5(loom_path, t(this_matrix))
write.table(metadata, metadata_path, sep = "\t", row.names = T, col.names = T, quote = F)
######
######
donor = "UCSFSKIN5_FORE8epi"
tmp = as.data.frame(fread(UCSFSKIN5_FORE8epi_path))
loom_path = paste0(output_dir, "/", donor, ".loom")
metadata_path = paste0(output_dir, "/", donor, "_metadata.tsv")
# load raw data
this_matrix = as.matrix(tmp[, -1])
rownames(this_matrix) = as.character(tmp[, 1])
this_matrix = geneid_to_genename(this_matrix, GTF_info[, "GeneID"], GTF_info[, "GeneName"], verbose = T)
# Metadata
donor_id = rep(donor, ncol(this_matrix))
experiment_id = rep("Cheng_Cho_CellReports_2018", ncol(this_matrix))
metadata = data.frame(row.names = colnames(this_matrix),
                      ExperimentID = experiment_id,
                      DonorID = paste0(experiment_id, "___", donor_id),
                      Age = rep(NA, ncol(this_matrix)),
                      Sex = rep(NA, ncol(this_matrix)),
                      Ethnicity1 = rep(NA, ncol(this_matrix)),
                      Ethnicity2 = rep(NA, ncol(this_matrix)),
                      EthnicityDetail = rep(NA, ncol(this_matrix)),
                      SubjectType = rep("Alive - Healthy", ncol(this_matrix)),
                      SampleID = paste0(experiment_id, "___", donor_id),
                      SampleType = rep("Surgical resection", ncol(this_matrix)),
                      SampledSiteCondition = rep("Healthy", ncol(this_matrix)),
                      Tissue = rep("Epidermis", ncol(this_matrix)),
                      BiologicalUnit = rep("Cells", ncol(this_matrix)),
                      LibraryPlatform = rep("10x_3'_v2", ncol(this_matrix)),
                      StrandSequence = rep("3'", ncol(this_matrix)),
                      SequencingPlatform = rep("lllumina HiSeq 2500/4000, Illumina NovaSeq 6000 S2", ncol(this_matrix)),
                      ReferenceGenome = rep("GRCh38", ncol(this_matrix)),
                      SampleStatus = rep(NA, ncol(this_matrix)),
                      SampleCultured = rep("No", ncol(this_matrix)),
                      AnatomicalRegionLevel1 = rep("Torso", ncol(this_matrix)),
                      AnatomicalRegionLevel2 = rep("Foreskin", ncol(this_matrix)),
                      AnatomicalRegionLevel3 = rep(NA, ncol(this_matrix)),
                      OriginalAnnotation = rep(NA, ncol(this_matrix)),
                      Y_X_ratio = colSums(this_matrix[which(rownames(this_matrix) %in% Y_Genes), ]) / colSums(this_matrix[which(rownames(this_matrix) %in% X_Genes), ]))
if(file.exists(loom_path)){
  unlink(loom_path)
}
save_h5(loom_path, t(this_matrix))
write.table(metadata, metadata_path, sep = "\t", row.names = T, col.names = T, quote = F)
######
######
donor = "UCSFSKIN6_FORE9epi"
tmp = as.data.frame(fread(UCSFSKIN6_FORE9epi_path))
loom_path = paste0(output_dir, "/", donor, ".loom")
metadata_path = paste0(output_dir, "/", donor, "_metadata.tsv")
# load raw data
this_matrix = as.matrix(tmp[, -1])
rownames(this_matrix) = as.character(tmp[, 1])
this_matrix = geneid_to_genename(this_matrix, GTF_info[, "GeneID"], GTF_info[, "GeneName"], verbose = T)
# Metadata
donor_id = rep(donor, ncol(this_matrix))
experiment_id = rep("Cheng_Cho_CellReports_2018", ncol(this_matrix))
metadata = data.frame(row.names = colnames(this_matrix),
                      ExperimentID = experiment_id,
                      DonorID = paste0(experiment_id, "___", donor_id),
                      Age = rep(NA, ncol(this_matrix)),
                      Sex = rep(NA, ncol(this_matrix)),
                      Ethnicity1 = rep(NA, ncol(this_matrix)),
                      Ethnicity2 = rep(NA, ncol(this_matrix)),
                      EthnicityDetail = rep(NA, ncol(this_matrix)),
                      SubjectType = rep("Alive - Healthy", ncol(this_matrix)),
                      SampleID = paste0(experiment_id, "___", donor_id),
                      SampleType = rep("Surgical resection", ncol(this_matrix)),
                      SampledSiteCondition = rep("Healthy", ncol(this_matrix)),
                      Tissue = rep("Epidermis", ncol(this_matrix)),
                      BiologicalUnit = rep("Cells", ncol(this_matrix)),
                      LibraryPlatform = rep("10x_3'_v2", ncol(this_matrix)),
                      StrandSequence = rep("3'", ncol(this_matrix)),
                      SequencingPlatform = rep("lllumina HiSeq 2500/4000, Illumina NovaSeq 6000 S2", ncol(this_matrix)),
                      ReferenceGenome = rep("GRCh38", ncol(this_matrix)),
                      SampleStatus = rep(NA, ncol(this_matrix)),
                      SampleCultured = rep("No", ncol(this_matrix)),
                      AnatomicalRegionLevel1 = rep("Torso", ncol(this_matrix)),
                      AnatomicalRegionLevel2 = rep("Foreskin", ncol(this_matrix)),
                      AnatomicalRegionLevel3 = rep(NA, ncol(this_matrix)),
                      OriginalAnnotation = rep(NA, ncol(this_matrix)),
                      Y_X_ratio = colSums(this_matrix[which(rownames(this_matrix) %in% Y_Genes), ]) / colSums(this_matrix[which(rownames(this_matrix) %in% X_Genes), ]))
if(file.exists(loom_path)){
  unlink(loom_path)
}
save_h5(loom_path, t(this_matrix))
write.table(metadata, metadata_path, sep = "\t", row.names = T, col.names = T, quote = F)
######
######
donor = "UCSFSKIN10_SCALP11epi"
tmp = as.data.frame(fread(UCSFSKIN10_SCALP11epi_path))
loom_path = paste0(output_dir, "/", donor, ".loom")
metadata_path = paste0(output_dir, "/", donor, "_metadata.tsv")
# load raw data
this_matrix = as.matrix(tmp[, -1])
rownames(this_matrix) = as.character(tmp[, 1])
this_matrix = geneid_to_genename(this_matrix, GTF_info[, "GeneID"], GTF_info[, "GeneName"], verbose = T)
# Metadata
donor_id = rep(donor, ncol(this_matrix))
experiment_id = rep("Cheng_Cho_CellReports_2018", ncol(this_matrix))
metadata = data.frame(row.names = colnames(this_matrix),
                      ExperimentID = experiment_id,
                      DonorID = paste0(experiment_id, "___", donor_id),
                      Age = rep(NA, ncol(this_matrix)),
                      Sex = rep(NA, ncol(this_matrix)),
                      Ethnicity1 = rep(NA, ncol(this_matrix)),
                      Ethnicity2 = rep(NA, ncol(this_matrix)),
                      EthnicityDetail = rep(NA, ncol(this_matrix)),
                      SubjectType = rep("Alive - Healthy", ncol(this_matrix)),
                      SampleID = paste0(experiment_id, "___", donor_id),
                      SampleType = rep("Surgical resection", ncol(this_matrix)),
                      SampledSiteCondition = rep("Healthy", ncol(this_matrix)),
                      Tissue = rep("Epidermis", ncol(this_matrix)),
                      BiologicalUnit = rep("Cells", ncol(this_matrix)),
                      LibraryPlatform = rep("10x_3'_v2", ncol(this_matrix)),
                      StrandSequence = rep("3'", ncol(this_matrix)),
                      SequencingPlatform = rep("lllumina HiSeq 2500/4000, Illumina NovaSeq 6000 S2", ncol(this_matrix)),
                      ReferenceGenome = rep("GRCh38", ncol(this_matrix)),
                      SampleStatus = rep(NA, ncol(this_matrix)),
                      SampleCultured = rep("No", ncol(this_matrix)),
                      AnatomicalRegionLevel1 = rep("Head", ncol(this_matrix)),
                      AnatomicalRegionLevel2 = rep("Scalp", ncol(this_matrix)),
                      AnatomicalRegionLevel3 = rep(NA, ncol(this_matrix)),
                      OriginalAnnotation = rep(NA, ncol(this_matrix)),
                      Y_X_ratio = colSums(this_matrix[which(rownames(this_matrix) %in% Y_Genes), ]) / colSums(this_matrix[which(rownames(this_matrix) %in% X_Genes), ]))
if(file.exists(loom_path)){
  unlink(loom_path)
}
save_h5(loom_path, t(this_matrix))
write.table(metadata, metadata_path, sep = "\t", row.names = T, col.names = T, quote = F)
######
######
donor = "UCSFSKIN11_SCALP26epi"
tmp = as.data.frame(fread(UCSFSKIN11_SCALP26epi_path))
loom_path = paste0(output_dir, "/", donor, ".loom")
metadata_path = paste0(output_dir, "/", donor, "_metadata.tsv")
# load raw data
this_matrix = as.matrix(tmp[, -1])
rownames(this_matrix) = as.character(tmp[, 1])
this_matrix = geneid_to_genename(this_matrix, GTF_info[, "GeneID"], GTF_info[, "GeneName"], verbose = T)
# Metadata
donor_id = rep(donor, ncol(this_matrix))
experiment_id = rep("Cheng_Cho_CellReports_2018", ncol(this_matrix))
metadata = data.frame(row.names = colnames(this_matrix),
                      ExperimentID = experiment_id,
                      DonorID = paste0(experiment_id, "___", donor_id),
                      Age = rep(NA, ncol(this_matrix)),
                      Sex = rep(NA, ncol(this_matrix)),
                      Ethnicity1 = rep(NA, ncol(this_matrix)),
                      Ethnicity2 = rep(NA, ncol(this_matrix)),
                      EthnicityDetail = rep(NA, ncol(this_matrix)),
                      SubjectType = rep("Alive - Healthy", ncol(this_matrix)),
                      SampleID = paste0(experiment_id, "___", donor_id),
                      SampleType = rep("Surgical resection", ncol(this_matrix)),
                      SampledSiteCondition = rep("Healthy", ncol(this_matrix)),
                      Tissue = rep("Epidermis", ncol(this_matrix)),
                      BiologicalUnit = rep("Cells", ncol(this_matrix)),
                      LibraryPlatform = rep("10x_3'_v2", ncol(this_matrix)),
                      StrandSequence = rep("3'", ncol(this_matrix)),
                      SequencingPlatform = rep("lllumina HiSeq 2500/4000, Illumina NovaSeq 6000 S2", ncol(this_matrix)),
                      ReferenceGenome = rep("GRCh38", ncol(this_matrix)),
                      SampleStatus = rep(NA, ncol(this_matrix)),
                      SampleCultured = rep("No", ncol(this_matrix)),
                      AnatomicalRegionLevel1 = rep("Head", ncol(this_matrix)),
                      AnatomicalRegionLevel2 = rep("Scalp", ncol(this_matrix)),
                      AnatomicalRegionLevel3 = rep(NA, ncol(this_matrix)),
                      OriginalAnnotation = rep(NA, ncol(this_matrix)),
                      Y_X_ratio = colSums(this_matrix[which(rownames(this_matrix) %in% Y_Genes), ]) / colSums(this_matrix[which(rownames(this_matrix) %in% X_Genes), ]))
if(file.exists(loom_path)){
  unlink(loom_path)
}
save_h5(loom_path, t(this_matrix))
write.table(metadata, metadata_path, sep = "\t", row.names = T, col.names = T, quote = F)
######
######
donor = "UCSFSKIN12_SCALP32epi"
tmp = as.data.frame(fread(UCSFSKIN12_SCALP32epi_path))
loom_path = paste0(output_dir, "/", donor, ".loom")
metadata_path = paste0(output_dir, "/", donor, "_metadata.tsv")
# load raw data
this_matrix = as.matrix(tmp[, -1])
rownames(this_matrix) = as.character(tmp[, 1])
this_matrix = geneid_to_genename(this_matrix, GTF_info[, "GeneID"], GTF_info[, "GeneName"], verbose = T)
# Metadata
donor_id = rep(donor, ncol(this_matrix))
experiment_id = rep("Cheng_Cho_CellReports_2018", ncol(this_matrix))
metadata = data.frame(row.names = colnames(this_matrix),
                      ExperimentID = experiment_id,
                      DonorID = paste0(experiment_id, "___", donor_id),
                      Age = rep(NA, ncol(this_matrix)),
                      Sex = rep(NA, ncol(this_matrix)),
                      Ethnicity1 = rep(NA, ncol(this_matrix)),
                      Ethnicity2 = rep(NA, ncol(this_matrix)),
                      EthnicityDetail = rep(NA, ncol(this_matrix)),
                      SubjectType = rep("Alive - Healthy", ncol(this_matrix)),
                      SampleID = paste0(experiment_id, "___", donor_id),
                      SampleType = rep("Surgical resection", ncol(this_matrix)),
                      SampledSiteCondition = rep("Healthy", ncol(this_matrix)),
                      Tissue = rep("Epidermis", ncol(this_matrix)),
                      BiologicalUnit = rep("Cells", ncol(this_matrix)),
                      LibraryPlatform = rep("10x_3'_v2", ncol(this_matrix)),
                      StrandSequence = rep("3'", ncol(this_matrix)),
                      SequencingPlatform = rep("lllumina HiSeq 2500/4000, Illumina NovaSeq 6000 S2", ncol(this_matrix)),
                      ReferenceGenome = rep("GRCh38", ncol(this_matrix)),
                      SampleStatus = rep(NA, ncol(this_matrix)),
                      SampleCultured = rep("No", ncol(this_matrix)),
                      AnatomicalRegionLevel1 = rep("Head", ncol(this_matrix)),
                      AnatomicalRegionLevel2 = rep("Scalp", ncol(this_matrix)),
                      AnatomicalRegionLevel3 = rep(NA, ncol(this_matrix)),
                      OriginalAnnotation = rep(NA, ncol(this_matrix)),
                      Y_X_ratio = colSums(this_matrix[which(rownames(this_matrix) %in% Y_Genes), ]) / colSums(this_matrix[which(rownames(this_matrix) %in% X_Genes), ]))
if(file.exists(loom_path)){
  unlink(loom_path)
}
save_h5(loom_path, t(this_matrix))
write.table(metadata, metadata_path, sep = "\t", row.names = T, col.names = T, quote = F)
######
healthy_loom_path = paste0(output_dir, "/healthy.loom")
healthy_metadata_path = paste0(output_dir, "/healthy_metadata.tsv")
#
loom_path_1 = "D:/Study/KI/Projects/###/MetaStudies/Data/Skin/Human/Cheng_Cho_CellReports_2018/convert_output/UCSFSKIN1_ABD4epi.loom"
metadata_path_1 = "D:/Study/KI/Projects/###/MetaStudies/Data/Skin/Human/Cheng_Cho_CellReports_2018/convert_output/UCSFSKIN1_ABD4epi_metadata.tsv"
loom_path_2 = "D:/Study/KI/Projects/###/MetaStudies/Data/Skin/Human/Cheng_Cho_CellReports_2018/convert_output/UCSFSKIN2_BRST41epi.loom"
metadata_path_2 = "D:/Study/KI/Projects/###/MetaStudies/Data/Skin/Human/Cheng_Cho_CellReports_2018/convert_output/UCSFSKIN2_BRST41epi_metadata.tsv"
loom_path_3 = "D:/Study/KI/Projects/###/MetaStudies/Data/Skin/Human/Cheng_Cho_CellReports_2018/convert_output/UCSFSKIN3_BRST53epi.loom"
metadata_path_3 = "D:/Study/KI/Projects/###/MetaStudies/Data/Skin/Human/Cheng_Cho_CellReports_2018/convert_output/UCSFSKIN3_BRST53epi_metadata.tsv"
loom_path_4 = "D:/Study/KI/Projects/###/MetaStudies/Data/Skin/Human/Cheng_Cho_CellReports_2018/convert_output/UCSFSKIN4_FORE12epi.loom"
metadata_path_4 = "D:/Study/KI/Projects/###/MetaStudies/Data/Skin/Human/Cheng_Cho_CellReports_2018/convert_output/UCSFSKIN4_FORE12epi_metadata.tsv"
loom_path_5 = "D:/Study/KI/Projects/###/MetaStudies/Data/Skin/Human/Cheng_Cho_CellReports_2018/convert_output/UCSFSKIN5_FORE8epi.loom"
metadata_path_5 = "D:/Study/KI/Projects/###/MetaStudies/Data/Skin/Human/Cheng_Cho_CellReports_2018/convert_output/UCSFSKIN5_FORE8epi_metadata.tsv"
loom_path_6 = "D:/Study/KI/Projects/###/MetaStudies/Data/Skin/Human/Cheng_Cho_CellReports_2018/convert_output/UCSFSKIN6_FORE9epi.loom"
metadata_path_6 = "D:/Study/KI/Projects/###/MetaStudies/Data/Skin/Human/Cheng_Cho_CellReports_2018/convert_output/UCSFSKIN6_FORE9epi_metadata.tsv"
loom_path_7 = "D:/Study/KI/Projects/###/MetaStudies/Data/Skin/Human/Cheng_Cho_CellReports_2018/convert_output/UCSFSKIN10_SCALP11epi.loom"
metadata_path_7 = "D:/Study/KI/Projects/###/MetaStudies/Data/Skin/Human/Cheng_Cho_CellReports_2018/convert_output/UCSFSKIN10_SCALP11epi_metadata.tsv"
loom_path_8 = "D:/Study/KI/Projects/###/MetaStudies/Data/Skin/Human/Cheng_Cho_CellReports_2018/convert_output/UCSFSKIN11_SCALP26epi.loom"
metadata_path_8 = "D:/Study/KI/Projects/###/MetaStudies/Data/Skin/Human/Cheng_Cho_CellReports_2018/convert_output/UCSFSKIN11_SCALP26epi_metadata.tsv"
loom_path_9 = "D:/Study/KI/Projects/###/MetaStudies/Data/Skin/Human/Cheng_Cho_CellReports_2018/convert_output/UCSFSKIN12_SCALP32epi.loom"
metadata_path_9 = "D:/Study/KI/Projects/###/MetaStudies/Data/Skin/Human/Cheng_Cho_CellReports_2018/convert_output/UCSFSKIN12_SCALP32epi_metadata.tsv"
py$adata_1 = load_data(loom_path_1, metadata_path_1)
py$adata_2 = load_data(loom_path_2, metadata_path_2)
py$adata_3 = load_data(loom_path_3, metadata_path_3)
py$adata_4 = load_data(loom_path_4, metadata_path_4)
py$adata_5 = load_data(loom_path_5, metadata_path_5)
py$adata_6 = load_data(loom_path_6, metadata_path_6)
py$adata_7 = load_data(loom_path_7, metadata_path_7)
py$adata_8 = load_data(loom_path_8, metadata_path_8)
py$adata_9 = load_data(loom_path_9, metadata_path_9)
######
py_run_string("adata = ad.concat({
              'UCSFSKIN1_ABD4epi': adata_1,
              'UCSFSKIN2_BRST41epi': adata_2, 
              'UCSFSKIN3_BRST53epi': adata_3, 
              'UCSFSKIN4_FORE12epi': adata_4, 
              'UCSFSKIN5_FORE8epi': adata_5, 
              'UCSFSKIN6_FORE9epi': adata_6,
              'UCSFSKIN10_SCALP11epi': adata_7,
              'UCSFSKIN11_SCALP26epi': adata_8,
              'UCSFSKIN12_SCALP32epi': adata_9
              }, axis=0, join='outer', index_unique='---')")
py$adata = unify_genename(py$adata, alias_gene_c_seurat)
py_run_string("adata.var_names_make_unique()")
py_run_string("sc.pp.filter_cells(adata, min_genes=1)")
py_run_string("sc.pp.filter_genes(adata, min_cells=1)")
# Sex Identification
tmp = aggregate(as.numeric(py$adata$obs["Y_X_ratio"][, 1]), by = list(py$adata$obs["DonorID"][, 1]), FUN = function(x){
  mean(x, na.rm = T)
})
mean_by_group = tmp$x
names(mean_by_group) = gsub("Cheng_Cho_CellReports_2018___", "", tmp$Group.1, fixed = T)
pdf(paste0(output_dir, '/DonorID_Sex_barplot.pdf'), width = 10, height = 5)
par(mar=c(11, 4, 4, 4))
barplot(mean_by_group, las=2)
dev.off()
###
metadata = py$adata$obs
metadata[py$adata$obs["DonorID"] == "Cheng_Cho_CellReports_2018___UCSFSKIN1_ABD4epi", "Sex"] = "F"
metadata[py$adata$obs["DonorID"] == "Cheng_Cho_CellReports_2018___UCSFSKIN10_SCALP11epi", "Sex"] = "M"
metadata[py$adata$obs["DonorID"] == "Cheng_Cho_CellReports_2018___UCSFSKIN11_SCALP26epi", "Sex"] = "F"
metadata[py$adata$obs["DonorID"] == "Cheng_Cho_CellReports_2018___UCSFSKIN12_SCALP32epi", "Sex"] = "F"
metadata[py$adata$obs["DonorID"] == "Cheng_Cho_CellReports_2018___UCSFSKIN2_BRST41epi", "Sex"] = "F"
metadata[py$adata$obs["DonorID"] == "Cheng_Cho_CellReports_2018___UCSFSKIN3_BRST53epi", "Sex"] = "F"
metadata[py$adata$obs["DonorID"] == "Cheng_Cho_CellReports_2018___UCSFSKIN4_FORE12epi", "Sex"] = "M"
metadata[py$adata$obs["DonorID"] == "Cheng_Cho_CellReports_2018___UCSFSKIN5_FORE8epi", "Sex"] = "M"
metadata[py$adata$obs["DonorID"] == "Cheng_Cho_CellReports_2018___UCSFSKIN6_FORE9epi", "Sex"] = "M"
if(file.exists(healthy_loom_path)){
  unlink(healthy_loom_path)
}
py$output_path = healthy_loom_path
py_run_string("loompy.create(output_path, adata.X.T, {'Gene': adata.var_names.values}, {'CellID': adata.obs_names.values})")
write.table(as.matrix(metadata), healthy_metadata_path, sep = "\t", row.names = T, col.names = T, quote = F)
