r_utilities = "D:/Study/KI/Projects/###/Analysis/code_library/utilities.r"
source(r_utilities)
######
X_gene_path = "D:/Study/KI/Projects/###/Analysis/useful_files/X_Gene_Homo_sapiens.GRCh38.103.rds"
Y_gene_path = "D:/Study/KI/Projects/###/Analysis/useful_files/Y_Gene_Homo_sapiens.GRCh38.103.rds"
alias_gene_c_seurat_path = "D:/Study/KI/Projects/###/Analysis/useful_files/hgnc_complete_set_c_seurat.rds"
alias_gene_c_seurat = readRDS(alias_gene_c_seurat_path)
######
# filtered
mtx_gz_path_filtered = "D:/Study/KI/Projects/###/MetaStudies/Data/Skin/Human/Boldo_Lyko_CommunicationsBiology_2020/GSE130973_matrix_filtered.mtx.gz"
barcodes_gz_path_filtered = "D:/Study/KI/Projects/###/MetaStudies/Data/Skin/Human/Boldo_Lyko_CommunicationsBiology_2020/GSE130973_barcodes_filtered.tsv.gz"
genes_gz_path_filtered = "D:/Study/KI/Projects/###/MetaStudies/Data/Skin/Human/Boldo_Lyko_CommunicationsBiology_2020/GSE130973_genes_filtered.tsv.gz"
# raw
#mtx_gz_path_raw = "D:/Study/KI/Projects/###/MetaStudies/Data/Skin/Human/Boldo_Lyko_CommunicationsBiology_2020/GSE130973_matrix_raw.mtx.gz"
#barcodes_gz_path_raw = "D:/Study/KI/Projects/###/MetaStudies/Data/Skin/Human/Boldo_Lyko_CommunicationsBiology_2020/GSE130973_barcodes_raw.tsv.gz"
#genes_gz_path_raw = "D:/Study/KI/Projects/###/MetaStudies/Data/Skin/Human/Boldo_Lyko_CommunicationsBiology_2020/GSE130973_genes_raw.tsv.gz"
# rds
rds_gz_path = "D:/Study/KI/Projects/###/MetaStudies/Data/Skin/Human/Boldo_Lyko_CommunicationsBiology_2020/GSE130973_seurat_analysis_lyko.rds.gz"
######
### Output
output_dir = "D:/Study/KI/Projects/###/MetaStudies/Data/Skin/Human/Boldo_Lyko_CommunicationsBiology_2020/convert_output"
dir.create(output_dir, showWarnings = F, recursive = T)
loom_path = paste0(output_dir, "/GSE130973.loom")
metadata_path = paste0(output_dir, "/GSE130973_metadata.tsv")
###
expression_matrix_filtered = ReadMtx(mtx = mtx_gz_path_filtered, cells = barcodes_gz_path_filtered, features = genes_gz_path_filtered) # 32738 genes X 16062 cells
#expression_matrix_raw = ReadMtx(mtx = gzfile(mtx_gz_path_raw), cells = gzfile(barcodes_gz_path_raw), features = gzfile(genes_gz_path_raw)) # 32738 genes X 3686400 cells
this_obj_0 = readRDS(gzcon(gzfile(rds_gz_path))) # 15457 cells, same amount as its paper (after filtering as "To remove possible cell doublets, we filtered out cells with more than 7500 expressed genes, and to remove potential apoptotic cells we discarded cells with more than 5% mitochondrial reads. The application of these filters resulted in a final dataset of 15,457 single-cell transcriptomes.")
this_obj = CreateSeuratObject(expression_matrix_filtered[, colnames(expression_matrix_filtered) %in% colnames(this_obj_0)], min.cells = 1)
this_obj = unify_genename_seurat(this_obj, alias_gene_c_seurat, verbose = T)
######
# metadata
#View(this_obj_0@meta.data)
cell_number = ncol(this_obj)
experiment_id = rep("Boldo_Lyko_CommunicationsBiology_2020", cell_number)
donor_id = as.character(this_obj_0@meta.data$subj)
subject_type = rep(NA, cell_number)
sample_id = paste0(experiment_id, "___", donor_id)
sample_type = rep("Biopsy", cell_number)
sampled_site_condition = rep("Healthy", cell_number)
tissue = rep(NA, cell_number)
biological_unit = rep("Cells", cell_number)
library_platform = rep("10x_3'_v2", cell_number)
strand_sequence = rep("3'", cell_number)
sequencing_platform = rep("Illumina HiSeq 4000", cell_number)
reference_genome = rep(NA, cell_number)
sample_status = rep("Fresh", cell_number) # immediately after resection, no longer than 1 h
sample_cultured = rep("No", cell_number)
anatomical_region_level1 = rep("Torso", cell_number)
anatomical_region_level2 = rep("Groins", cell_number)
anatomical_region_level3 = rep(NA, cell_number)
ethnicity_1 = rep(NA, cell_number)
ethnicity_2 = rep(NA, cell_number)
ethnicity_detail = rep(NA, cell_number)
age = rep(NA, cell_number)
sex = rep("M", cell_number)
###
#Gene expression signatures used for the definition of cell populations were: ACTA2, RGS5 and PDGFRB (pericytes, clusters #8 and #10)59; KRT5, KRT14, TP63, ITGB1, and ITGA6 (epidermal stem cells and other undifferentiated progenitors, clusters #7 and #15); KRT1, KRT10, SBSN, and KRTDAP (differentiated keratinocytes, cluster #5); PDGFRA, LUM, DCN, VIM, and COL1A2 (fibroblasts, clusters #1, #2, #3 and #9); AIF1, LYZ, HLA-DRA, CD68, and ITGAX (macrophages and dendritic cells (DC), clusters #0, #13 and #16)60; CD3D, CD3G, CD3E, and LCK (T cells, cluster #6)61; SELE, CLDN5, VWF, and CDH5 (vascular endothelial cells, cluster #4)62; PROX1, CLDN5, and LYVE1 (lymphatic endothelial cells, cluster #12)63; HBA1, HBA2, and HBB (erythrocytes, cluster #11)64 and PMEL, MLANA, TYRP1, and DCT (melanocytes, cluster #14)65.
original_annotation = as.character(this_obj_0@meta.data$integrated_snn_res.0.4)
original_annotation[original_annotation == "0"] = "macrophages and dendritic cells"
original_annotation[original_annotation == "1"] = "fibroblasts"
original_annotation[original_annotation == "2"] = "fibroblasts"
original_annotation[original_annotation == "3"] = "fibroblasts"
original_annotation[original_annotation == "4"] = "vascular endothelial cells"
original_annotation[original_annotation == "5"] = "differentiated keratinocytes"
original_annotation[original_annotation == "6"] = "T cells"
original_annotation[original_annotation == "7"] = "epidermal stem cells and other undifferentiated progenitors"
original_annotation[original_annotation == "8"] = "pericytes"
original_annotation[original_annotation == "9"] = "fibroblasts"
original_annotation[original_annotation == "10"] = "pericytes"
original_annotation[original_annotation == "11"] = "erythrocytes"
original_annotation[original_annotation == "12"] = "lymphatic endothelial cells"
original_annotation[original_annotation == "13"] = "macrophages and dendritic cells"
original_annotation[original_annotation == "14"] = "melanocytes"
original_annotation[original_annotation == "15"] = "epidermal stem cells and other undifferentiated progenitors"
original_annotation[original_annotation == "16"] = "macrophages and dendritic cells"
###
this_donor = "S1"
age[donor_id == this_donor] = 25
this_donor = "S2"
age[donor_id == this_donor] = 27
this_donor = "S3"
age[donor_id == this_donor] = 53
this_donor = "S4"
age[donor_id == this_donor] = 70
this_donor = "S5"
age[donor_id == this_donor] = 69
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
