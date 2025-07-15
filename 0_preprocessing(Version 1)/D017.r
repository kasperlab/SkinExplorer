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
py$loom_path = "D:/Study/KI/Projects/###/MetaStudies/Data/Skin/Human/Gao_Hu_CellDeathDisease_2021/GSE162183_Raw_gene_counts_matrix_LoomFile.loom"
###
output_dir = "D:/Study/KI/Projects/###/MetaStudies/Data/Skin/Human/Gao_Hu_CellDeathDisease_2021/convert_output"
dir.create(output_dir, showWarnings = F, recursive = T)
loom_path = paste0(output_dir, "/GSE162183_HC.loom")
metadata_path = paste0(output_dir, "/GSE162183_HC_metadata.tsv")
###
# HC
py_run_string("adata = sc.read_loom(loom_path)")
py$adata = unify_genename(py$adata, alias_gene_c_seurat, verbose = T)
py_run_string("adata.var_names_make_unique()")
py_run_string("sc.pp.filter_cells(adata, min_genes=1)")
py_run_string("sc.pp.filter_genes(adata, min_cells=1)")
this_obj = py$adata$obs
######
# There's UMAP and tSNE coordinates in the loom file, but the figure cannot map into the figure on its original paper...
py_run_string("sc.settings.verbosity = 3")
py_run_string("sc.logging.print_header()")
py_run_string("sc.settings.set_figure_params(dpi=300, dpi_save=300, facecolor='white', fontsize=10)")
py_run_string("sc.settings.autoshow = False")
py$output_dir = output_dir
py_run_string("figure_dir = output_dir + '/figures'")
py_run_string("sc.settings.figdir = figure_dir")
dir.create(py$figure_dir, showWarnings = F, recursive = T)
umap_coordinate = matrix(NA, nrow = py$adata$n_obs, ncol = 2, dimnames = list(py$adata$obs_names$values, c()))
umap_coordinate[, 1] = as.numeric(as.character(py$adata$obs[["_UMAP1"]]))
umap_coordinate[, 2] = as.numeric(as.character(py$adata$obs[["_UMAP2"]]))
py$tmp_umap_coordinate = umap_coordinate
py_run_string("adata.obsm['X_umap'] = tmp_umap_coordinate")
py_run_string("del tmp_umap_coordinate")
py_run_string("sc.pl.embedding(adata, 'X_umap', color=['Patient'], size=8, frameon=False, ncols=1, save='_Patient.pdf')")
######
# metadata
cell_number = py$adata$n_obs
experiment_id = rep("Gao_Hu_CellDeathDisease_2021", cell_number)
donor_id = as.character(py$adata$obs[["Patient"]])
subject_type = rep("Alive - Healthy", cell_number)
sample_id = paste0(experiment_id, "___", donor_id)
sample_type = rep("Surgical resection", cell_number)
sampled_site_condition = rep("Healthy", cell_number)
tissue = rep(NA, cell_number)
biological_unit = rep("Cells", cell_number)
library_platform = rep("10x_3'_v3", cell_number)
strand_sequence = rep("3'", cell_number)
sequencing_platform = rep("Illumina NovaSeq 6000", cell_number)
reference_genome = rep("GRCh38", cell_number)
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
original_annotation = rep(NA, cell_number)
######
this_donor = "Ctrl1"
age[donor_id == this_donor] = 32
sex[donor_id == this_donor] = "F"
anatomical_region_level1[donor_id == this_donor] = "Torso"
anatomical_region_level2[donor_id == this_donor] = "Abdominal"
anatomical_region_level3[donor_id == this_donor] = NA
this_donor = "Ctrl2"
age[donor_id == this_donor] = 23
sex[donor_id == this_donor] = "M"
anatomical_region_level1[donor_id == this_donor] = "Extremities"
anatomical_region_level2[donor_id == this_donor] = "Leg"
anatomical_region_level3[donor_id == this_donor] = NA
this_donor = "Ctrl3"
age[donor_id == this_donor] = 47
sex[donor_id == this_donor] = "M"
anatomical_region_level1[donor_id == this_donor] = "Extremities"
anatomical_region_level2[donor_id == this_donor] = "Leg"
anatomical_region_level3[donor_id == this_donor] = NA
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
py$HC_mask = donor_id %in% c("Ctrl1", "Ctrl2", "Ctrl3")
HC_obj = this_obj[py$HC_mask, ]
py_run_string("adata_HC = adata[HC_mask, :]")
######
py$output_path = loom_path
py_run_string("loompy.create(output_path, adata_HC.X.T, {'Gene': adata_HC.var_names.values}, {'CellID': adata_HC.obs_names.values})")
write.table(as.matrix(HC_obj), metadata_path, sep = "\t", row.names = T, col.names = T, quote = F)



