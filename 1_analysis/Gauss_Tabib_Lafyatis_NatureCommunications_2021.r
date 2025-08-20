r_utilities = "/home/haoy/projects/Analysis/code_library/utilities.r"
py_utilities = "/home/haoy/projects/Analysis/code_library/utilities.py"
source(r_utilities)
library(reticulate)
py_run_file(py_utilities, convert = F)
######
######
loom_path = "/home/haoy/projects/MetaStudiesAnalysis/Integration/Data/Tabib_Lafyatis_NatureCommunications_2021/convert_output/GSE138669_HC.loom"
metadata_path = "/home/haoy/projects/MetaStudiesAnalysis/Integration/Data/Tabib_Lafyatis_NatureCommunications_2021/convert_output/GSE138669_HC_metadata.tsv"
######
X_gene_path = "/home/haoy/projects/Analysis/useful_files/X_Gene_Homo_sapiens.GRCh38.103.rds"
Y_gene_path = "/home/haoy/projects/Analysis/useful_files/Y_Gene_Homo_sapiens.GRCh38.103.rds"
alias_gene_c_seurat_path = "/home/haoy/projects/Analysis/useful_files/hgnc_complete_set_c_seurat.rds"
CellCyele_gene_path = "/home/haoy/projects/MetaStudiesAnalysis/Integration/regev_lab_cell_cycle_genes.txt"
######
py$output_dir = "/home/haoy/projects/MetaStudiesAnalysis/Integration/Output/Tabib_Lafyatis_NatureCommunications_2021/ALL/v10_gauss"
py$output_dir_cache = "/home/haoy/projects/MetaStudiesAnalysis/Integration/Output/Tabib_Lafyatis_NatureCommunications_2021/ALL/v10_gauss"
dir.create(py$output_dir, showWarnings = F, recursive = T)
######
######
DEG_leiden_path_list = list()
knn_result_list = list()
knn_result_path_list = list()
cluster_name_list = list()
split_adata_name_list = list()
StudyID_list = list()
hvgene_list = list()
pca_cal_list = list()
mean_std_list = list()
other_cells_list = list()
DEG_leiden_list = list()
path_list = list()
key_list = list("divide" = list())
######
knn_result_path_list[["All"]] = paste0(py$output_dir, "/knn_result.rds")
# Settings
py_run_string("sc.settings.verbosity = 3")
py_run_string("sc.logging.print_header()")
py_run_string("sc.settings.set_figure_params(dpi=300, dpi_save=300, facecolor='white', fontsize=10)")
py_run_string("sc.settings.autoshow = False")
py_run_string("figure_dir = output_dir + '/figures'")
py_run_string("sc.settings.figdir = figure_dir")
dir.create(py$figure_dir, showWarnings = F, recursive = T)
# Clustering parameters
hvg_number = 2000
hvg_number_integrated = hvg_number
pca_dim = 50
pca_dim_integrated = pca_dim
random_state = 0
n_neighbors_clustering = 30
n_neighbors_visualization = 15
leiden_resolutions = c(0.01, 0.02, 0.03, 0.05, 0.08, 0.1, 0.12, 0.15, 0.2, 0.3, 0.5, 0.8, 1, 1.2, 1.5, 2, 3)
cell_cycle_regress_out = F
# DE parameters
top_DE_range = 50
# Other parameters
chunk_size = 500000
max_number = 5000
py_run_string("density_color = 'YlOrRd'")
# Other settings
do_DEG = F
use_altas_querying = F
use_cache = T
######
######
######
summary_list = list()
######
py$adata = load_data_wrap(loom_path, metadata_path, alias_unifying = readRDS(alias_gene_c_seurat_path))
######
py_run_string("adata.var_names_make_unique()")
py_run_string("adata.var['mt'] = adata.var_names.str.startswith('MT-')")
py_run_string("sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)")
######
py_run_string("merge_by = ['DonorID', 'LibraryPlatform']")
# Make a shorter DonorID
py$adata$obs["DonorID_short"] = sapply(as.character(py$adata$obs[["DonorID"]]), function(x){
  unlist(strsplit(x, split = "___", fixed = T))[2]
})
py_run_string("factors = ['DonorID_short', 'Age', 'Sex', 'Ethnicity1', 'LibraryPlatform']") # We use shorter DonorIDs in each separate analysis
######
py_run_string("annotation = ['OriginalAnnotation']")
py_run_string("adata.obs['merge_by'] = adata.obs[merge_by].agg('___'.join, axis=1).values") # Only use for record.
######
py_run_string("qc_metrics = ['n_genes_by_counts', 'total_counts', 'pct_counts_mt']")
py_run_string("sc.pl.violin(adata, qc_metrics, jitter=0.4, multi_panel=True, save='_before_QC.pdf')")
qc_state = "before_QC"
summary_list[[qc_state]] = list()
summary_list[[qc_state]][["dim"]] = c(py$adata$n_obs, py$adata$n_vars)
py_run_string("tmp_expressed_gene = adata.obs['n_genes_by_counts'].values.copy()")
py_run_string("tmp_expression = adata.obs['total_counts'].values.copy()")
summary_list[[qc_state]][["expressed_gene"]] = summary(py$tmp_expressed_gene)
summary_list[[qc_state]][["expression"]] = summary(py$tmp_expression)
py_run_string("del tmp_expressed_gene")
py_run_string("del tmp_expression")
for(ii in unlist(py$factors)){
  py_run_string(paste0("tmp_factor = adata.obs['", ii, "'].values.astype(str).copy()"))
  summary_list[[qc_state]][[ii]] = table(py$tmp_factor)
  if(length(summary_list[[qc_state]][[ii]]) == 0)(
    stop(Paste0("Factor ", ii, " consist with only NA!"))
  )
  py_run_string("del tmp_factor")
}
# filtering
# Didn't find the specific filtering parameters, we filter it using a normal one as other 10X datasets
filter =
  py$adata$obs["pct_counts_mt"] < 15 & # Use the same parameter as Dunlap_Rao_JCIInsight_2022
  py$adata$obs["n_genes_by_counts"] >= 200 & # Use the same parameter as Dunlap_Rao_JCIInsight_2022/Tabib_Lafyatis_JournalofInvestigativeDermatology_2017
  py$adata$obs["n_genes_by_counts"] <= 3500 # Use the same parameter as Dunlap_Rao_JCIInsight_2022
print(sum(filter))
py$adata = py$adata[filter]
py_run_string("sc.pp.filter_genes(adata, min_cells=1)")
py_run_string("sc.pl.violin(adata, qc_metrics, jitter=0.4, multi_panel=True, save='_after_QC.pdf')")
qc_state = "after_QC"
summary_list[[qc_state]] = list()
summary_list[[qc_state]][["dim"]] = c(py$adata$n_obs, py$adata$n_vars)
py_run_string("tmp_expressed_gene = adata.obs['n_genes_by_counts'].values.copy()")
py_run_string("tmp_expression = adata.obs['total_counts'].values.copy()")
summary_list[[qc_state]][["expressed_gene"]] = summary(py$tmp_expressed_gene)
summary_list[[qc_state]][["expression"]] = summary(py$tmp_expression)
py_run_string("del tmp_expressed_gene")
py_run_string("del tmp_expression")
for(ii in unlist(py$factors)){
  py_run_string(paste0("tmp_factor = adata.obs['", ii, "'].values.astype('str').copy()"))
  summary_list[[qc_state]][[ii]] = table(py$tmp_factor)
  if(length(summary_list[[qc_state]][[ii]]) == 0)(
    stop(Paste0("Factor ", ii, " consist with only NA!"))
  )
  py_run_string("del tmp_factor")
}
###
py$cell_sorting = py$adata$obs_names$values
saveRDS(py$cell_sorting, paste0(py$output_dir, "/cell_sorting.rds"))
print(summary_list)
######
# QC #
######
comparison_mat = matrix(ncol = 5, nrow = py$adata$n_obs,
                        dimnames = list(py$adata$obs_names$values,
                                        c("nUMI", "nFeatures", "MitoPrecentage", "XPrecentage", "YPrecentage")))
py$X_Genes = readRDS(X_gene_path)
py$Y_Genes = readRDS(Y_gene_path)
py_run_string("M_Genes = adata.var_names.str.startswith('MT-')")
py_run_string("tmp_numi = adata.X.sum(1)")
py_run_string("tmp_nfeatures = (adata.X > 0).sum(1)")
py_run_string("tmp_mito_pct = adata[:, M_Genes].X.sum(1) / tmp_numi")
py_run_string("tmp_x_pct = adata[:, np.intersect1d(X_Genes, adata.var_names.values)].X.sum(1) / tmp_numi")
py_run_string("tmp_y_pct = adata[:, np.intersect1d(Y_Genes, adata.var_names.values)].X.sum(1) / tmp_numi")
comparison_mat[, "nUMI"] = py$tmp_numi
comparison_mat[, "nFeatures"] = py$tmp_nfeatures
comparison_mat[, "MitoPrecentage"] = py$tmp_mito_pct
comparison_mat[, "XPrecentage"] = py$tmp_x_pct
comparison_mat[, "YPrecentage"] = py$tmp_y_pct
py_run_string("del tmp_numi")
py_run_string("del tmp_nfeatures")
py_run_string("del tmp_mito_pct")
py_run_string("del tmp_x_pct")
py_run_string("del tmp_y_pct")
comparison_df = cbind(as.data.frame(comparison_mat), cbind(as.character(py$adata$obs$DonorID_short), as.character(py$adata$obs$ExperimentID)))
colnames(comparison_df) = c(colnames(comparison_mat), "DonorID", "ExperimentID")
p = StackedVlnPlot(comparison_df, "DonorID", colnames(comparison_mat), "ExperimentID", ymaxs = c(200000, max(comparison_df$nFeatures), 0.2, 0.1, 0.001))
ggsave(paste0(py$figure_dir, "/qc.pdf"), p, width = 1.5 * length(unique(comparison_df[, "DonorID"])) + 3, height = 12, limitsize = FALSE)
saveRDS(comparison_df, paste0(py$output_dir, "/comparison_df.rds"))
write.table(comparison_df, paste0(py$output_dir, "/comparison_df.tsv"), sep = "\t", row.names = T, col.names = T, quote = F)
######
cell_number = py$adata$n_obs
############
# Cell cycle
############
# Following https://nbviewer.org/github/theislab/scanpy_usage/blob/master/180209_cell_cycle/cell_cycle.ipynb
CellCyele_gene = gsub(" ", "", toupper(as.matrix(read.table(CellCyele_gene_path))[, 1]))
CellCyele_gene[CellCyele_gene == "MLF1IP"] = "CENPU" # MLF1IP is coded by CENPU
py$s_genes = CellCyele_gene[1:43]
py$g2m_genes = CellCyele_gene[44:length(CellCyele_gene)]
py_run_string(paste0("sc.tl.score_genes_cell_cycle(adata, s_genes=s_genes, g2m_genes=g2m_genes, random_state=", random_state, ")"))
############
# Clustering
############
adata_imbalance_clustering_1_path = paste0(py$output_dir_cache, "/adata_imbalance_clustering_1.h5ad")
adata_imbalance_hvg_clustering_1_path = paste0(py$output_dir_cache, "/adata_imbalance_hvg_clustering_1.h5ad")
clustering_keys_1_path = paste0(py$output_dir_cache, "/clustering_keys_1.rds")
if(file.exists(adata_imbalance_clustering_1_path) && file.exists(adata_imbalance_hvg_clustering_1_path)){
  py_run_string(paste0("adata_imbalance = ad.read_h5ad('", adata_imbalance_clustering_1_path, "')"))
  py_run_string(paste0("adata_imbalance_hvg = ad.read_h5ad('", adata_imbalance_hvg_clustering_1_path, "')"))
  py_run_string("adata_imbalance.uns['log1p']['base'] = None")
  clustering_keys_1 = readRDS(clustering_keys_1_path)
}else{
  ############
  # Imbalanced
  ############
  # Clustering
  pca_cal_list[["All"]] = PCA_calculation_python(adata = py$adata, random_state = random_state, merge_by = unlist(py$merge_by), hvg_number = hvg_number, pca_dim = pca_dim, return_all = T)
  py$adata_imbalance_hvg = pca_cal_list[["All"]][["hvg_matrix"]]
  py$adata_imbalance = pca_cal_list[["All"]][["complete_matrix"]]
  py_run_string(paste0("sc.pl.pca_variance_ratio(adata_imbalance_hvg, n_pcs=", pca_dim, ", log=True, save='_imbalance.pdf')"))
  py_run_string(paste0("sc.pp.neighbors(adata_imbalance_hvg, n_neighbors=", n_neighbors_clustering, ", n_pcs=", pca_dim, ", use_rep='PCA_use', random_state=", random_state, ", key_added='clustering')"))
  # Leiden Clustering
  leiden_clustering_keys_1 = c()
  for(ii in leiden_resolutions){
    this_key = paste0("leiden_1_", gsub(pattern = ".", replacement = "", sprintf("%.2f", ii), fixed = T))
    py_run_string(paste0("sc.tl.leiden(adata_imbalance_hvg, resolution=", ii, ", key_added='", this_key, "', neighbors_key='clustering')"))
    py_run_string(paste0("adata_imbalance.obs['", this_key, "'] = adata_imbalance_hvg.obs['", this_key, "']"))
    leiden_clustering_keys_1 = c(leiden_clustering_keys_1, this_key)
  }
  # Visualization
  py_run_string(paste0("sc.pp.neighbors(adata_imbalance_hvg, n_neighbors=", n_neighbors_visualization, ", n_pcs=", pca_dim, ", use_rep='PCA_use', random_state=", random_state, ", key_added='visualization')"))
  py_run_string(paste0("sc.tl.umap(adata_imbalance_hvg, random_state=", random_state, ", neighbors_key='visualization')"))
  py_run_string("adata_imbalance.obsm['X_umap_imbalance'] = adata_imbalance_hvg.obsm['X_umap']")
  clustering_keys_1 = c(leiden_clustering_keys_1)
  ####################
  # Cache - Clustering
  ####################
  if(use_cache){
    py_run_string(paste0("adata_imbalance.write('", adata_imbalance_clustering_1_path, "')"))
    py_run_string(paste0("adata_imbalance_hvg.write('", adata_imbalance_hvg_clustering_1_path, "')"))
    saveRDS(clustering_keys_1, clustering_keys_1_path)
  }
}
py_run_string(paste0("sc.pl.pca_loadings(adata_imbalance_hvg, components = '", paste(seq(pca_dim), collapse = ", "), "', save='.pdf')"))
py_run_string(paste0("sc.pl.embedding(adata_imbalance, 'X_umap_imbalance', color=factors + ['phase'], legend_loc='right margin', legend_fontsize=6, size=", get_dot_size(cell_number), ", frameon=False, ncols=1, save='_factors.pdf')"))
py_run_string(paste0("sc.pl.embedding(adata_imbalance, 'X_umap_imbalance', color=['", paste(clustering_keys_1, collapse = "', '"), "'] + annotation, legend_loc='on data', legend_fontsize=6, size=", get_dot_size(cell_number), ", frameon=False, ncols=7, save='_clustering.pdf')"))
for(ii in c(clustering_keys_1, unlist(py$factors))){
  py_run_string(paste0("clustering_plot(adata_imbalance, '", ii, "', basis='X_umap_imbalance', size=", get_dot_size(cell_number), ", colorbar_loc=None, ncols=7, save='_detail_", ii, ".pdf')"))
}
############
# Signatures
############
markers_list = list()
markers_list[["Keratinocyte"]] = c("KRT14", "KRT15", "KRT5", "KRT1", "KRT10", "KRTDAP",
                                   "KRT2", "KRT6A", "KRT9", "KRT18", "KRT23", "KRT35",
                                   "MGST1", "NRARP", "POSTN", "TK1", "CENPU", "MKI67")
markers_list[["Fibroblast_MuralCell"]] = c("FBLN1", "COL1A1", "COL1A2", "COL3A1", "COL6A3", "DCN", "COCH", "SFRP2",
                                           "ACTA2", "RGS5", "RERGL", "TAGLN")
markers_list[["NeuralCrestderivedCell"]] = c("TYR", "DCT", "TYRP1", "MLANA",
                                             "PMEL", "SOX10", "SOX2", "MPZ")
markers_list[["EndothelialCell"]] = c("SELE", "SELP", "VWF", "PECAM1",
                                      "LYVE1", "CCL21", "TFF3", "MMRN1")
markers_list[["ImmuneCell"]] = c("TPSB2", "TPSAB1",
                                 "IL1B", "LYZ", "CD207", "CD1A",
                                 "CD3G", "CD3E", "CD3D", "TRBC1", "TRBC2", "TIGIT")
markers_list[["Plasma_Erythrocyte"]] = c("HBB", "HBA1", "HBM", "JCHAIN")
######
py$Keratinocyte_markers = markers_list[["Keratinocyte"]][markers_list[["Keratinocyte"]] %in% py$adata_imbalance$var_names$values]
py$Fibroblast_MuralCell_markers = markers_list[["Fibroblast_MuralCell"]][markers_list[["Fibroblast_MuralCell"]] %in% py$adata_imbalance$var_names$values]
py$NeuralCrestderivedCell_markers = markers_list[["NeuralCrestderivedCell"]][markers_list[["NeuralCrestderivedCell"]] %in% py$adata_imbalance$var_names$values]
py$EndothelialCell_markers = markers_list[["EndothelialCell"]][markers_list[["EndothelialCell"]] %in% py$adata_imbalance$var_names$values]
py$ImmuneCell_markers = markers_list[["ImmuneCell"]][markers_list[["ImmuneCell"]] %in% py$adata_imbalance$var_names$values]
py$Plasma_Erythrocyte_markers = markers_list[["Plasma_Erythrocyte"]][markers_list[["Plasma_Erythrocyte"]] %in% py$adata_imbalance$var_names$values]
######
py_run_string(paste0("sc.pl.embedding(adata_imbalance, 'X_umap_imbalance', color=Keratinocyte_markers, layer='normalized', legend_loc='on data', legend_fontsize=6, size=", get_dot_size(cell_number), ", ncols=6, frameon=False, save='_Signatures_Keratinocyte.pdf')"))
py_run_string(paste0("sc.pl.embedding(adata_imbalance, 'X_umap_imbalance', color=Fibroblast_MuralCell_markers, layer='normalized', legend_loc='on data', legend_fontsize=6, size=", get_dot_size(cell_number), ", ncols=4, frameon=False, save='_Signatures_Fibroblast_MuralCell.pdf')"))
py_run_string(paste0("sc.pl.embedding(adata_imbalance, 'X_umap_imbalance', color=NeuralCrestderivedCell_markers, layer='normalized', legend_loc='on data', legend_fontsize=6, size=", get_dot_size(cell_number), ", ncols=4, frameon=False, save='_Signatures_NeuralCrestderivedCell.pdf')"))
py_run_string(paste0("sc.pl.embedding(adata_imbalance, 'X_umap_imbalance', color=EndothelialCell_markers, layer='normalized', legend_loc='on data', legend_fontsize=6, size=", get_dot_size(cell_number), ", ncols=4, frameon=False, save='_Signatures_EndothelialCell.pdf')"))
py_run_string(paste0("sc.pl.embedding(adata_imbalance, 'X_umap_imbalance', color=ImmuneCell_markers, layer='normalized', legend_loc='on data', legend_fontsize=6, size=", get_dot_size(cell_number), ", ncols=6, frameon=False, save='_Signatures_ImmuneCell.pdf')"))
py_run_string(paste0("sc.pl.embedding(adata_imbalance, 'X_umap_imbalance', color=Plasma_Erythrocyte_markers, layer='normalized', legend_loc='on data', legend_fontsize=6, size=", get_dot_size(cell_number), ", ncols=4, frameon=False, save='_Signatures_Plasma_Erythrocyte.pdf')"))
######
##############
# DEG Analysis
##############
if(do_DEG){
  DEG_leiden_path_list[["All"]] = list()
  for(ii in c("leiden_1_015")){
    this_path = paste0(py$output_dir, "/DEG_", ii, ".rds")
    DEG_leiden_path_list[["All"]][[ii]] = this_path
    if(file.exists(this_path)){
      DEG_leiden_list[["All"]][[ii]] = readRDS(this_path)
    }else{
      DEG_leiden_list[["All"]][[ii]] = DEG_Analysis(py$adata_imbalance, ii, py$output_dir, prefix = paste0(ii, "_"), top_DE_range = top_DE_range)
      saveRDS(DEG_leiden_list[["All"]][[ii]], this_path)
    }
  }
}
#####################
# Cell identification
#####################
adata_imbalance_annotated_1_path = paste0(py$output_dir, "/adata_imbalance_annotated_1.h5ad")
if(file.exists(adata_imbalance_annotated_1_path)){
  py_run_string(paste0("adata_imbalance = ad.read_h5ad('", adata_imbalance_annotated_1_path, "')"))
  py_run_string("annotation = annotation if np.isin('SeparateAnnotation_1', annotation) else annotation + ['SeparateAnnotation_1']")
  ###
  py_run_string("tmp_leiden_all = adata_imbalance.obs['leiden_1'][cell_sorting].astype(str).values")
  leiden_all = py$tmp_leiden_all
  py_run_string("del tmp_leiden_all")
  py_run_string("tmp_SeparateAnnotation = adata_imbalance.obs['SeparateAnnotation_1'].astype(str).values")
  summary_list[[qc_state]][["SeparateAnnotation_1"]] = list("AllCells" = table(py$tmp_SeparateAnnotation))
  summary_list[[qc_state]][["leiden_1"]] = list("AllCells" = table(leiden_all))
  py_run_string("del tmp_SeparateAnnotation")
}else{
  py_run_string("tmp_leiden_all = adata_imbalance.obs['leiden_1_015'][cell_sorting].astype(str).values")
  leiden_all = py$tmp_leiden_all
  py_run_string("del tmp_leiden_all")
  cell_annotation = rep(NA, py$adata_imbalance$n_obs)
  names(cell_annotation) = py$cell_sorting
  leiden_1_015 = as.character(py$adata_imbalance$obs[['leiden_1_015']])
  ######
  #
  cell_annotation[leiden_1_015 == "0"] = "Fibroblast"
  cell_annotation[leiden_1_015 == "1"] = "Keratinocyte"
  cell_annotation[leiden_1_015 == "2"] = "Mural Cell"
  cell_annotation[leiden_1_015 == "3"] = "Vascular Endothelial Cell"
  cell_annotation[leiden_1_015 == "4"] = "Lymphocyte"
  cell_annotation[leiden_1_015 == "5"] = "Keratinocyte"
  cell_annotation[leiden_1_015 == "6"] = "Myeloid Cell"
  cell_annotation[leiden_1_015 == "7"] = "Mural Cell"
  cell_annotation[leiden_1_015 == "8"] = "Schwann Cell"
  cell_annotation[leiden_1_015 == "9"] = "Keratinocyte"
  cell_annotation[leiden_1_015 == "10"] = "Mast Cell"
  cell_annotation[leiden_1_015 == "11"] = "Keratinocyte"
  cell_annotation[leiden_1_015 == "12"] = "Melanocyte"
  cell_annotation[leiden_1_015 == "13"] = "Lymphatic Endothelial Cell"
  #
  print(sum(is.na(cell_annotation)))
  py$tmp_clustering_result = cell_annotation
  ###
  if(use_altas_querying){
    py_run_string("tmp_clustering_result = adata_imbalance.obs['Spread_mapping_representative_SeparateAnnotation_1_leiden_1_015'].astype(str).values")
  }
  ###
  py_run_string("adata_imbalance.obs['leiden_1'] = adata_imbalance.obs['leiden_1_015']")
  py_run_string("adata_imbalance.obs['SeparateAnnotation_1'] = tmp_clustering_result")
  py_run_string("del tmp_clustering_result")
  ###
  py_run_string("tmp_SeparateAnnotation = adata_imbalance.obs['SeparateAnnotation_1'].astype(str).values")
  summary_list[[qc_state]][["SeparateAnnotation_1"]] = list("AllCells" = table(py$tmp_SeparateAnnotation))
  summary_list[[qc_state]][["leiden_1"]] = list("AllCells" = table(leiden_all))
  py_run_string("del tmp_SeparateAnnotation")
  ######
  for(ii in unlist(py$annotation)){
    py_run_string(paste0("tmp_this_annotation = adata_imbalance.obs['", ii, "'].astype(str).values"))
    summary_list[[qc_state]][[ii]] = table(py$tmp_this_annotation)
    py_run_string("del tmp_this_annotation")
  }
  py_run_string("annotation = annotation if np.isin('SeparateAnnotation_1', annotation) else annotation + ['SeparateAnnotation_1']")
  ######################
  # Cache - Querying - 2
  ######################
  if(use_cache){
    py_run_string(paste0("adata_imbalance.write('", adata_imbalance_annotated_1_path, "')"))
  }
}
###
# Ploting will tickle .strings_to_categoricals() automatically, cause error when accessing adata attributes in r
py_run_string(paste0("sc.pl.embedding(adata_imbalance, 'X_umap_imbalance', color=['", paste(clustering_keys_1, collapse = "', '"), "'] + annotation, legend_loc='on data', legend_fontsize=6, size=", get_dot_size(cell_number), ", frameon=False, ncols=7, save='_annotation.pdf')"))
for(ii in c("phase", unlist(py$annotation))){
  py_run_string(paste0("clustering_plot(adata_imbalance, '", ii, "', basis='X_umap_imbalance', size=", get_dot_size(cell_number), ", colorbar_loc=None, ncols=5, save='_detail_", ii, ".pdf')"))
}
density_visualization("adata_imbalance", "umap_imbalance", "Sex", "imbalance.pdf", "density_color", "5", as.character(get_dot_size(cell_number)))
density_visualization("adata_imbalance", "umap_imbalance", "SeparateAnnotation_1", "imbalance.pdf", "density_color", "5", as.character(get_dot_size(cell_number)))
py_run_string("factor_distribution(adata_imbalance, factors + ['phase'], 'leiden_1', path=figure_dir + '/Distribution_clustering.pdf')")
py_run_string("factor_distribution_per_method(adata_imbalance, factors + ['phase'], 'leiden_1', precent='y', path=figure_dir + '/Distribution_clustering_divided.pdf')")
py_run_string("factor_distribution(adata_imbalance, factors + ['phase'], 'SeparateAnnotation_1', path=figure_dir + '/Distribution_SeparateAnnotation.pdf')")
py_run_string("factor_distribution_per_method(adata_imbalance, factors + ['phase'], 'SeparateAnnotation_1', precent='y', path=figure_dir + '/Distribution_SeparateAnnotation_divided.pdf')")
saveRDS(summary_list, paste0(py$output_dir, "/summary_list.rds"))
sink(paste0(py$output_dir, "/summary_list.txt"))
print(summary_list)
sink()
neighbor_list = list()
for(ii in c("clustering", "visualization")){
  this_neighbor_uns = py$adata_imbalance_hvg$uns[[ii]]
  neighbor_list[[ii]] = list("uns" = this_neighbor_uns,
                             "connectivities" = py$adata_imbalance_hvg$obsp[[this_neighbor_uns[["connectivities_key"]]]],
                             "distances" = py$adata_imbalance_hvg$obsp[[this_neighbor_uns[["distances_key"]]]])
}
saveRDS(neighbor_list, paste0(py$output_dir, "/all_neighbor_list.rds"))
############
py$loom_path_all = paste0(py$output_dir, "/all.loom")
py$loom_path_all_hvg = paste0(py$output_dir, "/all_hvg.loom")
metadata_path_all = paste0(py$output_dir, "/metadata_all.tsv")
if(file.exists(py$loom_path_all)){
  unlink(py$loom_path_all)
}
if(file.exists(py$loom_path_all_hvg)){
  unlink(py$loom_path_all_hvg)
}
py_run_string("loompy.create(loom_path_all, adata_imbalance.layers['raw_counts'].T, {'Gene': adata_imbalance.var_names.values}, {'CellID': adata_imbalance.obs_names.values})")
py_run_string("loompy.create(loom_path_all_hvg, adata_imbalance_hvg.layers['raw_counts'].T, {'Gene': adata_imbalance_hvg.var_names.values}, {'CellID': adata_imbalance_hvg.obs_names.values})")
py_run_string("tmp_metadata_all = adata_imbalance.obs.astype(str)")
py_run_string("del loom_path_all")
py_run_string("del loom_path_all_hvg")
write.table(py$tmp_metadata_all, metadata_path_all, sep = "\t", row.names = T, col.names = T, quote = F)
saveRDS(py$tmp_metadata_all, paste0(py$output_dir, "/metadata.rds"))
py_run_string("del tmp_metadata_all")
#############################
# End of Cell Type Clustering
#############################
py_run_string("annotation_1 = copy.deepcopy(annotation)")
annotation_1_path = paste0(py$output_dir, "/annotation_1.rds")
if(use_cache){
  saveRDS(py$annotation_1, annotation_1_path)
}
#####################
# Secondly Clustering
#####################
