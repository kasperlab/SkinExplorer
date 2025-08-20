r_utilities = "/home/haoy/projects/Analysis/code_library/utilities.r"
py_utilities = "/home/haoy/projects/Analysis/code_library/utilities.py"
py_visualization = "/home/haoy/projects/Analysis/code_library/visualization.py"
source(r_utilities)
library(reticulate)
py_run_file(py_utilities, convert = F)
py_run_file(py_visualization, convert = F)
######
######
h5ad_path = "/home/haoy/projects/MetaStudiesAnalysis/Integration/Data/Skin/Human/data_preparation_HCA_v1_0/Reynolds_Greenleaf_NatureGenetics_2023/GSE212447_HC.h5ad"
######
X_gene_path = "/home/haoy/projects/Analysis/useful_files/X_Gene_Homo_sapiens.GRCh38.103.rds"
Y_gene_path = "/home/haoy/projects/Analysis/useful_files/Y_Gene_Homo_sapiens.GRCh38.103.rds"
alias_gene_c_seurat_path = "/home/haoy/projects/Analysis/useful_files/hgnc_complete_set_c_seurat.rds"
CellCyele_gene_path = "/home/haoy/projects/MetaStudiesAnalysis/Integration/regev_lab_cell_cycle_genes.txt"
######
py$output_dir = "/home/haoy/projects/MetaStudiesAnalysis/Integration/Output/Reynolds_Greenleaf_NatureGenetics_2023/HC/v1_gauss"
dir.create(py$output_dir, showWarnings = F, recursive = T)
######
# Settings
py_run_string("sc.settings.verbosity = 3")
py_run_string("sc.logging.print_header()")
py_run_string("sc.settings.set_figure_params(dpi=600, dpi_save=600, facecolor='white', fontsize=10)")
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
UMAP_min_dist = 0.5
# DE parameters
top_DE_range = 50
DEG_embedding_highlight_max_number = 50
max_gene_number = 200
# Other parameters
py_run_string("density_color = 'YlOrRd'")
# Other settings
use_cache = T
######
######
######
summary_list = list()
######
py_run_string(paste0("adata = ad.read_h5ad('", h5ad_path, "')"))
py_run_string("print(adata.X)")
py_run_string("print(adata)")
py_run_string("print(adata.var)")
py_run_string("print(adata.obs[['sex_ontology_term', 'age_years']].value_counts())")
py_run_string("print(adata.obs['disease'].value_counts())")
py_run_string("print(adata.obs['sample_id'].value_counts())")
py_run_string("print(adata.obs['sample_id_running'].value_counts())")
py_run_string("print(adata.obs[['sample_collection_method', 'sample_preservation_method']].value_counts())")
######
py_run_string("adata.var['mt'] = adata.var['feature_name'].str.startswith('MT-')")
py_run_string("sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)")
######
py_run_string("merge_by = ['sample_id_running']")
# Make a shorter DonorID
py_run_string("adata.obs['DonorID_short'] = adata.obs['donor_id_running'].str.split(pat='___', regex=False).apply(lambda x: x[1])")
#
py_run_string("factors = ['DonorID_short', 'sample_id_running', 'age_years', 'disease', 'sample_collection_method', 'sample_preservation_method']") # We use shorter DonorIDs in each separate analysis
py_run_string("annotation = []")
#
py_run_string("adata.obs[factors + annotation] = adata.obs[factors + annotation].astype(str).apply(lambda x: np.where(x.str.lower().isin(['na', 'nan', 'unknown']), pd.NA, x))")
py_run_string("print(adata.obs[factors + annotation])")
#
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
# First, cells were removed if they had fewer than 200 genes expressed, fewer than 1,000 unique sequenced reads (unique molecular identifiers) or greater than 20% of counts corresponding to mitochondrial genes.
filter =
  py$adata$obs["pct_counts_mt"] <= 20 &
  py$adata$obs["n_genes_by_counts"] >= 200 &
  py$adata$obs["total_counts"] >= 1000
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
#
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
comparison_df = cbind(as.data.frame(comparison_mat), cbind(as.character(py$adata$obs$DonorID_short), as.character(py$adata$obs$experiment_id)))
colnames(comparison_df) = c(colnames(comparison_mat), "donor_id_running", "experiment_id")
p = StackedVlnPlot(comparison_df, "donor_id_running", colnames(comparison_mat), "experiment_id", ymaxs = c(200000, max(comparison_df$nFeatures), 0.2, 0.1, 0.001))
ggsave(paste0(py$figure_dir, "/qc.pdf"), p, width = 1.5 * length(unique(comparison_df[, "donor_id_running"])) + 3, height = 12, limitsize = FALSE)
saveRDS(comparison_df, paste0(py$output_dir, "/comparison_df.rds"))
write.table(comparison_df, paste0(py$output_dir, "/comparison_df.tsv"), sep = "\t", row.names = T, col.names = T, quote = F)
######
cell_number = py$adata$n_obs
py$dotsize = max(c(8, get_dot_size(cell_number)))
#
mapping_GeneID_GeneName = py$adata$var[["feature_name"]]
names(mapping_GeneID_GeneName) = py$adata$var_names$values
############
# Cell cycle
############
# Following https://nbviewer.org/github/theislab/scanpy_usage/blob/master/180209_cell_cycle/cell_cycle.ipynb
CellCyele_gene = gsub(" ", "", toupper(as.matrix(read.table(CellCyele_gene_path))[, 1]))
CellCyele_gene[CellCyele_gene == "MLF1IP"] = "CENPU" # MLF1IP is coded by CENPU
py$s_genes = names(mapping_GeneID_GeneName)[mapping_GeneID_GeneName %in% CellCyele_gene[1:43]]
py$g2m_genes = names(mapping_GeneID_GeneName)[mapping_GeneID_GeneName %in% CellCyele_gene[44:length(CellCyele_gene)]]
py_run_string(paste0("sc.tl.score_genes_cell_cycle(adata, s_genes=s_genes, g2m_genes=g2m_genes, random_state=", random_state, ")"))
############
# Clustering
############
adata_clustering_path = paste0(py$output_dir, "/adata_clustering.h5ad")
adata_clustering_hvg_path = paste0(py$output_dir, "/adata_clustering_hvg.h5ad")
clustering_keys_path = paste0(py$output_dir, "/clustering_keys.rds")
clustering_embedding = "PCA_L1"
visualization_embedding = "UMAP_L1"
if(file.exists(adata_clustering_path) && file.exists(adata_clustering_hvg_path)){
  py_run_string(paste0("adata_clustering = ad.read_h5ad('", adata_clustering_path, "')"))
  py_run_string(paste0("adata_clustering_hvg = ad.read_h5ad('", adata_clustering_hvg_path, "')"))
  tryCatch({
    py_run_string("del adata_clustering.uns")
  }, error = function(x){
    print(x)
  })
  clustering_keys = readRDS(clustering_keys_path)
}else{
  ############
  # Imbalanced
  ############
  # Clustering
  pca_cal_list = PCA_calculation_python(adata = py$adata, random_state = random_state, merge_by = unlist(py$merge_by), hvg_number = hvg_number, pca_dim = pca_dim, return_all = T)
  py$adata_clustering_hvg = pca_cal_list[["hvg_matrix"]]
  py$adata_clustering = pca_cal_list[["complete_matrix"]]
  py$adata_clustering$obsm[clustering_embedding] = pca_cal_list[["hvg_matrix"]]$obsm["PCA_use"]
  py_run_string("adata_clustering.X = adata_clustering.layers['raw_counts'].copy()")
  clustering_keys = pipeline_post_PCA("adata_clustering", clustering_embedding, visualization_embedding,
                                      n_neighbors_clustering, n_neighbors_visualization, leiden_resolutions,
                                      UMAP_min_dist=UMAP_min_dist,
                                      clustering_prefix = "_leiden_",
                                      random_state = random_state)
  ####################
  # Cache - Clustering
  ####################
  if(use_cache){
    py_run_string(paste0("adata_clustering.write('", adata_clustering_path, "')"))
    py_run_string(paste0("adata_clustering_hvg.write('", adata_clustering_hvg_path, "')"))
    saveRDS(clustering_keys, clustering_keys_path)
  }
}
py_run_string(paste0("sc.pl.pca_loadings(adata_clustering_hvg, components = '", paste(seq(pca_dim), collapse = ", "), "', save='.pdf')"))
py_run_string(paste0("sc.pl.embedding(adata_clustering, '", visualization_embedding, "', color=factors + ['phase'], legend_loc='right margin', legend_fontsize=6, size=dotsize, frameon=False, ncols=1, save='_factors.pdf')"))
py_run_string(paste0("sc.pl.embedding(adata_clustering, '", visualization_embedding, "', color=['", paste(clustering_keys, collapse = "', '"), "'] + annotation, legend_loc='on data', legend_fontsize=6, size=dotsize, frameon=False, ncols=7, save='_clustering.pdf')"))
for(ii in c(clustering_keys, unlist(py$factors))){
  py_run_string(paste0("clustering_plot(adata_clustering, '", ii, "', basis='", visualization_embedding, "', size=dotsize, colorbar_loc=None, ncols=7, save='_detail_", ii, ".pdf')"))
}
############
# Signatures
############
#
# Analysis
#
adata_str = "adata_clustering"
feature_name = eval(parse(text = paste0("py$", adata_str, "$var[['feature_name']]")))
#
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
#
py$Keratinocyte_markers = markers_list[["Keratinocyte"]][markers_list[["Keratinocyte"]] %in% feature_name]
py$Fibroblast_MuralCell_markers = markers_list[["Fibroblast_MuralCell"]][markers_list[["Fibroblast_MuralCell"]] %in% feature_name]
py$NeuralCrestderivedCell_markers = markers_list[["NeuralCrestderivedCell"]][markers_list[["NeuralCrestderivedCell"]] %in% feature_name]
py$EndothelialCell_markers = markers_list[["EndothelialCell"]][markers_list[["EndothelialCell"]] %in% feature_name]
py$ImmuneCell_markers = markers_list[["ImmuneCell"]][markers_list[["ImmuneCell"]] %in% feature_name]
py$Plasma_Erythrocyte_markers = markers_list[["Plasma_Erythrocyte"]][markers_list[["Plasma_Erythrocyte"]] %in% feature_name]
#
py_run_string(paste0("sc.pl.embedding(", adata_str, ", '", visualization_embedding, "', color=Keratinocyte_markers, gene_symbols='feature_name', layer='normalized', legend_loc='on data', legend_fontsize=6, size=8, cmap=gene_highlight_cmap, ncols=6, frameon=False, save='_", adata_str, "_Keratinocyte.pdf')"))
py_run_string(paste0("sc.pl.embedding(", adata_str, ", '", visualization_embedding, "', color=Fibroblast_MuralCell_markers, gene_symbols='feature_name', layer='normalized', legend_loc='on data', legend_fontsize=6, size=8, cmap=gene_highlight_cmap, ncols=4, frameon=False, save='_", adata_str, "_Fibroblast_MuralCell.pdf')"))
py_run_string(paste0("sc.pl.embedding(", adata_str, ", '", visualization_embedding, "', color=NeuralCrestderivedCell_markers, gene_symbols='feature_name', layer='normalized', legend_loc='on data', legend_fontsize=6, size=8, cmap=gene_highlight_cmap, ncols=4, frameon=False, save='_", adata_str, "_NeuralCrestderivedCell.pdf')"))
py_run_string(paste0("sc.pl.embedding(", adata_str, ", '", visualization_embedding, "', color=EndothelialCell_markers, gene_symbols='feature_name', layer='normalized', legend_loc='on data', legend_fontsize=6, size=8, cmap=gene_highlight_cmap, ncols=4, frameon=False, save='_", adata_str, "_EndothelialCell.pdf')"))
py_run_string(paste0("sc.pl.embedding(", adata_str, ", '", visualization_embedding, "', color=ImmuneCell_markers, gene_symbols='feature_name', layer='normalized', legend_loc='on data', legend_fontsize=6, size=8, cmap=gene_highlight_cmap, ncols=6, frameon=False, save='_", adata_str, "_ImmuneCell.pdf')"))
py_run_string(paste0("sc.pl.embedding(", adata_str, ", '", visualization_embedding, "', color=Plasma_Erythrocyte_markers, gene_symbols='feature_name', layer='normalized', legend_loc='on data', legend_fontsize=6, size=8, cmap=gene_highlight_cmap, ncols=3, frameon=False, save='_", adata_str, "_Plasma_Erythrocyte.pdf')"))
#
py_run_string("marker_list_main = ['KRT14', 'FBLN1', 'ACTA2', 'TYRP1', 'MPZ', 'VWF', 'LYVE1', 'TPSAB1', 'IL1B', 'TRBC2', 'HBA1', 'JCHAIN']")
#
py_run_string("marker_dict_KC_HF = {
  'Bulge (Related)': ['DIO2', 'LGR5', 'TIMP3', 'COMP', 'MGP', 'EPCAM', 'FBLN1', 'PTHLH'],
  'IRS/Cortex': ['KRT25', 'KRT27', 'KRT28', 'KRT71', 'TCHH', 'KRT35', 'KRT85', 'MT4'],
  'GL/ORS/CP': ['TK1', 'MCM3', 'HELLS', 'KRT5', 'KRT17', 'GJA1', 'KRT6B', 'CLDN4'],
  'uHF': ['PTN', 'SPRR1B', 'NCOA7', 'KRT1', 'COL17A1', 'POSTN'],
  'Gland Transition': ['APOE', 'KRT10', 'CLEC2B'],
  'Sebaceous Gland': ['MGST1', 'FASN', 'AWAT2', 'IL1R2'],
  'Sweat Gland/Duct': ['KRT19', 'AQP5', 'CA6', 'SNORC', 'CLDN10', 'WFDC2', 'SLPI', 'AZGP1', 'PIP', 'DCD', 'MUCL1', 'WFDC3', 'MMP7'],
  'Channel': ['KRT23', 'SEMA3C', 'FGF7', 'LAMB3', 'PLAUR', 'MYC', 'IER3', 'KRT17', 'CLCA2']
}")
py_run_string("marker_dict_KC_IFE = {
  'IFE Cornified': ['H19', 'LOR', 'FLG'],
  'IFE Granular': ['SPINK5', 'CALML5', 'KRT2'],
  'IFE Spinous': ['KRTDAP', 'KRT1', 'KRT10'],
  'IFE Basal': ['KRT14', 'KRT5', 'KRT15'],
  'IFE Subtype': ['KRT23', 'KRT77', 'KRT78', 'KRT80', 'CST6', 'CNFN', 'NCCRP1', 'SLURP1', 'SLC35F1', 'SLC24A3', 'SOX5', 'TRAPPC5'],
  'IRS/Cortex': ['KRT35', 'KRT85'],
  'IFE Cycle - G2M': ['GINS2'],
  'IFE Cycle - S': ['MKI67'],
  'Infundibulum': ['MOXD1', 'TM4SF1']
}")
py_run_string("marker_dict_FIB = {
  'Stressed (Pan)': ['CXCL1', 'CXCL2', 'CXCL3', 'HMOX1'],
  'FRC-like': ['APOE', 'C7', 'APOC1', 'IL33'],
  'FRC-like subtypes': ['MYOC', 'MEDAG', 'PTGS2', 'IL6', 'CCL19', 'IGFBP3'],
  'Papillary': ['WIF1', 'APCDD1', 'NKD2', 'COL23A1', 'COL13A1', 'NPTX2', 'COL6A5', 'HSPB3', 'CLEC2A'],
  'Reticular': ['SLPI', 'C1QTNF3'],
  'Reticular (stressed)': ['HSPB2'],
  'Interferon': ['ISG15', 'MX1', 'IFIT3', 'IFI44L', 'OAS1'],
  'HF-associated': ['ASPN', 'TNMD', 'PPP1R14A', 'TNN'],
  'HF-associated subtypes': ['INHBA', 'TFAP2A', 'WNT5A', 'COCH', 'POSTN', 'COL18A1', 'BGN', 'LOXL1', 'TNC', 'GPC3', 'COL11A1', 'DPEP1', 'PI16', 'CXCL12'],
  'Schwann-like (RAMP1+)': ['IGFBP2', 'FGFBP2', 'OLFML2A'],
  'Schwann-like (NGFR+)': ['ECRG4', 'EBF2'],
  'Schwann-like (NGFR+) subtypes': ['CDH19', 'ANGPTL7', 'ITGA6'],
  'Bridge (Reticular-Adipo) (stressed)': ['PTGS2', 'S100A2'],
  'MIX(KC)': ['KRT14', 'KRT5', 'KRT1', 'KRT10'],
  'MIX(LYM)': ['PTPRC', 'CD3D'],
  'MIX(MYE)': ['HLA-DRA', 'HLA-DPA1', 'CD74'],
  'MIX(MUR)': ['MCAM', 'RGS5']
}")
py_run_string("marker_dict_MUR = {
  'Pericyte': ['RERGL', 'RGS5'],
  'Pericyte_sub': ['NR4A2', 'RASD1', 'TCIM', 'HAS2', 'CYP26B1', 'CXCL12', 'FHL2', 'GGT5'],
  'Smooth Muscle Cell': ['DES', 'ACTG2'],
  'MIX(FIB)': ['LUM', 'VCAN'],
  'MIX(MYE)': ['CD74', 'HLA-DRA'],
  'MIX(LYM)': ['CD52', 'PTPRC'],
  'MIX(KC)': ['KRT5', 'KRT14']
}")
py_run_string("marker_dict_MYE = {
  'Monocyte': ['CD14'],
  'MDSC': ['S100A8', 'S100A9', 'TREM1'],
  'Macrophage': ['CD68', 'CD163', 'MARCO'],
  'Macrophage_sub': ['MSR1', 'GPNMB', 'FCGR3A', 'OLR1', 'STAT1', 'GBP1', 'CXCL2', 'CXCL3', 'MRC1', 'LYVE1'],
  'Dendritic Cell': ['CD1C', 'CLEC10A'],
  'DC_sub': ['IL23A', 'CCL17', 'CLEC9A', 'LAMP3', 'CCR7'],
  'Langerhans Cell': ['CD207', 'CD1A'],
  'Cycling': ['MKI67', 'PCLAF', 'CKS1B'],
  'Plasma': ['JCHAIN', 'IGKC'],
  'MIX(LYM)': ['CD3D', 'TRBC2'],
  'MIX(KC)': ['KRT10', 'KRT14'],
  'MIX(MUR)': ['TAGLN', 'MYL9'],
  'MIX(FIB)': ['COL1A2', 'COL6A2'],
  'MIX(MEL)': ['DCT', 'TYRP1']
}")
py_run_string("marker_dict_MAST = {
  'Mast - Subtype': ['SLC18A2', 'VWA5A', 'DLC1', 'NR4A1', 'IL1RL1', 'BIRC3', 'GCSAML'],
  'Mast - MIX(KC)': ['KRT14', 'KRT5', 'KRT1', 'KRT10'],
  'Mast - MIX(FIB)': ['COL1A1', 'COL1A2', 'COL3A1']
}")
py_run_string("marker_dict_LYM = {
  'Tc': ['CD8A', 'CD8B'],
  'NK': ['NKG7', 'GZMB', 'FCGR3A'],
  'ILC': ['GZMK', 'XCL1', 'XCL2', 'SPINK2'],
  'Th_sub': ['IL13', 'IL17A', 'CCL20', 'CCR6', 'PDE4D'],
  'Naive': ['CCR7', 'SELL'],
  'Treg': ['FOXP3', 'TIGIT', 'CTLA4'],
  'Plasma': ['JCHAIN', 'CD79A'],
  'Cycling': ['MKI67', 'STMN1', 'PCLAF', 'TYMS'],
  'ADIRF/IGFBP7/CCL2': ['ADIRF', 'IGFBP7', 'CCL2'],
  'MIX(KC)': ['KRT14', 'KRT5', 'KRT1', 'KRT10'],
  'MIX(MYE)': ['HLA-DRA', 'HLA-DPA1', 'CD74']
}")
py_run_string("marker_dict_LEC = {
  'Common': ['FLT4', 'NT5E', 'NUDT4'],
  'Lymphatic Valve Downstream': ['ADM', 'ANGPT2', 'SCG3', 'GJA4', 'CLDN11', 'PDLIM1', 'PROCR', 'RAMP3'],
  'Subcapsular Sinus Ceiling': ['NTS', 'PLAAT4', 'IL33'],
  'CXCL1/CXCL2/CXCL3/FABP5/FKBP1A': ['CXCL1', 'CXCL2', 'CXCL3', 'FABP5', 'FKBP1A'],
  'SNHG5/GAS5/SNHG29/SNHG6': ['SNHG5', 'GAS5', 'SNHG29', 'SNHG6'],
  'ADAMTS6/SLC41A1/NEO1/AGRN/CD200/SBSPON': ['ADAMTS6', 'SLC41A1', 'NEO1', 'AGRN', 'CD200', 'SBSPON'],
  'MIX(FIB)': ['FBLN1', 'COL1A1', 'COL1A2'],
  'MIX(MUR)': ['TAGLN', 'MYL9', 'ACTA2'],
  'MIX(MYE)': ['HLA-DRA', 'HLA-DPA1', 'CD74'],
  'MIX(KC)': ['KRT14', 'KRT5', 'KRT1', 'KRT10']
}")
py_run_string("marker_dict_VEC = {
  'Venous': ['VWF'],
  'Vein': ['SELE', 'ICAM1'],
  'Arterial': ['SEMA3G', 'GJA5', 'CXCL12', 'IGFBP3'],
  'MIX(MUR/FIB)': ['RGS5', 'MYL9', 'COL1A2', 'COL6A2'],
  'MIX(KC)': ['KRT14', 'KRT5', 'KRT1', 'KRT10'],
  'MIX(LYM)': ['CD3D', 'CD3E'],
  'MIX(MEL)': ['DCT', 'TYRP1', 'MLANA', 'PMEL']
}")
py_run_string("marker_dict_MEL = {
  'CLU/CRABP1/UQCR11': ['CLU', 'CRABP1', 'UQCR11'],
  'HLA-H': ['HLA-H'],
  'NDRG2/SNHG29/EEF1G/GAS5': ['NDRG2', 'SNHG29', 'EEF1G', 'GAS5'],
  'ADIRF/FXYD3/BSG/TTYH1/CADM3/CHI3L2/SUSD4': ['ADIRF', 'FXYD3', 'BSG', 'TTYH1', 'CADM3', 'CHI3L2', 'SUSD4'],
  'RETSAT/L1CAM/CCND1': ['RETSAT', 'L1CAM', 'CCND1'],
  'NME1/S100B': ['NME1', 'S100B'],
  'FOSB/NR4A2': ['FOSB', 'NR4A2'],
  'PAMR1/MYC': ['PAMR1', 'MYC'],
  'IFI44L/MX1/IFIT3/ISG15/PLSCR1/IFIT1/OAS1': ['IFI44L', 'MX1', 'IFIT3', 'ISG15', 'PLSCR1', 'IFIT1', 'OAS1'],
  'MIX(LYM)': ['CD3D', 'CD3E'],
  'MIX(FIB)': ['COL1A1', 'COL1A2', 'COL6A1'],
  'MIX(KC)': ['KRT14', 'KRT5', 'KRT1', 'KRT10']
}")
py_run_string("marker_dict_SC = {
  'Schwann Cell': ['S100B'],
  'SC - Remak': ['NRXN1', 'SCN7A'],
  'SC - Remak Subtypes (SOCS3-/EGR1-/S100A16-)': ['XKR4', 'FN1', 'ANK2', 'NRXN3', 'COL28A1', 'SPARCL1', 'CFD', 'TMOD2', 'SBSPON', 'CADM2'],
  'SC - Remak - MHC II': ['HLA-DRB1', 'CD74', 'HLA-DRB5'],
  'SOCS3/EGR1/S100A16': ['SOCS3', 'EGR1', 'S100A16'],
  'SC - Remak Subtypes (SOCS3+/EGR1+)': ['CCN5', 'IFI27', 'ABL2', 'PDLIM3'],
  'SC - Remak Repair': ['NGFR'],
  'SC - Nociceptive': ['BCHE', 'AQP1', 'PALMD'],
  'SC - Nociceptive (Extend)': ['PRSS23', 'PTN', 'IGFBP5', 'PTPRZ1'],
  'SC - Myelinated': ['PRX', 'GLDN', 'MLIP', 'NCMAP', 'BCAS1'],
  'SC - Myelinated (Stressed)': ['JUNB', 'ATF3', 'FOS', 'FOSB', 'HSPA1B', 'ZFP36'],
  'MPZ/SOX': ['MPZ', 'SOX2', 'SOX10'],
  'SC - MIX(MUR)': ['TAGLN', 'ACTA2'],
  'SC - MIX(KC/FIB)': ['DCN', 'CXCL14', 'APOD'],
  'SC - MIX(MEL)': ['DCT', 'KIT', 'ANXA2P2', 'QPCT']
}")
######
if(FALSE){
  py_run_string(paste0("sc.pl.embedding(", adata_str, ", '", visualization_embedding, "', color=[x for x in marker_list_main if np.isin(x, ", adata_str, ".var['feature_name'])], gene_symbols='feature_name', layer='normalized', legend_loc='on data', legend_fontsize=6, size=8, cmap=gene_highlight_cmap, ncols=4, frameon=False, add_outline=False, save='_", adata_str, "_marker_list_main.pdf')"))
  py_run_string("gc.collect()")
  py_run_string(paste0("sc.pl.embedding(", adata_str, ", '", visualization_embedding, "', color=[x for y in marker_dict_KC_HF.values() for x in y if np.isin(x, ", adata_str, ".var['feature_name'])], gene_symbols='feature_name', layer='normalized', legend_loc='on data', legend_fontsize=6, size=8, cmap=gene_highlight_cmap, ncols=5, frameon=False, add_outline=False, save='_", adata_str, "_marker_dict_KC_HF.pdf')"))
  py_run_string("gc.collect()")
  py_run_string(paste0("sc.pl.embedding(", adata_str, ", '", visualization_embedding, "', color=[x for y in marker_dict_KC_IFE.values() for x in y if np.isin(x, ", adata_str, ".var['feature_name'])], gene_symbols='feature_name', layer='normalized', legend_loc='on data', legend_fontsize=6, size=8, cmap=gene_highlight_cmap, ncols=5, frameon=False, add_outline=False, save='_", adata_str, "_marker_dict_KC_IFE.pdf')"))
  py_run_string("gc.collect()")
  py_run_string(paste0("sc.pl.embedding(", adata_str, ", '", visualization_embedding, "', color=[x for y in marker_dict_FIB.values() for x in y if np.isin(x, ", adata_str, ".var['feature_name'])], gene_symbols='feature_name', layer='normalized', legend_loc='on data', legend_fontsize=6, size=8, cmap=gene_highlight_cmap, ncols=5, frameon=False, add_outline=False, save='_", adata_str, "_marker_dict_FIB.pdf')"))
  py_run_string("gc.collect()")
  py_run_string(paste0("sc.pl.embedding(", adata_str, ", '", visualization_embedding, "', color=[x for y in marker_dict_MUR.values() for x in y if np.isin(x, ", adata_str, ".var['feature_name'])], gene_symbols='feature_name', layer='normalized', legend_loc='on data', legend_fontsize=6, size=8, cmap=gene_highlight_cmap, ncols=5, frameon=False, add_outline=False, save='_", adata_str, "_marker_dict_MUR.pdf')"))
  py_run_string("gc.collect()")
  py_run_string(paste0("sc.pl.embedding(", adata_str, ", '", visualization_embedding, "', color=[x for y in marker_dict_MYE.values() for x in y if np.isin(x, ", adata_str, ".var['feature_name'])], gene_symbols='feature_name', layer='normalized', legend_loc='on data', legend_fontsize=6, size=8, cmap=gene_highlight_cmap, ncols=5, frameon=False, add_outline=False, save='_", adata_str, "_marker_dict_MYE.pdf')"))
  py_run_string("gc.collect()")
  py_run_string(paste0("sc.pl.embedding(", adata_str, ", '", visualization_embedding, "', color=[x for y in marker_dict_MAST.values() for x in y if np.isin(x, ", adata_str, ".var['feature_name'])], gene_symbols='feature_name', layer='normalized', legend_loc='on data', legend_fontsize=6, size=8, cmap=gene_highlight_cmap, ncols=5, frameon=False, add_outline=False, save='_", adata_str, "_marker_dict_MAST.pdf')"))
  py_run_string("gc.collect()")
  py_run_string(paste0("sc.pl.embedding(", adata_str, ", '", visualization_embedding, "', color=[x for y in marker_dict_LYM.values() for x in y if np.isin(x, ", adata_str, ".var['feature_name'])], gene_symbols='feature_name', layer='normalized', legend_loc='on data', legend_fontsize=6, size=8, cmap=gene_highlight_cmap, ncols=5, frameon=False, add_outline=False, save='_", adata_str, "_marker_dict_LYM.pdf')"))
  py_run_string("gc.collect()")
  py_run_string(paste0("sc.pl.embedding(", adata_str, ", '", visualization_embedding, "', color=[x for y in marker_dict_LEC.values() for x in y if np.isin(x, ", adata_str, ".var['feature_name'])], gene_symbols='feature_name', layer='normalized', legend_loc='on data', legend_fontsize=6, size=8, cmap=gene_highlight_cmap, ncols=5, frameon=False, add_outline=False, save='_", adata_str, "_marker_dict_LEC.pdf')"))
  py_run_string("gc.collect()")
  py_run_string(paste0("sc.pl.embedding(", adata_str, ", '", visualization_embedding, "', color=[x for y in marker_dict_VEC.values() for x in y if np.isin(x, ", adata_str, ".var['feature_name'])], gene_symbols='feature_name', layer='normalized', legend_loc='on data', legend_fontsize=6, size=8, cmap=gene_highlight_cmap, ncols=5, frameon=False, add_outline=False, save='_", adata_str, "_marker_dict_VEC.pdf')"))
  py_run_string("gc.collect()")
  py_run_string(paste0("sc.pl.embedding(", adata_str, ", '", visualization_embedding, "', color=[x for y in marker_dict_MEL.values() for x in y if np.isin(x, ", adata_str, ".var['feature_name'])], gene_symbols='feature_name', layer='normalized', legend_loc='on data', legend_fontsize=6, size=8, cmap=gene_highlight_cmap, ncols=5, frameon=False, add_outline=False, save='_", adata_str, "_marker_dict_MEL.pdf')"))
  py_run_string("gc.collect()")
  py_run_string(paste0("sc.pl.embedding(", adata_str, ", '", visualization_embedding, "', color=[x for y in marker_dict_SC.values() for x in y if np.isin(x, ", adata_str, ".var['feature_name'])], gene_symbols='feature_name', layer='normalized', legend_loc='on data', legend_fontsize=6, size=8, cmap=gene_highlight_cmap, ncols=5, frameon=False, add_outline=False, save='_", adata_str, "_marker_dict_SC.pdf')"))
  py_run_string("gc.collect()")
  #
  # DEG analysis
  #
  DEG_key = "PCA_L1_leiden_020"
  py_run_string(paste0("sc.tl.dendrogram(", adata_str, ", groupby='", DEG_key, "', use_rep='", clustering_embedding, "')"))
  py_run_string(paste0("dot_fig = plot_dotplot_with_annotations(", adata_str, ", var_names=marker_list_main, groupby='", DEG_key, "', gene_symbols='feature_name', layer='normalized', dendrogram=True, n_genes=10, standard_scale = 'var', swap_axes = False, highlight_genes=None, highlight_params = {'facecolor': 'yellow', 'edgecolor': 'none', 'alpha': 0.5}, highlight_cluster_line = False, cluster_line_highlight_params = {'color': 'lightgreen', 'alpha': 0.2, 'zorder': -5},)"))
  py_run_string(paste0("dot_fig.savefig(f'{figure_dir}/dotplot_", adata_str, "_", DEG_key, "_marker_list_main_var.pdf', bbox_inches = 'tight')"))
  py_run_string(paste0("dot_fig = plot_dotplot_with_annotations(", adata_str, ", var_names={key: [x for x in value if np.isin(x, ", adata_str, ".var['feature_name'])] for key,value in marker_dict_KC_HF.items()}, groupby='", DEG_key, "', gene_symbols='feature_name', layer='normalized', dendrogram=True, n_genes=10, standard_scale = 'var', swap_axes = False, highlight_genes=None, highlight_params = {'facecolor': 'yellow', 'edgecolor': 'none', 'alpha': 0.5}, highlight_cluster_line = False, cluster_line_highlight_params = {'color': 'lightgreen', 'alpha': 0.2, 'zorder': -5},)"))
  py_run_string(paste0("dot_fig.savefig(f'{figure_dir}/dotplot_", adata_str, "_", DEG_key, "_marker_dict_KC_HF_var.pdf', bbox_inches = 'tight')"))
  py_run_string(paste0("dot_fig = plot_dotplot_with_annotations(", adata_str, ", var_names={key: [x for x in value if np.isin(x, ", adata_str, ".var['feature_name'])] for key,value in marker_dict_KC_IFE.items()}, groupby='", DEG_key, "', gene_symbols='feature_name', layer='normalized', dendrogram=True, n_genes=10, standard_scale = 'var', swap_axes = False, highlight_genes=None, highlight_params = {'facecolor': 'yellow', 'edgecolor': 'none', 'alpha': 0.5}, highlight_cluster_line = False, cluster_line_highlight_params = {'color': 'lightgreen', 'alpha': 0.2, 'zorder': -5},)"))
  py_run_string(paste0("dot_fig.savefig(f'{figure_dir}/dotplot_", adata_str, "_", DEG_key, "_marker_dict_KC_IFE_var.pdf', bbox_inches = 'tight')"))
  py_run_string(paste0("dot_fig = plot_dotplot_with_annotations(", adata_str, ", var_names={key: [x for x in value if np.isin(x, ", adata_str, ".var['feature_name'])] for key,value in marker_dict_FIB.items()}, groupby='", DEG_key, "', gene_symbols='feature_name', layer='normalized', dendrogram=True, n_genes=10, standard_scale = 'var', swap_axes = False, highlight_genes=None, highlight_params = {'facecolor': 'yellow', 'edgecolor': 'none', 'alpha': 0.5}, highlight_cluster_line = False, cluster_line_highlight_params = {'color': 'lightgreen', 'alpha': 0.2, 'zorder': -5},)"))
  py_run_string(paste0("dot_fig.savefig(f'{figure_dir}/dotplot_", adata_str, "_", DEG_key, "_marker_dict_FIB_var.pdf', bbox_inches = 'tight')"))
  py_run_string(paste0("dot_fig = plot_dotplot_with_annotations(", adata_str, ", var_names={key: [x for x in value if np.isin(x, ", adata_str, ".var['feature_name'])] for key,value in marker_dict_MUR.items()}, groupby='", DEG_key, "', gene_symbols='feature_name', layer='normalized', dendrogram=True, n_genes=10, standard_scale = 'var', swap_axes = False, highlight_genes=None, highlight_params = {'facecolor': 'yellow', 'edgecolor': 'none', 'alpha': 0.5}, highlight_cluster_line = False, cluster_line_highlight_params = {'color': 'lightgreen', 'alpha': 0.2, 'zorder': -5},)"))
  py_run_string(paste0("dot_fig.savefig(f'{figure_dir}/dotplot_", adata_str, "_", DEG_key, "_marker_dict_MUR_var.pdf', bbox_inches = 'tight')"))
  py_run_string(paste0("dot_fig = plot_dotplot_with_annotations(", adata_str, ", var_names={key: [x for x in value if np.isin(x, ", adata_str, ".var['feature_name'])] for key,value in marker_dict_MYE.items()}, groupby='", DEG_key, "', gene_symbols='feature_name', layer='normalized', dendrogram=True, n_genes=10, standard_scale = 'var', swap_axes = False, highlight_genes=None, highlight_params = {'facecolor': 'yellow', 'edgecolor': 'none', 'alpha': 0.5}, highlight_cluster_line = False, cluster_line_highlight_params = {'color': 'lightgreen', 'alpha': 0.2, 'zorder': -5},)"))
  py_run_string(paste0("dot_fig.savefig(f'{figure_dir}/dotplot_", adata_str, "_", DEG_key, "_marker_dict_MYE_var.pdf', bbox_inches = 'tight')"))
  py_run_string(paste0("dot_fig = plot_dotplot_with_annotations(", adata_str, ", var_names={key: [x for x in value if np.isin(x, ", adata_str, ".var['feature_name'])] for key,value in marker_dict_MAST.items()}, groupby='", DEG_key, "', gene_symbols='feature_name', layer='normalized', dendrogram=True, n_genes=10, standard_scale = 'var', swap_axes = False, highlight_genes=None, highlight_params = {'facecolor': 'yellow', 'edgecolor': 'none', 'alpha': 0.5}, highlight_cluster_line = False, cluster_line_highlight_params = {'color': 'lightgreen', 'alpha': 0.2, 'zorder': -5},)"))
  py_run_string(paste0("dot_fig.savefig(f'{figure_dir}/dotplot_", adata_str, "_", DEG_key, "_marker_dict_MAST_var.pdf', bbox_inches = 'tight')"))
  py_run_string(paste0("dot_fig = plot_dotplot_with_annotations(", adata_str, ", var_names={key: [x for x in value if np.isin(x, ", adata_str, ".var['feature_name'])] for key,value in marker_dict_LYM.items()}, groupby='", DEG_key, "', gene_symbols='feature_name', layer='normalized', dendrogram=True, n_genes=10, standard_scale = 'var', swap_axes = False, highlight_genes=None, highlight_params = {'facecolor': 'yellow', 'edgecolor': 'none', 'alpha': 0.5}, highlight_cluster_line = False, cluster_line_highlight_params = {'color': 'lightgreen', 'alpha': 0.2, 'zorder': -5},)"))
  py_run_string(paste0("dot_fig.savefig(f'{figure_dir}/dotplot_", adata_str, "_", DEG_key, "_marker_dict_LYM_var.pdf', bbox_inches = 'tight')"))
  py_run_string(paste0("dot_fig = plot_dotplot_with_annotations(", adata_str, ", var_names={key: [x for x in value if np.isin(x, ", adata_str, ".var['feature_name'])] for key,value in marker_dict_LEC.items()}, groupby='", DEG_key, "', gene_symbols='feature_name', layer='normalized', dendrogram=True, n_genes=10, standard_scale = 'var', swap_axes = False, highlight_genes=None, highlight_params = {'facecolor': 'yellow', 'edgecolor': 'none', 'alpha': 0.5}, highlight_cluster_line = False, cluster_line_highlight_params = {'color': 'lightgreen', 'alpha': 0.2, 'zorder': -5},)"))
  py_run_string(paste0("dot_fig.savefig(f'{figure_dir}/dotplot_", adata_str, "_", DEG_key, "_marker_dict_LEC_var.pdf', bbox_inches = 'tight')"))
  py_run_string(paste0("dot_fig = plot_dotplot_with_annotations(", adata_str, ", var_names={key: [x for x in value if np.isin(x, ", adata_str, ".var['feature_name'])] for key,value in marker_dict_VEC.items()}, groupby='", DEG_key, "', gene_symbols='feature_name', layer='normalized', dendrogram=True, n_genes=10, standard_scale = 'var', swap_axes = False, highlight_genes=None, highlight_params = {'facecolor': 'yellow', 'edgecolor': 'none', 'alpha': 0.5}, highlight_cluster_line = False, cluster_line_highlight_params = {'color': 'lightgreen', 'alpha': 0.2, 'zorder': -5},)"))
  py_run_string(paste0("dot_fig.savefig(f'{figure_dir}/dotplot_", adata_str, "_", DEG_key, "_marker_dict_VEC_var.pdf', bbox_inches = 'tight')"))
  py_run_string(paste0("dot_fig = plot_dotplot_with_annotations(", adata_str, ", var_names={key: [x for x in value if np.isin(x, ", adata_str, ".var['feature_name'])] for key,value in marker_dict_MEL.items()}, groupby='", DEG_key, "', gene_symbols='feature_name', layer='normalized', dendrogram=True, n_genes=10, standard_scale = 'var', swap_axes = False, highlight_genes=None, highlight_params = {'facecolor': 'yellow', 'edgecolor': 'none', 'alpha': 0.5}, highlight_cluster_line = False, cluster_line_highlight_params = {'color': 'lightgreen', 'alpha': 0.2, 'zorder': -5},)"))
  py_run_string(paste0("dot_fig.savefig(f'{figure_dir}/dotplot_", adata_str, "_", DEG_key, "_marker_dict_MEL_var.pdf', bbox_inches = 'tight')"))
  py_run_string(paste0("dot_fig = plot_dotplot_with_annotations(", adata_str, ", var_names={key: [x for x in value if np.isin(x, ", adata_str, ".var['feature_name'])] for key,value in marker_dict_SC.items()}, groupby='", DEG_key, "', gene_symbols='feature_name', layer='normalized', dendrogram=True, n_genes=10, standard_scale = 'var', swap_axes = False, highlight_genes=None, highlight_params = {'facecolor': 'yellow', 'edgecolor': 'none', 'alpha': 0.5}, highlight_cluster_line = False, cluster_line_highlight_params = {'color': 'lightgreen', 'alpha': 0.2, 'zorder': -5},)"))
  py_run_string(paste0("dot_fig.savefig(f'{figure_dir}/dotplot_", adata_str, "_", DEG_key, "_marker_dict_SC_var.pdf', bbox_inches = 'tight')"))
}
# Annotation
DEG_key = "PCA_L1_leiden_020"
for(used_table in c("top_p_by_score", "top_p_by_FC", "top_p_by_P")){
  this_path = paste0(py$output_dir, "/DEG_", adata_str, "_" , DEG_key, ".rds")
  if(file.exists(this_path)){
    DEG_result = readRDS(this_path)
  }else{
    DEG_result = DEG_Analysis(eval(parse(text = paste0("py$", adata_str))), DEG_key, py$output_dir, prefix = paste0(adata_str, "_" , DEG_key, "_"), top_DE_range = top_DE_range)
    saveRDS(DEG_result, this_path)
  }
  this_table = DEG_result[["rest"]][[used_table]]
  py_run_string(paste0("sc.tl.dendrogram(", adata_str, ", groupby='", DEG_key, "', use_rep='", clustering_embedding, "')"))
  this_clusters = unique(this_table$Cluster)
  top_DEGs = c()
  py_run_string("DEG_dict = {}")
  for(this_cluster in this_clusters){
    top_DEGs_this_cluster = mapping_GeneID_GeneName[head(this_table[this_table$Cluster == this_cluster, ], 10)$Gene]
    top_DEGs = c(top_DEGs, top_DEGs_this_cluster)
    py_run_string(paste0("DEG_dict['", this_cluster, "'] = ['", paste(top_DEGs_this_cluster, collapse = "', '"), "']"))
  }
  top_DEGs_ordered = c()
  for(this_cluster in eval(parse(text = paste0("py$", adata_str, "$uns[['dendrogram_", DEG_key, "']][['categories_ordered']]")))){
    top_DEGs_ordered = c(top_DEGs_ordered, mapping_GeneID_GeneName[head(this_table[this_table$Cluster == this_cluster, ], 10)$Gene])
  }
  if(length(top_DEGs) <= max_gene_number){
    py_run_string(paste0("dot_fig = plot_dotplot_with_annotations(", adata_str, ", var_names=DEG_dict, groupby='", DEG_key, "', gene_symbols='feature_name', layer='normalized', dendrogram=True, n_genes=10, standard_scale = 'var', swap_axes = False, highlight_genes = DEG_dict, highlight_params = {'facecolor': 'yellow', 'edgecolor': 'none', 'alpha': 0.5}, highlight_cluster_line = False, cluster_line_highlight_params = {'color': 'lightgreen', 'alpha': 0.2, 'zorder': -5},)"))
    py_run_string(paste0("dot_fig.savefig(f'{figure_dir}/dotplot_", adata_str, "_", DEG_key, "_", used_table, ".pdf', bbox_inches = 'tight')"))
  }else{
    for(ii in 1:ceiling(length(top_DEGs) / max_gene_number)){
      this_index = seq(max_gene_number) + max_gene_number * (ii - 1)
      this_index = this_index[this_index <= length(top_DEGs)]
      py$tmp_genes = top_DEGs[this_index]
      py$tmp_genes_ordered = top_DEGs_ordered[this_index]
      py_run_string(paste0("dot_fig = plot_dotplot_with_annotations(", adata_str, ", var_names=tmp_genes_ordered, groupby='", DEG_key, "', gene_symbols='feature_name', layer='normalized', dendrogram=True, n_genes=10, standard_scale = 'var', swap_axes = False, highlight_genes = DEG_dict, highlight_params = {'facecolor': 'yellow', 'edgecolor': 'none', 'alpha': 0.5}, highlight_cluster_line = False, cluster_line_highlight_params = {'color': 'lightgreen', 'alpha': 0.2, 'zorder': -5},)"))
      py_run_string(paste0("dot_fig.savefig(f'{figure_dir}/dotplot_", adata_str, "_", DEG_key, "_", ii, "_", used_table, ".pdf', bbox_inches = 'tight')"))
      py_run_string("del tmp_genes")
    }
  }
}
#####################
# Cell identification
#####################
adata_annotated_1_path = paste0(py$output_dir, "/adata_annotated_1.h5ad")
if(file.exists(adata_annotated_1_path)){
  py_run_string(paste0("adata_annotated = ad.read_h5ad('", adata_annotated_1_path, "')"))
  py_run_string("annotation = annotation if np.isin('individual_annotation_L1', annotation) else annotation + ['individual_annotation_L1']")
  ###
  py_run_string("tmp_SeparateAnnotation = adata_annotated.obs['individual_annotation_L1'].astype(str).values")
  summary_list[[qc_state]][["individual_annotation_L1"]] = list("AllCells" = table(py$tmp_SeparateAnnotation))
  py_run_string("del tmp_SeparateAnnotation")
}else{
  py_run_string(paste0("adata_annotated = adata_clustering.copy()"))
  cell_annotation = rep(NA, py$adata_annotated$n_obs)
  names(cell_annotation) = py$adata_annotated$obs_names$values
  PCA_L1_leiden_020 = as.character(py$adata_annotated$obs[['PCA_L1_leiden_020']])
  PCA_L1_leiden_300 = as.character(py$adata_annotated$obs[['PCA_L1_leiden_300']])
  #
  cell_annotation[PCA_L1_leiden_020 == "0"] = "Keratinocyte"
  cell_annotation[PCA_L1_leiden_020 == "1"] = "Fibroblast"
  cell_annotation[PCA_L1_leiden_020 == "2"] = "Keratinocyte"
  cell_annotation[PCA_L1_leiden_020 == "3"] = "Mural Cell"
  cell_annotation[PCA_L1_leiden_020 == "4"] = "Lymphocyte"
  cell_annotation[PCA_L1_leiden_020 == "5"] = "Vascular Endothelial Cell"
  cell_annotation[PCA_L1_leiden_020 == "6"] = "Myeloid Cell"
  cell_annotation[PCA_L1_leiden_020 == "7"] = "Mast Cell"
  cell_annotation[PCA_L1_leiden_020 == "8"] = "Keratinocyte"
  cell_annotation[PCA_L1_leiden_020 == "9"] = "Keratinocyte + Fibroblast + Myeloid Cell + Mast Cell"
  cell_annotation[PCA_L1_leiden_020 == "10"] = "Keratinocyte"
  cell_annotation[PCA_L1_leiden_020 == "11"] = "Melanocyte"
  cell_annotation[PCA_L1_leiden_020 == "12"] = "Keratinocyte"
  cell_annotation[PCA_L1_leiden_020 == "13"] = "Lymphatic Endothelial Cell"
  cell_annotation[PCA_L1_leiden_020 == "14"] = "Plasma"
  cell_annotation[PCA_L1_leiden_020 == "15"] = "Myeloid Cell"
  cell_annotation[PCA_L1_leiden_020 == "16"] = "Keratinocyte"
  cell_annotation[PCA_L1_leiden_020 == "17"] = "Schwann Cell"
  #
  cell_annotation[PCA_L1_leiden_300 == "42"] = "Fibroblast + Myeloid Cell"
  cell_annotation[PCA_L1_leiden_300 == "51"] = "Keratinocyte + Mural Cell"
  cell_annotation[PCA_L1_leiden_020 == "3" & PCA_L1_leiden_300 == "44"] = "Keratinocyte + Mural Cell"
  #
  print(sum(is.na(cell_annotation)))
  py$tmp_clustering_result = cell_annotation
  py_run_string("adata_annotated.obs['individual_annotation_L1'] = tmp_clustering_result")
  py_run_string("del tmp_clustering_result")
  ###
  py_run_string("tmp_SeparateAnnotation = adata_annotated.obs['individual_annotation_L1'].astype(str).values")
  summary_list[[qc_state]][["individual_annotation_L1"]] = list("AllCells" = table(py$tmp_SeparateAnnotation))
  py_run_string("del tmp_SeparateAnnotation")
  ######
  for(ii in unlist(py$annotation)){
    py_run_string(paste0("tmp_this_annotation = adata_annotated.obs['", ii, "'].astype(str).values"))
    summary_list[[qc_state]][[ii]] = table(py$tmp_this_annotation)
    py_run_string("del tmp_this_annotation")
  }
  py_run_string("annotation = annotation if np.isin('individual_annotation_L1', annotation) else annotation + ['individual_annotation_L1']")
  ######################
  # Cache - Querying - 2
  ######################
  if(use_cache){
    py_run_string(paste0("adata_annotated.write('", adata_annotated_1_path, "')"))
  }
}
###
adata_str = "adata_annotated"
# Ploting will tickle .strings_to_categoricals() automatically, cause error when accessing adata attributes in r
py_run_string(paste0("sc.pl.embedding(", adata_str, ", '", visualization_embedding, "', color=['", paste(clustering_keys, collapse = "', '"), "'] + annotation, legend_loc='on data', legend_fontsize=6, size=dotsize, frameon=False, ncols=7, save='_", adata_str, "_annotation.pdf')"))
for(ii in c("phase", unlist(py$annotation))){
  py_run_string(paste0("clustering_plot(", adata_str, ", '", ii, "', basis='", visualization_embedding, "', size=dotsize, ncols=5, save='_detail_", adata_str, "_", ii, ".pdf')"))
}
density_visualization(adata_str, visualization_embedding, "sex_ontology_term", ".pdf", "density_color", "5", as.character(py$dotsize))
density_visualization(adata_str, visualization_embedding, "individual_annotation_L1", ".pdf", "density_color", "5", as.character(py$dotsize))
py_run_string(paste0("factor_distribution(", adata_str, ", factors + ['phase'], 'individual_annotation_L1', path=figure_dir + '/Distribution_", adata_str, "_SeparateAnnotation.pdf')"))
py_run_string(paste0("factor_distribution_per_method(", adata_str, ", factors + ['phase'], 'individual_annotation_L1', precent='y', path=figure_dir + '/Distribution_", adata_str, "_SeparateAnnotation_divided.pdf')"))
saveRDS(summary_list, paste0(py$output_dir, "/summary_list.rds"))
sink(paste0(py$output_dir, "/summary_list.txt"))
print(summary_list)
sink()

