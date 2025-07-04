py$output_dir = "/home/haoy/projects/MetaStudiesAnalysis/Integration/Output/Integration/v70_gauss"
py$this_output_dir = paste0(py$output_dir, "/IndividualAnalyses")
dir.create(py$this_output_dir, showWarnings = F, recursive = T)
py_run_string("this_figure_dir = this_output_dir + '/figures'")
py_run_string("sc.settings.figdir = this_figure_dir")
dir.create(py$this_figure_dir, showWarnings = F, recursive = T)
###
py_run_string("from matplotlib import font_manager as fm, rcParams")
py_run_string("from matplotlib.gridspec import GridSpec")
py_run_string("fm.fontManager.addfont('/home/haoy/projects/Analysis/fonts/ARIAL.TTF')")
py_run_string("fm.fontManager.addfont('/home/haoy/projects/Analysis/fonts/ARIALBD.TTF')")
py_run_string("plt.rcParams['font.size'] = 6")
py_run_string("plt.rcParams['font.family'] = 'Arial'")
###
na_value = c("NA", "nan", "unknown", "True")
py$classify_dict = list(
  "Keratinocyte" = "Keratinocyte",
  "Melanocyte" = "Melanocyte",
  "Schwann Cell" = "Schwann Cell",
  "Fibroblast" = "Fibroblast",
  "Mural Cell" = "Mural Cell",
  "Vascular Endothelial Cell" = "Vascular Endothelial Cell",
  "Lymphatic Endothelial Cell" = "Lymphatic Endothelial Cell",
  "Erythrocyte" = "Erythrocyte",
  "Lymphocyte" = "Lymphocyte",
  "Myeloid Cell" = "Myeloid Cell",
  "Mast Cell" = "Mast Cell",
  "Plasma" = "Plasma",
  "nan" = "Unknown"
)
py$rename_dict = list(
  "Keratinocyte" = "Keratinocytes",
  "Melanocyte" = "Melanocytes",
  "Schwann Cell" = "Schwann Cells",
  "Fibroblast" = "Fibroblasts",
  "Mural Cell" = "Mural Cells",
  "Vascular Endothelial Cell" = "Vascular Endothelial Cells",
  "Lymphatic Endothelial Cell" = "Lymphatic Endothelial Cells",
  "Erythrocyte" = "Erythrocytes",
  "Lymphocyte" = "Lymphocytes",
  "Myeloid Cell" = "Myeloid Cells",
  "Mast Cell" = "Mast Cells",
  "Plasma" = "Plasma Cells",
  "Unknown" = "Unknown"
)
py$index_name_dict = list(
  "Keratinocytes" = "KC",
  "Melanocytes" = "MEL",
  "Schwann Cells" = "SC",
  "Fibroblasts" = "FIB",
  "Mural Cells" = "MUR",
  "Vascular Endothelial Cells" = "VEC",
  "Lymphatic Endothelial Cells" = "LEC",
  "Erythrocytes" = "ERY",
  "Lymphocytes" = "LYM",
  "Myeloid Cells" = "MYE",
  "Mast Cells" = "MAST",
  "Plasma Cells" = "PC",
  "Unknown" = "NA"
)
py_run_string("color_dict_CT_FULL={'KC: Keratinocytes': colorsys.hsv_to_rgb(0, 0.9, 0.5),
                                   'MEL: Melanocytes': colorsys.hsv_to_rgb(0.1, 0.7, 0.7),
                                   'SC: Schwann Cells': colorsys.hsv_to_rgb(0.2, 0.5, 0.9),
                                   'FIB: Fibroblasts': colorsys.hsv_to_rgb(0.3, 0.9, 0.5),
                                   'MUR: Mural Cells': colorsys.hsv_to_rgb(0.4, 0.5, 0.9),
                                   'VEC: Vascular Endothelial Cells': colorsys.hsv_to_rgb(0.5, 0.9, 0.5),
                                   'LEC: Lymphatic Endothelial Cells': colorsys.hsv_to_rgb(0.6, 0.5, 0.9),
                                   'ERY: Erythrocytes': colorsys.hsv_to_rgb(0.85, 1, 1),
                                   'LYM: Lymphocytes': colorsys.hsv_to_rgb(0.7, 0.9, 0.5),
                                   'MYE: Myeloid Cells': colorsys.hsv_to_rgb(0.8, 0.7, 0.7),
                                   'MAST: Mast Cells': colorsys.hsv_to_rgb(0.9, 0.5, 0.9),
                                   'PC: Plasma Cells': colorsys.hsv_to_rgb(0.75, 1, 1),
                                   'NA: Unknown': colorsys.hsv_to_rgb(0, 0, 0.9)
}")
py_run_string("marker_dict = {
  'KC': ['KRT14', 'KRT15', 'KRT5', 'KRT1', 'KRT10', 'KRTDAP', 'KRT2', 'KRT6A', 'KRT9', 'KRT18', 'KRT23', 'KRT35', 'MGST1', 'NRARP', 'POSTN', 'TK1', 'CENPU', 'MKI67'],
  'MEL': ['TYRP1', 'TYR', 'DCT', 'MLANA', 'PMEL'],
  'SC': ['MPZ', 'SOX10', 'SOX2'],
  'FIB': ['FBLN1', 'COL1A1', 'COL1A2', 'COL3A1', 'COL6A3', 'DCN', 'COCH', 'SFRP2'],
  'MUR': ['ACTA2', 'RGS5', 'RERGL', 'TAGLN'],
  'VEC': ['VWF', 'SELE', 'SELP', 'PECAM1'],
  'LEC': ['LYVE1', 'CCL21', 'TFF3', 'MMRN1'],
  'ERY': ['HBA1', 'HBB', 'HBM'],
  'LYM': ['TRBC2', 'CD3G', 'CD3E', 'CD3D', 'TRBC1', 'TIGIT'],
  'MYE': ['IL1B', 'LYZ', 'CD207', 'CD1A'],
  'MAST': ['TPSAB1', 'TPSB2'],
  'PC': ['JCHAIN'],
}")
#
dataset_information_df = data.frame(row.names = c("Gaydosik_Fuschiotti_ClinicalCancerResearch_2019", "Liu_Landen_CellStemCell_2024", "Reynolds_Haniffa_Science_2021", "Wang_Atwood_NatureCommunications_2020",
                                                  "Takahashi_Lowry_JournalofInvestigativeDermatology_2020", "Belote_Torres_NatureCellBiology_2021", "Boldo_Lyko_CommunicationsBiology_2020", "Cheng_Cho_CellReports_2018",
                                                  "Ji_Khavari_Cell_2020", "Tabib_Lafyatis_NatureCommunications_2021", "He_Guttman_JournalofAllergyandClinicalImmunology_2020", "Zou_Liu_DevelopmentalCell_2020",
                                                  "Wiedemann_Andersen_CellReports_2023", "Alkon_Stingl_JournalofAllergyandClinicalImmunology_2022", "Dunlap_Rao_JCIInsight_2022", "Boothby_Rosenblum_Nature_2021",
                                                  "Gao_Hu_CellDeathDisease_2021", "Hughes_Shalek_Immunity_2020", "Kim_Nagao_NatureMedicine_2020", "Rindler_Brunner_MolecularCancer_2021",
                                                  "Rojahn_Brunner_JournalofAllergyandClinicalImmunology_2020", "Theocharidis_Bhasin_NatureCommunications_2022", "Vorstandlechner_Mildner_NatureCommunications_2021", "Xu_Chen_Nature_2022"))
dataset_information_df["Gaydosik_Fuschiotti_ClinicalCancerResearch_2019", c("StudyID_short", "PMID", "Title")] = c("Gaydosik 2019", "31010835", "Single-Cell Lymphocyte Heterogeneity in Advanced Cutaneous T-cell Lymphoma Skin Tumors")
dataset_information_df["Liu_Landen_CellStemCell_2024", c("StudyID_short", "PMID", "Title")] = c("Liu 2024", "39729995", "Spatiotemporal single-cell roadmap of human skin wound healing")
dataset_information_df["Reynolds_Haniffa_Science_2021", c("StudyID_short", "PMID", "Title")] = c("Reynolds 2021", "33479125", "Developmental cell programs are co-opted in inflammatory skin disease")
dataset_information_df["Wang_Atwood_NatureCommunications_2020", c("StudyID_short", "PMID", "Title")] = c("Wang 2020", "32843640", "Single cell transcriptomics of human epidermis identifies basal stem cell transition states")
dataset_information_df["Takahashi_Lowry_JournalofInvestigativeDermatology_2020", c("StudyID_short", "PMID", "Title")] = c("Takahashi 2020", "31676413", "Defining Transcriptional Signatures of Human Hair Follicle Cell States")
dataset_information_df["Belote_Torres_NatureCellBiology_2021", c("StudyID_short", "PMID", "Title")] = c("Belote 2021", "34475532", "Human melanocyte development and melanoma dedifferentiation at single-cell resolution")
dataset_information_df["Boldo_Lyko_CommunicationsBiology_2020", c("StudyID_short", "PMID", "Title")] = c("Sole-Boldo 2020", "32327715", "Single-cell transcriptomes of the human skin reveal age-related loss of fibroblast priming")
dataset_information_df["Cheng_Cho_CellReports_2018", c("StudyID_short", "PMID", "Title")] = c("Cheng 2018", "30355494", "Transcriptional Programming of Normal and Inflamed Human Epidermis at Single-Cell Resolution")
dataset_information_df["Ji_Khavari_Cell_2020", c("StudyID_short", "PMID", "Title")] = c("Ji 2020", "32579974", "Multimodal Analysis of Composition and Spatial Architecture in Human Squamous Cell Carcinoma")
dataset_information_df["Tabib_Lafyatis_NatureCommunications_2021", c("StudyID_short", "PMID", "Title")] = c("Tabib 2021", "34282151", "Myofibroblast transcriptome indicates SFRP2hi fibroblast progenitors in systemic sclerosis skin")
dataset_information_df["He_Guttman_JournalofAllergyandClinicalImmunology_2020", c("StudyID_short", "PMID", "Title")] = c("He 2020", "32035984", "Single-cell transcriptome analysis of human skin identifies novel fibroblast subpopulation and enrichment of immune subsets in atopic dermatitis")
dataset_information_df["Zou_Liu_DevelopmentalCell_2020", c("StudyID_short", "PMID", "Title")] = c("Zou 2020", "33238152", "A Single-Cell Transcriptomic Atlas of Human Skin Aging")
dataset_information_df["Wiedemann_Andersen_CellReports_2023", c("StudyID_short", "PMID", "Title")] = c("Wiedemann 2023", "36732947", "Differential cell composition and split epidermal differentiation in human palm, sole, and hip skin")
dataset_information_df["Alkon_Stingl_JournalofAllergyandClinicalImmunology_2022", c("StudyID_short", "PMID", "Title")] = c("Alkon 2022", "34363841", "Single-cell analysis reveals innate lymphoid cell lineage infidelity in atopic dermatitis")
dataset_information_df["Dunlap_Rao_JCIInsight_2022", c("StudyID_short", "PMID", "Title")] = c("Dunlap 2022", "35290245", "Single-cell transcriptomics reveals distinct effector profiles of infiltrating T cells in lupus skin and kidney")
dataset_information_df["Boothby_Rosenblum_Nature_2021", c("StudyID_short", "PMID", "Title")] = c("Boothby 2021", "34707292", "Early-life inflammation primes a T helper 2 cell-fibroblast niche in skin")
dataset_information_df["Gao_Hu_CellDeathDisease_2021", c("StudyID_short", "PMID", "Title")] = c("Gao 2021", "33958582", "Single cell transcriptional zonation of human psoriasis skin identifies an alternative immunoregulatory axis conducted by skin resident cells")
dataset_information_df["Hughes_Shalek_Immunity_2020", c("StudyID_short", "PMID", "Title")] = c("Hughes 2020", "33053333", "Second-Strand Synthesis-Based Massively Parallel scRNA-Seq Reveals Cellular States and Molecular Features of Human Inflammatory Skin Pathologies")
dataset_information_df["Kim_Nagao_NatureMedicine_2020", c("StudyID_short", "PMID", "Title")] = c("Kim 2020", "31959990", "Targeted therapy guided by single-cell transcriptomic analysis in drug-induced hypersensitivity syndrome: a case report")
dataset_information_df["Rindler_Brunner_MolecularCancer_2021", c("StudyID_short", "PMID", "Title")] = c("Rindler 2021", "34583709", "Single-cell RNA sequencing reveals markers of disease progression in primary cutaneous T-cell lymphoma")
dataset_information_df["Rojahn_Brunner_JournalofAllergyandClinicalImmunology_2020", c("StudyID_short", "PMID", "Title")] = c("Rojahn 2020", "32344053", "Single-cell transcriptomics combined with interstitial fluid proteomics defines cell type-specific immune regulation in atopic dermatitis")
dataset_information_df["Theocharidis_Bhasin_NatureCommunications_2022", c("StudyID_short", "PMID", "Title")] = c("Theocharidis 2022", "35013299", "Single cell transcriptomic landscape of diabetic foot ulcers")
dataset_information_df["Vorstandlechner_Mildner_NatureCommunications_2021", c("StudyID_short", "PMID", "Title")] = c("Vorstandlechner 2021", "34716325", "The serine proteases dipeptidyl-peptidase 4 and urokinase are key molecules in human and mouse scar formation")
dataset_information_df["Xu_Chen_Nature_2022", c("StudyID_short", "PMID", "Title")] = c("Xu 2022", "34912121", "Anatomically distinct fibroblast subsets determine skin autoimmune patterns")
#
#
#
for(ii in 1:24){
  adata_str = paste0("adata_", ii)
  this_StudyID = rownames(dataset_information_df)[ii]
  this_StudyID_short = dataset_information_df[ii, "StudyID_short"]
  this_PMID = dataset_information_df[ii, "PMID"]
  this_Title = dataset_information_df[ii, "Title"]
  #
  if(this_StudyID == "Reynolds_Haniffa_Science_2021"){
    py_run_string(paste0(adata_str, " = ad.read_h5ad('", raw_input_paths[["D025"]], "')"))
    py_run_string(paste0("individual_barcode_0 = ", adata_str, ".obs_names.to_series().astype(str).str.split('___', expand=True)[1].str.split('_', n=1, expand=True)[1].values"))
    py_run_string(paste0("individual_barcode_1 = ", adata_str, ".obs_names.to_series().astype(str).str.split('___', expand=True)[1].str.split('_', n=1, expand=True)[0].values"))
    py_run_string(paste0(adata_str, ".obs['individual_barcode'] = individual_barcode_0 + '-' + individual_barcode_1"))
    py_run_string(paste0(adata_str, "_annotation = ad.read_h5ad('", individual_annotation_result_paths[["D025"]], "')[", adata_str, ".obs['individual_barcode']]"))
    py_run_string(paste0(adata_str, "_hvg = ad.read_h5ad('", individual_hvg_clustering_paths[["D025"]], "')[", adata_str, ".obs['individual_barcode']]"))
    py_run_string(paste0(adata_str, ".obsm = {'UMAP': ", adata_str, "_hvg.obsm['X_umap'].copy()}"))
    py_run_string(paste0(adata_str, ".layers = {'raw_counts': ", adata_str, ".X.copy()}"))
    py_run_string(paste0("print(", adata_str, ".layers['raw_counts'])"))
    py_run_string(paste0("sc.pp.normalize_total(", adata_str, ", target_sum=1e4)"))
    py_run_string(paste0("sc.pp.log1p(", adata_str, ")"))
    py_run_string(paste0(adata_str, ".layers['normalized'] = ", adata_str, ".X.copy()"))
    for(jj in leiden_resolutions){
      py_run_string(paste0(adata_str, ".obs['", "individual_leiden_L1_", gsub(pattern = ".", replacement = "", sprintf("%.2f", jj), fixed = T), "'] = ", adata_str, "_annotation.obs['leiden_1_", gsub(pattern = ".", replacement = "", sprintf("%.2f", jj), fixed = T), "'].values"))
    }
    py_run_string(paste0(adata_str, ".obs['individual_annotation_L1'] = ", adata_str, "_annotation.obs['SeparateAnnotation_1'].values"))
  }else{
    py_run_string(paste0(adata_str, " = adata[adata.obs['study_id'] == '", this_StudyID, "'].copy()"))
    py_run_string(paste0(adata_str, ".obsm = {'UMAP': ", adata_str, ".obsm['individual_UMAP_L1'].copy()}"))
  }
  # Take out genes that have no expression in this dataset.
  py_run_string(paste0("sc.pp.filter_genes(", adata_str, ", min_cells=1)"))
  #
  py_run_string(paste0("running_marker_dict = {this_key: np.array(this_value)[np.isin(this_value, ", adata_str, ".var['feature_name'].values)].tolist() for this_key, this_value in marker_dict.items()}"))
  #
  visualization_embedding = "UMAP"
  cell_number = eval(parse(text = paste0("py$", adata_str, "$n_obs")))
  py$dotsize = max(c(8, get_dot_size(cell_number)))
  ######
  py$figure_size = c(8.27, 11.69)
  SW = 210
  SH = 297
  BT = 20
  BB = 17
  BL = 20
  BR = 20
  CW = SW - BL - BR
  CH = SH - BT - BB
  #
  # SF_1
  #
  py_run_string("fig = plt.figure(figsize=figure_size)")
  py_run_string(paste0("fig.text(", 1 - 5/SW, ", ", 1 - 5/SH, ", '", this_StudyID_short, "', fontstyle='normal', fontfamily='Arial', fontsize=18, fontweight='bold', ha='right', va='top', color='black')"))
  py_run_string(paste0("gs = GridSpec(", as.integer(CH), ", ", as.integer(CW), ", figure=fig, left=", BL/SW, ", bottom=", BB/SH, ", right=", 1 - BR/SW, ", top=", 1 - BT/SH, ", wspace=0, hspace=0, width_ratios=None, height_ratios=None)"))
  #
  py_run_string(paste0("fig.text(10/", SW, ", 1 - 20/", SH, ", 'A', fontstyle='normal', fontfamily='Arial', fontsize=18, fontweight='bold', ha='center', va='top', color='black')"))
  py_run_string(paste0("fig.text(", 25/SW, ", ", 1 - 30/SH, ", 'PMID: ", this_PMID, "', fontsize=10, fontweight='normal', ha='left', va='center', color='black')"))
  py_run_string(paste0("fig.text(90/", SW, ", 1 - 20/", SH, ", 'B', fontstyle='normal', fontfamily='Arial', fontsize=18, fontweight='bold', ha='center', va='top', color='black')"))
  py_run_string("ax_b = fig.add_subplot(gs[0:60, 70:130])")
  py_run_string(paste0(adata_str, ".obs['AnnotationName'] = ", adata_str, ".obs['individual_annotation_L1'].map(classify_dict).astype('str').map(classify_dict).map({index: f'{name}' for index, name in rename_dict.items()})"))
  py_run_string(paste0(adata_str, ".obs['AnnotationABBR'] = ", adata_str, ".obs['AnnotationName'].map({index: f'{name}' for index, name in index_name_dict.items()})"))
  py_run_string(paste0(adata_str, ".obs['AnnotationFULL'] = ", adata_str, ".obs['AnnotationName'].map({index: f'{name}: {index}' for index, name in index_name_dict.items()})"))
  py_run_string(paste0(adata_str, ".obs['AnnotationABBR_I'] = ", adata_str, ".obs['AnnotationABBR'].map({abbr: f'{index}:{abbr}' for index, abbr in enumerate(index_name_dict.values())})"))
  py_run_string("order_AnnotationFULL = np.array([f'{name}: {index}' for index, name in index_name_dict.items()])")
  py_run_string("order_AnnotationABBR_I = np.array([f'{index}:{abbr}' for index, abbr in enumerate(index_name_dict.values())])")
  #
  py_run_string(paste0(adata_str, ".obs['AnnotationFULL'] = ", adata_str, ".obs['AnnotationFULL'].astype('category').cat.reorder_categories(order_AnnotationFULL[np.isin(order_AnnotationFULL, ", adata_str, ".obs['AnnotationFULL'])])"))
  py_run_string(paste0(adata_str, ".obs['AnnotationABBR_I'] = ", adata_str, ".obs['AnnotationABBR_I'].astype('category').cat.reorder_categories(order_AnnotationABBR_I[np.isin(order_AnnotationABBR_I, ", adata_str, ".obs['AnnotationABBR_I'])])"))
  #
  py_run_string(paste0("tmp_adata = ad.AnnData(obs=", adata_str, ".obs.copy(), obsm={'", visualization_embedding, "': ", adata_str, ".obsm['", visualization_embedding, "'].copy()})"))
  py_run_string("np.random.seed(0)")
  py_run_string("tmp_adata_random = tmp_adata[np.random.permutation(list(range(tmp_adata.n_obs))), :]")
  py_run_string(paste0("sc.pl.embedding(tmp_adata_random,
                                      basis='", visualization_embedding, "',
                                      color=['AnnotationABBR'],
                                      legend_loc='on data',
                                      legend_fontoutline=0.5,
                                      legend_fontsize=6,
                                      show=False,
                                      ncols=1,
                                      size=0,
                                      frameon=False,
                                      add_outline=False,
                                      ax=ax_b
  )"))
  py_run_string(paste0("sc.pl.embedding(tmp_adata_random,
                                      basis='", visualization_embedding, "',
                                      color=['AnnotationFULL'],
                                      palette=color_dict_CT_FULL,
                                      legend_loc='right margin',
                                      legend_fontsize=8,
                                      show=False,
                                      size=dotsize,
                                      add_outline=False,
                                      ax=ax_b
  )"))
  py_run_string("ax_b.title.set_visible(False)")
  py_run_string("ax_b.get_legend().set_bbox_to_anchor((1, 0.5))")
  py_run_string("del tmp_adata")
  py_run_string("gc.collect()")
  #
  col_number = 4
  WI = 10
  HI = 10
  start_LC = 0
  start_TC = 72
  WH = 35
  this_BL = 0
  py_run_string(paste0("fig.text(10/", SW, ", 1 - ", start_TC + BT - 5, "/", SH, ", 'C', fontstyle='normal', fontfamily='Arial', fontsize=18, fontweight='bold', ha='left', va='center', color='black')"))
  for(MG_I in (seq(length(py$running_marker_dict)) - 1)){
    this_CT = names(py$running_marker_dict)[MG_I + 1]
    this_row = MG_I %/% col_number
    this_col = MG_I %% col_number
    this_LC = start_LC + this_col * (WH + WI)
    this_TC = start_TC + this_row * (WH + HI)
    if(length(py$running_marker_dict[[this_CT]]) > 0){
      py_run_string(paste0("ax = fig.add_subplot(gs[", as.integer(this_TC), ":", as.integer(this_TC + WH), ", ", as.integer(this_BL + this_LC), ":", as.integer(this_BL + this_LC + WH), "])"))
      py_run_string(paste0("this_gene = running_marker_dict['", this_CT, "'][0]"))
      py_run_string(paste0("sc.pl.embedding(", adata_str, ", '", visualization_embedding, "', color=[this_gene], gene_symbols='feature_name', layer='normalized', cmap=gene_highlight_cmap, legend_loc='on data', legend_fontsize=6, size=dotsize, frameon=False, ax=ax)"))
      py_run_string(paste0("ax.set_title('", this_CT, " - ' + this_gene, fontstyle='normal', fontfamily='Arial', fontsize=8)"))
      py_run_string(paste0("ax.collections[0].colorbar.ax.tick_params(labelsize=6)"))
    }else{
      this_fontsize = 8
      py_run_string(paste0("fig.text(", this_LC + 0.5 * WH + BL + this_BL, " / ", SW, ", 1 - ", this_TC + 0.5 * WH + BT - 0.25 * this_fontsize, "/",  SH,", 'No ", this_CT, " marker detected', fontstyle='normal', fontfamily='Arial', fontsize=", this_fontsize, ", fontweight='normal', ha='center', va='center', color='black')"))
    }
  }
  #
  py_run_string(paste0("fig.text(10/", SW, ", ", 1 - (210 + BT)/SH, ", 'D', fontstyle='normal', fontfamily='Arial', fontsize=18, fontweight='bold', ha='left', va='center', color='black')"))
  py_run_string("ax_d = fig.add_subplot(gs[205:250, :])")
  py_run_string("this_dotplot_var_names = running_marker_dict")
  py_run_string(paste0("this_dotplot = sc.pl.DotPlot(", adata_str, ", var_names=this_dotplot_var_names, groupby='AnnotationABBR_I', gene_symbols='feature_name', layer='normalized', standard_scale='var', var_group_rotation=90, ax=ax_d).style(largest_dot=40, cmap='Reds').legend(width=2.5)"))
  py_run_string("this_dotplot_ax = this_dotplot.get_axes()['mainplot_ax']")
  py_run_string("[this_dotplot_ax.axhspan(y, y + 1, facecolor='lightgray', alpha=0.3, zorder=0.5) for y in range(1, len(this_dotplot_ax.get_yticks()), 2)]")
  py_run_string("[this_dotplot_ax.axvline(i, ymin=0, ls=':', lw=0.5, c='dimgray') for i in np.cumsum([len(genes) for genes in this_dotplot_var_names.values()])[:-1]]")
  #  
  py_run_string(paste0("fig.text(", 1 - 5/SW, ", ", 5/SH, ", 'Dataset ID: ", as.integer(ii), "      Page 1 of 4', fontsize=10, fontweight='normal', ha='right', va='bottom', color='black')"))
  #
  py_run_string(paste0("fig.savefig(f'{this_figure_dir}/", adata_str, "_SF_1.pdf')"))
  py_run_string("gc.collect()")
  py_run_string("plt.close('all')")
  #
  # SF_2
  #
  py_run_string("fig = plt.figure(figsize=figure_size)")
  py_run_string(paste0("fig.text(", 1 - 5/SW, ", ", 1 - 5/SH, ", '", this_StudyID_short, "', fontstyle='normal', fontfamily='Arial', fontsize=18, fontweight='bold', ha='right', va='top', color='black')"))
  py_run_string(paste0("gs = GridSpec(", as.integer(CH), ", ", as.integer(CW), ", figure=fig, left=", BL/SW, ", bottom=", BB/SH, ", right=", 1 - BR/SW, ", top=", 1 - BT/SH, ", wspace=0, hspace=0, width_ratios=None, height_ratios=None)"))
  #
  py_run_string(paste0("tmp_adata = ad.AnnData(obs=", adata_str, ".obs.copy(), obsm={'", visualization_embedding, "': ", adata_str, ".obsm['", visualization_embedding, "'].copy()})"))
  py_run_string("np.random.seed(0)")
  py_run_string("tmp_adata_random = tmp_adata[np.random.permutation(list(range(tmp_adata.n_obs))), :].copy()")
  py_run_string("del tmp_adata")
  py_run_string("gc.collect()")
  #
  col_number = 4
  WI = 6
  HI = 10
  start_LC = 0
  start_TC = 0
  WH = 38
  this_BL = 0
  #
  py_run_string(paste0("fig.text(10/", SW, ", 1 - ", start_TC + BT - 3, "/", SH, ", 'E', fontstyle='normal', fontfamily='Arial', fontsize=18, fontweight='bold', ha='left', va='center', color='black')"))
  #
  for(RES_I in (seq(12) - 1)){
    this_RES = leiden_resolutions[RES_I + 1]
    this_row = RES_I %/% col_number
    this_col = RES_I %% col_number
    this_LC = start_LC + this_col * (WH + WI)
    this_TC = start_TC + this_row * (WH + HI)
    #
    py_run_string(paste0("ax = fig.add_subplot(gs[", as.integer(this_TC), ":", as.integer(this_TC + WH), ", ", as.integer(this_BL + this_LC), ":", as.integer(this_BL + this_LC + WH), "])"))
    this_key = paste0("individual_leiden_L1_", gsub(pattern = ".", replacement = "", sprintf("%.2f", this_RES), fixed = T))
    py_run_string(paste0("sc.pl.embedding(tmp_adata_random, '", visualization_embedding, "', color=['", this_key, "'], legend_loc='on data', legend_fontoutline=0.5, legend_fontsize=6, size=dotsize/2, frameon=False, ax=ax)"))
    py_run_string(paste0("ax.set_title('Leiden Res ", sprintf("%.2f", this_RES), "', fontstyle='normal', fontfamily='Arial', fontsize=8)"))
  }
  col_number = 3
  WI = 7
  HI = 11
  start_LC = 0
  start_TC = this_TC + WH + HI
  WH = 52
  this_BL = 0
  #
  py_run_string(paste0("fig.text(10/", SW, ", 1 - ", start_TC + BT - 3, "/", SH, ", 'F', fontstyle='normal', fontfamily='Arial', fontsize=18, fontweight='bold', ha='left', va='center', color='black')"))
  #
  for(RES_I in (seq(from=1 + 12, to=length(leiden_resolutions)) - 1 - 12)){
    this_RES = leiden_resolutions[RES_I + 1 + 12]
    this_row = RES_I %/% col_number
    this_col = RES_I %% col_number
    this_LC = start_LC + this_col * (WH + WI)
    this_TC = start_TC + this_row * (WH + HI)
    #
    py_run_string(paste0("ax = fig.add_subplot(gs[", as.integer(this_TC), ":", as.integer(this_TC + WH), ", ", as.integer(this_BL + this_LC), ":", as.integer(this_BL + this_LC + WH), "])"))
    this_key = paste0("individual_leiden_L1_", gsub(pattern = ".", replacement = "", sprintf("%.2f", this_RES), fixed = T))
    py_run_string(paste0("sc.pl.embedding(tmp_adata_random, '", visualization_embedding, "', color=['", this_key, "'], legend_loc='on data', legend_fontoutline=0.5, legend_fontsize=6, size=dotsize/2, frameon=False, ax=ax)"))
    py_run_string(paste0("ax.set_title('Leiden Res ", sprintf("%.2f", this_RES), "', fontstyle='normal', fontfamily='Arial', fontsize=8)"))
  }
  #
  if("author_cell_type_1" %in% colnames(py$tmp_adata_random$obs)){
    this_row = (RES_I + 1) %/% col_number
    this_col = (RES_I + 1) %% col_number
    this_LC = start_LC + this_col * (WH + WI)
    this_TC = start_TC + this_row * (WH + HI)
    #
    this_OA = as.character(py$tmp_adata_random$obs[['author_cell_type_1']])
    this_OA[tolower(this_OA) %in% tolower(na_value)] = NA
    py_run_string("del tmp_adata_random.obs['author_cell_type_1']")
    py$tmp_adata_random$obs[["author_cell_type_1"]] = as.factor(this_OA)
    if(sum(!is.na(this_OA)) > 0){
      py_run_string(paste0("ax = fig.add_subplot(gs[", as.integer(this_TC), ":", as.integer(this_TC + WH), ", ", as.integer(this_BL + this_LC), ":", as.integer(this_BL + this_LC + WH), "])"))
      py_run_string(paste0("sc.pl.embedding(tmp_adata_random, '", visualization_embedding, "', color=['author_cell_type_1'], legend_loc='on data', legend_fontoutline=0.5, legend_fontsize=6, size=dotsize/2, frameon=False, ax=ax)"))
      py_run_string("ax.set_title('Original Annotation', fontstyle='normal', fontfamily='Arial', fontsize=8)")
    }
    rm(this_OA)
  }
  py_run_string("del tmp_adata_random")
  py_run_string("gc.collect()")
  #
  #
  #
  py_run_string(paste0("fig.text(", 1 - 5/SW, ", ", 5/SH, ", 'Dataset ID: ", as.integer(ii), "      Page 2 of 4', fontsize=10, fontweight='normal', ha='right', va='bottom', color='black')"))
  #
  py_run_string(paste0("fig.savefig(f'{this_figure_dir}/", adata_str, "_SF_2.pdf')"))
  py_run_string("gc.collect()")
  py_run_string("plt.close('all')")
  #
  # SF_3_1
  #
  py_run_string("fig = plt.figure(figsize=figure_size)")
  py_run_string(paste0("fig.text(", 1 - 5/SW, ", ", 1 - 5/SH, ", '", this_StudyID_short, "', fontstyle='normal', fontfamily='Arial', fontsize=18, fontweight='bold', ha='right', va='top', color='black')"))
  py_run_string(paste0("gs = GridSpec(", as.integer(CH), ", ", as.integer(CW), ", figure=fig, left=", BL/SW, ", bottom=", BB/SH, ", right=", 1 - BR/SW, ", top=", 1 - BT/SH, ", wspace=0, hspace=0, width_ratios=None, height_ratios=None)"))
  ######
  gene_names = eval(parse(text = paste0("py$", adata_str, "$var[['feature_name']]")))
  #
  this_CT = "KC"
  this_CT_title = "Keratinocyte"
  plot_genes = py$marker_dict[[this_CT]]
  main_genes = py$running_marker_dict[[this_CT]][1]
  #
  col_number = 6
  WI = 4
  HI = 10
  start_LC = 0
  start_TC = 10
  WH = 25
  #
  py_run_string(paste0("fig.text(10/", SW, ", 1 - ", start_TC + BT - 3, "/", SH, ", 'G', fontstyle='normal', fontfamily='Arial', fontsize=18, fontweight='bold', ha='left', va='center', color='black')"))
  #
  has_incomplete_row = length(plot_genes) %% col_number > 0
  expected_rows = length(plot_genes) %/% col_number + as.integer(has_incomplete_row)
  expected_EBC = start_TC + (WH + HI) * expected_rows
  py_run_string(paste0("fig.text(", 17/SW, ", 1 - ", 0.5 * (start_TC + expected_EBC) + BT - 5, "/",  SH,", '", this_CT_title, "', rotation='vertical', fontstyle='normal', fontfamily='Arial', fontsize=10, fontweight='bold', ha='center', va='center', color='black')"))
  if(sum(gene_names %in% plot_genes) > 0){
    for(MG_I in (seq(length(plot_genes)) - 1)){
      this_row = MG_I %/% col_number
      this_col = MG_I %% col_number
      this_LC = start_LC + this_col * (WH + WI)
      this_TC = start_TC + this_row * (WH + HI)
      this_gene = plot_genes[MG_I + 1]
      if(this_gene %in% gene_names){
        py_run_string(paste0("ax = fig.add_subplot(gs[", as.integer(this_TC), ":", as.integer(this_TC + WH), ", ", as.integer(this_LC), ":", as.integer(this_LC + WH), "])"))
        py_run_string(paste0("sc.pl.embedding(", adata_str, ", '", visualization_embedding, "', color=['", this_gene, "'], gene_symbols='feature_name', layer='normalized', cmap=gene_highlight_cmap, legend_loc='on data', legend_fontsize=6, size=dotsize, frameon=False, ax=ax)"))
        if(this_gene %in% main_genes){
          py_run_string(paste0("ax.set_title('", this_gene, "', fontstyle='normal', fontfamily='Arial', fontsize=8, fontweight='bold', ha='center', va='center', color='red')"))
        }else{
          py_run_string(paste0("ax.set_title('", this_gene, "', fontstyle='normal', fontfamily='Arial', fontsize=8, fontweight='normal', ha='center', va='center', color='black')"))
        }
        py_run_string(paste0("ax.collections[0].colorbar.ax.tick_params(labelsize=6)"))
      }else{
        this_fontsize = 8
        py_run_string(paste0("fig.text(", this_LC + BL + 0.5 * WH, " / ", SW, ", 1 - ", this_TC + BT + 0.5 * WH - 0.25 * this_fontsize, "/",  SH,", 'No ", this_gene, "', fontstyle='normal', fontfamily='Arial', fontsize=", this_fontsize, ", fontweight='normal', ha='center', va='center', color='black')"))
      }
    }
  }else{
    this_fontsize = 18
    py_run_string(paste0("fig.text(", 0.5 * (SW + BL), " / ", SW, ", 1 - ", 0.5 * (start_TC + expected_EBC) + BT - 0.25 * this_fontsize, "/",  SH,", 'No ", this_CT_title, " marker detected', fontstyle='normal', fontfamily='Arial', fontsize=", this_fontsize, ", fontweight='bold', ha='center', va='center', color='black')"))
  }
  #
  #
  #
  this_CT = c("MAST", "MYE", "LYM")
  this_CT_title = "Immune Cell"
  plot_genes = c(py$marker_dict[["MAST"]], py$marker_dict[["MYE"]], py$marker_dict[["LYM"]])
  main_genes = c(py$running_marker_dict[["MAST"]][1], py$running_marker_dict[["MYE"]][1], py$running_marker_dict[["LYM"]][1])
  col_number = 6
  WI = 4
  HI = 10
  start_LC = 0
  start_TC = expected_EBC + 20
  WH = 25
  this_BL = 0
  #
  py_run_string(paste0("fig.text(10/", SW, ", 1 - ", start_TC + BT - 3, "/", SH, ", 'H', fontstyle='normal', fontfamily='Arial', fontsize=18, fontweight='bold', ha='left', va='center', color='black')"))
  #
  has_incomplete_row = length(plot_genes) %% col_number > 0
  expected_rows = length(plot_genes) %/% col_number + as.integer(has_incomplete_row)
  expected_EBC = start_TC + (WH + HI) * expected_rows
  py_run_string(paste0("fig.text(", 17/SW, ", 1 - ", 0.5 * (start_TC + expected_EBC) + BT - 5, "/",  SH,", '", this_CT_title, "', rotation='vertical', fontstyle='normal', fontfamily='Arial', fontsize=10, fontweight='bold', ha='center', va='center', color='black')"))
  if(sum(gene_names %in% plot_genes) > 0){
    for(MG_I in (seq(length(plot_genes)) - 1)){
      this_row = MG_I %/% col_number
      this_col = MG_I %% col_number
      this_LC = start_LC + this_col * (WH + WI)
      this_TC = start_TC + this_row * (WH + HI)
      this_gene = plot_genes[MG_I + 1]
      if(this_gene %in% gene_names){
        py_run_string(paste0("ax = fig.add_subplot(gs[", as.integer(this_TC), ":", as.integer(this_TC + WH), ", ", as.integer(this_BL + this_LC), ":", as.integer(this_BL + this_LC + WH), "])"))
        py_run_string(paste0("sc.pl.embedding(", adata_str, ", '", visualization_embedding, "', color=['", this_gene, "'], gene_symbols='feature_name', layer='normalized', cmap=gene_highlight_cmap, legend_loc='on data', legend_fontsize=6, size=dotsize, frameon=False, ax=ax)"))
        if(this_gene %in% main_genes){
          py_run_string(paste0("ax.set_title('", this_gene, "', fontstyle='normal', fontfamily='Arial', fontsize=8, fontweight='bold', ha='center', va='center', color='red')"))
        }else{
          py_run_string(paste0("ax.set_title('", this_gene, "', fontstyle='normal', fontfamily='Arial', fontsize=8, fontweight='normal', ha='center', va='center', color='black')"))
        }
        py_run_string(paste0("ax.collections[0].colorbar.ax.tick_params(labelsize=6)"))
      }else{
        this_fontsize = 8
        py_run_string(paste0("fig.text(", this_LC + BL + 0.5 * WH, " / ", SW, ", 1 - ", this_TC + BT + 0.5 * WH - 0.25 * this_fontsize, "/",  SH,", 'No ", this_gene, "', fontstyle='normal', fontfamily='Arial', fontsize=", this_fontsize, ", fontweight='normal', ha='center', va='center', color='black')"))
      }
    }
  }else{
    this_fontsize = 18
    py_run_string(paste0("fig.text(", 0.5 * (SW + BL), " / ", SW, ", 1 - ", 0.5 * (start_TC + expected_EBC) + BT - 0.25 * this_fontsize, "/",  SH,", 'No ", this_CT_title, " marker detected', fontstyle='normal', fontfamily='Arial', fontsize=", this_fontsize, ", fontweight='bold', ha='center', va='center', color='black')"))
  }
  #
  #
  #
  this_CT = c("ERY", "PC")
  this_CT_title = "ERY & PC"
  plot_genes = c(py$marker_dict[["ERY"]], py$marker_dict[["PC"]])
  main_genes = c(py$running_marker_dict[["ERY"]][1], py$running_marker_dict[["PC"]][1])
  col_number = 4
  WI = 4
  HI = 10
  start_LC = 29
  start_TC = expected_EBC + 20
  WH = 25
  #
  py_run_string(paste0("fig.text(10/", SW, ", 1 - ", start_TC + BT - 3, "/", SH, ", 'I', fontstyle='normal', fontfamily='Arial', fontsize=18, fontweight='bold', ha='left', va='center', color='black')"))
  #
  has_incomplete_row = length(plot_genes) %% col_number > 0
  expected_rows = length(plot_genes) %/% col_number + as.integer(has_incomplete_row)
  expected_EBC = start_TC + (WH + HI) * expected_rows
  py_run_string(paste0("fig.text(", 17/SW, ", 1 - ", 0.5 * (start_TC + expected_EBC) + BT - 5, "/",  SH,", '", this_CT_title, "', rotation='vertical', fontstyle='normal', fontfamily='Arial', fontsize=10, fontweight='bold', ha='center', va='center', color='black')"))
  if(sum(gene_names %in% plot_genes) > 0){
    for(MG_I in (seq(length(plot_genes)) - 1)){
      this_row = MG_I %/% col_number
      this_col = MG_I %% col_number
      this_LC = start_LC + this_col * (WH + WI)
      this_TC = start_TC + this_row * (WH + HI)
      this_gene = plot_genes[MG_I + 1]
      if(this_gene %in% gene_names){
        py_run_string(paste0("ax = fig.add_subplot(gs[", as.integer(this_TC), ":", as.integer(this_TC + WH), ", ", as.integer(this_LC), ":", as.integer(this_LC + WH), "])"))
        py_run_string(paste0("sc.pl.embedding(", adata_str, ", '", visualization_embedding, "', color=['", this_gene, "'], gene_symbols='feature_name', layer='normalized', cmap=gene_highlight_cmap, legend_loc='on data', legend_fontsize=6, size=dotsize, frameon=False, ax=ax)"))
        if(this_gene %in% main_genes){
          py_run_string(paste0("ax.set_title('", this_gene, "', fontstyle='normal', fontfamily='Arial', fontsize=8, fontweight='bold', ha='center', va='center', color='red')"))
        }else{
          py_run_string(paste0("ax.set_title('", this_gene, "', fontstyle='normal', fontfamily='Arial', fontsize=8, fontweight='normal', ha='center', va='center', color='black')"))
        }
        py_run_string(paste0("ax.collections[0].colorbar.ax.tick_params(labelsize=6)"))
      }else{
        this_fontsize = 8
        py_run_string(paste0("fig.text(", this_LC + BL + 0.5 * WH, " / ", SW, ", 1 - ", this_TC + BT + 0.5 * WH - 0.25 * this_fontsize, "/",  SH,", 'No ", this_gene, "', fontstyle='normal', fontfamily='Arial', fontsize=", this_fontsize, ", fontweight='normal', ha='center', va='center', color='black')"))
      }
    }
  }else{
    this_fontsize = 18
    py_run_string(paste0("fig.text(", 0.5 * (SW + BL), " / ", SW, ", 1 - ", 0.5 * (start_TC + expected_EBC) + BT - 0.25 * this_fontsize, "/",  SH,", 'No ", this_CT_title, " marker detected', fontstyle='normal', fontfamily='Arial', fontsize=", this_fontsize, ", fontweight='bold', ha='center', va='center', color='black')"))
  }
  #
  py_run_string(paste0("fig.text(", 1 - 5/SW, ", ", 5/SH, ", 'Dataset ID: ", as.integer(ii), "      Page 3 of 4', fontsize=10, fontweight='normal', ha='right', va='bottom', color='black')"))
  #
  py_run_string(paste0("fig.savefig(f'{this_figure_dir}/", adata_str, "_SF_3_1.pdf')"))
  py_run_string("gc.collect()")
  py_run_string("plt.close('all')")
  #
  # SF_3_2
  #
  py_run_string("fig = plt.figure(figsize=figure_size)")
  py_run_string(paste0("fig.text(", 1 - 5/SW, ", ", 1 - 5/SH, ", '", this_StudyID_short, "', fontstyle='normal', fontfamily='Arial', fontsize=18, fontweight='bold', ha='right', va='top', color='black')"))
  py_run_string(paste0("gs = GridSpec(", as.integer(CH), ", ", as.integer(CW), ", figure=fig, left=", BL/SW, ", bottom=", BB/SH, ", right=", 1 - BR/SW, ", top=", 1 - BT/SH, ", wspace=0, hspace=0, width_ratios=None, height_ratios=None)"))
  ######
  gene_names = eval(parse(text = paste0("py$", adata_str, "$var[['feature_name']]")))
  #
  this_CT = c("FIB", "MUR")
  this_CT_title = "Fibroblast & Mural"
  plot_genes = c(py$marker_dict[["FIB"]], py$marker_dict[["MUR"]])
  main_genes = c(py$running_marker_dict[["FIB"]][1], py$running_marker_dict[["MUR"]][1])
  #
  col_number = 4
  WI = 4
  HI = 10
  start_LC = 29
  start_TC = 0
  WH = 25
  #
  py_run_string(paste0("fig.text(10/", SW, ", 1 - ", start_TC + BT - 3, "/", SH, ", 'J', fontstyle='normal', fontfamily='Arial', fontsize=18, fontweight='bold', ha='left', va='center', color='black')"))
  #
  has_incomplete_row = length(plot_genes) %% col_number > 0
  expected_rows = length(plot_genes) %/% col_number + as.integer(has_incomplete_row)
  expected_EBC = start_TC + (WH + HI) * expected_rows
  py_run_string(paste0("fig.text(", 17/SW, ", 1 - ", 0.5 * (start_TC + expected_EBC) + BT - 5, "/",  SH,", '", this_CT_title, "', rotation='vertical', fontstyle='normal', fontfamily='Arial', fontsize=10, fontweight='bold', ha='center', va='center', color='black')"))
  if(sum(gene_names %in% plot_genes) > 0){
    for(MG_I in (seq(length(plot_genes)) - 1)){
      this_row = MG_I %/% col_number
      this_col = MG_I %% col_number
      this_LC = start_LC + this_col * (WH + WI)
      this_TC = start_TC + this_row * (WH + HI)
      this_gene = plot_genes[MG_I + 1]
      if(this_gene %in% gene_names){
        py_run_string(paste0("ax = fig.add_subplot(gs[", as.integer(this_TC), ":", as.integer(this_TC + WH), ", ", as.integer(this_LC), ":", as.integer(this_LC + WH), "])"))
        py_run_string(paste0("sc.pl.embedding(", adata_str, ", '", visualization_embedding, "', color=['", this_gene, "'], gene_symbols='feature_name', layer='normalized', cmap=gene_highlight_cmap, legend_loc='on data', legend_fontsize=6, size=dotsize, frameon=False, ax=ax)"))
        if(this_gene %in% main_genes){
          py_run_string(paste0("ax.set_title('", this_gene, "', fontstyle='normal', fontfamily='Arial', fontsize=8, fontweight='bold', ha='center', va='center', color='red')"))
        }else{
          py_run_string(paste0("ax.set_title('", this_gene, "', fontstyle='normal', fontfamily='Arial', fontsize=8, fontweight='normal', ha='center', va='center', color='black')"))
        }
        py_run_string(paste0("ax.collections[0].colorbar.ax.tick_params(labelsize=6)"))
      }else{
        this_fontsize = 8
        py_run_string(paste0("fig.text(", this_LC + BL + 0.5 * WH, " / ", SW, ", 1 - ", this_TC + BT + 0.5 * WH - 0.25 * this_fontsize, "/",  SH,", 'No ", this_gene, "', fontstyle='normal', fontfamily='Arial', fontsize=", this_fontsize, ", fontweight='normal', ha='center', va='center', color='black')"))
      }
    }
  }else{
    this_fontsize = 18
    py_run_string(paste0("fig.text(", 0.5 * (SW + BL), " / ", SW, ", 1 - ", 0.5 * (start_TC + expected_EBC) + BT - 0.25 * this_fontsize, "/",  SH,", 'No ", this_CT_title, " marker detected', fontstyle='normal', fontfamily='Arial', fontsize=", this_fontsize, ", fontweight='bold', ha='center', va='center', color='black')"))
  }
  #
  #
  #
  this_CT = c("MEL", "SC")
  this_CT_title = "Neural crest-derived"
  plot_genes = c(py$marker_dict[["MEL"]], py$marker_dict[["SC"]])
  main_genes = c(py$running_marker_dict[["MEL"]][1], py$running_marker_dict[["SC"]][1])
  col_number = 4
  WI = 4
  HI = 10
  start_LC = 29
  start_TC = expected_EBC + 10
  WH = 25
  this_BL = 0
  #
  py_run_string(paste0("fig.text(10/", SW, ", 1 - ", start_TC + BT - 3, "/", SH, ", 'K', fontstyle='normal', fontfamily='Arial', fontsize=18, fontweight='bold', ha='left', va='center', color='black')"))
  #
  has_incomplete_row = length(plot_genes) %% col_number > 0
  expected_rows = length(plot_genes) %/% col_number + as.integer(has_incomplete_row)
  expected_EBC = start_TC + (WH + HI) * expected_rows
  py_run_string(paste0("fig.text(", 17/SW, ", 1 - ", 0.5 * (start_TC + expected_EBC) + BT - 5, "/",  SH,", '", this_CT_title, "', rotation='vertical', fontstyle='normal', fontfamily='Arial', fontsize=10, fontweight='bold', ha='center', va='center', color='black')"))
  if(sum(gene_names %in% plot_genes) > 0){
    for(MG_I in (seq(length(plot_genes)) - 1)){
      this_row = MG_I %/% col_number
      this_col = MG_I %% col_number
      this_LC = start_LC + this_col * (WH + WI)
      this_TC = start_TC + this_row * (WH + HI)
      this_gene = plot_genes[MG_I + 1]
      if(this_gene %in% gene_names){
        py_run_string(paste0("ax = fig.add_subplot(gs[", as.integer(this_TC), ":", as.integer(this_TC + WH), ", ", as.integer(this_BL + this_LC), ":", as.integer(this_BL + this_LC + WH), "])"))
        py_run_string(paste0("sc.pl.embedding(", adata_str, ", '", visualization_embedding, "', color=['", this_gene, "'], gene_symbols='feature_name', layer='normalized', cmap=gene_highlight_cmap, legend_loc='on data', legend_fontsize=6, size=dotsize, frameon=False, ax=ax)"))
        if(this_gene %in% main_genes){
          py_run_string(paste0("ax.set_title('", this_gene, "', fontstyle='normal', fontfamily='Arial', fontsize=8, fontweight='bold', ha='center', va='center', color='red')"))
        }else{
          py_run_string(paste0("ax.set_title('", this_gene, "', fontstyle='normal', fontfamily='Arial', fontsize=8, fontweight='normal', ha='center', va='center', color='black')"))
        }
        py_run_string(paste0("ax.collections[0].colorbar.ax.tick_params(labelsize=6)"))
      }else{
        this_fontsize = 8
        py_run_string(paste0("fig.text(", this_LC + BL + 0.5 * WH, " / ", SW, ", 1 - ", this_TC + BT + 0.5 * WH - 0.25 * this_fontsize, "/",  SH,", 'No ", this_gene, "', fontstyle='normal', fontfamily='Arial', fontsize=", this_fontsize, ", fontweight='normal', ha='center', va='center', color='black')"))
      }
    }
  }else{
    this_fontsize = 18
    py_run_string(paste0("fig.text(", 0.5 * (SW + BL), " / ", SW, ", 1 - ", 0.5 * (start_TC + expected_EBC) + BT - 0.25 * this_fontsize, "/",  SH,", 'No ", this_CT_title, " marker detected', fontstyle='normal', fontfamily='Arial', fontsize=", this_fontsize, ", fontweight='bold', ha='center', va='center', color='black')"))
  }
  #
  #
  #
  this_CT = c("VEC", "LEC")
  this_CT_title = "Endothelial Cell"
  plot_genes = c(py$marker_dict[["VEC"]], py$marker_dict[["LEC"]])
  main_genes = c(py$running_marker_dict[["VEC"]][1], py$running_marker_dict[["LEC"]][1])
  col_number = 4
  WI = 4
  HI = 10
  start_LC = 29
  start_TC = expected_EBC + 10
  WH = 25
  #
  py_run_string(paste0("fig.text(10/", SW, ", 1 - ", start_TC + BT - 3, "/", SH, ", 'L', fontstyle='normal', fontfamily='Arial', fontsize=18, fontweight='bold', ha='left', va='center', color='black')"))
  #
  has_incomplete_row = length(plot_genes) %% col_number > 0
  expected_rows = length(plot_genes) %/% col_number + as.integer(has_incomplete_row)
  expected_EBC = start_TC + (WH + HI) * expected_rows
  py_run_string(paste0("fig.text(", 17/SW, ", 1 - ", 0.5 * (start_TC + expected_EBC) + BT - 5, "/",  SH,", '", this_CT_title, "', rotation='vertical', fontstyle='normal', fontfamily='Arial', fontsize=10, fontweight='bold', ha='center', va='center', color='black')"))
  if(sum(gene_names %in% plot_genes) > 0){
    for(MG_I in (seq(length(plot_genes)) - 1)){
      this_row = MG_I %/% col_number
      this_col = MG_I %% col_number
      this_LC = start_LC + this_col * (WH + WI)
      this_TC = start_TC + this_row * (WH + HI)
      this_gene = plot_genes[MG_I + 1]
      if(this_gene %in% gene_names){
        py_run_string(paste0("ax = fig.add_subplot(gs[", as.integer(this_TC), ":", as.integer(this_TC + WH), ", ", as.integer(this_LC), ":", as.integer(this_LC + WH), "])"))
        py_run_string(paste0("sc.pl.embedding(", adata_str, ", '", visualization_embedding, "', color=['", this_gene, "'], gene_symbols='feature_name', layer='normalized', cmap=gene_highlight_cmap, legend_loc='on data', legend_fontsize=6, size=dotsize, frameon=False, ax=ax)"))
        if(this_gene %in% main_genes){
          py_run_string(paste0("ax.set_title('", this_gene, "', fontstyle='normal', fontfamily='Arial', fontsize=8, fontweight='bold', ha='center', va='center', color='red')"))
        }else{
          py_run_string(paste0("ax.set_title('", this_gene, "', fontstyle='normal', fontfamily='Arial', fontsize=8, fontweight='normal', ha='center', va='center', color='black')"))
        }
        py_run_string(paste0("ax.collections[0].colorbar.ax.tick_params(labelsize=6)"))
      }else{
        this_fontsize = 8
        py_run_string(paste0("fig.text(", this_LC + BL + 0.5 * WH, " / ", SW, ", 1 - ", this_TC + BT + 0.5 * WH - 0.25 * this_fontsize, "/",  SH,", 'No ", this_gene, "', fontstyle='normal', fontfamily='Arial', fontsize=", this_fontsize, ", fontweight='normal', ha='center', va='center', color='black')"))
      }
    }
  }else{
    this_fontsize = 18
    py_run_string(paste0("fig.text(", 0.5 * (SW + BL), " / ", SW, ", 1 - ", 0.5 * (start_TC + expected_EBC) + BT - 0.25 * this_fontsize, "/",  SH,", 'No ", this_CT_title, " marker detected', fontstyle='normal', fontfamily='Arial', fontsize=", this_fontsize, ", fontweight='bold', ha='center', va='center', color='black')"))
  }
  #
  py_run_string(paste0("fig.text(", 1 - 5/SW, ", ", 5/SH, ", 'Dataset ID: ", as.integer(ii), "      Page 4 of 4', fontsize=10, fontweight='normal', ha='right', va='bottom', color='black')"))
  #
  py_run_string(paste0("fig.savefig(f'{this_figure_dir}/", adata_str, "_SF_3_2.pdf')"))
  py_run_string("gc.collect()")
  py_run_string("plt.close('all')")
  #
  #
  #
  print(ii)
}
#
#
#
py$figure_size = c(8.27, 11.69)
SW = 210
SH = 297
BT = 20
BB = 17
BL = 20
BR = 20
CW = SW - BL - BR
CH = SH - BT - BB
#
py_run_string("fig = plt.figure(figsize=figure_size)")
py_run_string(paste0("gs = GridSpec(", as.integer(CH), ", ", as.integer(CW), ", figure=fig, left=", BL/SW, ", bottom=", BB/SH, ", right=", 1 - BR/SW, ", top=", 1 - BT/SH, ", wspace=0, hspace=0, width_ratios=None, height_ratios=None)"))
#
col_number = 6
WI = 4
HI = 10
start_LC = 0
start_TC = 0
WH = 25
this_BL = 0
#
for(ii in 1:24){
  adata_str = paste0("adata_", ii)
  visualization_embedding = "UMAP"
  cell_number = eval(parse(text = paste0("py$", adata_str, "$n_obs")))
  py$dotsize = max(c(8, get_dot_size(cell_number))) / col_number
  ######
  actual_index = ii - 1
  this_row = actual_index %/% col_number
  this_col = actual_index %% col_number
  this_LC = start_LC + this_col * (WH + WI)
  this_TC = start_TC + this_row * (WH + HI)
  #
  py_run_string(paste0("ax = fig.add_subplot(gs[", as.integer(this_TC), ":", as.integer(this_TC + WH), ", ", as.integer(this_BL + this_LC), ":", as.integer(this_BL + this_LC + WH), "])"))
  py_run_string(paste0(adata_str, ".obs['AnnotationName'] = ", adata_str, ".obs['individual_annotation_L1'].map(classify_dict).astype('str').map(classify_dict).map({index: f'{name}' for index, name in rename_dict.items()})"))
  py_run_string(paste0(adata_str, ".obs['AnnotationFULL'] = ", adata_str, ".obs['AnnotationName'].map({index: f'{name}: {index}' for index, name in index_name_dict.items()})"))
  py_run_string("order_AnnotationFULL = np.array([f'{name}: {index}' for index, name in index_name_dict.items()])")
  #
  py_run_string(paste0(adata_str, ".obs['AnnotationFULL'] = ", adata_str, ".obs['AnnotationFULL'].astype('category').cat.reorder_categories(order_AnnotationFULL[np.isin(order_AnnotationFULL, ", adata_str, ".obs['AnnotationFULL'])])"))
  #
  py_run_string(paste0("tmp_adata = ad.AnnData(obs=", adata_str, ".obs.copy(), obsm={'", visualization_embedding, "': ", adata_str, ".obsm['", visualization_embedding, "'].copy()})"))
  py_run_string("np.random.seed(0)")
  py_run_string("tmp_adata_random = tmp_adata[np.random.permutation(list(range(tmp_adata.n_obs))), :]")
  py_run_string("del tmp_adata")
  py_run_string("gc.collect()")
  #
  py_run_string(paste0("sc.pl.embedding(tmp_adata_random,
                                        basis='", visualization_embedding, "',
                                        color=['AnnotationFULL'],
                                        palette=color_dict_CT_FULL,
                                        legend_loc='none',
                                        legend_fontoutline=0.5,
                                        legend_fontsize=6,
                                        show=False,
                                        size=dotsize,
                                        frameon=False,
                                        add_outline=False,
                                        ax=ax
  )")) 
  py_run_string("ax.title.set_visible(False)")
  #
  #
  #
  py_run_string("del tmp_adata_random")
  py_run_string("gc.collect()")
  #
  print(ii)
}
py_run_string(paste0("fig.savefig(f'{this_figure_dir}/ALL_UMAP_4x6.pdf')"))
py_run_string("plt.close('all')")
py_run_string("gc.collect()")
#
metadata_fields = c('experiment_id', 'donor_id_running', 'sex_ontology_term', 'ethnicity_1', 'age_range', 'sample_id_running', 'anatomical_region', 'skin_tissue', 'sample_preservation_method', 'assay_ontology_term', 'author_cell_type_1', 'author_cell_type_2', 'author_cell_type_3')
metadata_status_matrix = matrix(nrow = 0, ncol = length(metadata_fields), dimnames = list(c(), metadata_fields))
for(ii in 1:24){
  adata_str = paste0("adata_", ii)
  this_metadata_status = rep(NA, length(metadata_fields))
  names(this_metadata_status) = metadata_fields
  for(jj in metadata_fields){
    py_run_string(paste0("tmp_this_metadata = ", adata_str, ".obs['", jj, "'].astype(str).copy()"))
    py_run_string(paste0("tmp_this_cell_number = ", adata_str, ".n_obs"))
    this_metadata_status[jj] = sum(py$tmp_this_metadata != "nan") / py$tmp_this_cell_number
    py_run_string("del tmp_this_metadata")
    py_run_string("del tmp_this_cell_number")
  }
  metadata_status_matrix = rbind(metadata_status_matrix, this_metadata_status)
  print(ii)
}
rownames(metadata_status_matrix) = dataset_information_df[, "StudyID_short"]
print(metadata_status_matrix)
write.table(metadata_status_matrix, paste0(py$this_figure_dir, "/metadata_status_matrix.txt"), sep = "\t", quote = F, row.names = TRUE, col.names = TRUE)
#
#
#
GeneName_key = "feature_name"
adata_str_in = adata_str
adata_str_out = adata_str

GeneID_to_GeneName = function(adata_str_in, adata_str_out, GeneName_key, verbose=T){
  GeneName = eval(parse(text = paste0("py$", adata_str_in, "$var[['", GeneName_key, "']]")))
  py$tmp_unique_index = which(!duplicated(GeneName)) - 1
  py_run_string(paste0("tmp_X = ", adata_str_in, ".X[:, tmp_unique_index].copy()"))
  GeneName_duplicated = unique(GeneName[duplicated(GeneName)])
  GeneName_unique = GeneName[!duplicated(GeneName)]
  for(ii in GeneName_duplicated){
    py$tmp_index_old = which(GeneName == ii) - 1
    py$tmp_index_new = which(GeneName_unique == ii) - 1
    py_run_string(paste0("tmp_X[:, np.array(tmp_index_new, int)] = ", adata_str_in, ".X[:, np.array(tmp_index_old, int)].sum(1)"))
    if(verbose){
      print(paste0(which(GeneName_duplicated == ii), " - ", length(GeneName_duplicated)))
    }
  }
  py$tmp_var_names = GeneName_unique
  py_run_string(paste0(adata_str_out, " = ", adata_str_in, "[:, tmp_unique_index].copy()"))
  py_run_string(paste0(adata_str_out, ".var_names = tmp_var_names"))
  py_run_string(paste0(adata_str_out, ".X = tmp_X.copy()"))
  py_run_string("del tmp_unique_index")
  py_run_string("del tmp_var_names")
  py_run_string("del tmp_X")
  py_run_string("gc.collect()")
}
#
#
#
py_run_string("new_colnames = {'study_id': 'Study ID',
                               'experiment_id': 'Experiment ID',
                               'donor_id_running_short': 'Donor ID',
                               'sex_ontology_term': 'Sex',
                               'ethnicity_1': 'Ethnicity 1',
                               'ethnicity_2': 'Ethnicity 2',
                               'ethnicity_details': 'Ethnicity details',
                               'self_reported_ethnicity_ontology_term_id': 'Self reported ethnicity (ontology term ID)',
                               'genotype': 'Genotype',
                               'age_years': 'Age (years)',
                               'age_range': 'Age (range)',
                               'development_stage_ontology_term': 'Development stage (ontology term)',
                               'development_stage_ontology_term_id': 'Development stage (ontology term ID)',
                               'disease': 'Disease',
                               'disease_ontology_term_id': 'Disease (ontology term ID)',
                               'disease_status': 'Disease status',
                               'treatment_status': 'Treatment status',
                               'bmi': 'BMI',
                               'skin_tone': 'Skin tone',
                               'manner_of_death': 'Manner of death',
                               'sample_id_running_short': 'Sample ID',
                               'sample_source': 'Sample source',
                               'sample_collection_method': 'Sample collection method',
                               'sampled_site_condition': 'Sampled site condition',
                               'anatomical_region_level_1': 'Anatomical region level 1',
                               'anatomical_region_level_2': 'Anatomical region level 2',
                               'anatomical_region_level_3': 'Anatomical region level 3',
                               'anatomical_region': 'Anatomical region',
                               'tissue_ontology_term': 'Tissue (ontology term)',
                               'tissue_ontology_term_id': 'Tissue (ontology term ID)',
                               'skin_tissue': 'Skin tissue',
                               'sample_collection_year': 'Sample collection year',
                               'library_id': 'Library ID',
                               'library_id_repository': 'Library (ID repository)',
                               'library_preparation_batch': 'Library preparation batch',
                               'library_sequencing_run': 'Library sequencing run',
                               'sample_preservation_method': 'Sample preservation method',
                               'dissociation_protocol': 'Dissociation protocol',
                               'cell_enrichment': 'Cell enrichment',
                               'cell_viability_percentage': 'Cell viability percentage',
                               'cell_number_loaded': 'Cell number loaded',
                               'assay_ontology_term': 'Assay (ontology term)',
                               'assay_ontology_term_id': 'Assay (ontology term ID)',
                               'sequenced_fragment': 'Sequenced fragment',
                               'sequencing_platform': 'Sequencing platform',
                               'reference_genome': 'Reference genome',
                               'alignment_software': 'Alignment software',
                               'author_cell_type_1': 'Author cell type 1',
                               'author_cell_type_2': 'Author cell type 2',
                               'author_cell_type_3': 'Author cell type 3',
                               'gene_expression': 'Total RNA counts',
                               'expressed_genes': 'Number of detected genes',
                               'individual_leiden_L1_001': 'Leiden L1 001',
                               'individual_leiden_L1_002': 'Leiden L1 002',
                               'individual_leiden_L1_003': 'Leiden L1 003',
                               'individual_leiden_L1_005': 'Leiden L1 005',
                               'individual_leiden_L1_008': 'Leiden L1 008',
                               'individual_leiden_L1_010': 'Leiden L1 010',
                               'individual_leiden_L1_012': 'Leiden L1 012',
                               'individual_leiden_L1_015': 'Leiden L1 015',
                               'individual_leiden_L1_020': 'Leiden L1 020',
                               'individual_leiden_L1_030': 'Leiden L1 030',
                               'individual_leiden_L1_050': 'Leiden L1 050',
                               'individual_leiden_L1_080': 'Leiden L1 080',
                               'individual_leiden_L1_100': 'Leiden L1 100',
                               'individual_leiden_L1_120': 'Leiden L1 120',
                               'individual_leiden_L1_150': 'Leiden L1 150',
                               'individual_leiden_L1_200': 'Leiden L1 200',
                               'individual_leiden_L1_300': 'Leiden L1 300',
                               'AnnotationName': 'Cell type',
                               'AnnotationABBR': 'Cell type (abbr)'}")
py_run_string("new_colnames_keeping_list = [x for x in new_colnames.keys()]")
#
data_output_dir = paste0(py$this_output_dir, "/data_h5ad")
dir.create(data_output_dir, showWarnings = F, recursive = T)
py_run_string("adata_dict = {}")
for(ii in 1:24){
  adata_str = paste0("adata_", ii)
  py_run_string(paste0("tmp_adata = ", adata_str, "[:, ", adata_str, ".var_names.str.startswith('ENSG')].copy()"))
  py_run_string("tmp_adata.X = tmp_adata.layers['raw_counts'].copy()")
  GeneID_to_GeneName("tmp_adata", "tmp_adata", "feature_name")
  py_run_string("tmp_adata.layers['raw_counts'] = tmp_adata.X.copy()")
  py_run_string("sc.pp.normalize_total(tmp_adata, target_sum=1e4)")
  py_run_string("sc.pp.log1p(tmp_adata)")
  py_run_string("del tmp_adata.X")
  py_run_string("tmp_adata.obs['donor_id_running_short'] = tmp_adata.obs['donor_id_running'].str.split('___').str[1].to_list()")
  py_run_string("tmp_adata.obs['sample_id_running_short'] = tmp_adata.obs['sample_id_running'].str.split('___').str[1].to_list()")
  py_run_string("tmp_adata.obs['gene_expression'] = tmp_adata.layers['raw_counts'].sum(1)")
  py_run_string("tmp_adata.obs['expressed_genes'] = np.sum(tmp_adata.layers['raw_counts'] > 0, axis=1)")
  # Remove all columns that are not in the dict
  py_run_string("tmp_adata.obs = tmp_adata.obs.loc[:, new_colnames_keeping_list].copy()")
  py_run_string("tmp_adata.obs.rename(new_colnames, axis=1, inplace=True)")
  py_run_string(paste0("adata_dict['", ii, "'] = tmp_adata"))
  print(ii)
}
py_run_string("adata_output = ad.concat(adata_dict, axis=0, join='outer')")
py_run_string("del adata_dict")
#
hgnc_mapping_genename_ensemblID = readRDS("/home/haoy/projects/Analysis/useful_files/hgnc_mapping_genename_ensemblID.rds")
GRCh38v103_mapping_genename_ensemblID = readRDS("/home/haoy/projects/Analysis/useful_files/GRCh38v103_mapping_genename_ensemblID.rds")
#
GeneName = py$adata_output$var_names$values
GeneName_to_GeneID = mapping_GeneName_GeneID(GeneName, hgnc_mapping_genename_ensemblID)
GeneName_to_GeneID = mapping_GeneName_GeneID(GeneName_to_GeneID, GRCh38v103_mapping_genename_ensemblID)
names(GeneName_to_GeneID) = GeneName
GeneID = GeneName_to_GeneID[GeneName]
GeneID[is.na(GeneID)] = GeneName[is.na(GeneID)]
py$adata_output$var[["feature_id"]] = GeneID
py_run_string("gc.collect()")
#
for(ii in 1:24){
  adata_str = paste0("adata_", ii, "_output")
  this_StudyID = rownames(dataset_information_df)[ii]
  py_run_string(paste0(adata_str, " = adata_output[adata_output.obs['Study ID'] == '", this_StudyID, "'].copy()"))
  py_run_string(paste0(adata_str, ".write('", paste0(data_output_dir, "/", paste0("adata_", ii), ".h5ad"), "')"))
}
#
#
linking_str = "___linking___"
#
this_metadata_field = "AnnotationABBR"
this_metadata_categories = c("KC", "MEL", "SC", "FIB", "MUR", "VEC", "LEC", "ERY", "LYM", "MYE", "MAST", "PC", "NA")
cell_type_CN_matrix = matrix(nrow = 24, ncol = length(this_metadata_categories), dimnames = list(dataset_information_df[, "StudyID_short"], this_metadata_categories))
cell_type_DN_matrix = matrix(nrow = 24, ncol = length(this_metadata_categories), dimnames = list(dataset_information_df[, "StudyID_short"], this_metadata_categories))
for(ii in 1:24){
  adata_str = paste0("adata_", ii)
  #
  this_metadata = eval(parse(text = paste0("py$", adata_str, "$obs[['", this_metadata_field, "']]")))
  this_CN_count = table(this_metadata)[this_metadata_categories]
  this_CN_count[is.na(this_CN_count)] = 0
  #
  donor_id_running = eval(parse(text = paste0("py$", adata_str, "$obs[['donor_id_running']]")))
  this_DN_combination = unique(paste0(this_metadata, linking_str, donor_id_running))
  this_DN_unique = sapply(this_DN_combination, function(x){
    unlist(strsplit(x, linking_str))[1]
  })
  this_DN_count = table(this_DN_unique)[this_metadata_categories]
  this_DN_count[is.na(this_DN_count)] = 0
  #
  cell_type_CN_matrix[ii, ] = this_CN_count
  cell_type_DN_matrix[ii, ] = this_DN_count
}
write.table(cell_type_CN_matrix, paste0(py$this_figure_dir, "/cell_type_CN_matrix.txt"), sep = "\t", quote = F, row.names = TRUE, col.names = TRUE)
write.table(cell_type_DN_matrix, paste0(py$this_figure_dir, "/cell_type_DN_matrix.txt"), sep = "\t", quote = F, row.names = TRUE, col.names = TRUE)
#
this_metadata_field = "anatomical_region_level_1"
this_metadata_categories = c("Head", "Extremities", "Torso", "unknown")
bodysites_CN_matrix = matrix(nrow = 24, ncol = length(this_metadata_categories), dimnames = list(dataset_information_df[, "StudyID_short"], this_metadata_categories))
bodysites_DN_matrix = matrix(nrow = 24, ncol = length(this_metadata_categories), dimnames = list(dataset_information_df[, "StudyID_short"], this_metadata_categories))
for(ii in 1:24){
  adata_str = paste0("adata_", ii)
  #
  this_metadata = eval(parse(text = paste0("py$", adata_str, "$obs[['", this_metadata_field, "']]")))
  this_CN_count = table(this_metadata)[this_metadata_categories]
  this_CN_count[is.na(this_CN_count)] = 0
  #
  donor_id_running = eval(parse(text = paste0("py$", adata_str, "$obs[['donor_id_running']]")))
  this_DN_combination = unique(paste0(this_metadata, linking_str, donor_id_running))
  this_DN_unique = sapply(this_DN_combination, function(x){
    unlist(strsplit(x, linking_str))[1]
  })
  this_DN_count = table(this_DN_unique)[this_metadata_categories]
  this_DN_count[is.na(this_DN_count)] = 0
  #
  bodysites_CN_matrix[ii, ] = this_CN_count
  bodysites_DN_matrix[ii, ] = this_DN_count
}
write.table(bodysites_CN_matrix, paste0(py$this_figure_dir, "/bodysites_CN_matrix.txt"), sep = "\t", quote = F, row.names = TRUE, col.names = TRUE)
write.table(bodysites_DN_matrix, paste0(py$this_figure_dir, "/bodysites_DN_matrix.txt"), sep = "\t", quote = F, row.names = TRUE, col.names = TRUE)
#
this_metadata_field = "age_range"
this_metadata_categories = c("0-14", "15-19", "20-29", "30-39", "40-49", "50-59", "60-69", "70-79", "80-89", "90-99", "unknown")
age_range_CN_matrix = matrix(nrow = 24, ncol = length(this_metadata_categories), dimnames = list(dataset_information_df[, "StudyID_short"], this_metadata_categories))
age_range_DN_matrix = matrix(nrow = 24, ncol = length(this_metadata_categories), dimnames = list(dataset_information_df[, "StudyID_short"], this_metadata_categories))
for(ii in 1:24){
  adata_str = paste0("adata_", ii)
  #
  this_metadata = eval(parse(text = paste0("py$", adata_str, "$obs[['", this_metadata_field, "']]")))
  this_CN_count = table(this_metadata)[this_metadata_categories]
  this_CN_count[is.na(this_CN_count)] = 0
  #
  donor_id_running = eval(parse(text = paste0("py$", adata_str, "$obs[['donor_id_running']]")))
  this_DN_combination = unique(paste0(this_metadata, linking_str, donor_id_running))
  this_DN_unique = sapply(this_DN_combination, function(x){
    unlist(strsplit(x, linking_str))[1]
  })
  this_DN_count = table(this_DN_unique)[this_metadata_categories]
  this_DN_count[is.na(this_DN_count)] = 0
  #
  age_range_CN_matrix[ii, ] = this_CN_count
  age_range_DN_matrix[ii, ] = this_DN_count
}
write.table(age_range_CN_matrix, paste0(py$this_figure_dir, "/age_range_CN_matrix.txt"), sep = "\t", quote = F, row.names = TRUE, col.names = TRUE)
write.table(age_range_DN_matrix, paste0(py$this_figure_dir, "/age_range_DN_matrix.txt"), sep = "\t", quote = F, row.names = TRUE, col.names = TRUE)
#
this_metadata_field = "sex_ontology_term"
this_metadata_categories = c("female", "male", "unknown")
sex_ontology_term_CN_matrix = matrix(nrow = 24, ncol = length(this_metadata_categories), dimnames = list(dataset_information_df[, "StudyID_short"], this_metadata_categories))
sex_ontology_term_DN_matrix = matrix(nrow = 24, ncol = length(this_metadata_categories), dimnames = list(dataset_information_df[, "StudyID_short"], this_metadata_categories))
for(ii in 1:24){
  adata_str = paste0("adata_", ii)
  #
  this_metadata = eval(parse(text = paste0("py$", adata_str, "$obs[['", this_metadata_field, "']]")))
  this_CN_count = table(this_metadata)[this_metadata_categories]
  this_CN_count[is.na(this_CN_count)] = 0
  #
  donor_id_running = eval(parse(text = paste0("py$", adata_str, "$obs[['donor_id_running']]")))
  this_DN_combination = unique(paste0(this_metadata, linking_str, donor_id_running))
  this_DN_unique = sapply(this_DN_combination, function(x){
    unlist(strsplit(x, linking_str))[1]
  })
  this_DN_count = table(this_DN_unique)[this_metadata_categories]
  this_DN_count[is.na(this_DN_count)] = 0
  #
  sex_ontology_term_CN_matrix[ii, ] = this_CN_count
  sex_ontology_term_DN_matrix[ii, ] = this_DN_count
}
write.table(sex_ontology_term_CN_matrix, paste0(py$this_figure_dir, "/sex_ontology_term_CN_matrix.txt"), sep = "\t", quote = F, row.names = TRUE, col.names = TRUE)
write.table(sex_ontology_term_DN_matrix, paste0(py$this_figure_dir, "/sex_ontology_term_DN_matrix.txt"), sep = "\t", quote = F, row.names = TRUE, col.names = TRUE)
#
this_metadata_field = "Anatomical region"
this_metadata_categories = names(table(py$adata_output$obs[['Anatomical region']]))
anatomical_region_CN_matrix = matrix(nrow = 24, ncol = length(this_metadata_categories), dimnames = list(dataset_information_df[, "StudyID_short"], this_metadata_categories))
anatomical_region_DN_matrix = matrix(nrow = 24, ncol = length(this_metadata_categories), dimnames = list(dataset_information_df[, "StudyID_short"], this_metadata_categories))
for(ii in 1:24){
  adata_str = paste0("adata_", ii, "_output")
  #
  this_metadata = eval(parse(text = paste0("py$", adata_str, "$obs[['", this_metadata_field, "']]")))
  this_CN_count = table(this_metadata)[this_metadata_categories]
  this_CN_count[is.na(this_CN_count)] = 0
  #
  donor_id = eval(parse(text = paste0("py$", adata_str, "$obs[['Donor ID']]")))
  experiment_id = eval(parse(text = paste0("py$", adata_str, "$obs[['Experiment ID']]")))
  this_DN_combination = unique(paste0(this_metadata, linking_str, experiment_id, "___", donor_id))
  this_DN_unique = sapply(this_DN_combination, function(x){
    unlist(strsplit(x, linking_str))[1]
  })
  this_DN_count = table(this_DN_unique)[this_metadata_categories]
  this_DN_count[is.na(this_DN_count)] = 0
  #
  anatomical_region_CN_matrix[ii, ] = this_CN_count
  anatomical_region_DN_matrix[ii, ] = this_DN_count
}
write.table(anatomical_region_CN_matrix, paste0(py$this_figure_dir, "/anatomical_region_CN_matrix.txt"), sep = "\t", quote = F, row.names = TRUE, col.names = TRUE)
write.table(anatomical_region_DN_matrix, paste0(py$this_figure_dir, "/anatomical_region_DN_matrix.txt"), sep = "\t", quote = F, row.names = TRUE, col.names = TRUE)
#
donor_number = c()
for(ii in 1:24){
  adata_str = paste0("adata_", ii)
  donor_number = c(donor_number, length(table(eval(parse(text = paste0("py$", adata_str, "$obs[['donor_id_running']]"))))))
}
names(donor_number) = dataset_information_df[, "StudyID_short"]

