r_utilities = "/home/haoy/projects/Analysis/code_library/utilities.r"
py_utilities = "/home/haoy/projects/Analysis/code_library/utilities.py"
py_visualization = "/home/haoy/projects/Analysis/code_library/visualization.py"
source(r_utilities)
library(reticulate)
py_run_file(py_utilities, convert = F)
py_run_file(py_visualization, convert = F)
######
hgnc_mapping_ensemblID_genename = readRDS("/home/haoy/projects/Analysis/useful_files/hgnc_mapping_ensemblID_genename_seurat.rds")
######
######
######
py$output_dir = "/home/haoy/projects/MetaStudiesAnalysis/Integration/Output/Integration/v70_gauss"
py$output_dir_cache = "/home/haoy/projects/MetaStudiesAnalysis/Integration/Output/Integration/v70_gauss"
dir.create(py$output_dir, showWarnings = F, recursive = T)
######
raw_input_paths = list()
######
raw_input_paths[["D001"]] = "/home/haoy/projects/MetaStudiesAnalysis/Integration/Data/Skin/Human/data_preparation_HCA_v1_0/SS3E_Human_Skin_20K/SS3E_Human_Skin_20K_ALL.h5ad"
raw_input_paths[["D002"]] = "/home/haoy/projects/MetaStudiesAnalysis/Integration/Data/Skin/Human/data_preparation_HCA_v1_0/Liu_Landen_CellStemCell_2024/GSE241132_Ning.h5ad"
raw_input_paths[["D003"]] = "/home/haoy/projects/MetaStudiesAnalysis/Integration/Data/Skin/Human/data_preparation_HCA_v1_0/Gaydosik_Fuschiotti_ClinicalCancerResearch_2019/GSE128531_HC.h5ad"
raw_input_paths[["D004"]] = "/home/haoy/projects/MetaStudiesAnalysis/Integration/Data/Skin/Human/data_preparation_HCA_v1_0/Wang_Atwood_NatureCommunications_2020/GSE147482_ALL.h5ad"
raw_input_paths[["D005"]] = "/home/haoy/projects/MetaStudiesAnalysis/Integration/Data/Skin/Human/data_preparation_HCA_v1_0/Takahashi_Lowry_JournalofInvestigativeDermatology_2020/GSE129611_ALL.h5ad"
raw_input_paths[["D006"]] = "/home/haoy/projects/MetaStudiesAnalysis/Integration/Data/Skin/Human/data_preparation_HCA_v1_0/Belote_Torres_NatureCellBiology_2021/GSE151091_Belote.h5ad"
raw_input_paths[["D007"]] = "/home/haoy/projects/MetaStudiesAnalysis/Integration/Data/Skin/Human/data_preparation_HCA_v1_0/Boldo_Lyko_CommunicationsBiology_2020/GSE130973_ALL.h5ad"
raw_input_paths[["D008"]] = "/home/haoy/projects/MetaStudiesAnalysis/Integration/Data/Skin/Human/data_preparation_HCA_v1_0/Cheng_Cho_CellReports_2018/Cheng_Cho_CellReports_2018_HC.h5ad"
raw_input_paths[["D009"]] = "/home/haoy/projects/MetaStudiesAnalysis/Integration/Data/Skin/Human/data_preparation_HCA_v1_0/Ji_Khavari_Cell_2020/GSE144236_ADJACENT.h5ad"
raw_input_paths[["D010"]] = "/home/haoy/projects/MetaStudiesAnalysis/Integration/Data/Skin/Human/data_preparation_HCA_v1_0/Tabib_Lafyatis_NatureCommunications_2021/GSE138669_HC.h5ad"
raw_input_paths[["D011"]] = "/home/haoy/projects/MetaStudiesAnalysis/Integration/Data/Skin/Human/data_preparation_HCA_v1_0/He_Guttman_JournalofAllergyandClinicalImmunology_2020/GSE147424_HC.h5ad"
raw_input_paths[["D012"]] = "/home/haoy/projects/MetaStudiesAnalysis/Integration/Data/Skin/Human/data_preparation_HCA_v1_0/Zou_Liu_DevelopmentalCell_2020/HRA000395_ALL.h5ad"
raw_input_paths[["D013"]] = "/home/haoy/projects/MetaStudiesAnalysis/Integration/Data/Skin/Human/data_preparation_HCA_v1_0/Wiedemann_Andersen_CellReports_2023/GSE202352_ALL.h5ad"
raw_input_paths[["D014"]] = "/home/haoy/projects/MetaStudiesAnalysis/Integration/Data/Skin/Human/data_preparation_HCA_v1_0/Alkon_Stingl_JournalofAllergyandClinicalImmunology_2022/GSE180885_HC.h5ad"
raw_input_paths[["D015"]] = "/home/haoy/projects/MetaStudiesAnalysis/Integration/Data/Skin/Human/data_preparation_HCA_v1_0/Dunlap_Rao_JCIInsight_2022/GSE186476_HC.h5ad"
raw_input_paths[["D016"]] = "/home/haoy/projects/MetaStudiesAnalysis/Integration/Data/Skin/Human/data_preparation_HCA_v1_0/Boothby_Rosenblum_Nature_2021/GSE183031_HC.h5ad"
raw_input_paths[["D017"]] = "/home/haoy/projects/MetaStudiesAnalysis/Integration/Data/Skin/Human/data_preparation_HCA_v1_0/Gao_Hu_CellDeathDisease_2021/GSE162183_HC.h5ad"
raw_input_paths[["D018"]] = "/home/haoy/projects/MetaStudiesAnalysis/Integration/Data/Skin/Human/data_preparation_HCA_v1_0/Hughes_Shalek_Immunity_2020/GSE150672_HC.h5ad"
raw_input_paths[["D019"]] = "/home/haoy/projects/MetaStudiesAnalysis/Integration/Data/Skin/Human/data_preparation_HCA_v1_0/Kim_Nagao_NatureMedicine_2020/GSE132802_HC.h5ad"
raw_input_paths[["D020"]] = "/home/haoy/projects/MetaStudiesAnalysis/Integration/Data/Skin/Human/data_preparation_HCA_v1_0/Rindler_Brunner_MolecularCancer_2021/GSE173205_HC.h5ad"
raw_input_paths[["D021"]] = "/home/haoy/projects/MetaStudiesAnalysis/Integration/Data/Skin/Human/data_preparation_HCA_v1_0/Rojahn_Brunner_JournalofAllergyandClinicalImmunology_2020/GSE153760_HC.h5ad"
raw_input_paths[["D022"]] = "/home/haoy/projects/MetaStudiesAnalysis/Integration/Data/Skin/Human/data_preparation_HCA_v1_0/Theocharidis_Bhasin_NatureCommunications_2022/GSE165816_SKIN_HC.h5ad"
raw_input_paths[["D023"]] = "/home/haoy/projects/MetaStudiesAnalysis/Integration/Data/Skin/Human/data_preparation_HCA_v1_0/Vorstandlechner_Mildner_NatureCommunications_2021/GSE156326_HC.h5ad"
raw_input_paths[["D024"]] = "/home/haoy/projects/MetaStudiesAnalysis/Integration/Data/Skin/Human/data_preparation_HCA_v1_0/Xu_Chen_Nature_2022/OMIX691_HC.h5ad"
raw_input_paths[["D025"]] = "/home/haoy/projects/MetaStudiesAnalysis/Integration/Data/Skin/Human/data_preparation_HCA_v1_0/Reynolds_Haniffa_Science_2021/Reynolds_Haniffa_Science_2021_10X_POSTPARTUM_NORMAL.h5ad"
raw_input_paths[["D026"]] = "/home/haoy/projects/MetaStudiesAnalysis/Integration/Data/Skin/Human/data_preparation_HCA_v1_0/Reynolds_Greenleaf_NatureGenetics_2023/GSE212447_HC.h5ad"
#
individual_annotation_result_paths = list()
individual_hvg_clustering_paths = list()
######
individual_annotation_result_paths[["D001"]] = paste0("/home/haoy/projects/MetaStudiesAnalysis/Integration/Output/SS3_HumanSkin_20K/ALL/v17_gauss/adata_imbalance_annotated_2_20240814_SampleID_Fixed.h5ad")
individual_hvg_clustering_paths[["D001"]] = paste0("/home/haoy/projects/MetaStudiesAnalysis/Integration/Output/SS3_HumanSkin_20K/ALL/v17_gauss/adata_imbalance_hvg_clustering_1.h5ad")
individual_annotation_result_paths[["D002"]] = paste0("/home/haoy/projects/MetaStudiesAnalysis/Integration/Output/Ning/ALL/v10_gauss/adata_imbalance_annotated_1.h5ad")
individual_hvg_clustering_paths[["D002"]] = paste0("/home/haoy/projects/MetaStudiesAnalysis/Integration/Output/Ning/ALL/v10_gauss/adata_imbalance_hvg_clustering_1.h5ad")
individual_annotation_result_paths[["D003"]] = paste0("/home/haoy/projects/MetaStudiesAnalysis/Integration/Output/Gaydosik_Fuschiotti_ClinicalCancerResearch_2019/ALL/v10_gauss/adata_imbalance_annotated_1.h5ad")
individual_hvg_clustering_paths[["D003"]] = paste0("/home/haoy/projects/MetaStudiesAnalysis/Integration/Output/Gaydosik_Fuschiotti_ClinicalCancerResearch_2019/ALL/v10_gauss/adata_imbalance_hvg_clustering_1.h5ad")
individual_annotation_result_paths[["D004"]] = paste0("/home/haoy/projects/MetaStudiesAnalysis/Integration/Output/Wang_Atwood_NatureCommunications_2020/ALL/v10_gauss/adata_imbalance_annotated_1.h5ad")
individual_hvg_clustering_paths[["D004"]] = paste0("/home/haoy/projects/MetaStudiesAnalysis/Integration/Output/Wang_Atwood_NatureCommunications_2020/ALL/v10_gauss/adata_imbalance_hvg_clustering_1.h5ad")
individual_annotation_result_paths[["D005"]] = paste0("/home/haoy/projects/MetaStudiesAnalysis/Integration/Output/Takahashi_Lowry_JournalofInvestigativeDermatology_2020/ALL/v10_gauss/adata_imbalance_annotated_1.h5ad")
individual_hvg_clustering_paths[["D005"]] = paste0("/home/haoy/projects/MetaStudiesAnalysis/Integration/Output/Takahashi_Lowry_JournalofInvestigativeDermatology_2020/ALL/v10_gauss/adata_imbalance_hvg_clustering_1.h5ad")
individual_annotation_result_paths[["D006"]] = paste0("/home/haoy/projects/MetaStudiesAnalysis/Integration/Output/Belote_Torres_NatureCellBiology_2021/ALL/v10_gauss/adata_imbalance_annotated_1.h5ad")
individual_hvg_clustering_paths[["D006"]] = paste0("/home/haoy/projects/MetaStudiesAnalysis/Integration/Output/Belote_Torres_NatureCellBiology_2021/ALL/v10_gauss/adata_imbalance_hvg_clustering_1.h5ad")
individual_annotation_result_paths[["D007"]] = paste0("/home/haoy/projects/MetaStudiesAnalysis/Integration/Output/Boldo_Lyko_CommunicationsBiology_2020/ALL/v10_gauss/adata_imbalance_annotated_1.h5ad")
individual_hvg_clustering_paths[["D007"]] = paste0("/home/haoy/projects/MetaStudiesAnalysis/Integration/Output/Boldo_Lyko_CommunicationsBiology_2020/ALL/v10_gauss/adata_imbalance_hvg_clustering_1.h5ad")
individual_annotation_result_paths[["D008"]] = paste0("/home/haoy/projects/MetaStudiesAnalysis/Integration/Output/Cheng_Cho_CellReports_2018/ALL/v10_gauss/adata_imbalance_annotated_1.h5ad")
individual_hvg_clustering_paths[["D008"]] = paste0("/home/haoy/projects/MetaStudiesAnalysis/Integration/Output/Cheng_Cho_CellReports_2018/ALL/v10_gauss/adata_imbalance_hvg_clustering_1.h5ad")
individual_annotation_result_paths[["D009"]] = paste0("/home/haoy/projects/MetaStudiesAnalysis/Integration/Output/Ji_Khavari_Cell_2020/ADJACENT/v10_gauss/adata_imbalance_annotated_1.h5ad")
individual_hvg_clustering_paths[["D009"]] = paste0("/home/haoy/projects/MetaStudiesAnalysis/Integration/Output/Ji_Khavari_Cell_2020/ADJACENT/v10_gauss/adata_imbalance_hvg_clustering_1.h5ad")
individual_annotation_result_paths[["D010"]] = paste0("/home/haoy/projects/MetaStudiesAnalysis/Integration/Output/Tabib_Lafyatis_NatureCommunications_2021/ALL/v10_gauss/adata_imbalance_annotated_1.h5ad")
individual_hvg_clustering_paths[["D010"]] = paste0("/home/haoy/projects/MetaStudiesAnalysis/Integration/Output/Tabib_Lafyatis_NatureCommunications_2021/ALL/v10_gauss/adata_imbalance_hvg_clustering_1.h5ad")
individual_annotation_result_paths[["D011"]] = paste0("/home/haoy/projects/MetaStudiesAnalysis/Integration/Output/He_Guttman_JournalofAllergyandClinicalImmunology_2020/ALL/v10_gauss/adata_imbalance_annotated_1.h5ad")
individual_hvg_clustering_paths[["D011"]] = paste0("/home/haoy/projects/MetaStudiesAnalysis/Integration/Output/He_Guttman_JournalofAllergyandClinicalImmunology_2020/ALL/v10_gauss/adata_imbalance_hvg_clustering_1.h5ad")
individual_annotation_result_paths[["D012"]] = paste0("/home/haoy/projects/MetaStudiesAnalysis/Integration/Output/Zou_Liu_DevelopmentalCell_2020/ALL/v10_gauss/adata_imbalance_annotated_1.h5ad")
individual_hvg_clustering_paths[["D012"]] = paste0("/home/haoy/projects/MetaStudiesAnalysis/Integration/Output/Zou_Liu_DevelopmentalCell_2020/ALL/v10_gauss/adata_imbalance_hvg_clustering_1.h5ad")
individual_annotation_result_paths[["D013"]] = paste0("/home/haoy/projects/MetaStudiesAnalysis/Integration/Output/Wiedemann_Andersen_CellReports_2023/ALL/v11_gauss/adata_imbalance_annotated_1.h5ad")
individual_hvg_clustering_paths[["D013"]] = paste0("/home/haoy/projects/MetaStudiesAnalysis/Integration/Output/Wiedemann_Andersen_CellReports_2023/ALL/v11_gauss/adata_imbalance_hvg_clustering_1.h5ad")
individual_annotation_result_paths[["D014"]] = paste0("/home/haoy/projects/MetaStudiesAnalysis/Integration/Output/Alkon_Stingl_JournalofAllergyandClinicalImmunology_2022/ALL/v10_gauss/adata_imbalance_annotated_1.h5ad")
individual_hvg_clustering_paths[["D014"]] = paste0("/home/haoy/projects/MetaStudiesAnalysis/Integration/Output/Alkon_Stingl_JournalofAllergyandClinicalImmunology_2022/ALL/v10_gauss/adata_imbalance_hvg_clustering_1.h5ad")
individual_annotation_result_paths[["D015"]] = paste0("/home/haoy/projects/MetaStudiesAnalysis/Integration/Output/Dunlap_Rao_JCIInsight_2022/ALL/v11_gauss/adata_imbalance_annotated_1.h5ad")
individual_hvg_clustering_paths[["D015"]] = paste0("/home/haoy/projects/MetaStudiesAnalysis/Integration/Output/Dunlap_Rao_JCIInsight_2022/ALL/v11_gauss/adata_imbalance_hvg_clustering_1.h5ad")
individual_annotation_result_paths[["D016"]] = paste0("/home/haoy/projects/MetaStudiesAnalysis/Integration/Output/Boothby_Rosenblum_Nature_2021/ALL/v10_gauss/adata_imbalance_annotated_1.h5ad")
individual_hvg_clustering_paths[["D016"]] = paste0("/home/haoy/projects/MetaStudiesAnalysis/Integration/Output/Boothby_Rosenblum_Nature_2021/ALL/v10_gauss/adata_imbalance_hvg_clustering_1.h5ad")
individual_annotation_result_paths[["D017"]] = paste0("/home/haoy/projects/MetaStudiesAnalysis/Integration/Output/Gao_Hu_CellDeathDisease_2021/ALL/v10_gauss/adata_imbalance_annotated_1.h5ad")
individual_hvg_clustering_paths[["D017"]] = paste0("/home/haoy/projects/MetaStudiesAnalysis/Integration/Output/Gao_Hu_CellDeathDisease_2021/ALL/v10_gauss/adata_imbalance_hvg_clustering_1.h5ad")
individual_annotation_result_paths[["D018"]] = paste0("/home/haoy/projects/MetaStudiesAnalysis/Integration/Output/Hughes_Shalek_Immunity_2020/HC/v10_gauss/adata_imbalance_annotated_1.h5ad")
individual_hvg_clustering_paths[["D018"]] = paste0("/home/haoy/projects/MetaStudiesAnalysis/Integration/Output/Hughes_Shalek_Immunity_2020/HC/v10_gauss/adata_imbalance_hvg_clustering_1.h5ad")
individual_annotation_result_paths[["D019"]] = paste0("/home/haoy/projects/MetaStudiesAnalysis/Integration/Output/Kim_Nagao_NatureMedicine_2020/ALL/v10_gauss/adata_imbalance_annotated_1.h5ad")
individual_hvg_clustering_paths[["D019"]] = paste0("/home/haoy/projects/MetaStudiesAnalysis/Integration/Output/Kim_Nagao_NatureMedicine_2020/ALL/v10_gauss/adata_imbalance_hvg_clustering_1.h5ad")
individual_annotation_result_paths[["D020"]] = paste0("/home/haoy/projects/MetaStudiesAnalysis/Integration/Output/Rindler_Brunner_MolecularCancer_2021/ALL/v10_gauss/adata_imbalance_annotated_1.h5ad")
individual_hvg_clustering_paths[["D020"]] = paste0("/home/haoy/projects/MetaStudiesAnalysis/Integration/Output/Rindler_Brunner_MolecularCancer_2021/ALL/v10_gauss/adata_imbalance_hvg_clustering_1.h5ad")
individual_annotation_result_paths[["D021"]] = paste0("/home/haoy/projects/MetaStudiesAnalysis/Integration/Output/Rojahn_Brunner_JournalofAllergyandClinicalImmunology_2020/ALL/v10_gauss/adata_imbalance_annotated_1.h5ad")
individual_hvg_clustering_paths[["D021"]] = paste0("/home/haoy/projects/MetaStudiesAnalysis/Integration/Output/Rojahn_Brunner_JournalofAllergyandClinicalImmunology_2020/ALL/v10_gauss/adata_imbalance_hvg_clustering_1.h5ad")
individual_annotation_result_paths[["D022"]] = paste0("/home/haoy/projects/MetaStudiesAnalysis/Integration/Output/Theocharidis_Bhasin_NatureCommunications_2022/ALL/v11_gauss/adata_imbalance_annotated_1.h5ad")
individual_hvg_clustering_paths[["D022"]] = paste0("/home/haoy/projects/MetaStudiesAnalysis/Integration/Output/Theocharidis_Bhasin_NatureCommunications_2022/ALL/v11_gauss/adata_imbalance_hvg_clustering_1.h5ad")
individual_annotation_result_paths[["D023"]] = paste0("/home/haoy/projects/MetaStudiesAnalysis/Integration/Output/Vorstandlechner_Mildner_NatureCommunications_2021/ALL/v10_gauss/adata_imbalance_annotated_1.h5ad")
individual_hvg_clustering_paths[["D023"]] = paste0("/home/haoy/projects/MetaStudiesAnalysis/Integration/Output/Vorstandlechner_Mildner_NatureCommunications_2021/ALL/v10_gauss/adata_imbalance_hvg_clustering_1.h5ad")
individual_annotation_result_paths[["D024"]] = paste0("/home/haoy/projects/MetaStudiesAnalysis/Integration/Output/Xu_Chen_Nature_2022/ALL/v10_gauss/adata_imbalance_annotated_1.h5ad")
individual_hvg_clustering_paths[["D024"]] = paste0("/home/haoy/projects/MetaStudiesAnalysis/Integration/Output/Xu_Chen_Nature_2022/ALL/v10_gauss/adata_imbalance_hvg_clustering_1.h5ad")
individual_annotation_result_paths[["D025"]] = paste0("/home/haoy/projects/MetaStudiesAnalysis/Integration/Output/Reynolds_Haniffa_Science_2021/10X_POSTPARTUM_NORMAL/v10_gauss/adata_imbalance_annotated_1.h5ad")
individual_hvg_clustering_paths[["D025"]] = paste0("/home/haoy/projects/MetaStudiesAnalysis/Integration/Output/Reynolds_Haniffa_Science_2021/10X_POSTPARTUM_NORMAL/v10_gauss/adata_imbalance_hvg_clustering_1.h5ad")
individual_annotation_result_paths[["D026"]] = paste0("/home/haoy/projects/MetaStudiesAnalysis/Integration/Output/Reynolds_Greenleaf_NatureGenetics_2023/HC/v1_gauss/adata_annotated_1.h5ad")
individual_hvg_clustering_paths[["D026"]] = paste0("/home/haoy/projects/MetaStudiesAnalysis/Integration/Output/Reynolds_Greenleaf_NatureGenetics_2023/HC/v1_gauss/adata_clustering_hvg.h5ad")
#for(ii in seq(24)){
#  individual_result_paths[[sprintf("D%.3d", ii)]] = paste0("/home/haoy/projects/MetaStudiesAnalysis/Integration/Output/Integration/v67_gauss/adata_", ii + 1, ".h5ad")
#  individual_hvg_result_paths[[sprintf("D%.3d", ii)]] = paste0("/home/haoy/projects/MetaStudiesAnalysis/Integration/Output/Integration/v67_gauss/adata_", ii + 1, "_hvg.h5ad")
#}
# Clustering parameters
hvg_number = 2000
pca_dim = 50
random_state = 0
n_neighbors_clustering = 30
n_neighbors_visualization = 15
leiden_resolutions = c(0.01, 0.02, 0.03, 0.05, 0.08, 0.1, 0.12, 0.15, 0.2, 0.3, 0.5, 0.8, 1, 1.2, 1.5, 2, 3)
cell_cycle_regress_out = F
# DE parameters
top_DE_range = 50
# Other parameters
max_number = 4000
py_run_string("density_color = 'YlOrRd'")
# Other settings
do_DEG = F
use_cache = T
######
X_gene_path = "/home/haoy/projects/Analysis/useful_files/X_GeneName_Homo_sapiens.GRCh38.103.rds"
Y_gene_path = "/home/haoy/projects/Analysis/useful_files/Y_GeneName_Homo_sapiens.GRCh38.103.rds"
py$X_Genes = readRDS(X_gene_path)
py$Y_Genes = readRDS(Y_gene_path)
alias_gene_c_seurat_path = "/home/haoy/projects/Analysis/useful_files/hgnc_complete_set_c_seurat.rds"
CellCyele_gene_path = "/home/haoy/projects/MetaStudiesAnalysis/Integration/regev_lab_cell_cycle_genes.txt"
alias_gene_c_seurat = readRDS(alias_gene_c_seurat_path)
######
# Following https://nbviewer.org/github/theislab/scanpy_usage/blob/master/180209_cell_cycle/cell_cycle.ipynb
CellCyele_gene = gsub(" ", "", toupper(as.matrix(read.table(CellCyele_gene_path))[, 1]))
CellCyele_gene[CellCyele_gene == "MLF1IP"] = "CENPU" # MLF1IP is coded by CENPU
py$s_genes = CellCyele_gene[1:43]
py$g2m_genes = CellCyele_gene[44:length(CellCyele_gene)]
# Settings
py_run_string("sc.settings.verbosity = 3")
py_run_string("sc.logging.print_header()")
py_run_string("sc.settings.set_figure_params(dpi=600, dpi_save=600, facecolor='white', fontsize=10)")
py_run_string("sc.settings.autoshow = False")
py_run_string("figure_dir = output_dir + '/figures'")
py_run_string("sc.settings.figdir = figure_dir")
dir.create(py$figure_dir, showWarnings = F, recursive = T)
######
summary_list = list("Selection" = list())
