r_utilities = "D:/Study/KI/Projects/###/Analysis/code_library/utilities.r"
py_utilities = "D:/Study/KI/Projects/###/Analysis/code_library/utilities.py"
source(r_utilities)
library(reticulate)
py_run_file(py_utilities, convert = F)
######
X_GeneName_path = "D:/Study/KI/Projects/###/Analysis/useful_files/X_GeneName_Homo_sapiens.GRCh38.103.rds"
Y_GeneName_path = "D:/Study/KI/Projects/###/Analysis/useful_files/Y_GeneName_Homo_sapiens.GRCh38.103.rds"
#
hgnc_mapping_alias_genename = readRDS("D:/Study/KI/Projects/###/Analysis/useful_files/hgnc_mapping_alias_genename_seurat.rds")
hgnc_mapping_genename_ensemblID = readRDS("D:/Study/KI/Projects/###/Analysis/useful_files/hgnc_mapping_genename_ensemblID.rds")
GRCh38v103_mapping_genename_ensemblID = readRDS("D:/Study/KI/Projects/###/Analysis/useful_files/GRCh38v103_mapping_genename_ensemblID.rds")
GRCh37_mapping_genename_ensemblID = readRDS("D:/Study/KI/Projects/###/Analysis/useful_files/GRCh37_mapping_genename_ensemblID.rds")
IntegratedAnnotation_2_metadata = readRDS("D:/Study/KI/Projects/###/MetaStudies/Analysis/Integration/0_preprocessing/IntegratedAnnotation_2_metadata.rds")
IndividualAnalysis_metadata_list = readRDS("D:/Study/KI/Projects/###/MetaStudies/Analysis/Integration/0_preprocessing/IndividualAnalysis_metadata_list.rds")
#
# IO - Setting
#
py$output_dir = "D:/Study/KI/Projects/###/MetaStudies/Data/Skin/Human/data_preparation_HCA_v1_0/Theocharidis_Bhasin_NatureCommunications_2022"
py_run_string("os.makedirs(output_dir, exist_ok=True)")
h5ad_ALL_file = "GSE165816_ALL.h5ad"
metadata_ALL_file = "GSE165816_ALL_metadata.tsv"
h5ad_SKIN_file = "GSE165816_SKIN.h5ad"
metadata_SKIN_file = "GSE165816_SKIN_metadata.tsv"
h5ad_SKIN_HC_file = "GSE165816_SKIN_HC.h5ad"
metadata_SKIN_HC_file = "GSE165816_SKIN_HC_metadata.tsv"
######
input_folder = "D:/Study/KI/Projects/###/MetaStudies/Data/Skin/Human/Theocharidis_Bhasin_NatureCommunications_2022/GSE165816_RAW"
lib_id_c = c()
lib_id_repo_c = c()
py_run_string("GEX_dict = {}")
for(this_file in list.files(input_folder, pattern = "counts.csv.gz")){
  this_path = paste0(input_folder, "/", this_file)
  str_fragments = unlist(strsplit(gsub("counts.csv.gz", "", this_file, fixed = T), "_", fixed = T))
  this_lib_repo = str_fragments[1]
  this_lib = str_fragments[2]
  tmp_df = as.data.frame(fread(this_path))
  py$tmp_df = tmp_df[, -1]
  py$tmp_var_df = data.frame(row.names = dplyr::case_match(
    as.character(tmp_df[, 1]),
    "1-Mar" ~ "MTARC1", # MARC1
    "1-Mar.1" ~ "MARCHF1", # MARCH1
    "2-Mar" ~ "MTARC2", # MARC2
    "2-Mar.1" ~ "MARCHF2", # MARCH2
    "3-Mar" ~ "MARCHF3",
    "4-Mar" ~ "MARCHF4",
    "5-Mar" ~ "MARCHF5",
    "6-Mar" ~ "MARCHF6",
    "7-Mar" ~ "MARCHF7",
    "8-Mar" ~ "MARCHF8",
    "9-Mar" ~ "MARCHF9",
    "10-Mar" ~ "MARCHF10",
    "11-Mar" ~ "MARCHF11",
    "1-Sep" ~ "SEPTIN1",
    "2-Sep" ~ "SEPTIN2",
    "3-Sep" ~ "SEPTIN3",
    "4-Sep" ~ "SEPTIN4",
    "5-Sep" ~ "SEPTIN5",
    "6-Sep" ~ "SEPTIN6",
    "7-Sep" ~ "SEPTIN7",
    "8-Sep" ~ "SEPTIN8",
    "9-Sep" ~ "SEPTIN9",
    "10-Sep" ~ "SEPTIN10",
    "11-Sep" ~ "SEPTIN11",
    "12-Sep" ~ "SEPTIN12",
    "13-Sep" ~ "SEPTIN13",
    "14-Sep" ~ "SEPTIN14",
    .default = as.character(tmp_df[, 1])
  ))
  py_run_string(paste0("GEX_dict['", this_lib, "'] = ad.AnnData(X=scipy.sparse.csr_matrix(tmp_df.values.transpose()).astype(np.float32), obs=pd.DataFrame(index='", this_lib, "_' + tmp_df.columns.astype(str).values), var=tmp_var_df)"))
  lib_id_c = c(lib_id_c, rep(this_lib, eval(parse(text = paste0("py$GEX_dict$'", this_lib, "'$n_obs")))))
  lib_id_repo_c = c(lib_id_repo_c, rep(this_lib_repo, eval(parse(text = paste0("py$GEX_dict$'", this_lib, "'$n_obs")))))
  print(this_lib)
  rm(tmp_df)
  py_run_string("del tmp_df")
  py_run_string("del tmp_var_df")
  py_run_string("gc.collect()")
  gc()
}
#
py_run_string("adata_raw = ad.concat(GEX_dict, axis=0, join='outer')")
py_run_string("del GEX_dict")
GeneName = featureMapMerge("adata_raw", "adata_genename", hgnc_mapping_alias_genename)
py_run_string("del adata_raw")
GeneName_to_GeneID = mapping_GeneName_GeneID(GeneName, hgnc_mapping_genename_ensemblID)
GeneName_to_GeneID = mapping_GeneName_GeneID(GeneName_to_GeneID, GRCh38v103_mapping_genename_ensemblID)
names(GeneName_to_GeneID) = GeneName
GeneID = featureMapMerge("adata_genename", "adata", GeneName_to_GeneID)
py_run_string("del adata_genename")
py$var_data_frame = data.frame(row.names = GeneID, feature_is_filtered=F, feature_biotype="gene", feature_reference="NCBITaxon:9606", feature_name=names(GeneID))
py_run_string(paste0("var_data_frame.loc[var_data_frame.index.str.startswith('ERCC-'), 'feature_biotype'] = 'spike-in'"))
py_run_string(paste0("adata.var = var_data_frame.copy()"))
py_run_string("del var_data_frame")
py_run_string("gc.collect()")
#
# Metadata
#
HCA_metadata_list = load_HCA_metadata_xlsx("D:/Study/KI/Projects/###/MetaStudies/Data/Skin/Human/HCA_metadata/Theocharidis_2022_HCA_tier 1_metadata.xlsx")
HCA_metadata_list_original = HCA_metadata_list # Keep the original for checking later
rownames(HCA_metadata_list$DN) = HCA_metadata_list$DN[, "donor_id"]
HCA_metadata_list$DN[is.na(HCA_metadata_list$DN)] = "unknown"
rownames(HCA_metadata_list$SP) = HCA_metadata_list$SP[, "library_id"]
HCA_metadata_list$SP[is.na(HCA_metadata_list$SP)] = "unknown"
#
donor_id_c = HCA_metadata_list$SP[lib_id_c, "donor_id"]
# Sex Identification
X_GeneName = readRDS(X_GeneName_path)
Y_GeneName = readRDS(Y_GeneName_path)
adata_X = py$adata$X
dimnames(adata_X) = list(py$adata$obs_names$values, py$adata$var[["feature_name"]])
ALL_expression = rowSums(adata_X)
Y_expression = rowSums(adata_X[, which(colnames(adata_X) %in% Y_GeneName)])
X_expression = rowSums(adata_X[, which(colnames(adata_X) %in% X_GeneName)])
Y_ratio = Y_expression / ALL_expression
X_ratio = X_expression / ALL_expression
Y_ratio_aggregate = aggregate(Y_ratio, by = list(donor_id_c), FUN = function(x){
  return(mean(x, na.rm = T))
})
X_ratio_aggregate = aggregate(X_ratio, by = list(donor_id_c), FUN = function(x){
  return(mean(x, na.rm = T))
})
Y_X_ratio = Y_ratio_aggregate$x / X_ratio_aggregate$x
names(Y_X_ratio) = Y_ratio_aggregate$Group.1
pdf(paste0(py$output_dir, '/DonorID_Sex_barplot.pdf'), width = 10, height = 5)
par(mar=c(11, 4, 4, 4))
barplot(Y_X_ratio, las=2)
dev.off()
#
### Study level
py$adata$obs["study_id"] = "Theocharidis_Bhasin_NatureCommunications_2022"
### Experiment level
py$adata$obs["experiment_id"] = as.character(py$adata$obs[["study_id"]])
py$adata$obs["institute"] = HCA_metadata_list$SP[lib_id_c, "institute"]
py$adata$obs["author_batch_notes"] = HCA_metadata_list$SP[lib_id_c, "author_batch_notes"]
### Donor background level
py$adata$obs["organism"] = "Homo sapiens"
py$adata$obs["organism_ontology_term_id"] = HCA_metadata_list$DN[donor_id_c, "organism_ontology_term_id"]
py$adata$obs["donor_id"] = donor_id_c
py$adata$obs["donor_id_running"] = paste0(py$adata$obs[["experiment_id"]], "___", donor_id_c)
py$adata$obs["sex_ontology_term"] = HCA_metadata_list$DN[donor_id_c, "sex_ontology_term"]
py$adata$obs["sex_ontology_term_id"] = HCA_metadata_list$DN[donor_id_c, "sex_ontology_term_id"]
py$adata$obs["ethnicity_1"] = dplyr::case_match(
  lib_id_c,
  "G6" ~ "European Ancestry",
  "G11" ~ "European Ancestry",
  "G21" ~ "European Ancestry",
  "G4" ~ "European Ancestry",
  "G23" ~ "European Ancestry",
  "G12" ~ "European Ancestry",
  "G22" ~ "European Ancestry",
  "G2" ~ "European Ancestry",
  "G13" ~ "European Ancestry",
  "G25" ~ "European Ancestry",
  "G15" ~ "European Ancestry",
  "G17" ~ "European Ancestry",
  "G27" ~ "European Ancestry",
  "G9" ~ "European Ancestry",
  "G35" ~ "European Ancestry",
  "G37" ~ "European Ancestry",
  "G33" ~ "European Ancestry",
  "G34" ~ "European Ancestry",
  "G39" ~ "European Ancestry",
  "G7" ~ "African Ancestry",
  "G8" ~ "African Ancestry",
  "G49" ~ "European Ancestry",
  "G42" ~ "African Ancestry",
  "G45" ~ "European Ancestry",
  "G2A" ~ "European Ancestry",
  "G4A" ~ "European Ancestry",
  "G19" ~ "European Ancestry",
  "G46" ~ "European Ancestry",
  "G48" ~ "European Ancestry",
  "G47" ~ "European Ancestry",
  "G38" ~ "European Ancestry",
  "G41" ~ "European Ancestry",
  "G1A" ~ "European Ancestry",
  "G3A" ~ "European Ancestry",
  "G3" ~ "European Ancestry",
  "G5" ~ "European Ancestry",
  "G1" ~ "European Ancestry",
  "G20" ~ "European Ancestry",
  "G16" ~ "European Ancestry",
  "G18" ~ "European Ancestry",
  "G26" ~ "European Ancestry",
  "G28" ~ "European Ancestry",
  "G29" ~ "European Ancestry",
  "G30" ~ "European Ancestry",
  "G31" ~ "European Ancestry",
  "G32" ~ "European Ancestry",
  "G36" ~ "European Ancestry",
  "G50" ~ "European Ancestry",
  "G10" ~ "European Ancestry",
  "G24" ~ "European Ancestry",
  "G43" ~ "European Ancestry",
  "G14" ~ "European Ancestry",
  "G44" ~ "European Ancestry",
  "G40" ~ "European Ancestry",
  .default = lib_id_c
)
py$adata$obs["ethnicity_2"] = "unknown"
py$adata$obs["ethnicity_details"] = dplyr::case_match(
  lib_id_c,
  "G6" ~ "White",
  "G11" ~ "White",
  "G21" ~ "White",
  "G4" ~ "White",
  "G23" ~ "White",
  "G12" ~ "White",
  "G22" ~ "White",
  "G2" ~ "White",
  "G13" ~ "White",
  "G25" ~ "White",
  "G15" ~ "White",
  "G17" ~ "White",
  "G27" ~ "White",
  "G9" ~ "White",
  "G35" ~ "White",
  "G37" ~ "White",
  "G33" ~ "White",
  "G34" ~ "White",
  "G39" ~ "White",
  "G7" ~ "African American",
  "G8" ~ "African American",
  "G49" ~ "White",
  "G42" ~ "African American",
  "G45" ~ "White",
  "G2A" ~ "White",
  "G4A" ~ "White",
  "G19" ~ "White",
  "G46" ~ "White",
  "G48" ~ "White",
  "G47" ~ "White",
  "G38" ~ "White",
  "G41" ~ "White",
  "G1A" ~ "White",
  "G3A" ~ "White",
  "G3" ~ "White",
  "G5" ~ "White",
  "G1" ~ "White",
  "G20" ~ "White",
  "G16" ~ "White",
  "G18" ~ "White",
  "G26" ~ "White",
  "G28" ~ "White",
  "G29" ~ "White",
  "G30" ~ "White",
  "G31" ~ "White",
  "G32" ~ "White",
  "G36" ~ "White",
  "G50" ~ "White",
  "G10" ~ "White",
  "G24" ~ "White",
  "G43" ~ "White",
  "G14" ~ "White",
  "G44" ~ "White",
  "G40" ~ "White",
  .default = lib_id_c
)
py$adata$obs["self_reported_ethnicity_ontology_term_id"] = "unknown"
py$adata$obs["genotype"] = "unknown"
### Donor status level
py$adata$obs["age_years"] = dplyr::case_match(
  lib_id_c,
  "G6" ~ "54",
  "G11" ~ "54",
  "G21" ~ "54",
  "G4" ~ "63",
  "G23" ~ "63",
  "G12" ~ "63",
  "G22" ~ "63",
  "G2" ~ "35",
  "G13" ~ "35",
  "G25" ~ "35",
  "G15" ~ "37",
  "G17" ~ "37",
  "G27" ~ "37",
  "G9" ~ "34",
  "G35" ~ "34",
  "G37" ~ "34",
  "G33" ~ "60",
  "G34" ~ "60",
  "G39" ~ "51",
  "G7" ~ "53",
  "G8" ~ "53",
  "G49" ~ "81",
  "G42" ~ "76",
  "G45" ~ "63",
  "G2A" ~ "58",
  "G4A" ~ "58",
  "G19" ~ "58",
  "G46" ~ "68",
  "G48" ~ "68",
  "G47" ~ "68",
  "G38" ~ "71",
  "G41" ~ "81",
  "G1A" ~ "52",
  "G3A" ~ "52",
  "G3" ~ "52",
  "G5" ~ "52",
  "G1" ~ "56",
  "G20" ~ "56",
  "G16" ~ "62",
  "G18" ~ "62",
  "G26" ~ "62",
  "G28" ~ "74",
  "G29" ~ "74",
  "G30" ~ "74",
  "G31" ~ "75",
  "G32" ~ "75",
  "G36" ~ "75",
  "G50" ~ "66",
  "G10" ~ "59",
  "G24" ~ "59",
  "G43" ~ "68",
  "G14" ~ "44",
  "G44" ~ "51",
  "G40" ~ "34",
  .default = lib_id_c
)
age_years_character = py$adata$obs[["age_years"]]
age_years_numeric = as.numeric(age_years_character)
development_stage_ontology_term_id = rep(NA, length(age_years_character))
development_stage_ontology_term_id[age_years_numeric >= 0 & age_years_numeric <= 14] = "HsapDv:0000264"
development_stage_ontology_term_id[age_years_numeric >= 15 & age_years_numeric <= 19] = "HsapDv:0000268"
development_stage_ontology_term_id[age_years_numeric >= 20 & age_years_numeric <= 29] = "HsapDv:0000237"
development_stage_ontology_term_id[age_years_numeric >= 30 & age_years_numeric <= 39] = "HsapDv:0000238"
development_stage_ontology_term_id[age_years_numeric >= 40 & age_years_numeric <= 49] = "HsapDv:0000239"
development_stage_ontology_term_id[age_years_numeric >= 50 & age_years_numeric <= 59] = "HsapDv:0000240"
development_stage_ontology_term_id[age_years_numeric >= 60 & age_years_numeric <= 69] = "HsapDv:0000241"
development_stage_ontology_term_id[age_years_numeric >= 70 & age_years_numeric <= 79] = "HsapDv:0000242"
development_stage_ontology_term_id[age_years_numeric >= 80 & age_years_numeric <= 89] = "HsapDv:0000243"
development_stage_ontology_term_id[age_years_numeric >= 90 & age_years_numeric <= 99] = "HsapDv:0000244"
development_stage_ontology_term_id[age_years_numeric >= 100 & age_years_numeric <= 109] = "HsapDv:0000245"
development_stage_ontology_term_id[is.na(development_stage_ontology_term_id)] = "unknown"
development_stage = dplyr::case_match(
  development_stage_ontology_term_id,
  "HsapDv:0000264" ~ "pediatric stage",
  "HsapDv:0000268" ~ "15-19 year-old",
  "HsapDv:0000237" ~ "third decade stage",
  "HsapDv:0000238" ~ "fourth decade stage",
  "HsapDv:0000239" ~ "fifth decade stage",
  "HsapDv:0000240" ~ "sixth decade stage",
  "HsapDv:0000241" ~ "seventh decade stage",
  "HsapDv:0000242" ~ "eighth decade stage",
  "HsapDv:0000243" ~ "ninth decade stage",
  "HsapDv:0000244" ~ "tenth decade stage",
  "HsapDv:0000245" ~ "eleventh decade stage",
  .default = development_stage_ontology_term_id
)
age_range_c = dplyr::case_match(
  development_stage_ontology_term_id,
  "HsapDv:0000264" ~ "0-14",
  "HsapDv:0000268" ~ "15-19",
  "HsapDv:0000237" ~ "20-29",
  "HsapDv:0000238" ~ "30-39",
  "HsapDv:0000239" ~ "40-49",
  "HsapDv:0000240" ~ "50-59",
  "HsapDv:0000241" ~ "60-69",
  "HsapDv:0000242" ~ "70-79",
  "HsapDv:0000243" ~ "80-89",
  "HsapDv:0000244" ~ "90-99",
  "HsapDv:0000245" ~ "100-109",
  .default = development_stage_ontology_term_id
)
py$adata$obs["age_range"] = age_range_c
py$adata$obs["development_stage_ontology_term"] = development_stage
py$adata$obs["development_stage_ontology_term_id"] = HCA_metadata_list$SP[lib_id_c, "development_stage_ontology_term_id"]
py$adata$obs["disease"] = dplyr::case_match(
  lib_id_c,
  "G6" ~ "diabetes mellitus",
  "G11" ~ "diabetes mellitus",
  "G21" ~ "diabetes mellitus",
  "G4" ~ "diabetes mellitus",
  "G23" ~ "diabetes mellitus",
  "G12" ~ "diabetes mellitus",
  "G22" ~ "diabetes mellitus",
  "G2" ~ "diabetes mellitus",
  "G13" ~ "diabetes mellitus",
  "G25" ~ "diabetes mellitus",
  "G15" ~ "diabetes mellitus",
  "G17" ~ "diabetes mellitus",
  "G27" ~ "diabetes mellitus",
  "G9" ~ "diabetes mellitus",
  "G35" ~ "diabetes mellitus",
  "G37" ~ "diabetes mellitus",
  "G33" ~ "diabetes mellitus",
  "G34" ~ "diabetes mellitus",
  "G39" ~ "diabetes mellitus",
  "G7" ~ "diabetes mellitus",
  "G8" ~ "diabetes mellitus",
  "G49" ~ "diabetes mellitus",
  "G42" ~ "diabetes mellitus",
  "G45" ~ "diabetes mellitus",
  "G2A" ~ "diabetes mellitus",
  "G4A" ~ "diabetes mellitus",
  "G19" ~ "diabetes mellitus",
  "G46" ~ "diabetes mellitus",
  "G48" ~ "diabetes mellitus",
  "G47" ~ "diabetes mellitus",
  "G38" ~ "diabetes mellitus",
  "G41" ~ "diabetes mellitus",
  "G1A" ~ "diabetes mellitus",
  "G3A" ~ "diabetes mellitus",
  "G3" ~ "diabetes mellitus",
  "G5" ~ "diabetes mellitus",
  "G1" ~ "normal",
  "G20" ~ "normal",
  "G16" ~ "prediabetes syndrome",
  "G18" ~ "prediabetes syndrome",
  "G26" ~ "prediabetes syndrome",
  "G28" ~ "normal",
  "G29" ~ "normal",
  "G30" ~ "normal",
  "G31" ~ "Hypertension",
  "G32" ~ "Hypertension",
  "G36" ~ "Hypertension",
  "G50" ~ "Hypertension",
  "G10" ~ "normal",
  "G24" ~ "normal",
  "G43" ~ "Hypertension",
  "G14" ~ "hypothyroidism",
  "G44" ~ "Asthma",
  "G40" ~ "normal",
  .default = lib_id_c
)
py$adata$obs["disease_ontology_term_id"] = HCA_metadata_list$SP[lib_id_c, "disease_ontology_term_id"]
disease_c = as.character(py$adata$obs[["disease"]])
py$adata$obs["disease_status"] = dplyr::case_match(
  disease_c,
  "Asthma" ~ "unknown",
  "hypothyroidism" ~ "unknown",
  "Hypertension" ~ "unknown",
  "prediabetes syndrome" ~ "unknown",
  "diabetes mellitus" ~ "unknown",
  "normal" ~ "not applicable",
  .default = disease_c
)
py$adata$obs["treatment_status"] = dplyr::case_match(
  lib_id_c,
  "G6" ~ "heparin,losartan potassium,hydrochlorothiazide,atorvastatin,insulin SC,vancomycin,ceftriaxone,metronidazole,clopidogrel",
  "G11" ~ "heparin,losartan potassium,hydrochlorothiazide,atorvastatin,insulin SC,vancomycin,ceftriaxone,metronidazole,clopidogrel",
  "G21" ~ "heparin,losartan potassium,hydrochlorothiazide,atorvastatin,insulin SC,vancomycin,ceftriaxone,metronidazole,clopidogrel",
  "G4" ~ "duraglutide,insulin glargine,metformin,metoprolol succinate,minocycline,olmesartan-hydrochlorothiazide,simvastatin",
  "G23" ~ "duraglutide,insulin glargine,metformin,metoprolol succinate,minocycline,olmesartan-hydrochlorothiazide,simvastatin",
  "G12" ~ "duraglutide,insulin glargine,metformin,metoprolol succinate,minocycline,olmesartan-hydrochlorothiazide,simvastatin",
  "G22" ~ "duraglutide,insulin glargine,metformin,metoprolol succinate,minocycline,olmesartan-hydrochlorothiazide,simvastatin",
  "G2" ~ "amlodipine,fenofibrate,glipizide,hydroxyzine pamoate,lisinopril,metformin,omeprazole,tramadol",
  "G13" ~ "amlodipine,fenofibrate,glipizide,hydroxyzine pamoate,lisinopril,metformin,omeprazole,tramadol",
  "G25" ~ "amlodipine,fenofibrate,glipizide,hydroxyzine pamoate,lisinopril,metformin,omeprazole,tramadol",
  "G15" ~ "albuterol sulfate,clonazepam,duloxetine,hydroxyzine pamoate,lamotrigine,levonorgestrel,lidocaine,meloxicam,nystatin,olanzepine,prazosin,pregabalin,sumatriptan succinate",
  "G17" ~ "albuterol sulfate,clonazepam,duloxetine,hydroxyzine pamoate,lamotrigine,levonorgestrel,lidocaine,meloxicam,nystatin,olanzepine,prazosin,pregabalin,sumatriptan succinate",
  "G27" ~ "albuterol sulfate,clonazepam,duloxetine,hydroxyzine pamoate,lamotrigine,levonorgestrel,lidocaine,meloxicam,nystatin,olanzepine,prazosin,pregabalin,sumatriptan succinate",
  "G9" ~ "amoxicilllin,erythromycin,insulin lispro,lisinopril,naloxone,ondansetron,oxycodone,polyethylene glycol,torsemide",
  "G35" ~ "amoxicilllin,erythromycin,insulin lispro,lisinopril,naloxone,ondansetron,oxycodone,polyethylene glycol,torsemide",
  "G37" ~ "amoxicilllin,erythromycin,insulin lispro,lisinopril,naloxone,ondansetron,oxycodone,polyethylene glycol,torsemide",
  "G33" ~ "metformin,gemfibrozol,HCTZ,lisinopril,omega 3,jardiance,Lantus",
  "G34" ~ "metformin,gemfibrozol,HCTZ,lisinopril,omega 3,jardiance,Lantus",
  "G39" ~ "amoxicillin,insulin glargine,insulin glulisine,levothyroxine",
  "G7" ~ "pregabalin,metformin,alprazolam,lisinopril,Benadryl,Victoza,Lantus,hydrocodone,iron,B12,melatonin,hydrocodone-acetaminophen,liraglutide",
  "G8" ~ "pregabalin,metformin,alprazolam,lisinopril,Benadryl,Victoza,Lantus,hydrocodone,iron,B12,melatonin,hydrocodone-acetaminophen,liraglutide",
  "G49" ~ "amlodipinen,BD ultra-fine NDL,calcitriol,febuxostat,hydralazine,insulin aspart,insulin glargine,isosorbide mononitrate,rosuvastatin,testosterone,torsemide,valsartan",
  "G42" ~ "torsemide,rosuvastatin,apixaban,budesonide,carvedilol,ciprofloxacin,clindamycin,clopidogrel,insulin,lisinopril,metformin,omeprazole,ondansetron,pantoprazole",
  "G45" ~ "Neurontin,metformin,glipizide,lisinopril,metoprolol",
  "G2A" ~ "albuterol sulfate,insulin glargine,insulin lispro,lisinopril,metformin,aripiprazole,atorvastatin,semaglutide,venlafaxine",
  "G4A" ~ "albuterol sulfate,insulin glargine,insulin lispro,lisinopril,metformin,aripiprazole,atorvastatin,semaglutide,venlafaxine",
  "G19" ~ "albuterol sulfate,insulin glargine,insulin lispro,lisinopril,metformin,aripiprazole,atorvastatin,semaglutide,venlafaxine",
  "G46" ~ "atenolol,gabapentin,icosapent ethyl,insulin glargine,losartan,metformin,omeprazole,prep,resmed vpap,EERS,sertaline,simvastatin,tadalafil,tamsulosin,trazodone",
  "G48" ~ "atenolol,gabapentin,icosapent ethyl,insulin glargine,losartan,metformin,omeprazole,prep,resmed vpap,EERS,sertaline,simvastatin,tadalafil,tamsulosin,trazodone",
  "G47" ~ "atenolol,gabapentin,icosapent ethyl,insulin glargine,losartan,metformin,omeprazole,prep,resmed vpap,EERS,sertaline,simvastatin,tadalafil,tamsulosin,trazodone",
  "G38" ~ "Sodium chloride,heparin,glipizide,lisinopril,insulin,vancomycin",
  "G41" ~ "allopurinol,lipitor,diltiazem,doxycycline,finasteride,gabapentin,levoxyl,tradjenta,metoprolol,oxybutynin,torsemide,coumadin",
  "G1A" ~ "citalopram,clopidogrel,cyclobenzaprine,duloxetine,folic acid,gabapentin,meloxicam,pentosan polysulfate sodium,ranitidine HCl,rosuvastatin",
  "G3A" ~ "citalopram,clopidogrel,cyclobenzaprine,duloxetine,folic acid,gabapentin,meloxicam,pentosan polysulfate sodium,ranitidine HCl,rosuvastatin",
  "G3" ~ "Albuterol sulfate,alendronate,atorvastatin,clonazepam,diazepam,dulaglutide,duloxetine,furosemide,gabapentin,hydrochlorothiazide,levothyroxine,lisinopril,metformin,omeprazole,tramadol,zonisamide",
  "G5" ~ "Albuterol sulfate,alendronate,atorvastatin,clonazepam,diazepam,dulaglutide,duloxetine,furosemide,gabapentin,hydrochlorothiazide,levothyroxine,lisinopril,metformin,omeprazole,tramadol,zonisamide",
  "G1" ~ "levonorgestrel,tretinoin,norethindrone acetate,clobetasol,minocycline,celecoxib",
  "G20" ~ "levonorgestrel,tretinoin,norethindrone acetate,clobetasol,minocycline,celecoxib",
  "G16" ~ "cefadroxil,aspirin",
  "G18" ~ "cefadroxil,aspirin",
  "G26" ~ "cefadroxil,aspirin",
  "G28" ~ "No",
  "G29" ~ "No",
  "G30" ~ "No",
  "G31" ~ "amlodipine,atenolol,Lipitor",
  "G32" ~ "amlodipine,atenolol,Lipitor",
  "G36" ~ "amlodipine,atenolol,Lipitor",
  "G50" ~ "Levothyroxine,fluoxetine,buspirone",
  "G10" ~ "No",
  "G24" ~ "No",
  "G43" ~ "diclofenac sodium,amlodipine,losartan,naproxen,sertaline",
  "G14" ~ "Levoxyl,meclizine,sertraline",
  "G44" ~ "hydrochlorothiazide",
  "G40" ~ "No",
  .default = lib_id_c
)
#
py$adata$obs["smoking_status"] = "unknown"
py$adata$obs["smoking_history"] = "unknown"
py$adata$obs["bmi"] = dplyr::case_match(
  lib_id_c,
  "G6" ~ "28.4470635484056",
  "G11" ~ "28.4470635484056",
  "G21" ~ "28.4470635484056",
  "G4" ~ "40.73528039602",
  "G23" ~ "40.73528039602",
  "G12" ~ "40.73528039602",
  "G22" ~ "40.73528039602",
  "G2" ~ "34.8615974857911",
  "G13" ~ "34.8615974857911",
  "G25" ~ "34.8615974857911",
  "G15" ~ "45.4788944583558",
  "G17" ~ "45.4788944583558",
  "G27" ~ "45.4788944583558",
  "G9" ~ "71.2244448420232",
  "G35" ~ "71.2244448420232",
  "G37" ~ "71.2244448420232",
  "G33" ~ "31.8041864604543",
  "G34" ~ "31.8041864604543",
  "G39" ~ "40.8858238235705",
  "G7" ~ "36.793487003764",
  "G8" ~ "36.793487003764",
  "G49" ~ "30.0636984859662",
  "G42" ~ "27.7961202546646",
  "G45" ~ "33.4656015825637",
  "G2A" ~ "25.9813658228334",
  "G4A" ~ "25.9813658228334",
  "G19" ~ "25.9813658228334",
  "G46" ~ "34.4018432360667",
  "G48" ~ "34.4018432360667",
  "G47" ~ "34.4018432360667",
  "G38" ~ "43.8820876509314",
  "G41" ~ "25.7975821394854",
  "G1A" ~ "22.1471093794421",
  "G3A" ~ "22.1471093794421",
  "G3" ~ "40.1427811123724",
  "G5" ~ "40.1427811123724",
  "G1" ~ "25.6906964686289",
  "G20" ~ "25.6906964686289",
  "G16" ~ "31.6218683443778",
  "G18" ~ "31.6218683443778",
  "G26" ~ "31.6218683443778",
  "G28" ~ "25.9945590961862",
  "G29" ~ "25.9945590961862",
  "G30" ~ "25.9945590961862",
  "G31" ~ "33.2319866830414",
  "G32" ~ "33.2319866830414",
  "G36" ~ "33.2319866830414",
  "G50" ~ "24.3659310665817",
  "G10" ~ "30.2093893294063",
  "G24" ~ "30.2093893294063",
  "G43" ~ "26.6291447668637",
  "G14" ~ "21.8824422706597",
  "G44" ~ "28.184521476585",
  "G40" ~ "unknown",
  .default = lib_id_c
)
py$adata$obs["systemic_disorder"] = "unknown"
py$adata$obs["autoimmune_component"] = "unknown"
py$adata$obs["menopause_status"] = "unknown"
py$adata$obs["puperty_status"] = "unknown"
py$adata$obs["sun_exposure"] = "unknown"
py$adata$obs["skin_tone"] = "unknown"
py$adata$obs["skin_care"] = "unknown"
py$adata$obs["manner_of_death"] = HCA_metadata_list$DN[donor_id_c, "manner_of_death"]
### Sample level
py$adata$obs["sample_id"] = HCA_metadata_list$SP[lib_id_c, "sample_id"]
py$adata$obs["sample_id_running"] = paste0(py$adata$obs[["experiment_id"]], "___", lib_id_c)
py$adata$obs["sample_source"] = HCA_metadata_list$SP[lib_id_c, "sample_source"]
py$adata$obs["sample_collection_method"] = HCA_metadata_list$SP[lib_id_c, "sample_collection_method"]
py$adata$obs["sampled_site_condition"] = HCA_metadata_list$SP[lib_id_c, "sampled_site_condition"]
py$adata$obs["anatomical_region_level_1"] = dplyr::case_match(
  lib_id_c,
  "G6" ~ "Extremities",
  "G11" ~ "Extremities",
  "G21" ~ "not applicable",
  "G4" ~ "Extremities",
  "G23" ~ "Extremities",
  "G12" ~ "Extremities",
  "G22" ~ "not applicable",
  "G2" ~ "Extremities",
  "G13" ~ "Extremities",
  "G25" ~ "not applicable",
  "G15" ~ "Extremities",
  "G17" ~ "Extremities",
  "G27" ~ "not applicable",
  "G9" ~ "Extremities",
  "G35" ~ "Extremities",
  "G37" ~ "not applicable",
  "G33" ~ "Extremities",
  "G34" ~ "Extremities",
  "G39" ~ "Extremities",
  "G7" ~ "Extremities",
  "G8" ~ "Extremities",
  "G49" ~ "Extremities",
  "G42" ~ "Extremities",
  "G45" ~ "Extremities",
  "G2A" ~ "Extremities",
  "G4A" ~ "Extremities",
  "G19" ~ "not applicable",
  "G46" ~ "Extremities",
  "G48" ~ "Extremities",
  "G47" ~ "not applicable",
  "G38" ~ "Extremities",
  "G41" ~ "Extremities",
  "G1A" ~ "Extremities",
  "G3A" ~ "Extremities",
  "G3" ~ "Extremities",
  "G5" ~ "Extremities",
  "G1" ~ "Extremities",
  "G20" ~ "not applicable",
  "G16" ~ "Extremities",
  "G18" ~ "Extremities",
  "G26" ~ "not applicable",
  "G28" ~ "Extremities",
  "G29" ~ "Extremities",
  "G30" ~ "not applicable",
  "G31" ~ "Extremities",
  "G32" ~ "Extremities",
  "G36" ~ "Extremities",
  "G50" ~ "Extremities",
  "G10" ~ "Extremities",
  "G24" ~ "Extremities",
  "G43" ~ "Extremities",
  "G14" ~ "Extremities",
  "G44" ~ "Extremities",
  "G40" ~ "Extremities",
  .default = lib_id_c
)
py$adata$obs["anatomical_region_level_2"] = dplyr::case_match(
  lib_id_c,
  "G6" ~ "Foot",
  "G11" ~ "Arm",
  "G21" ~ "not applicable",
  "G4" ~ "Foot",
  "G23" ~ "Foot",
  "G12" ~ "Arm",
  "G22" ~ "not applicable",
  "G2" ~ "Foot",
  "G13" ~ "Arm",
  "G25" ~ "not applicable",
  "G15" ~ "Foot",
  "G17" ~ "Arm",
  "G27" ~ "not applicable",
  "G9" ~ "Foot",
  "G35" ~ "Arm",
  "G37" ~ "not applicable",
  "G33" ~ "Foot",
  "G34" ~ "Foot",
  "G39" ~ "Foot",
  "G7" ~ "Foot",
  "G8" ~ "Foot",
  "G49" ~ "Foot",
  "G42" ~ "Foot",
  "G45" ~ "Foot",
  "G2A" ~ "Foot",
  "G4A" ~ "Arm",
  "G19" ~ "not applicable",
  "G46" ~ "Foot",
  "G48" ~ "Arm",
  "G47" ~ "not applicable",
  "G38" ~ "Foot",
  "G41" ~ "Foot",
  "G1A" ~ "Foot",
  "G3A" ~ "Foot",
  "G3" ~ "Foot",
  "G5" ~ "Foot",
  "G1" ~ "Arm",
  "G20" ~ "not applicable",
  "G16" ~ "Foot",
  "G18" ~ "Arm",
  "G26" ~ "not applicable",
  "G28" ~ "Foot",
  "G29" ~ "Arm",
  "G30" ~ "not applicable",
  "G31" ~ "Foot",
  "G32" ~ "Foot",
  "G36" ~ "Arm",
  "G50" ~ "Foot",
  "G10" ~ "Foot",
  "G24" ~ "Foot",
  "G43" ~ "Foot",
  "G14" ~ "Foot",
  "G44" ~ "Foot",
  "G40" ~ "Foot",
  .default = lib_id_c
)
py$adata$obs["anatomical_region_level_3"] = dplyr::case_match(
  lib_id_c,
  "G6" ~ "Dorsal foot",
  "G11" ~ "Forearm",
  "G21" ~ "not applicable",
  "G4" ~ "Dorsal foot",
  "G23" ~ "Dorsal foot",
  "G12" ~ "Forearm",
  "G22" ~ "not applicable",
  "G2" ~ "Dorsal foot",
  "G13" ~ "Forearm",
  "G25" ~ "not applicable",
  "G15" ~ "Dorsal foot",
  "G17" ~ "Forearm",
  "G27" ~ "not applicable",
  "G9" ~ "Dorsal foot",
  "G35" ~ "Forearm",
  "G37" ~ "not applicable",
  "G33" ~ "Dorsal foot",
  "G34" ~ "Dorsal foot",
  "G39" ~ "Dorsal foot",
  "G7" ~ "Dorsal foot",
  "G8" ~ "Dorsal foot",
  "G49" ~ "Dorsal foot",
  "G42" ~ "Dorsal foot",
  "G45" ~ "Dorsal foot",
  "G2A" ~ "Dorsal foot",
  "G4A" ~ "Forearm",
  "G19" ~ "not applicable",
  "G46" ~ "Dorsal foot",
  "G48" ~ "Forearm",
  "G47" ~ "not applicable",
  "G38" ~ "Dorsal foot",
  "G41" ~ "Dorsal foot",
  "G1A" ~ "Dorsal foot",
  "G3A" ~ "Dorsal foot",
  "G3" ~ "Dorsal foot",
  "G5" ~ "Dorsal foot",
  "G1" ~ "Forearm",
  "G20" ~ "not applicable",
  "G16" ~ "Dorsal foot",
  "G18" ~ "Forearm",
  "G26" ~ "not applicable",
  "G28" ~ "Dorsal foot",
  "G29" ~ "Forearm",
  "G30" ~ "not applicable",
  "G31" ~ "Dorsal foot",
  "G32" ~ "Dorsal foot",
  "G36" ~ "Forearm",
  "G50" ~ "Dorsal foot",
  "G10" ~ "Dorsal foot",
  "G24" ~ "Dorsal foot",
  "G43" ~ "Dorsal foot",
  "G14" ~ "Dorsal foot",
  "G44" ~ "Dorsal foot",
  "G40" ~ "Dorsal foot",
  .default = lib_id_c
)
metadata_AR = as.character(py$adata$obs[["anatomical_region_level_3"]])
metadata_AR[metadata_AR %in% c("unknown", "not applicable")] = NA
metadata_AR[is.na(metadata_AR)] = as.character(py$adata$obs[["anatomical_region_level_2"]])[is.na(metadata_AR)]
metadata_AR[metadata_AR %in% c("unknown", "not applicable")] = NA
metadata_AR[is.na(metadata_AR)] = as.character(py$adata$obs[["anatomical_region_level_1"]])[is.na(metadata_AR)]
metadata_AR[metadata_AR %in% c("unknown", "not applicable")] = NA
metadata_AR[is.na(metadata_AR)] = "unknown"
py$adata$obs["anatomical_region"] = metadata_AR
py$adata$obs["tissue_type"] = HCA_metadata_list$SP[lib_id_c, "tissue_type"]
py$adata$obs["tissue_free_text"] = HCA_metadata_list$SP[lib_id_c, "tissue_free_text"]
py$adata$obs["tissue_ontology_term"] = HCA_metadata_list$SP[lib_id_c, "tissue_ontology_term"]
py$adata$obs["tissue_ontology_term_id"] = HCA_metadata_list$SP[lib_id_c, "tissue_ontology_term_id"]
py$adata$obs["skin_tissue"] = dplyr::case_match(
  lib_id_c,
  "G6" ~ "unknown",
  "G11" ~ "unknown",
  "G21" ~ "No",
  "G4" ~ "unknown",
  "G23" ~ "unknown",
  "G12" ~ "unknown",
  "G22" ~ "No",
  "G2" ~ "unknown",
  "G13" ~ "unknown",
  "G25" ~ "No",
  "G15" ~ "unknown",
  "G17" ~ "unknown",
  "G27" ~ "No",
  "G9" ~ "unknown",
  "G35" ~ "unknown",
  "G37" ~ "No",
  "G33" ~ "unknown",
  "G34" ~ "unknown",
  "G39" ~ "unknown",
  "G7" ~ "unknown",
  "G8" ~ "unknown",
  "G49" ~ "unknown",
  "G42" ~ "unknown",
  "G45" ~ "unknown",
  "G2A" ~ "unknown",
  "G4A" ~ "unknown",
  "G19" ~ "No",
  "G46" ~ "unknown",
  "G48" ~ "unknown",
  "G47" ~ "No",
  "G38" ~ "unknown",
  "G41" ~ "unknown",
  "G1A" ~ "unknown",
  "G3A" ~ "unknown",
  "G3" ~ "unknown",
  "G5" ~ "unknown",
  "G1" ~ "unknown",
  "G20" ~ "No",
  "G16" ~ "unknown",
  "G18" ~ "unknown",
  "G26" ~ "No",
  "G28" ~ "unknown",
  "G29" ~ "unknown",
  "G30" ~ "No",
  "G31" ~ "unknown",
  "G32" ~ "unknown",
  "G36" ~ "unknown",
  "G50" ~ "unknown",
  "G10" ~ "unknown",
  "G24" ~ "unknown",
  "G43" ~ "unknown",
  "G14" ~ "unknown",
  "G44" ~ "unknown",
  "G40" ~ "unknown",
  .default = lib_id_c
)
py$adata$obs["sample_collection_site"] = HCA_metadata_list$SP[lib_id_c, "sample_collection_site"]
py$adata$obs["sample_cultured"] = "No"
### Sequencing level
py$adata$obs["sequencing_type"] = "GEX"
py$adata$obs["sample_collection_year"] = HCA_metadata_list$SP[lib_id_c, "sample_collection_year"]
py$adata$obs["sample_collection_relative_time_point"] = HCA_metadata_list$SP[lib_id_c, "sample_collection_relative_time_point"]
py$adata$obs["library_id"] = lib_id_c
py$adata$obs["library_id_repository"] = lib_id_repo_c
py$adata$obs["library_preparation_batch"] = HCA_metadata_list$SP[lib_id_c, "library_preparation_batch"]
py$adata$obs["library_sequencing_run"] = HCA_metadata_list$SP[lib_id_c, "library_sequencing_run"]
py$adata$obs["sample_preservation_method"] = HCA_metadata_list$SP[lib_id_c, "sample_preservation_method"]
py$adata$obs["protocol_url"] = "unknown"
py$adata$obs["dissociation_protocol"] = HCA_metadata_list$SP[lib_id_c, "dissociation_protocol"]
py$adata$obs["cell_enrichment"] = HCA_metadata_list$SP[lib_id_c, "cell_enrichment"]
py$adata$obs["cell_viability_percentage"] = HCA_metadata_list$SP[lib_id_c, "cell_viability_percentage"]
py$adata$obs["cell_number_loaded"] = HCA_metadata_list$SP[lib_id_c, "cell_number_loaded"]
py$adata$obs["suspension_type"] = HCA_metadata_list$SP[lib_id_c, "suspension_type"]
py$adata$obs["assay_ontology_term"] = "10x 3' v2"
py$adata$obs["assay_ontology_term_id"] = "EFO:0009899"
py$adata$obs["sequenced_fragment"] = "3 prime tag"
py$adata$obs["sequencing_platform"] = "Illumina NovaSeq 6000"
### Preprocessing level
py$adata$obs["reference_genome"] = "GRCh38"
py$adata$obs["alignment_software"] = "cell ranger"
py$adata$obs["intron_inclusion"] = "unknown"
### Analysis Level
py$adata$obs["is_primary_data"] = FALSE
py$adata$obs["author_cell_type_1"] = "unknown"
py$adata$obs["author_cell_type_2"] = "unknown"
py$adata$obs["author_cell_type_3"] = "unknown"
py$adata$obs["author_cell_type_4"] = "unknown"
py$adata$obs["author_cell_type_5"] = "unknown"
py$adata$obs["author_cell_type_6"] = "unknown"
### Previous Barcode
this_dataset_individual = "Theocharidis 2022"
this_dataset_integrated = "Theocharidis_Bhasin_NatureCommunications_2022"
individual_barcode_used = rownames(IndividualAnalysis_metadata_list[[this_dataset_individual]])
names(individual_barcode_used) = individual_barcode_used
individual_barcode_all = individual_barcode_used[py$adata$obs_names$values]
individual_barcode_all[is.na(individual_barcode_all)] = "unknown"
py$adata$obs["individual_barcode"] = individual_barcode_all
if(sum(individual_barcode_all != "unknown") != length(individual_barcode_used)){
  stop("Individual barcodes didn't match well!")
}else{
  print("Pass: Individual")
}
###
integrated_barcode_used = rownames(IntegratedAnnotation_2_metadata)[sapply(rownames(IntegratedAnnotation_2_metadata), function(x){
  str_fragments = unlist(strsplit(x, "---", fixed = T))
  return(str_fragments[length(str_fragments)])
}) == this_dataset_integrated]
names(integrated_barcode_used) = sapply(integrated_barcode_used, function(x){
  str_fragments = unlist(strsplit(x, "---", fixed = T))
  return(paste(str_fragments[-length(str_fragments)], collapse = "---"))
})
integrated_barcode_all = integrated_barcode_used[py$adata$obs_names$values]
integrated_barcode_all[is.na(integrated_barcode_all)] = "unknown"
py$adata$obs["integrated_barcode"] = integrated_barcode_all
if(sum(integrated_barcode_all != "unknown") != length(integrated_barcode_used)){
  stop("Integrated barcodes didn't match well!")
}else{
  print("Pass: Integrated")
}
#
# Embeddings
#
#
# Put unique identifier for this experiment
#
py_run_string("adata.obs_names = adata.obs['experiment_id'].astype(str) + '___' + adata.obs_names.astype(str)")
#
# Check before writing
#
print(py$adata)
print(py$adata$obs)
print(py$adata$var)
#
# Write
#
py_run_string("sc.pp.filter_cells(adata, min_genes=1)")
py_run_string(paste0("adata.write(output_dir + '/", h5ad_ALL_file, "')"))
py_run_string(paste0("adata.obs.to_csv(output_dir + '/", metadata_ALL_file, "', sep='\t', index=True)"))
py_run_string("adata_SKIN = adata[adata.obs['skin_tissue'] != 'No']")
py_run_string(paste0("adata_SKIN.write(output_dir + '/", h5ad_SKIN_file, "')"))
py_run_string(paste0("adata_SKIN.obs.to_csv(output_dir + '/", metadata_SKIN_file, "', sep='\t', index=True)"))
py_run_string("adata_SKIN_HC = adata_SKIN[adata_SKIN.obs['sampled_site_condition'] == 'healthy']")
py_run_string(paste0("adata_SKIN_HC.write(output_dir + '/", h5ad_SKIN_HC_file, "')"))
py_run_string(paste0("adata_SKIN_HC.obs.to_csv(output_dir + '/", metadata_SKIN_HC_file, "', sep='\t', index=True)"))
py_run_string("print(output_dir)")
#
metadata_matrix = as.matrix(py$adata$obs)
write.table(x = HCA_metadata_usage(HCA_metadata=HCA_metadata_list_original$DT, ADATA_metadata=metadata_matrix,
                                   metadata_rowname_column="assay_ontology_term_id"),
            file = paste0(py$output_dir, "/", "HCA_metadata_usage_DT.txt"), sep = "\t", col.names=T, row.names=F, quote=T)
write.table(x = HCA_metadata_usage(HCA_metadata=HCA_metadata_list_original$DN, ADATA_metadata=metadata_matrix,
                                   metadata_rowname_column="donor_id"),
            file = paste0(py$output_dir, "/", "HCA_metadata_usage_DN.txt"), sep = "\t", col.names=T, row.names=F, quote=T)
write.table(x = HCA_metadata_usage(HCA_metadata=HCA_metadata_list_original$SP, ADATA_metadata=metadata_matrix,
                                   metadata_rowname_column="sample_id"),
            file = paste0(py$output_dir, "/", "HCA_metadata_usage_SP.txt"), sep = "\t", col.names=T, row.names=F, quote=T)

