r_utilities = "/home/haoy/projects/Analysis/code_library/utilities.r"
source(r_utilities)
###
# Edited Seurat(4.0.2)'s buggable function
ReadMtx = function (mtx, cells, features, cell.column = 1, feature.column = 2, 
                    skip.cell = 0, skip.feature = 0, unique.features = TRUE, 
                    strip.suffix = FALSE) 
{
  all.files <- list(`expression matrix` = mtx, `barcode list` = cells, 
                    `feature list` = features)
  cell.barcodes <- read.table(file = all.files[["barcode list"]], 
                              header = FALSE, sep = "\t", row.names = NULL, skip = skip.cell)
  feature.names <- read.table(file = all.files[["feature list"]], 
                              header = FALSE, sep = "\t", row.names = NULL, skip = skip.feature)
  bcols <- ncol(x = cell.barcodes)
  if (bcols < cell.column) {
    stop("cell.column was set to ", cell.column, " but ", 
         cells, " only has ", bcols, " columns.", 
         " Try setting the cell.column argument to a value <= to ", 
         bcols, ".")
  }
  cell.names <- cell.barcodes[, cell.column]
  if (all(grepl(pattern = "\\-1$", x = cell.names)) & 
      strip.suffix) {
    cell.names <- as.vector(x = as.character(x = sapply(X = cell.names, 
                                                        FUN = ExtractField, field = 1, delim = "-")))
  }
  fcols <- ncol(x = feature.names)
  if (fcols < feature.column) {
    stop("feature.column was set to ", feature.column, 
         " but ", features, " only has ", fcols, 
         " column(s).", " Try setting the feature.column argument to a value <= to ", 
         fcols, ".")
  }
  if (any(is.na(x = feature.names[, feature.column]))) {
    na.features <- which(x = is.na(x = feature.names[, feature.column]))
    replacement.column <- ifelse(test = feature.column == 
                                   2, yes = 1, no = 2)
    if (replacement.column > fcols) {
      stop("Some features names are NA in column ", 
           feature.column, ". Try specifiying a different column.", 
           call. = FALSE)
    }
    else {
      warning("Some features names are NA in column ", 
              feature.column, ". Replacing NA names with ID from column ", 
              replacement.column, ".", call. = FALSE)
    }
    feature.names[na.features, feature.column] <- feature.names[na.features, 
                                                                replacement.column]
  }
  feature.names <- feature.names[, feature.column]
  if (unique.features) {
    feature.names <- make.unique(names = feature.names)
  }
  data <- readMM(file = all.files[["expression matrix"]])
  if (length(x = cell.names) != ncol(x = data)) {
    stop("Matrix has ", ncol(data), " columns but found ", 
         length(cell.names), " barcodes. ", ifelse(test = length(x = cell.names) > 
                                                     ncol(x = data), yes = "Try increasing `skip.cell`. ", 
                                                   no = ""), call. = FALSE)
  }
  if (length(x = feature.names) != nrow(x = data)) {
    stop("Matrix has ", ncol(data), " rows but found ", 
         length(feature.names), " features. ", ifelse(test = length(x = feature.names) > 
                                                        nrow(x = data), yes = "Try increasing `skip.feature`. ", 
                                                      no = ""), call. = FALSE)
  }
  colnames(x = data) <- cell.names
  rownames(x = data) <- feature.names
  data <- as(data, Class = "dgCMatrix")
  return(data)
}
ReadMtx_1 = function(data_dir){
  return(ReadMtx(mtx = paste0(data_dir, "/matrix.mtx.gz"), cells = paste0(data_dir, "/barcodes.tsv.gz"), features = paste0(data_dir, "/features.tsv.gz")))
}
######
output_dir = "/home/haoy/projects/MetaStudiesAnalysis/Integration/Data/Zou_Liu_DevelopmentalCell_2020/convert_output"
dir.create(output_dir, showWarnings = F, recursive = T)
loom_path = paste0(output_dir, "/Zou_Liu_DevelopmentalCell_2020.loom")
metadata_path = paste0(output_dir, "/Zou_Liu_DevelopmentalCell_2020_metadata.tsv")
######
gene_bc_list = list()
cell_id_list = list()
donor_id = c()
age = c()
for(ii in 118996:119004){
  this_sample = paste0("HRS", ii)
  this_path = paste0("/home/haoy/projects/MetaStudiesAnalysis/Integration/Data/Zou_Liu_DevelopmentalCell_2020/", this_sample, "/outs/filtered_feature_bc_matrix")
  gene_bc_list[[this_sample]] = ReadMtx_1(data_dir = this_path)
  cell_id_list[[this_sample]] = paste0(this_sample, "_", colnames(gene_bc_list[[this_sample]]))
  if(ii == 118996){
    donor_id = c(donor_id, rep("Skin-Y-18", length(cell_id_list[[this_sample]])))
    age = c(age, rep(18, length(cell_id_list[[this_sample]])))
  }else if(ii == 118997){
    donor_id = c(donor_id, rep("Skin-Y-22", length(cell_id_list[[this_sample]])))
    age = c(age, rep(22, length(cell_id_list[[this_sample]])))
  }else if(ii == 118998){
    donor_id = c(donor_id, rep("Skin-Y-23", length(cell_id_list[[this_sample]])))
    age = c(age, rep(23, length(cell_id_list[[this_sample]])))
  }else if(ii == 118999){
    donor_id = c(donor_id, rep("Skin-M-44", length(cell_id_list[[this_sample]])))
    age = c(age, rep(44, length(cell_id_list[[this_sample]])))
  }else if(ii == 119000){
    donor_id = c(donor_id, rep("Skin-M-47", length(cell_id_list[[this_sample]])))
    age = c(age, rep(47, length(cell_id_list[[this_sample]])))
  }else if(ii == 119001){
    donor_id = c(donor_id, rep("Skin-M-48", length(cell_id_list[[this_sample]])))
    age = c(age, rep(48, length(cell_id_list[[this_sample]])))
  }else if(ii == 119002){
    donor_id = c(donor_id, rep("Skin-O-70", length(cell_id_list[[this_sample]])))
    age = c(age, rep(70, length(cell_id_list[[this_sample]])))
  }else if(ii == 119003){
    donor_id = c(donor_id, rep("Skin-O-73", length(cell_id_list[[this_sample]])))
    age = c(age, rep(73, length(cell_id_list[[this_sample]])))
  }else if(ii == 119004){
    donor_id = c(donor_id, rep("Skin-O-76", length(cell_id_list[[this_sample]])))
    age = c(age, rep(76, length(cell_id_list[[this_sample]])))
  }else{
    stop("Invalid ii")
  }
}
this_matrix = merge_matrices(gene_bc_list, is_sparse = T)
rm(gene_bc_list)
this_obj = CreateSeuratObject(this_matrix, min.cells = 1)
rm(this_matrix)
######
######
experiment_id = rep("Zou_Liu_DevelopmentalCell_2020", ncol(this_obj))
technology = rep("10X", ncol(this_obj))
isolation_site = rep("Droplet", ncol(this_obj))
body_site = rep("Eyelids", ncol(this_obj)) # Hair Transplantation Receivers'micrografts
tissue = rep(NA, ncol(this_obj))
sex = rep("F", ncol(this_obj))
original_annotation = rep(NA, ncol(this_obj))
condition = rep("Healthy", ncol(this_obj))
sample_status = rep(NA, ncol(this_obj))
this_obj[["DonorID"]] = paste0(experiment_id, "___", donor_id)
this_obj[["Tissue"]] = tissue
this_obj[["Site"]] = body_site
this_obj[["Age"]] = age
this_obj[["Sex"]] = sex
this_obj[["OriginalAnnotation"]] = original_annotation
this_obj[["ExperimentID"]] = experiment_id
this_obj[["Technology"]] = technology
this_obj[["IsolationSite"]] = isolation_site
this_obj[["Condition"]] = condition
this_obj[["SampleStatus"]] = sample_status
this_obj[["percent.mt"]] = PercentageFeatureSet(this_obj, pattern = "^MT-")
filter_cell =
  this_obj$nFeature_RNA >= 500 &
  this_obj$nFeature_RNA <= 5000 &
  this_obj[["percent.mt"]] <= 10
filtered_obj = this_obj[, filter_cell]
filter_gene = rowSums(filtered_obj) > 0
filtered_obj = filtered_obj[filter_gene, ]
if(file.exists(loom_path)){
  unlink(loom_path)
}
save_h5(loom_path, t(as_matrix(filtered_obj@assays$RNA@counts)))
write.table(filtered_obj@meta.data, metadata_path, sep = "\t", row.names = T, col.names = T, quote = F)
