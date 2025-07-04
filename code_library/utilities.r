library(Seurat)
library(data.table)
library(Matrix)
library(rhdf5)
library(hdf5r)
library(ggplot2)
library(dplyr)
library(R.utils)
library(patchwork)
library(ggrepel)
library(ggvenn)
library(openxlsx)
library(openxlsx2)

options(max.print=1000)
message("max.print is set as 1000!")

extract_gtf_information = function(gtf_path = NULL, rds_path = NULL, txt_path = NULL){
  if(is.null(gtf_path) & is.null(rds_path) & is.null(txt_path)){
    stop("One of 'gtf_path', 'rds_path' or 'txt_path' should be input!")
  }
  status = "loading"
  if(!is.null(rds_path)){
    if(file.exists(rds_path)){
      gtf_information = readRDS(rds_path)
      status = "loaded"
    }else{
      if(!is.null(txt_path)){
        if(file.exists(txt_path)){
          gtf_information = read.table(txt_path) ###
          status = "loaded"
        }
      }
    }
  }
  if(status == "loading"){
    gtf_file = fread(gtf_path, sep = "\t")
    gtf_information = t(apply(as.matrix(as.data.frame(gtf_file[, c(1, 9)])), 1, function(x){
      ChrID = x[1]
      tmp = unlist(strsplit(x[2], 'gene_id \"', fixed = T))[2]
      GeneID = unlist(strsplit(tmp, '\"', fixed = T))[1]
      tmp = unlist(strsplit(x[2], 'transcript_id \"', fixed = T))[2]
      TranscriptID  = unlist(strsplit(tmp, '\"', fixed = T))[1]
      tmp = unlist(strsplit(x[2], 'gene_name \"', fixed = T))[2]
      GeneName = unlist(strsplit(tmp, '\"', fixed = T))[1]
      return(c(ChrID, GeneID, TranscriptID, GeneName))
    }))
    rownames(gtf_information) = apply(gtf_information, 1, paste, collapse = "___")
    colnames(gtf_information) = c("ChrID", "GeneID", "TranscriptID", "GeneName")
    gtf_information = gtf_information[!is.na(gtf_information[, "TranscriptID"]), ]
    gtf_information = gtf_information[!duplicated(rownames(gtf_information)), ]
    rownames(gtf_information) = c()
    if(!is.null(rds_path)){
      saveRDS(gtf_information, rds_path)
    }
    if(!is.null(txt_path)){
      write.table(gtf_information, txt_path, sep = "\t", quote = F)
    }
  }
  return(gtf_information)
}
geneid_to_genename = function(gene_bc_mat, geneid, correspond_genename, verbose = F){
  # Will keep geneID if its correspond geneName is NA, use GeneName___GeneID for duplicated genes.
  ID_Name = unique(paste0(geneid, "__geneid_to_genename__", correspond_genename))
  ID_Name_matrix = sapply(ID_Name, function(x){
    return(unlist(strsplit(x, "__geneid_to_genename__", fixed = T)))
  })
  ID_Name_mapping = ID_Name_matrix[2, ]
  names(ID_Name_mapping) = ID_Name_matrix[1, ]
  ID_Name_mapping[ID_Name_mapping == "NA"] = names(ID_Name_mapping)[ID_Name_mapping == "NA"]
  ID_Name_mapping[duplicated(ID_Name_mapping)] = paste(ID_Name_mapping[duplicated(ID_Name_mapping)], names(ID_Name_mapping)[duplicated(ID_Name_mapping)], sep = "___")
  if(verbose){
    print(head(ID_Name_mapping))
    print(gene_bc_mat[1:5, 1:5])
  }
  target_genename = ID_Name_mapping[rownames(gene_bc_mat)]
  target_genename[is.na(target_genename)] = rownames(gene_bc_mat)[is.na(target_genename)]
  rownames(gene_bc_mat) = target_genename
  if(verbose){
    print(gene_bc_mat[1:5, 1:5])
  }
  return(gene_bc_mat)
}

if(packageVersion("Seurat") >= "4.0.2" & packageVersion("Seurat") <= "4.0.4"){
  # Edited Seurat(4.0.2 - 4.0.4)'s buggable function
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
    #data <- as(data, Class = "dgCMatrix")
    data <- as(data, Class = "CsparseMatrix")
    return(data)
  }
}

ReadMtx_1 = function(data_dir, prefix="", suffix="", matrix="matrix", barcodes="barcodes", features="features", ...){
  return(ReadMtx(mtx = paste0(data_dir, "/", prefix, matrix, suffix, ".mtx.gz"), cells = paste0(data_dir, "/", prefix, barcodes, suffix, ".tsv.gz"), features = paste0(data_dir, "/", prefix, features, suffix, ".tsv.gz"), ...))
}

ReadMtx_py = function(output_str, data_dir, prefix="", suffix="", matrix="matrix", barcodes="barcodes", features="features", feature_column_name = c("gene_ids", "feature_name", "feature_types"), is_compressed_file=T){
  py$feature_column_name = feature_column_name
  if(is_compressed_file){
    mtx_suffix = ".mtx.gz"
    tsv_suffix = ".tsv.gz"
  }else{
    mtx_suffix = ".mtx"
    tsv_suffix = ".tsv"
  }
  py_run_string(paste0(output_str, " = read_10x_mtx(mtx_path='", paste0(data_dir, "/", prefix, matrix, suffix, mtx_suffix), "', barcode_path='", paste0(data_dir, "/", prefix, barcodes, suffix, tsv_suffix), "', feature_path='", paste0(data_dir, "/", prefix, features, suffix, tsv_suffix), "', feature_column_name=feature_column_name)"))
  py_run_string("del feature_column_name")
  py_run_string("gc.collect()")
}

load_HCA_metadata_xlsx = function(HCA_metadata_xlsx_path){
  HCA_metadata_wb = openxlsx2::wb_load(HCA_metadata_xlsx_path)
  HCA_metadata_wb_DT_0 = HCA_metadata_wb %>% openxlsx2::wb_to_df(sheet = "Tier 1 Dataset Metadata")
  HCA_metadata_wb_DN_0 = HCA_metadata_wb %>% openxlsx2::wb_to_df(sheet = "Tier 1 Donor Metadata")
  HCA_metadata_wb_SP_0 = HCA_metadata_wb %>% openxlsx2::wb_to_df(sheet = "Tier 1 Sample Metadata")
  HCA_metadata_wb_CT_0 = HCA_metadata_wb %>% openxlsx2::wb_to_df(sheet = "Tier 1 Celltype Metadata")
  #
  HCA_metadata_wb_DT_1 = HCA_metadata_wb_DT_0[1 + which(HCA_metadata_wb_DT_0[, 1] == "FILL OUT INFORMATION BELOW THIS ROW"): nrow(HCA_metadata_wb_DT_0), ]
  HCA_metadata_wb_DT = HCA_metadata_wb_DT_1[rowSums(!is.na(HCA_metadata_wb_DT_1)) > 0, ]
  colnames(HCA_metadata_wb_DT) = gsub(" ", "", colnames(HCA_metadata_wb_DT), fixed = T)
  HCA_metadata_wb_DN_1 = HCA_metadata_wb_DN_0[1 + which(HCA_metadata_wb_DN_0[, 1] == "FILL OUT INFORMATION BELOW THIS ROW"): nrow(HCA_metadata_wb_DN_0), ]
  HCA_metadata_wb_DN = HCA_metadata_wb_DN_1[rowSums(!is.na(HCA_metadata_wb_DN_1)) > 0, ]
  colnames(HCA_metadata_wb_DN) = gsub(" ", "", colnames(HCA_metadata_wb_DN), fixed = T)
  HCA_metadata_wb_SP_1 = HCA_metadata_wb_SP_0[1 + which(HCA_metadata_wb_SP_0[, 1] == "FILL OUT INFORMATION BELOW THIS ROW"): nrow(HCA_metadata_wb_SP_0), ]
  HCA_metadata_wb_SP = HCA_metadata_wb_SP_1[rowSums(!is.na(HCA_metadata_wb_SP_1)) > 0, ]
  colnames(HCA_metadata_wb_SP) = gsub(" ", "", colnames(HCA_metadata_wb_SP), fixed = T)
  HCA_metadata_wb_CT_1 = HCA_metadata_wb_CT_0[1 + which(HCA_metadata_wb_CT_0[, 1] == "FILL OUT INFORMATION BELOW THIS ROW"): nrow(HCA_metadata_wb_CT_0), ]
  HCA_metadata_wb_CT = HCA_metadata_wb_CT_1[rowSums(!is.na(HCA_metadata_wb_CT_1)) > 0, ]
  colnames(HCA_metadata_wb_CT) = gsub(" ", "", colnames(HCA_metadata_wb_CT), fixed = T)
  message("Spaces (' ') in column names will be removed!")
  return(list(DT=HCA_metadata_wb_DT, DN=HCA_metadata_wb_DN, SP=HCA_metadata_wb_SP, CT=HCA_metadata_wb_CT))
}

HCA_metadata_usage = function(HCA_metadata, ADATA_metadata, metadata_rowname_column){
  HCA_metadata_matrix = as.matrix(HCA_metadata)
  rownames(HCA_metadata_matrix) = HCA_metadata_matrix[, metadata_rowname_column]
  HCA_metadata_matrix[is.na(HCA_metadata_matrix)] = "EMPTY"
  ###
  common_columns = colnames(HCA_metadata_matrix)[colnames(HCA_metadata_matrix) %in% colnames(ADATA_metadata)]
  not_used_columns = colnames(HCA_metadata_matrix)[!colnames(HCA_metadata_matrix) %in% colnames(ADATA_metadata)]
  HCA_metadata_CC = HCA_metadata_matrix[, common_columns, drop=F]
  #
  metadata_column_connector = "___METADATACOLUMN___"
  ADATA_metadata_unique_text = unique(apply(ADATA_metadata[, common_columns, drop=F], 1, function(x){
    return(paste(x, collapse = metadata_column_connector))
  }))
  #
  ADATA_metadata_CC = t(as.data.frame(x = strsplit(ADATA_metadata_unique_text, metadata_column_connector, fixed = T), row.names = common_columns))
  rownames(ADATA_metadata_CC) = ADATA_metadata_CC[, metadata_rowname_column]
  #
  common_rows = rownames(HCA_metadata_CC)[rownames(HCA_metadata_CC) %in% rownames(ADATA_metadata_CC)]
  not_used_rows = rownames(HCA_metadata_CC)[!rownames(HCA_metadata_CC) %in% rownames(ADATA_metadata_CC)]
  not_included_rows = rownames(ADATA_metadata_CC)[!rownames(ADATA_metadata_CC) %in% rownames(HCA_metadata_CC)]
  #
  HCA_metadata_latest = matrix(NA,
                               nrow = length(common_rows) + length(not_used_rows) + length(not_included_rows),
                               ncol = ncol(HCA_metadata_matrix),
                               dimnames = list(c(common_rows, not_used_rows, not_included_rows), colnames(HCA_metadata_matrix)))
  HCA_metadata_latest[c(common_rows, not_used_rows), ] = HCA_metadata_matrix[c(common_rows, not_used_rows), , drop=F]
  #
  HCA_metadata_latest[not_used_rows, ] = paste0("(", HCA_metadata_latest[not_used_rows, , drop=F], ")[RowNotUsed]")
  HCA_metadata_latest[common_rows, not_used_columns] = paste0("(", HCA_metadata_latest[common_rows, not_used_columns, drop=F], ")[ColNotUsed]")
  #
  HCA_metadata_CC_CR = HCA_metadata_CC[common_rows, ]
  ###
  metadata_latest_CC_CR = ADATA_metadata_CC[common_rows, ]
  metadata_latest_CC_CR_DM = metadata_latest_CC_CR != HCA_metadata_CC_CR
  metadata_latest_CC_CR[metadata_latest_CC_CR_DM] = paste0(metadata_latest_CC_CR[metadata_latest_CC_CR_DM], "(", HCA_metadata_CC_CR[metadata_latest_CC_CR_DM], ")")
  HCA_metadata_latest[common_rows, common_columns] = metadata_latest_CC_CR
  HCA_metadata_latest[not_included_rows, ] = "(MISSING)"
  HCA_metadata_latest[not_included_rows, common_columns] = paste0(ADATA_metadata_CC[not_included_rows, , drop=F], "(MISSING)")
  return(HCA_metadata_latest)
}

as_matrix = function(input_mat, chunk_size = 1024, verbose = F){
  tmp_mat = matrix(nrow = nrow(input_mat), ncol = ncol(input_mat), dimnames = dimnames(input_mat))
  for(ii in seq(ceiling(ncol(input_mat) / chunk_size))){
    writing_index = seq((ii - 1) * chunk_size + 1, min(c(ii * chunk_size, ncol(input_mat))))
    tmp_mat[, writing_index] = as.matrix(input_mat[, writing_index])
    if(verbose){
      print(paste0(min(writing_index), " - ", max(writing_index)))
    }
  }
  return(tmp_mat)
}
save_h5 = function(output_path, bc_gene_mat){
  if(file.exists(output_path)){
    unlink(output_path)
  }
  h5createFile(output_path)
  h5createGroup(output_path, "col_attrs")
  h5write(rownames(bc_gene_mat), output_path,"col_attrs/CellID")
  h5createGroup(output_path, "row_attrs")
  h5write(colnames(bc_gene_mat), output_path,"row_attrs/Gene")
  h5createGroup(output_path, "layers")
  h5createGroup(output_path, "col_graphs")
  h5createGroup(output_path, "row_graphs")
  h5createDataset(file = output_path,
                  dataset = "matrix",
                  dims = dim(bc_gene_mat),
                  storage.mode = "double",
                  chunk=c(1, ncol(bc_gene_mat)))
  h5write(bc_gene_mat, file=output_path, name="matrix")
}
merge_2_matrices = function(mat1, mat2){
  cat("Input1: ", dim(mat1), "\n")
  cat("Input2: ", dim(mat2), "\n")
  rownames_merged = union(rownames(mat1), rownames(mat2))
  mat_merged = matrix(0, nrow = length(rownames_merged), ncol = ncol(mat1) + ncol(mat2), dimnames = list(rownames_merged, c(colnames(mat1), colnames(mat2))))
  mat_merged[rownames(mat1), 1:ncol(mat1)] = mat1
  mat_merged[rownames(mat2), (ncol(mat1) + 1): (ncol(mat1) + ncol(mat2))] = mat2
  cat("Output: ", dim(mat_merged), "\n")
  return(mat_merged)
}
merge_2_sparse_matrices = function(mat1, mat2){
  cat("Input1: ", dim(mat1), "\n")
  cat("Input2: ", dim(mat2), "\n")
  rownames_merged = union(rownames(mat1), rownames(mat2))
  tmp_mat1_0 = as(mat1, "dgCMatrix")
  tmp_mat1_0@Dim[1] = length(rownames_merged)
  rownames(tmp_mat1_0) = c(rownames(tmp_mat1_0), rownames_merged[!rownames_merged %in% rownames(tmp_mat1_0)])
  tmp_mat1 = tmp_mat1_0[rownames_merged, ]
  tmp_mat2_0 = as(mat2, "dgCMatrix")
  tmp_mat2_0@Dim[1] = length(rownames_merged)
  rownames(tmp_mat2_0) = c(rownames(tmp_mat2_0), rownames_merged[!rownames_merged %in% rownames(tmp_mat2_0)])
  tmp_mat2 = tmp_mat2_0[rownames_merged, ]
  mat_merged = Matrix::sparseMatrix(i = integer(0), j = integer(0), x = numeric(0), dims = c(length(rownames_merged), ncol(mat1) + ncol(mat2)), dimnames = list(rownames_merged, c(colnames(mat1), colnames(mat2))))
  mat_merged@i = c(tmp_mat1@i, tmp_mat2@i)
  mat_merged@x = c(tmp_mat1@x, tmp_mat2@x)
  mat_merged@p = c(tmp_mat1@p, tmp_mat2@p[-1] + max(tmp_mat1@p))
  cat("Output: ", dim(mat_merged), "\n")
  return(mat_merged)
}
merge_matrices = function(..., is_sparse = T){
  if(is_sparse){
    merge_function = merge_2_sparse_matrices
  }else{
    merge_function = merge_2_matrices
  }
  if(length(...) < 2){
    stop("Should input more than 2 matrices!")
  }else{
    output_mat = ...[[1]]
    for(ii in 2:length(...)){
      output_mat = merge_function(output_mat, ...[[ii]])
    }
    return(output_mat)
  }
}

load_data = function(loom_path, metadata_path, seurat_genename_adjustion = T){
  # Only cells in metadata can be read
  py$this_loom_path = loom_path
  py_run_string("tmp_adata = sc.read_loom(this_loom_path)")
  if(length(grep("_", py$tmp_adata$var_names$values)) != 0 & seurat_genename_adjustion){
    message("Adjustion for Seurat: Replacing underscores ('_') with dashes ('-')")
    py$tmp_adata$var_names = gsub("_", "-", py$tmp_adata$var_names$values, fixed = T)
  }
  metadata_fread_df = as.data.frame(fread(metadata_path, sep = "\t"))
  metadata = metadata_fread_df[, -1]
  rownames(metadata) = metadata_fread_df[, 1]
  py$this_cell = rownames(metadata)
  py_run_string("tmp_adata = tmp_adata[this_cell, :]")
  numeric_columns = c("n_genes", "n_genes_by_counts", "total_counts", "total_counts_mt", "pct_counts_mt",
                      "Ambiguity", "Exon", "Intergenic", "Intron", "Unmapped", "nReads", "rRNA", "pct_exon",
                      "nCount_RNA", "nFeature_RNA")
  for(ii in colnames(metadata)){
    this_metadata = metadata[, ii]
    this_metadata[is.na(this_metadata) | tolower(as.character(this_metadata)) %in% c("na", "nan")] = NA
    if(ii %in% numeric_columns){
      py$tmp_adata$obs[ii] = as.numeric(this_metadata)
    }else{
      py$tmp_adata$obs[ii] = as.factor(this_metadata)
    }
  }
  py_run_string("del this_loom_path")
  py_run_string("del this_cell")
  tmp_adata = py$tmp_adata
  py_run_string("del tmp_adata")
  return(tmp_adata)
}

seurat_hvg_selection = function(input_counts, metadata, hvg_number = 2000, min_cell_number = 10, all_hvgs_from_batches = F, verbose = F){
  seurat_obj = CreateSeuratObject(CreateAssayObject(input_counts))
  if(!is.null(metadata)){
    seurat_obj@meta.data = metadata
    per_sample_list = list(seurat_obj)
    for(ii in colnames(metadata)){
      per_sample_list = unlist(lapply(per_sample_list, function(x){
        return(SplitObject(x, split.by = ii))
      }), recursive = F)
    }
    batch_number = length(per_sample_list)
    removed_batch_number = 0
    removed_batch_min = Inf
    removed_batch_max = -Inf
    for(ii in names(per_sample_list)){
      if(ncol(per_sample_list[[ii]]) < min_cell_number){
        removed_batch_min = min(c(removed_batch_min, ncol(per_sample_list[[ii]])))
        removed_batch_max = max(c(removed_batch_max, ncol(per_sample_list[[ii]])))
        per_sample_list[[ii]] = NULL
        removed_batch_number = removed_batch_number + 1
      }
    }
    if(removed_batch_number > 0){
      message(paste0("Removed (", removed_batch_number, " / ", batch_number, ") batches of cells from ", removed_batch_min, " to ", removed_batch_max, "."))
      if(removed_batch_number == batch_number){
        return(NULL)
      }
    }
    if(all_hvgs_from_batches){
      per_sample_list = lapply(per_sample_list, function(x){
        x = NormalizeData(x, verbose = verbose)
        return(FindVariableFeatures(x, selection.method = "mvp", verbose = verbose))
      })
      return(unique(unlist(lapply(per_sample_list, VariableFeatures))))
    }else{
      per_sample_list = lapply(per_sample_list, function(x){
        return(FindVariableFeatures(x, selection.method = "vst", nfeatures = hvg_number, verbose = verbose))
      })
      return(SelectIntegrationFeatures(object.list = per_sample_list, nfeatures = hvg_number))
    }
  }else{
    return(VariableFeatures(FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = hvg_number, verbose = verbose)))
  }
}

PCA_calculation_python = function(adata, random_state = 0, merge_by = NULL, gene_used = NULL, hvg = NULL, hvg_number = 2000, min_cell_number = 10, all_hvgs_from_batches = F, pca_dim = 50, verbose = F, keep = F, return_all = F, factors_regress_out = NULL){
  return_list = list()
  py$tmp_adata = adata$copy()
  py_run_string("tmp_adata.layers['raw_counts'] = tmp_adata.X.copy()")
  py_run_string("sc.pp.normalize_total(tmp_adata, target_sum=1e4)")
  py_run_string("sc.pp.log1p(tmp_adata)")
  py_run_string("tmp_adata.layers['normalized'] = tmp_adata.X.copy()")
  if(!is.null(gene_used)){
    py$tmp_gene_used = gene_used
    py_run_string("tmp_adata = tmp_adata[:, tmp_gene_used].copy()")
    py_run_string("del tmp_gene_used")
  }
  if(hvg_number == "all"){
    py_run_string("tmp_adata.var['highly_variable'] = True")
    py_run_string("tmp_adata_hvg = tmp_adata.copy()")
  }else{
    if(!is.null(hvg)){
      py$tmp_hvg = hvg
    }else{
      if(all_hvgs_from_batches){
        input_counts = t(as(py$tmp_adata$layers['normalized'], "CsparseMatrix"))
      }else{
        input_counts = t(as(py$tmp_adata$layers['raw_counts'], "CsparseMatrix"))
      }
      dimnames(input_counts) = list(py$tmp_adata$var_names$values, py$tmp_adata$obs_names$values)
      if(!is.null(merge_by)){
        py$tmp_merge_by = merge_by
        py_run_string("tmp_selected_metadata = tmp_adata.obs[tmp_merge_by].astype(str).copy()")
        tmp_selected_metadata = py$tmp_selected_metadata
        py_run_string("del tmp_selected_metadata")
        selected_metadata = as.data.frame(apply(tmp_selected_metadata, length(dim(tmp_selected_metadata)), as.character))
        rownames(selected_metadata) = py$tmp_adata$obs_names$values
        colnames(selected_metadata) = merge_by
        py_run_string("del tmp_merge_by")
      }else{
        selected_metadata = NULL
      }
      py$tmp_hvg = seurat_hvg_selection(input_counts = input_counts, metadata = selected_metadata, hvg_number = hvg_number, min_cell_number = min_cell_number, all_hvgs_from_batches = all_hvgs_from_batches, verbose = verbose)
      rm(input_counts)
    }
    if(is.null(py$tmp_hvg) | py$tmp_adata$n_obs - 1 < pca_dim){
      py_run_string("del tmp_adata")
      py_run_string("del tmp_hvg")
      return(NULL)
    }
    py_run_string("tmp_adata.var['highly_variable'] = False")
    py_run_string("tmp_adata.var.loc[tmp_hvg, ['highly_variable']] = True")
    py_run_string("tmp_adata_hvg = tmp_adata[:, tmp_hvg]")
    py_run_string("del tmp_hvg")
  }
  if(keep){
    py_run_string("tmp_adata_hvg.obsm['raw_counts'] = tmp_adata.layers['raw_counts']")
  }
  if(return_all){
    return_list[["complete_matrix"]] = py$tmp_adata
  }
  py_run_string("del tmp_adata")
  if(!is.null(factors_regress_out)){
    py_run_string(paste0("sc.pp.regress_out(tmp_adata_hvg, ['", paste(factors_regress_out, collapse = "', '"), "'])"))
  }
  py_run_string("sc.pp.scale(tmp_adata_hvg)")
  py_run_string("tmp_adata_hvg.layers['scaled'] = tmp_adata_hvg.X.copy()")
  py_run_string(paste0("sc.pp.pca(tmp_adata_hvg, n_comps=", pca_dim, ", svd_solver='arpack', random_state=", random_state, ", zero_center=False, use_highly_variable=False)"))
  if(!is.null(merge_by)){
    merge_by_str = paste(merge_by, collapse = "', '")
    py_run_string(paste0("tmp_adata_hvg.obs[['", merge_by_str, "']] = tmp_adata_hvg.obs[['", merge_by_str, "']].astype(str).astype('object')"))
    py_run_string(paste0("sc.external.pp.harmony_integrate(tmp_adata_hvg, key=['", merge_by_str, "'], basis='X_pca', adjusted_basis='PCA_use', random_state=", random_state, ")"))
  }else{
    py_run_string("tmp_adata_hvg.obsm['PCA_use'] = tmp_adata_hvg.obsm['X_pca']")
  }
  return_list[["hvg_matrix"]] = py$tmp_adata_hvg
  py_run_string("del tmp_adata_hvg")
  if(return_all){
    return(return_list)
  }else{
    return(return_list[["hvg_matrix"]])
  }
}

pipeline_post_PCA = function(adata_str, clustering_embedding, visualization_embedding, n_neighbors_clustering, n_neighbors_visualization, leiden_resolutions, UMAP_min_dist = 0.5, clustering_prefix = "_leiden_", random_state = 0){
  # Visualization
  py_run_string(paste0("pca_dim = ", adata_str, ".obsm['", clustering_embedding, "'].shape[1]"))
  py_run_string(paste0("sc.pp.neighbors(", adata_str, ", use_rep='", clustering_embedding, "', n_neighbors=", n_neighbors_visualization, ", n_pcs=pca_dim, random_state=", random_state, ", key_added='visualization')"))
  py_run_string(paste0("sc.tl.umap(", adata_str, ", random_state=", random_state, ", neighbors_key='visualization', min_dist=", UMAP_min_dist, ")"))
  py_run_string(paste0(adata_str, ".obsm['", visualization_embedding, "'] = ", adata_str, ".obsm['X_umap'].copy()"))
  py_run_string(paste0("del ", adata_str, ".obsm['X_umap']"))
  # Clustering
  py_run_string(paste0("sc.pp.neighbors(", adata_str, ", use_rep='", clustering_embedding, "', n_neighbors=", n_neighbors_clustering, ", n_pcs=pca_dim, random_state=", random_state, ", key_added='clustering')"))
  clustering_keys = c()
  for(ii in leiden_resolutions){
    this_key = paste0(clustering_embedding, clustering_prefix, gsub(pattern = ".", replacement = "", sprintf("%.2f", ii), fixed = T))
    py_run_string(paste0("sc.tl.leiden(", adata_str, ", resolution=", ii, ", key_added='", this_key, "', neighbors_key='clustering')"))
    clustering_keys = c(clustering_keys, this_key)
  }
  return(clustering_keys)
}

matrix_clipping = function(adata){
  py$tmp_adata = adata$copy()
  message("Cells and genes with no expression, ERCC spikes and rRNA in are removed from the matrix.")
  py_run_string("tmp_adata = tmp_adata[:, ~tmp_adata.var_names.str.startswith('ERCC-')]")
  py_run_string("tmp_adata = tmp_adata[:, ~tmp_adata.var_names.str.startswith('rRNA-')]")
  py_run_string("sc.pp.filter_cells(tmp_adata, min_genes=1)")
  py_run_string("sc.pp.filter_genes(tmp_adata, min_cells=1)")
  adata = py$tmp_adata
  py_run_string("del tmp_adata")
  return(adata)
}

unify_genename_seurat = function(seurat_obj, alias_gene_c_seurat, verbose=F){
  message("Adjustion for Seurat: Replacing underscores ('_') with dashes ('-')")
  genename_original = rownames(seurat_obj)
  genename_original_seurat = gsub("_", "-", genename_original, fixed = T)
  genename = alias_gene_c_seurat[genename_original_seurat]
  genename[is.na(genename)] = genename_original_seurat[is.na(genename)]
  unique_index = which(!duplicated(genename))
  this_sp = as(seurat_obj[["RNA"]]$counts[unique_index, ], "RsparseMatrix")
  genename_duplicated = unique(genename[duplicated(genename)])
  genename_unique = genename[!duplicated(genename)]
  for(ii in genename_duplicated){
    this_index_old = which(genename == ii)
    this_index_new = which(genename_unique == ii)
    this_sp[this_index_new, ] = colSums(seurat_obj[["RNA"]]$counts[this_index_old, ])
    if(verbose){
      print(paste0(which(genename_duplicated == ii), " - ", length(genename_duplicated)))
    }
  }
  rownames(this_sp) = genename_unique
  seurat_obj_updated = CreateSeuratObject(as(this_sp, "CsparseMatrix"))
  seurat_obj_updated@meta.data = seurat_obj@meta.data
  return(seurat_obj_updated)
}

unify_genename = function(adata, alias_gene_c, verbose=F){
  py$tmp_adata = adata
  message("Adjustion for Seurat: Replacing underscores ('_') with dashes ('-')")
  py_run_string("tmp_var_names = tmp_adata.var_names.values.copy()")
  genename_original_seurat = gsub("_", "-", py$tmp_var_names, fixed = T)
  genename = alias_gene_c[genename_original_seurat]
  genename[is.na(genename)] = genename_original_seurat[is.na(genename)]
  py$tmp_unique_index = which(!duplicated(genename)) - 1
  py_run_string("this_sp = tmp_adata.X[:, tmp_unique_index].copy()")
  genename_duplicated = unique(genename[duplicated(genename)])
  genename_unique = genename[!duplicated(genename)]
  for(ii in genename_duplicated){
    py$tmp_this_index_old = which(genename == ii) - 1
    py$tmp_this_index_new = which(genename_unique == ii) - 1
    py_run_string("this_sp[:, np.array(tmp_this_index_new, int)] = tmp_adata.X[:, np.array(tmp_this_index_old, int)].sum(1)")
    if(verbose){
      print(paste0(which(genename_duplicated == ii), " - ", length(genename_duplicated)))
    }
  }
  py$tmp_var_names = genename_unique
  py_run_string("tmp_adata = tmp_adata[:, tmp_unique_index]")
  py_run_string("tmp_adata.var_names = tmp_var_names")
  py_run_string("tmp_adata.X = this_sp")
  py_run_string("del tmp_unique_index")
  py_run_string("del tmp_var_names")
  py_run_string("del this_sp")
  adata = py$tmp_adata
  return(adata)
}

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

load_data_wrap = function(loom_path, metadata_path, max_number = NA, sample_method = "geosketch", random_state = 0, merge_by = NULL, hvg_number = 2000, min_cell_number = 10, pca_dim = 50, keep = F, seurat_genename_adjustion = T, alias_unifying=NULL, use_embedding = NULL){
  adata = load_data(loom_path, metadata_path, seurat_genename_adjustion)
  adata = matrix_clipping(adata)
  if(!is.null(alias_unifying)){
    adata = unify_genename(adata, alias_unifying)
  }
  return(adata)
}

DEG_Analysis = function(adata, clustering, output_dir, layer = "normalized", prefix = "", top_DE_range = 50, calculate_every_cluster = F){
  py$tmp_adata = py$ad$AnnData(X=adata$layers[[layer]], obs=adata$obs, var=adata$var)$copy()
  return_list = list()
  return_list[["rest"]] = list()
  py_run_string(paste0("sc.tl.rank_genes_groups(tmp_adata, '", clustering, "', method='wilcoxon', layer=None, use_raw=False)"))
  py_run_string("tmp_DEG_result = tmp_adata.uns['rank_genes_groups']")
  py_run_string("tmp_ClusterID = tmp_DEG_result['names'].dtype.names")
  py_run_string("tmp_DEG_table_0 = np.array([[tmp_DEG_result[x][y] for y in tmp_ClusterID] for x in ['names', 'pvals_adj', 'logfoldchanges', 'scores']])")
  py_run_string("tmp_DEG_table_0_add = np.repeat([[tmp_ClusterID]], tmp_DEG_table_0.shape[2], axis=0).transpose([1, 2, 0])")
  py_run_string("tmp_DEG_table_1 = np.vstack([tmp_DEG_table_0_add, tmp_DEG_table_0])")
  py_run_string("tmp_DEG_table = pd.DataFrame(tmp_DEG_table_1.reshape((tmp_DEG_table_1.shape[0], tmp_DEG_table_1.shape[1] * tmp_DEG_table_1.shape[2])), index=['Cluster', 'Gene', 'P_adj', 'Log2FC', 'Scores']).T")
  py_run_string("tmp_DEG_table_filtered = tmp_DEG_table[tmp_DEG_table['P_adj'] <= 0.95]") # Only filter out genes that probably make no sense (even in another side of 0.05) for a more clear vision
  py_run_string(paste0("tmp_DEG_top_p_by_FC = tmp_DEG_table_filtered[tmp_DEG_table_filtered['Log2FC'] > 0].groupby('Cluster', sort=False).apply(lambda x: x.sort_values(['Log2FC'], ascending=False).head(", top_DE_range, ")).reset_index(drop=True)"))
  py_run_string(paste0("tmp_DEG_top_n_by_FC = tmp_DEG_table_filtered[tmp_DEG_table_filtered['Log2FC'] < 0].groupby('Cluster', sort=False).apply(lambda x: x.sort_values(['Log2FC'], ascending=True).head(", top_DE_range, ")).reset_index(drop=True)"))
  py_run_string(paste0("tmp_DEG_top_p_by_score = tmp_DEG_table_filtered[tmp_DEG_table_filtered['Log2FC'] > 0].groupby('Cluster', sort=False).apply(lambda x: x.sort_values(['Scores'], ascending=False).head(", top_DE_range, ")).reset_index(drop=True)"))
  py_run_string(paste0("tmp_DEG_top_n_by_score = tmp_DEG_table_filtered[tmp_DEG_table_filtered['Log2FC'] < 0].groupby('Cluster', sort=False).apply(lambda x: x.sort_values(['Scores'], ascending=True).head(", top_DE_range, ")).reset_index(drop=True)"))
  py_run_string(paste0("tmp_DEG_top_p_by_P = tmp_DEG_table_filtered[tmp_DEG_table_filtered['Log2FC'] > 0].groupby('Cluster', sort=False).apply(lambda x: x.sort_values(['P_adj'], ascending=True).head(", top_DE_range, ")).reset_index(drop=True)"))
  py_run_string(paste0("tmp_DEG_top_n_by_P = tmp_DEG_table_filtered[tmp_DEG_table_filtered['Log2FC'] < 0].groupby('Cluster', sort=False).apply(lambda x: x.sort_values(['P_adj'], ascending=True).head(", top_DE_range, ")).reset_index(drop=True)"))
  DEG_table_output = matrix(nrow = 0, ncol = 5)
  DEG_top_p_by_FC_output = matrix(nrow = 0, ncol = 5)
  DEG_top_n_by_FC_output = matrix(nrow = 0, ncol = 5)
  DEG_top_p_by_score_output = matrix(nrow = 0, ncol = 5)
  DEG_top_n_by_score_output = matrix(nrow = 0, ncol = 5)
  DEG_top_p_by_P_output = matrix(nrow = 0, ncol = 5)
  DEG_top_n_by_P_output = matrix(nrow = 0, ncol = 5)
  for(ii in unlist(py$tmp_ClusterID)){
    DEG_table_output = rbind(DEG_table_output, colnames(py$tmp_DEG_table))
    DEG_table_output = rbind(DEG_table_output, as.matrix(py$tmp_DEG_table[py$tmp_DEG_table$Cluster == ii, ]))
    DEG_top_p_by_FC_output = rbind(DEG_top_p_by_FC_output, colnames(py$tmp_DEG_top_p_by_FC))
    DEG_top_p_by_FC_output = rbind(DEG_top_p_by_FC_output, as.matrix(py$tmp_DEG_top_p_by_FC[py$tmp_DEG_top_p_by_FC$Cluster == ii, ]))
    DEG_top_n_by_FC_output = rbind(DEG_top_n_by_FC_output, colnames(py$tmp_DEG_top_n_by_FC))
    DEG_top_n_by_FC_output = rbind(DEG_top_n_by_FC_output, as.matrix(py$tmp_DEG_top_n_by_FC[py$tmp_DEG_top_n_by_FC$Cluster == ii, ]))
    DEG_top_p_by_score_output = rbind(DEG_top_p_by_score_output, colnames(py$tmp_DEG_top_p_by_score))
    DEG_top_p_by_score_output = rbind(DEG_top_p_by_score_output, as.matrix(py$tmp_DEG_top_p_by_score[py$tmp_DEG_top_p_by_score$Cluster == ii, ]))
    DEG_top_n_by_score_output = rbind(DEG_top_n_by_score_output, colnames(py$tmp_DEG_top_n_by_score))
    DEG_top_n_by_score_output = rbind(DEG_top_n_by_score_output, as.matrix(py$tmp_DEG_top_n_by_score[py$tmp_DEG_top_n_by_score$Cluster == ii, ]))
    DEG_top_p_by_P_output = rbind(DEG_top_p_by_P_output, colnames(py$tmp_DEG_top_p_by_P))
    DEG_top_p_by_P_output = rbind(DEG_top_p_by_P_output, as.matrix(py$tmp_DEG_top_p_by_P[py$tmp_DEG_top_p_by_P$Cluster == ii, ]))
    DEG_top_n_by_P_output = rbind(DEG_top_n_by_P_output, colnames(py$tmp_DEG_top_n_by_P))
    DEG_top_n_by_P_output = rbind(DEG_top_n_by_P_output, as.matrix(py$tmp_DEG_top_n_by_P[py$tmp_DEG_top_n_by_P$Cluster == ii, ]))
  }
  write.table(DEG_table_output, paste0(output_dir, "/", prefix, "DEG_table.tsv"), row.names = F, col.names = F, quote = F, sep = "\t")
  write.table(DEG_top_p_by_FC_output, paste0(output_dir, "/", prefix, "DEG_P_by_FC_top_", top_DE_range, ".tsv"), row.names = F, col.names = F, quote = F, sep = "\t")
  write.table(DEG_top_n_by_FC_output, paste0(output_dir, "/", prefix, "DEG_N_by_FC_top_", top_DE_range, ".tsv"), row.names = F, col.names = F, quote = F, sep = "\t")
  write.table(DEG_top_p_by_score_output, paste0(output_dir, "/", prefix, "DEG_P_by_score_top_", top_DE_range, ".tsv"), row.names = F, col.names = F, quote = F, sep = "\t")
  write.table(DEG_top_n_by_score_output, paste0(output_dir, "/", prefix, "DEG_N_by_score_top_", top_DE_range, ".tsv"), row.names = F, col.names = F, quote = F, sep = "\t")
  write.table(DEG_top_p_by_P_output, paste0(output_dir, "/", prefix, "DEG_P_by_P_top_", top_DE_range, ".tsv"), row.names = F, col.names = F, quote = F, sep = "\t")
  write.table(DEG_top_n_by_P_output, paste0(output_dir, "/", prefix, "DEG_N_by_P_top_", top_DE_range, ".tsv"), row.names = F, col.names = F, quote = F, sep = "\t")
  py_run_string("del tmp_DEG_table_0")
  py_run_string("del tmp_DEG_table_0_add")
  py_run_string("del tmp_DEG_table_1")
  return_list[["rest"]][["table"]] = py$tmp_DEG_table
  return_list[["rest"]][["top_p_by_FC"]] = py$tmp_DEG_top_p_by_FC
  return_list[["rest"]][["top_n_by_FC"]] = py$tmp_DEG_top_n_by_FC
  return_list[["rest"]][["top_p_by_score"]] = py$tmp_DEG_top_p_by_score
  return_list[["rest"]][["top_n_by_score"]] = py$tmp_DEG_top_n_by_score
  return_list[["rest"]][["top_p_by_P"]] = py$tmp_DEG_top_p_by_P
  return_list[["rest"]][["top_n_by_P"]] = py$tmp_DEG_top_n_by_P
  py_run_string("del tmp_DEG_result")
  py_run_string("del tmp_DEG_table")
  py_run_string("del tmp_DEG_table_filtered")
  py_run_string("del tmp_DEG_top_p_by_FC")
  py_run_string("del tmp_DEG_top_n_by_FC")
  py_run_string("del tmp_DEG_top_p_by_score")
  py_run_string("del tmp_DEG_top_n_by_score")
  py_run_string("del tmp_DEG_top_p_by_P")
  py_run_string("del tmp_DEG_top_n_by_P")
  if(calculate_every_cluster){
    return_list[["every_cluster"]] = list()
    py_run_string("tmp_DEG_table_list = []")
    for(ii in py$tmp_ClusterID){
      ii_name = gsub(" ", "", ii)
      py_run_string(paste0("sc.tl.rank_genes_groups(tmp_adata, '", clustering, "', reference='", ii, "', method='wilcoxon', layer=None, use_raw=False)"))
      py_run_string(paste0("tmp_DEG_result_", ii_name, " = tmp_adata.uns['rank_genes_groups']"))
      py_run_string(paste0("tmp_ClusterID_", ii_name, " = tmp_DEG_result_", ii_name, "['names'].dtype.names"))
      py_run_string(paste0("tmp_DEG_table_0_", ii_name, " = np.array([pd.DataFrame(np.array([tmp_DEG_result_", ii_name, "[x][y] for x in ['names', 'pvals_adj', 'logfoldchanges', 'scores']]).transpose(), index=tmp_DEG_result_", ii_name, "['names'][y]).reindex(tmp_adata.var_names.values) for y in tmp_ClusterID_", ii_name, "]).transpose([2, 0, 1])"))
      py_run_string(paste0("tmp_DEG_table_0_add_0_", ii_name, " = np.repeat([[tmp_ClusterID_", ii_name, "]], tmp_DEG_table_0_", ii_name, ".shape[2], axis=0).transpose([1, 2, 0])"))
      py_run_string(paste0("tmp_DEG_table_0_add_1_", ii_name, " = np.repeat([[np.repeat('", ii, "', len(tmp_ClusterID_", ii_name, "))]], tmp_DEG_table_0_", ii_name, ".shape[2], axis=0).transpose([1, 2, 0])"))
      py_run_string(paste0("tmp_DEG_table_1_", ii_name, " = np.vstack([tmp_DEG_table_0_add_0_", ii_name, ", tmp_DEG_table_0_add_1_", ii_name, ", tmp_DEG_table_0_", ii_name, "])"))
      py_run_string(paste0("tmp_DEG_table_", ii_name, " = pd.DataFrame(tmp_DEG_table_1_", ii_name, ".reshape((tmp_DEG_table_1_", ii_name, ".shape[0], tmp_DEG_table_1_", ii_name, ".shape[1] * tmp_DEG_table_1_", ii_name, ".shape[2])), index=['Cluster', 'Reference', 'Gene', 'P_adj', 'Log2FC', 'Scores']).T"))
      py_run_string(paste0("tmp_DEG_table_list.append(tmp_DEG_table_", ii_name, ")"))
      py_run_string(paste0("del tmp_DEG_result_", ii_name))
      py_run_string(paste0("del tmp_ClusterID_", ii_name))
      py_run_string(paste0("del tmp_DEG_table_0_", ii_name))
      py_run_string(paste0("del tmp_DEG_table_0_add_0_", ii_name))
      py_run_string(paste0("del tmp_DEG_table_0_add_1_", ii_name))
      py_run_string(paste0("del tmp_DEG_table_1_", ii_name))
      py_run_string(paste0("del tmp_DEG_table_", ii_name))
    }
    py_run_string("tmp_DEG_table_every_cluster = pd.concat(tmp_DEG_table_list, 0)")
    py_run_string("del tmp_DEG_table_list")
    return_list[["every_cluster"]][["table"]] = py$tmp_DEG_table_every_cluster
    py_run_string("del tmp_DEG_table_every_cluster")
  }
  py_run_string("del tmp_ClusterID")
  py_run_string("del tmp_adata")
  py_run_string("gc.collect()")
  return(return_list)
}

density_visualization = function(str_adata, str_basis, str_groupby, str_save, str_color_map, str_ncols="4", str_dotsize="8", show=FALSE){
  # Somehow tl.embedding_density can cause LinAlgError: singular matrix
  # Somehow tl.embedding_density can cause error if there's only 1 or 2 cells in a category
  if(show){
    str_show = "True"
  }else{
    str_show = "False"
  }
  tryCatch(expr = {
    py_run_string(paste0("tmp_elements, tmp_counts = np.unique(", str_adata, ".obs['", str_groupby, "'].astype(str).values, return_counts=True)"))
    py_run_string("element_cannot_calculate = tmp_elements[tmp_counts <= 2]")
    print("Only have one cells:")
    print(py$element_cannot_calculate)
    print("Removed from the figure")
    py_run_string(paste0("tmp_adata = ", str_adata, "[np.logical_not(np.isin(", str_adata, ".obs['", str_groupby, "'], element_cannot_calculate))]"))
    py_run_string(paste0("basis_lower = '", str_basis, "'.lower()"))
    py_run_string("x_basis_lower = f'X_{basis_lower}'")
    py_run_string(paste0("tmp_adata.obsm[x_basis_lower] = tmp_adata.obsm['", str_basis, "'] if x_basis_lower not in tmp_adata.obsm_keys() else tmp_adata.obsm[x_basis_lower]"))
    py_run_string(paste0("sc.tl.embedding_density(tmp_adata, basis=basis_lower, groupby='", str_groupby, "')"))
    py_run_string(paste0("sc.pl.embedding_density(tmp_adata, basis=basis_lower, key=f'{basis_lower}_density_", str_groupby, "', save='", str_save, "', fg_dotsize=", str_dotsize, ", bg_dotsize=", str_dotsize, ", frameon=False, ncols=", str_ncols, ", color_map=", str_color_map, ", show=", str_show, ")"))
    py_run_string("del tmp_elements")
    py_run_string("del tmp_counts")
    py_run_string("del tmp_adata")
  }, error = function(e){
    message(e)
  })
}

modify_vlnplot = function(input_df, x, y, fill) {
  p_0 = eval(parse(text = paste0("ggplot(input_df, aes(x = ", x, ", y = ", y, ", fill = ", fill, "))")))
  p = p_0 + geom_violin() + geom_jitter(size = 0.05) + xlab("") + ylab(y) + theme_classic() + 
    theme(legend.position = "none", 
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(), 
          axis.title.y = element_text(size = rel(1), angle = 0), 
          axis.text.y = element_text(size = rel(1)),
          plot.margin = unit(c(-0.75, 0, -0.75, 0), "cm"))
  return(p)
}

## extract the max value of the y axis
extract_max = function(p){
  ymax = max(ggplot_build(p)$layout$panel_scales_y[[1]]$range$range)
  return(ceiling(ymax))
}


## main function
StackedVlnPlot = function(input_df, x, ys, fill, ymaxs = NULL) {
  plot_list = purrr::map(ys, function(y) modify_vlnplot(input_df = input_df, x = x, y = y, fill = fill))
  # Add back x-axis title to bottom plot. patchwork is going to support this?
  plot_list[[length(plot_list)]] = plot_list[[length(plot_list)]] + theme(axis.text.x=element_text(), axis.ticks.x = element_line())
  # change the y-axis tick to only max value 
  if(is.null(ymaxs)){
    ymaxs = purrr::map_dbl(plot_list, extract_max)
  }
  plot_list = purrr::map2(plot_list, ymaxs, function(x,y) x + ylim(0, y))
  p = patchwork::wrap_plots(plotlist = plot_list, ncol = 1)
  return(p)
}

get_dot_size = function(cell_number){
  if(cell_number <= 200){
    return(100)
  }else if(cell_number > 200 & cell_number <= 300){
    return(90)
  }else if(cell_number > 300 & cell_number <= 400){
    return(85)
  }else if(cell_number > 400 & cell_number <= 500){
    return(80)
  }else if(cell_number > 500 & cell_number <= 600){
    return(78)
  }else if(cell_number > 600 & cell_number <= 700){
    return(75)
  }else if(cell_number > 700 & cell_number <= 800){
    return(72)
  }else if(cell_number > 800 & cell_number <= 900){
    return(70)
  }else if(cell_number > 900 & cell_number <= 1000){
    return(68)
  }else if(cell_number > 1000 & cell_number <= 1100){
    return(65)
  }else if(cell_number > 1100 & cell_number <= 1200){
    return(62)
  }else if(cell_number > 1200 & cell_number <= 1300){
    return(60)
  }else if(cell_number > 1300 & cell_number <= 1400){
    return(58)
  }else if(cell_number > 1400 & cell_number <= 1500){
    return(55)
  }else if(cell_number > 1500 & cell_number <= 1600){
    return(54)
  }else if(cell_number > 1600 & cell_number <= 1700){
    return(53)
  }else if(cell_number > 1700 & cell_number <= 1800){
    return(52)
  }else if(cell_number > 1800 & cell_number <= 1900){
    return(51)
  }else if(cell_number > 1900 & cell_number <= 2000){
    return(50)
  }else if(cell_number > 2000 & cell_number <= 200000){
    return(100000/cell_number)
  }else if(cell_number > 200000){
    return(0.5)
  }else{
    stop("cell_number is out of range!")
  }
}
#
featureMapMerge = function(input_adata_str, output_adata_str, mapping_c, RC_layer=NULL){
  # Map and merge (e.g. gene id to gene name)
  if(input_adata_str == output_adata_str){
    stop("input_adata_str and output_adata_str should be different!")
  }
  GeneID = eval(parse(text = paste0("py$", input_adata_str, "$var_names$values")))
  GeneID_to_GeneName = mapping_c[GeneID]
  GeneID_to_GeneName[is.na(GeneID_to_GeneName)] = GeneID[is.na(GeneID_to_GeneName)]
  #
  mask_unique_GeneName = !duplicated(GeneID_to_GeneName)
  GeneName_unique = GeneID_to_GeneName[mask_unique_GeneName]
  names(GeneName_unique) = GeneID[mask_unique_GeneName]
  duplicated_GeneName = GeneID_to_GeneName[duplicated(GeneID_to_GeneName)]
  unique_DGeneName = duplicated_GeneName[!duplicated(duplicated_GeneName)]
  py$mask_DGeneName_on_UGeneName = GeneName_unique %in% unique_DGeneName
  #
  py$mask_unique_GeneName = mask_unique_GeneName
  py_run_string(paste0(output_adata_str, " = ", input_adata_str, "[:, mask_unique_GeneName].copy()"))
  if(is.null(RC_layer)){
    message("RC (raw count) layer is not specified, will copy .X to .layers['raw_counts']")
    py_run_string(paste0(output_adata_str, ".layers['raw_counts'] = ", output_adata_str, ".X.copy()"))
    RCL_running = "raw_counts"
    this_running_code = ".X.sum(1)"
  }else{
    RCL_running = RC_layer
    this_running_code = paste0(".layers['", RCL_running, "'].sum(1)")
  }
  if(length(unique_DGeneName) > 0){
    RC_unique_DGeneName = sapply(unique_DGeneName, function(x){
      py_run_string(paste0("tmp_RC = ", input_adata_str, "[:, ", input_adata_str, ".var_names.isin(['", x, "'])]", this_running_code))
      return(as.numeric(py$tmp_RC))
    })
    py_run_string("del tmp_RC")
    py$RC_unique_DGeneName = RC_unique_DGeneName
    py_run_string(paste0(output_adata_str, ".layers['", RCL_running, "'][:, mask_DGeneName_on_UGeneName] = RC_unique_DGeneName.copy()"))    
    py_run_string("del RC_unique_DGeneName")
  }
  if(is.null(RC_layer)){
    py_run_string(paste0(output_adata_str, ".X = ", output_adata_str, ".layers['raw_counts'].copy()"))
    RCL_running = "raw_counts"
  }
  py_run_string("gc.collect()")
  return(GeneName_unique)
}
#
mapping_GeneID_GeneName = function(input_adata_str, output_adata_str, mapping_GeneID_GeneName_c, RC_layer, feature_reference="NCBITaxon:9606"){
  # Map and merge GeneID to GeneName
  GeneName_unique = featureMapMerge(input_adata_str, output_adata_str, mapping_GeneID_GeneName_c, RC_layer=RC_layer)
  py$var_data_frame = data.frame(feature_is_filtered=F, feature_biotype="gene", feature_reference=feature_reference, feature_name=GeneName_unique)
  row.names(py$var_data_frame) = GeneName_unique
  py_run_string(paste0(output_adata_str, ".var = var_data_frame.copy()"))
  py_run_string(paste0(output_adata_str, ".var.loc[", output_adata_str, ".var_names.str.startswith('ERCC-'), 'feature_biotype'] = 'spike-in'"))
  py_run_string("del var_data_frame")
  py_run_string("gc.collect()")
}
#
mapping_GeneName_GeneID = function(GeneName, mapping_GeneName_GeneID_c){
  GeneName_to_GeneID = mapping_GeneName_GeneID_c[GeneName]
  GeneName_to_GeneID[is.na(GeneName_to_GeneID)] = GeneName[is.na(GeneName_to_GeneID)]
  return(GeneName_to_GeneID)
}