py_run_string("adata_dict = {}")
for(ii in seq(24)){
  this_dataset_index = sprintf("D%.3d", ii)
  print(this_dataset_index)
  py_run_string(paste0("tmp_adata_raw = ad.read_h5ad('", raw_input_paths[[this_dataset_index]], "')"))
  py_run_string(paste0("tmp_adata_raw = tmp_adata_raw[tmp_adata_raw.obs['individual_barcode'] != 'unknown']"))
  py_run_string(paste0("tmp_adata_individual = ad.read_h5ad('", individual_annotation_result_paths[[this_dataset_index]], "')"))
  py_run_string(paste0("tmp_adata_individual_hvg = ad.read_h5ad('", individual_hvg_clustering_paths[[this_dataset_index]], "')"))
  ###
  ###
  ###
  if(this_dataset_index == "D013"){
    py_run_string(paste0("tmp_adata_raw.obs['individual_barcode'] = tmp_adata_raw.obs_names.to_series().astype(str).str.split('___', expand=True)[1].values"))
    py_run_string(paste0("print(tmp_adata_raw.obs['individual_barcode'])"))
    py_run_string(paste0("print(tmp_adata_individual.obs_names)"))
  }
  if(this_dataset_index == "D015"){
    py_run_string(paste0("tmp_adata_raw.obs['individual_barcode'] = tmp_adata_raw.obs_names.to_series().astype(str).str.split('___', expand=True)[1].str.split('_', n=1, expand=True)[1].values"))
    py_run_string(paste0("print(tmp_adata_raw.obs['individual_barcode'])"))
    py_run_string(paste0("print(tmp_adata_individual.obs_names)"))
  }
  ###
  ###
  ###
  for(jj in leiden_resolutions){
    py_run_string(paste0("tmp_adata_raw.obs['", "individual_leiden_L1_", gsub(pattern = ".", replacement = "", sprintf("%.2f", jj), fixed = T), "'] = tmp_adata_individual[tmp_adata_raw.obs['individual_barcode'].values.tolist(), :].obs['leiden_1_", gsub(pattern = ".", replacement = "", sprintf("%.2f", jj), fixed = T), "'].values"))
  }
  py_run_string(paste0("tmp_adata_raw.obs['individual_annotation_L1'] = tmp_adata_individual[tmp_adata_raw.obs['individual_barcode'].values.tolist(), :].obs['SeparateAnnotation_1'].values"))
  ###
  ###
  ###
  py_run_string(paste0("tmp_adata_raw.obsm['individual_PCA_L1'] = tmp_adata_individual_hvg[tmp_adata_raw.obs['individual_barcode'].values.tolist(), :].obsm['PCA_use']"))
  py_run_string(paste0("tmp_adata_raw.obsm['individual_UMAP_L1'] = tmp_adata_individual_hvg[tmp_adata_raw.obs['individual_barcode'].values.tolist(), :].obsm['X_umap']"))
  py$tmp_adata_raw = sampling(adata = py$tmp_adata_raw, max_number = max_number, random_state = random_state, sample_method = "geosketch", keep = T, use_embedding = "individual_PCA_L1")
  py_run_string("print(tmp_adata_raw.X)")
  py_run_string(paste0("adata_dict['", this_dataset_index, "'] = tmp_adata_raw"))
}
py_run_string("adata = ad.concat(adata_dict, axis=0, join='outer')")
py_run_string("adata.obs['individual_annotation_L2'] = adata.obs['individual_annotation_L2'].astype(str)")
py_run_string("adata.obs.loc[adata.obs['individual_annotation_L2'] == 'nan', 'individual_annotation_L2'] = 'unknown'")
py_run_string(paste0("adata_obs_to_string(adata).write('", paste0(py$output_dir, "/adata_merged_individual_analyses.h5ad"), "')"))
#
py_run_string("print(adata.X)")
py_run_string("adata.layers['raw_counts'] = adata.X.copy()")
py_run_string("sc.pp.normalize_total(adata, target_sum=1e4)")
py_run_string("sc.pp.log1p(adata)")
py_run_string("adata.layers['normalized'] = adata.X.copy()")
py_run_string("adata.X = adata.layers['raw_counts'].copy()")
#
GeneID = py$adata$var_names$values
py$var_data_frame = data.frame(row.names = GeneID, feature_is_filtered=F, feature_biotype="gene", feature_reference="NCBITaxon:9606", feature_name=mapping_GeneName_GeneID(GeneID, hgnc_mapping_ensemblID_genename))
py_run_string(paste0("var_data_frame.loc[var_data_frame.index.str.startswith('ERCC-'), 'feature_biotype'] = 'spike-in'"))
py_run_string(paste0("adata.var = var_data_frame.copy()"))
py_run_string("del var_data_frame")
py_run_string("gc.collect()")
#
py_run_string(paste0("adata_obs_to_string(adata).write('", paste0(py$output_dir, "/adata_merged_I_20250624.h5ad"), "')"))
#
