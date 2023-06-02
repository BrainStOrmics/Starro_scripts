SpatialRNA1 <- function( counts, nUMI = NULL,  require_int = TRUE) {
  #counts <- check_counts(counts, 'SpatialRNA', require_int = require_int)
  if(is.null(nUMI)) {
    nUMI = colSums(counts)
    names(nUMI) <- colnames(counts)
  } else {
    check_UMI(nUMI, 'SpatialRNA', require_int = require_int)
  }
  barcodes <- intersect(names(nUMI), colnames(counts))
  if(length(barcodes) == 0)
    stop('SpatialRNA: coords, counts, and nUMI do not share any barcode names. Please ensure that rownames(coords) matches colnames(counts) and names(nUMI)')
  if(length(barcodes) < max(length(nUMI),dim(counts)[2]))
    warning('SpatialRNA: some barcodes in nUMI, coords, or counts were not mutually shared. Such barcodes were removed.')
  if(sum(nUMI[barcodes] != colSums(counts[,barcodes])) > 0)
    warning('SpatialRNA: nUMI does not match colSums of counts. If this is unintended, please correct this discrepancy. If this is intended, there is no problem.')
  new("SpatialRNA",  counts = counts[,barcodes], nUMI = nUMI[barcodes])
}

Reference <- function(counts, cell_types, nUMI = NULL, require_int = TRUE, n_max_cells = 10000, min_UMI = 50) {
  #counts <- check_counts(counts, 'Reference', require_2d = T, require_int = require_int)
  if(is.null(nUMI)) {
    nUMI = colSums(counts)
    names(nUMI) <- colnames(counts)
  } else {
    check_UMI(nUMI, 'Reference', require_2d = T, require_int = require_int, min_UMI = min_UMI)
  }
  check_cell_types(cell_types)
  barcodes <- intersect(intersect(names(nUMI), names(cell_types)), colnames(counts))
  if(length(barcodes) == 0)
    stop('Reference: cell_types, counts, and nUMI do not share any barcode names. Please ensure that names(cell_types) matches colnames(counts) and names(nUMI)')
  if(length(barcodes) < max(length(nUMI),length(cell_types),dim(counts)[2]))
    warning('Reference: some barcodes in nUMI, cell_types, or counts were not mutually shared. Such barcodes were removed.')
  barcodes <- names(which(nUMI[barcodes] >= min_UMI))
  if(length(barcodes) < 1)
    stop('Reference: no barcodes were included with nUMI at least min_UMI. Please lower the parameter min_UMI or ensure that cells have sufficient UMI counts.')
  if(sum(nUMI[barcodes] != colSums(counts[,barcodes])) > 0)
    warning('Reference: nUMI does not match colSums of counts. If this is unintended, please correct this discrepancy. If this is intended, there is no problem.')
  missing_cell_types <- names(which(table(cell_types[barcodes]) == 0))
  if(length(missing_cell_types) > 0)
    warning(paste('Reference: missing cell types with no occurences: ',paste(missing_cell_types,collapse=', ')))
  reference <- new("Reference", cell_types = cell_types[barcodes], counts = counts[,barcodes], nUMI = nUMI[barcodes])
  cur_count <- max(table(reference@cell_types))
  if(cur_count > n_max_cells) {
    warning(paste0('Reference: number of cells per cell type is ', cur_count, ', larger than maximum allowable of ', n_max_cells,
                   '. Downsampling number of cells to: ', n_max_cells))
    reference <- create_downsampled_data(reference, n_samples = n_max_cells)
  }
  return(reference)
}

makeRCTD=function(sc_obj,st_count){
  celltype=sc_obj$annotation##修改
  celltype=gsub("/","_",celltype)
  celltype=as.factor(celltype)
  names(celltype)=rownames(sc_obj@meta.data)
  counts=as(sc_obj@assays$RNA@counts,"dgCMatrix")
  reference=Reference(counts, celltype, nUMI = NULL, require_int = TRUE, n_max_cells = 10000, min_UMI = 50) 
  #st_count=st_count*10
  #st_count <- as(st_count,"matrix")
  st_count <- as(st_count,"dgCMatrix")
  puck= SpatialRNA1(st_count)
  myRCTD <- create.RCTD(puck, reference, max_cores =32, test_mode = FALSE,CELL_MIN_INSTANCE = 10)
  myRCTD <- run.RCTD(myRCTD, doublet_mode = 'doublet')
  return(myRCTD)
}

make_exp=function(sc_obj){
  sc_meta<-sc_obj@meta.data
  sc_data<-as.matrix(sc_obj@assays[["RNA"]]@counts)
  mouse_sc <- createGiottoObject(raw_exprs = sc_data)
  mouse_sc <- normalizeGiotto(gobject = mouse_sc, scalefactor = 6000, verbose = T)
  mouse_sc <- calculateHVG(gobject = mouse_sc)
  gene_metadata = fDataDT(mouse_sc)
  featgenes = gene_metadata[hvg == 'yes']$gene_ID
  mouse_sc <- runPCA(gobject = mouse_sc, genes_to_use = featgenes, scale_unit = F)
  signPCA(mouse_sc, genes_to_use = featgenes, scale_unit = F)
  #######calculate Sig for deconvolution, This step use DEG function implemented in Giotto
  mouse_sc@cell_metadata$leiden_clus <- as.character(sc_meta$annotation)##修改
  scran_markers_subclusters = findMarkers_one_vs_all(gobject = mouse_sc,
                                                     method = 'scran',
                                                     expression_values = 'normalized',
                                                     cluster_column = 'leiden_clus')
  Sig_scran <- unique(scran_markers_subclusters$genes[which(scran_markers_subclusters$ranking <= 200)])
  ########Calculate median expression value of signature genes in each cell type
  norm_exp<-2^(mouse_sc@norm_expr)-1
  id<-mouse_sc@cell_metadata$leiden_clus
  ExprSubset<-norm_exp[Sig_scran,]
  Sig_exp<-NULL
  for (i in unique(id)){
    Sig_exp<-cbind(Sig_exp,(apply(ExprSubset,1,function(y) mean(y[which(id==i)]))))
  }
  colnames(Sig_exp)<-unique(id)
  return(Sig_exp)
}
#intersect(rownames(st_count),rownames(Sig_exp))
makeSpatialDWLS=function(sc_obj,st_count){
  Sig_exp=make_exp(sc_obj)
  spatial_exp<-as.matrix(st_count)
  ST_data <- createGiottoObject(raw_exprs = spatial_exp#spatial_locs = spatial_loc,
  )
  ST_data <- normalizeGiotto(gobject = ST_data)
  ST_data <- calculateHVG(gobject = ST_data,logbase = 0.5, 
                          zscore_threshold = 1,  difference_in_cov = 0.01)
  gene_metadata = fDataDT(ST_data)
  featgenes = gene_metadata[hvg == 'yes']$gene_ID
  ST_data <- runPCA(gobject = ST_data, genes_to_use = featgenes, scale_unit = F)
  signPCA(ST_data, genes_to_use = featgenes, scale_unit = F)
  ST_data <- createNearestNetwork(gobject = ST_data, dimensions_to_use = 1:10, k = 10)
  ST_data <- doLeidenCluster(gobject = ST_data, resolution = 0.4, n_iterations = 100)
  ST_data <- runDWLSDeconv(gobject = ST_data, sign_matrix = Sig_exp)
  plot_data <- as.data.frame(ST_data@spatial_enrichment$DWLS)
  rownames(plot_data)=plot_data$cell_ID
  plot_data=plot_data[,-1]
  return(plot_data)
}