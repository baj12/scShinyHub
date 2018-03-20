# here we define reactive values/variables
scaterReads <- reactive({
  if(DEBUG)cat(file=stderr(), "scaterReads\n")
  gbm = gbm()
  gbm_log = gbm_log()
  fd = featureDataReact()
  if( is.null(gbm) | is.null(gbm_log))
    return(NULL)
  
  counts = as.matrix(exprs( gbm))
  
  anno = pData(gbm)
  anno$sample_id = anno$barcode
  anno$fixed = "red"
  # anno$individual= "NA1"
  # anno$replicate = "r1"
  # anno$well = "A01"
  # anno$batch = "b1"
  pheno_data <- new("AnnotatedDataFrame", anno)
  rownames(pheno_data) <- pheno_data$sample_id
  
  reads = as.matrix(counts)
  rownames(reads) = make.unique(fd[rownames(reads),"Associated.Gene.Name"])
  rownames(reads)[is.na(rownames(reads)) ] = "na"
  reads = SingleCellExperiment(
    assays = list(counts = reads), 
    colData = anno
  )
  ercc <- rownames(reads)[grepl("ERCC-", rownames(reads))]
  
  mt = rownames(reads)[grepl("^MT",rownames(reads))]
  
  reads <- scater::calculateQCMetrics(
    reads
  )
  filter_by_expr_features <- (reads$total_features > 200)
  reads$use <- (
    # sufficient features (genes)
    filter_by_expr_features 
    # sufficient molecules counted
    #filter_by_total_counts &
    # sufficient endogenous RNA
    #filter_by_ERCC &
    # remove cells with unusual number of reads in MT genes
    #filter_by_MT
  )
  reads
})

# declare function as heavy
myHeavyCalculations = list(c("scaterReads", "scaterReads"))
