prioritized_genes = reactive({
  tsne.data = tsne.data()
  gbm = gbm()
  if(is.null(tsne) | is.null(gbm)){
    if(DEBUG)cat(file=stderr(), "tsne.data: NULL\n")
    return(NULL)
  }
  set.seed(seed = seed)
  prioritize_top_genes(gbm, tsne.data$dbCluster, "sseq",
                       logscale = FALSE, 
                       min_mean=0.5, 
                       p_cutoff=0.05,
                       order_by='pvalue')
})
