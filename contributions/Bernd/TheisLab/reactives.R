# require(ggplot2)
# 
# gbmPheatmap <- function(gbm, genes_to_plot, cells_to_plot, n_genes = 5, colour = NULL,
#                         limits = c(-3, 3)) {
#   if (DEBUG) cat(file = stderr(), "gbmPheatmap\n")
#   if (DEBUGSAVE) {
#     save(file = "~/scShinyHubDebug/gbmPheatmap.RData", list = c(ls(), ls(envir = globalenv())))
#   }
#   # load(file="~/scShinyHubDebug/gbmPheatmap.RData")
#   if (!is.list(genes_to_plot)) {
#     cat("Plotting one gene set instead of multiple cluster-specific gene sets\n")
#     gene_indices <- sapply(genes_to_plot, function(x) get_gene_index(
#         gbm,
#         x
#       ))
#     gene_annotation <- NULL
#   } else {
#     if ("significant" %in% names(genes_to_plot[[1]])) {
#       gene_indices <- unlist(lapply(genes_to_plot, function(x) with(
#           x,
#           head(ix[significant], n_genes)
#         )))
#       gene_grouping <- unlist(lapply(
#         names(genes_to_plot),
#         function(nm) rep(nm, with(
#             genes_to_plot[[nm]],
#             length(head(ix[significant], n_genes))
#           ))
#       ))
#     }
#     else {
#       gene_indices <- unlist(lapply(genes_to_plot, function(x) x$ix[1:n_genes]))
#       gene_grouping <- rep(names(genes_to_plot), each = n_genes)
#     }
#     gene_annotation <- data.frame(ClusterID = as.factor(gene_grouping))
#   }
#   cell_indices <- unlist(lapply(cells_to_plot, function(x) x$ix))
#   value <- t(scale(t(as.matrix(exprs(gbm))[gene_indices, cell_indices])))
#   value[value < limits[1]] <- limits[1]
#   value[value > limits[2]] <- limits[2]
#   rownames(value) <- make.unique(fData(gbm)$symbol[gene_indices])
#   cell_grouping <- unlist(lapply(1:length(cells_to_plot), function(x) {
#     rep(names(cells_to_plot)[x], length(cells_to_plot[[x]]$barcode))
#   }))
#   cell_annotation <- data.frame(ClusterID = as.factor(cell_grouping))
#   rownames(cell_annotation) <- colnames(value)
#   if (!is.null(gene_annotation)) {
#     rownames(gene_annotation) <- rownames(value)
#   }
#   if (is.null(colour)) {
#     anno_colors <- NULL
#   } else {
#     names(colour) <- names(cells_to_plot)
#     anno_colors <- list(ClusterID = colour)
#   }
#   list(
#     mat = value, cluster_rows = FALSE, cluster_cols = FALSE,
#     show_colnames = FALSE, annotation_row = gene_annotation,
#     annotation_col = cell_annotation, annotation_names_row = FALSE,
#     annotation_names_col = FALSE, annotation_colors = anno_colors
#   )
# }
# 
# 
# crHeatImage_func <- function(gbm, projections, prioritized_genes) {
#   example_K <- length(levels(projections$dbCluster))
#   # max is 30 for number of clusters
#   if (example_K > 12) {
#     example_Cols <- rep(rev(brewer.pal(12, "Set3")), 3)[1:example_K]
#   } else {
#     example_Cols <- rev(brewer.pal(12, "Set3")) # customize plotting colors
#   }
#   cells_to_plot <- order_cell_by_clusters(gbm, as.numeric(as.character(projections$dbCluster)))
# 
#   example_col <- example_Cols[1:example_K]
# 
#   # For high-res displays, this will be greater than 1
#   # if (!is.null(session)){
#   #   pixelratio <- session$clientData$pixelratio
#   #   width <- session$clientData$output_plot_width
#   #   height <- session$clientData$output_plot_height
#   # }else{
#   pixelratio <- NULL
#   width <- NULL
#   height <- NULL
#   # }
#   if (is.null(pixelratio)) pixelratio <- 1
#   if (is.null(width)) {
#     width <- 96 * 7
#   } # 7x7 inch output
#   if (is.null(height)) {
#     height <- 96 * 12
#   }
# 
#   # px to inch conversion
#   myPNGwidth <- width / 96
#   myPNGheight <- height / 96
# 
#   # outfile <- paste0(tempdir(),'/crHeatImage.svg')
#   outfile <- paste0(tempdir(), "/crHeatImage.png")
#   if (DEBUG) cat(file = stderr(), paste("output file: ", outfile, "\n"))
#   if (DEBUG) cat(file = stderr(), paste("output file normalized: ", normalizePath(outfile, mustWork = FALSE), "\n"))
#   logGB <- log_gene_bc_matrix(gbm)
# 
#   retVal <- tryCatch({
#     gbmPheatmap(
#       gbm = logGB,
#       genes_to_plot = prioritized_genes,
#       cells_to_plot = cells_to_plot,
#       n_genes = 10,
#       colour = example_col,
#       limits = c(-3, 3)
#     )
#   },
#   error = function(cond) {
#     cat(file = stderr(), "crHeatImage: problem, in gbm_pheatmap?\n")
#     if (!is.null(getDefaultReactiveDomain())) {
#       removeNotification(id = "crHeatMap")
#     }
#     return(NULL)
#   },
#   warning = function(cond) {
#     cat(file = stderr(), "crHeatImage: problem, in gbm_pheatmap?\n")
#     if (!is.null(getDefaultReactiveDomain())) {
#       removeNotification(id = "crHeatMap")
#     }
#     return(NULL)
#   }
#   )
# }
# 
# crHeatImage <- reactive({
#   if (DEBUG) cat(file = stderr(), "output$crHeatImage\n")
#   gbm <- gbm()
#   # projections = projections()
#   projections <- projections()
#   prioritized_genes <- prioritized_genes()
#   if (is.null(gbm) | is.null(gbm) | is.null(prioritized_genes)) {
#     if (DEBUG) cat(file = stderr(), "output$crHeatImage:NULL\n")
#     return(NULL)
#   }
#   if (!is.null(getDefaultReactiveDomain())) {
#     showNotification("cell ranger heat map", id = "crHeatMap", duration = NULL)
#   }
#   if (DEBUGSAVE) {
#     save(file = "~/scShinyHubDebug/crHeatImage.RData", list = c(ls(), ls(envir = globalenv())))
#   }
#   # load(file='~/scShinyHubDebug/crHeatImage.RData')
#   retVal <- crHeatImage_func(gbm, projections, prioritized_genes)
#   if (!is.null(getDefaultReactiveDomain())) {
#     removeNotification(id = "crHeatMap")
#   }
# 
#   return(retVal)
# })
# 
# 
# prioritized_genes_func <- function(gbm, projections, seed) {
#   set.seed(seed = seed)
#   retVal <- tryCatch({
#     prioritize_top_genes(gbm, as.numeric(as.character(projections$dbCluster)), "sseq",
#       logscale = FALSE,
#       min_mean = 0.5,
#       p_cutoff = 0.05,
#       order_by = "pvalue"
#     )
#   },
#   error = function(cond) {
#     cat(file = stderr(), "prioritized_genes.Rdata: problem, are id and sample column persent?\n")
#     if (!is.null(getDefaultReactiveDomain())) {
#       removeNotification(id = "crpriotGenes")
#     }
#     return(NULL)
#   },
#   warning = function(cond) {
#     return("prioritize warning")
#   }
#   )
#   if (!is.null(getDefaultReactiveDomain())) {
#     removeNotification(id = "crpriotGenes")
#   }
# 
#   return(retVal)
# }
# 
# crPrioGenesTable <- reactive({
#   require("forcats")
#   require("tidyverse")
#   if (DEBUG) cat(file = stderr(), "output$crPrioGenes\n")
#   prioritized_genes <- prioritized_genes()
#   # cl5 = input$cluster5
#   # if (is.null(prioritized_genes) | is.null(cl5)) {
#   if (is.null(prioritized_genes)) {
#     return(NULL)
#   }
# 
# 
#   if (DEBUGSAVE) {
#     save(file = "~/scShinyHubDebug/crPrioGenes.RData", list = c(ls(), ls(envir = globalenv())))
#   }
#   # load(file="~/scShinyHubDebug/crPrioGenes.RData")
# 
#   dt <- data.frame()
#   for (listIter in 1:length(prioritized_genes)) {
#     prioritized_genes[[listIter]]$cluster <- listIter
#     dt <- rbind(dt, prioritized_genes[[listIter]])
#   }
#   rownames(dt) <- make.unique(as.character(dt$gene_name), sep = "___")
#   dt$cluster <- factor(dt$cluster)
# 
#   # move cluster column to second position
#   cnames <- colnames(dt)
#   clNr <- which(cnames == "cluster")
#   sigCol <- which(cnames == "significant")
#   adjCol <- which(cnames == "p_adj")
#   dt <- dt[, c(1, clNr, sigCol, adjCol, c(1:length(cnames))[-c(1, clNr, sigCol, adjCol)])]
# 
#   return(dt)
# })
# 
# 
# prioritized_genes <- reactive({
#   projections <- projections()
#   gbm <- gbm()
#   if (is.null(projections) | is.null(gbm)) {
#     if (DEBUG) cat(file = stderr(), "projections: NULL\n")
#     return(NULL)
#   }
#   if (!is.null(getDefaultReactiveDomain())) {
#     showNotification("prioritizing genes", id = "crpriotGenes", duration = NULL)
#   }
# 
#   if (DEBUGSAVE) {
#     save(file = "~/scShinyHubDebug/prioritized_genes.Rdata", list = c(ls(), ls(envir = globalenv())))
#   }
#   # load(file='~/scShinyHubDebug/prioritized_genes.Rdata')
#   retVal <- prioritized_genes_func(gbm, projections, seed)
#   set.seed(seed = seed)
# 
#   return(retVal)
# })
# 
# # myHeavyCalculations = list(c("prioritized_genes", "prioritized_genes"), c("crHeatImage", "crHeatImage"))
