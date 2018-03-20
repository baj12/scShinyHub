# here we define reactive values/variables
# e.g.
# pca = reactive({
#   if(DEBUG)cat(file=stderr(), "pca\n")
#   gbm_log = gbm_log()
#   if(is.null(gbm_log)){
#    if(DEBUG)cat(file=stderr(), "pca:NULL\n")
#    return(NULL)

#   }
#   if(DEBUG)cat(file=stderr(), "pca:done\n")
#   return(run_pca(gbm_log))
# })

# declare heavy calculations
myHeavyCalculations = NULL
