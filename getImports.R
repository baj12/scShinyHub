fps = dir(path = ".", pattern = "*.R$", recursive = TRUE, full.names = TRUE )

pattern = ".* ([[:alpha:]]*)\\(.*"

functionList = data.frame("filename"=character(), "func"=as.character(), "matchFun"=character(), stringsAsFactors = F)
options(stringsAsFactors = FALSE)
rn = 0
for (fp in fps) {
  print (fp)
  getImports(fp)
  
}


uniqFuncs = unique(functionList[functionList$filename=="./server.R" & !functionList$matchFun == "", ])


getImports <- function(fp) {
  require(dplyr)
  options(stringsAsFactors = FALSE)
  pattern = ".* ([[:alpha:]]*)\\(.*"
  functionList = data.frame("filename"=character(), "func"=as.character(), "matchFun"=character(), stringsAsFactors = F)
  for (ln in readLines(fp)){
    if (length(grep(pattern, ln)) > 0) {
      mt = gsub (pattern, "\\1", ln)
      matFun = tryCatch(match.fun(mt), error = function(e){})
      if (nchar(mt) > 0) {
        # print(c(fp, mt, environmentName(environment(matFun))))
        functionList = rbind(functionList, c(fp, mt, environmentName(environment(matFun))))
      }
    }
  }
  colnames(functionList) = c("filename", "func", "matchFun")
  uniqFuncs = unique(functionList[!functionList$matchFun == "", ])
  uniqFuncs %>% dplyr::group_by(matchFun) %>% 
    dplyr::mutate(allFunc = paste0(func, collapse = " "))   %>% 
    select(filename, matchFun, allFunc) %>% 
    unique %>% 
    print
  print(uniqFuncs)
}

