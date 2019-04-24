dir1='~/Rstudio/shHubgit/'
dir2='~/Rstudio/scPackageTest/inst/app/'

type='*.R$'

fps1l = dir(path = dir1, pattern = type, recursive = T, full.names = T)
fps1s = dir(path = dir1, pattern = type, recursive = T, full.names = F)
fps2l = dir(path = dir2, pattern = type, recursive = T, full.names = T)
fps2s = dir(path = dir2, pattern = type, recursive = T, full.names = F)
library(diffr)
for (fp in fps1s) {
  if (fp %in% fps2s) {
    idx1 = which (fp == fps1s)
    idx2 = which (fp == fps2s)
    print(diffr(fps1l[idx1], fps2l[idx2]))
  }
}
