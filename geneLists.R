# download from https://www.proteinatlas.org/about/download protienatlas.tsv.zip
# uncompress
library(stringi)

rm(list = ls())


dt = read.delim(file = "proteinatlas.tsv", sep = "\t", stringsAsFactors = FALSE)
load('Jim2.Rds')
rownames(featuredata)
rownames(dt) = dt$Ensembl
rownames(featuredata)[!rownames(featuredata) %in% rownames(dt)]
sum(!rownames(featuredata) %in% rownames(dt))
featureIdx = rownames(featuredata) %in% rownames(dt)
dt[rownames(featuredata)[featureIdx], "featureNames"] = featuredata[featureIdx,"Associated.Gene.Name"]
dat = dt[!is.na(dt$featureNames),]


proteinClasses = unique(stri_trim(unlist(strsplit(dat$Protein.class, ","))))
proteinClassList = list()
for(clIdx in 1:length(proteinClasses)){
  proteinClassList[[proteinClasses[clIdx]]] = dat[grep(proteinClasses[clIdx], dat$Protein.class), "Ensembl"]
}

hasAntibody = list("has Antibody" = dat[dat$Antibody == "", "Ensembl"])

locClasses = unique(stri_trim(unlist(strsplit(dat$Subcellular.location, "<br>"))))
subCellularLocationList = list()
for(clIdx in 1:length(locClasses)){
  subCellularLocationList[[locClasses[clIdx]]] = dat[grep(locClasses[clIdx], dat$Subcellular.location), "Ensembl"]
}

geneLists = list(
  "protein Classes" = proteinClassList,
  "Antibody" = hasAntibody,
  "subcellular Location" = subCellularLocationList
)

save(file = "geneLists.RData", list = c("geneLists"))
