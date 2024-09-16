# install dependencies
if (!require(BiocManager, quietly = TRUE)) install.packages("BiocManager")

if (!require(orthogene, quietly = TRUE)) BiocManager::install("orthogene")

if (!require(ComplexHeatmap, quietly = TRUE)) BiocManager::install("ComplexHeatmap")

if (!require(scales, quietly = TRUE)) install.packages("scales")

if (!require(dplyr, quietly = TRUE)) install.packages("dplyr")

if (!require(gplots, quietly = TRUE)) install.packages("gplots")

if (!require(ggpubr, quietly = TRUE)) install.packages("ggpubr")

if (!require(remotes, quietly = TRUE)) install.packages("remotes")

if (!require(magrittr, quietly = TRUE)) install.packages("magrittr")

# install EVqualityMS

library(remotes)

remotes::install_github("ruma1974/EVqualityMS@master")

library(EVqualityMS)

#- test 1 - data in package

data(DFmouseTest)

data(DFspecie)

data(DFurine)

#- test 2 - EV and contaminant marker plot

res=EVqualityMS::HeatmapEVmarkers(DFurine,fontSizeRow = 20)
res$Heatmap

#- save and change dimensions
png("filename.png", res = 300, width = 4500, height = 2000)
print(res$Heatmap)
dev.off()

#- test 3 - ISEV marker plot

res=EVqualityMS::heatmapEVqualityMISEV( DFurine,Fac = NULL,color = NULL,CatColL = NULL,IncludeAlbumine = TRUE,plotCatDescription = TRUE)

#- test 4 - mapping mouse gene symbols

Fac=append(rep("Normal",9),rep("HFD",9))

res=EVqualityMS::HeatmapEVmarkers(DFmouseTest,Fac=Fac,specie="Mus musculus",fontSizeRow = 20)
res$Heatmap

#- test 5 - scatter plot of quality metrics for reference data

EVqualityMS::scatterPlot()

#- test 6 - scatter plot of quality metrics for all reference data and new data

EVqualityMS::scatterPlot(x = DFurine, Name = "Urine", src = "all")

#- test 7 - scatter plot of quality metrics for selected reference data and new data

res=EVqualityMS::scatterPlot(x = DFurine, Name = "Urine", src = "plasma",palette =c("grey","yellow"))

#- test 8 list common species

EVqualityMS::DFspecie

#- test 9 MaxQuant import

filename_txt="path/proteinGroups.txt"
DF=EVqualityMS::ImportMaxQuantProteinGroups(filename_txt, data_s = "iBAQ")

res=EVqualityMS::heatmapEVqualityMISEV( DF,Fac = NULL,color = NULL,CatColL = NULL,IncludeAlbumine = TRUE,plotCatDescription = TRUE)
res
