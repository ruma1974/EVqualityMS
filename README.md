# install dependencies

BiocManager::install("orthogene")

install.packages("scales")

install.packages("dplyr")

install.packages("ComplexHeatmap")

install.packages("gplots")

install.packages("ggpubr")

# install EVqualityMS

library(remotes)

remotes::install_github("ruma1974/EVqualityMS@master")

library(EVqualityMS)

#- test 1 - data in package
data(DFmouseTest)
data(DFspecie)
data(DFurine)

#- test 2 - EV and contaminant marker plot

EVqualityMS::HeatmapEVmarkers(DFurine,fontSizeRow = 20)

#- test 3 - ISEV marker plot

res=EVqualityMS::heatmapEVqualityISEV( DFurine,Fac = NULL,color = NULL,CatColL = NULL,IncludeAlbumine = TRUE,plotCatDescription = TRUE)

#- test 4 - mapping mouse gene symbols

Fac=append(rep("Normal",9),rep("HFD",9))
EVqualityMS::HeatmapEVmarkers(DFmouseTest,Fac=Fac,specie="Mus musculus",fontSizeRow = 20)

#- test 5 - scatter plot of quality metrics for reference data
exosomeRM::scatterPlot()

#- test 6 - scatter plot of quality metrics for all reference data and new data
EVqualityMS::scatterPlot(x = DFurine, Name = "Urine", src = "all")

#- test 7 - scatter plot of quality metrics for selected reference data and new data

res=EVqualityMS::scatterPlot(x = DFurine, Name = "Urine", src = "plasma",palette =c("grey","yellow"))

#- test 7 list common species

EVqualityMS::DFspecie