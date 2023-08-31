# Windows
The code below have been tested on several Windows systems using R4.3.1, Rtools (rtools43-5550-5548) and RStudio-2023.06.2-561. If you allready have the dependency packages install then you can jump directly to the section "Install EVqualityMS"

# Linux
On linux the below installation instructions should work directly. Drop a note if this is not the case.

# Install dependencies

BiocManager::install("orthogene")

BiocManager::install("ComplexHeatmap")

install.packages("scales")

install.packages("dplyr")

install.packages("gplots")

install.packages("ggpubr")

install.packages("magrittr")

install.packages("remotes")

# Install EVqualityMS

library(remotes)

remotes::install_github("ruma1974/EVqualityMS@master")

library(EVqualityMS)

library(magrittr)

#- Test 1 - data in package

data(DFmouseTest)

data(DFspecie)

data(DFurine)

#- Test 2 - EV and contaminant marker plot

EVqualityMS::HeatmapEVmarkers(DFurine,fontSizeRow = 20)

#- Test 3 - ISEV marker plot

res=EVqualityMS::heatmapEVqualityISEV( DFurine,Fac = NULL,color = NULL,CatColL = NULL,IncludeAlbumine = TRUE,plotCatDescription = TRUE)

#- Test 4 - mapping mouse gene symbols

Fac=append(rep("Normal",9),rep("HFD",9))
EVqualityMS::HeatmapEVmarkers(DFmouseTest,Fac=Fac,specie="Mus musculus",fontSizeRow = 20)

#- Test 5 - scatter plot of quality metrics for reference data

exosomeRM::scatterPlot()

#- Test 6 - scatter plot of quality metrics for all reference data and new data

EVqualityMS::scatterPlot(x = DFurine, Name = "Urine", src = "all")

#- Test 7 - scatter plot of quality metrics for selected reference data and new data

res=EVqualityMS::scatterPlot(x = DFurine, Name = "Urine", src = "plasma",palette =c("grey","yellow"))

#- Test 7 list common species

EVqualityMS::DFspecie
