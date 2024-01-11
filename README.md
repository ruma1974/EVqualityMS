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

# Use
library(EVqualityMS)

library(magrittr)

#- Test 1 - data in package

#- Mouse test data set. First columns are quantitative values and last column is gene symboles

data(DFmouseTest)

#- Human urinary EV test data set. First columns are quantitative values and last column is gene symboles
data(DFurine)

#- Data frame with main species supported
data(DFspecie)

# Test of functions

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

#- Test 8 list common species

EVqualityMS::DFspecie

#- Test 9 EV subpopulations
EVqualityMS::heatmapEVsub( DF=DFurine,Fac = NULL,color = NULL,CatColL = NULL,IncludeAlbumine = TRUE)

#- Test 10 EV CD marker coverage
EVqualityMS::heatmapEVsubCD( DFurine,Fac = NULL,color = NULL,CatColL = NULL,IncludeAlbumine = TRUE)


# FAQ

1. Do I need to install the packages each time I start R?
   Answer: No
2. EVqualityMS::HeatmapEVmarkers(DFspecie) gives an error
   Answer: DFspecie is a data frame containing information about supported species and is not test data.
3. EVqualityMS::HeatmapEVmarkers(DFmouseTest) gives a wired looking heatmap
   Answer: Data from non human species will require that the parameter specie="Mus musculus" is set in the function call
 

# please cite one of below articles

1. Extracellular Vesicles in Diffuse Large B Cell Lymphoma: Characterization and Diagnostic Potential.
Matthiesen R, Gameiro P, Henriques A, Bodo C, Moraes MCS, Costa-Silva B, Cabeçadas J, Gomes da Silva M, Beck HC, Carvalho AS.
Int J Mol Sci. 2022 Nov 1;23(21):13327. doi: 10.3390/ijms232113327.
PMID: 36362114 

2. Proteomic Landscape of Extracellular Vesicles for Diffuse Large B-Cell Lymphoma Subtyping.
Carvalho AS, Baeta H, Henriques AFA, Ejtehadifar M, Tranfield EM, Sousa AL, Farinho A, Silva BC, Cabeçadas J, Gameiro P, Silva MGD, Beck HC, Matthiesen R.
Int J Mol Sci. 2021 Oct 12;22(20):11004. doi: 10.3390/ijms222011004.

3. Is the Proteome of Bronchoalveolar Lavage Extracellular Vesicles a Marker of Advanced Lung Cancer?
Carvalho AS, Moraes MCS, Hyun Na C, Fierro-Monti I, Henriques A, Zahedi S, Bodo C, Tranfield EM, Sousa AL, Farinho A, Rodrigues LV, Pinto P, Bárbara C, Mota L, Abreu TT, Semedo J, Seixas S, Kumar P, Costa-Silva B, Pandey A, Matthiesen R.
Cancers (Basel). 2020 Nov 20;12(11):3450. doi: 10.3390/cancers12113450.
PMID: 33233545
