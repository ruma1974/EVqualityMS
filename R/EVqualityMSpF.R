# _"/media/rune/Jochime2/RscripterCool/Rprojects/rmpackges/EVqualityMSp.R"


calculate_coverage_percentages <- function(sL1, sL2,caseSensitive=TRUE) {
  if (length(sL1) == 0 || length(sL2) == 0) {
    stop("Both vectors must contain elements.")
  }
  
# Convert vectors to lowercase to make the comparison case-insensitive
if (caseSensitive==FALSE){
  sL1 <- tolower(sL1)
  sL2 <- tolower(sL2)	
}
# Initialize an empty list to store coverage percentages for each term in vector1
coverage_percentages <- list()
# Loop through each term in vector1
sLu=unique(sL2)
for (term in sLu) {
# Calculate the number of times the term appears in vector2
term_count_in_vector1 <- sum(sL1 == term)
term_count_in_vector2 <- sum(sL2 == term)
# Calculate the percentage of coverage for the term
coverage_percentage <- ( term_count_in_vector1/term_count_in_vector2) * 100
# Store the coverage percentage in the list
coverage_percentages[[term]] <- coverage_percentage
}  
return(coverage_percentages)
}



#' Calculate summary stat from DF
#' @param DF (DF - first column is values and second column is grouping vector)
#' @return list
#' @keywords stat
#' @examples
#' ---
#' @export
statFunc <- function(DF){
#
#
names(DF)=c("numbers","groups")
res <- doBy::summaryBy(numbers ~ groups, data = DF, FUN = function(x) {
  c(
    mean = mean(x,na.rm=TRUE),
    median = median(x,na.rm=TRUE),
    IQR = IQR(x,na.rm=TRUE),
    min = min(x,na.rm=TRUE),
    max = max(x,na.rm=TRUE),
    SD = sd(x,na.rm=TRUE)
  )
})

return(res) # invisible(res)
}





#' Quality control of EV data
#' @param x (data frame with a GN column)
#' @param Fac (Fac - should have length ncol(DF)-1)
#' @param color=NULL (vector with ten colors)
#' @param plotHeat (bool)
#' @param IncludeClinical (bool)
#' @param IncludeAlbumine (bool)
#' @param fontSizeCol (int)
#' @param fontSizeRow (int)
#' @param specie (string)
#' @return list
#' @keywords Quality_control
#' @examples
#' ---
#' @export
HeatmapEVmarkers <- function(x,Fac=NULL,color=NULL,plotHeat=TRUE,IncludeClinical=TRUE,IncludeAlbumine=TRUE,fontSizeCol=12,fontSizeRow=12,specie="Homo sapiens"){
# 
#
if (is.null(color)){color=c("grey60","grey30","gray15","steelblue3","navy","seagreen1","seagreen4","palegreen","violetred4","red")}
#if (length(color)!=10){stop("Length of color must be 10")}
if (!is.data.frame(x)){stop("x must be a data frame with a column named GN")}
# median of sample groups
DF=x
if (!is.null(Fac)){
idx=which(names(DF)=="GN")
mat=DF[,-idx]
FUN = function(x) {median(x, na.rm = T)}
mat=t(apply(mat, 1, function(x) {tapply(x, Fac, FUN = FUN)}))
x=as.data.frame(mat)
x$GN=DF[,idx]
}
# map to humans if not allready
if (specie!="Homo sapiens"){
gnL=x$GN
#gnL
res <- orthogene::convert_orthologs(gene_df = gnL, input_species = specie, output_species = "Homo sapiens", method = "gprofiler")
x$GN[res$index]=rownames(res)
}

# rank scale data

for (i in 1:ncol(x)){
if (names(x)[i]!="GN"){
xL=as.numeric(x[,i])
xL[is.na(xL)]=0
xL=scales::rescale(xL,to=c(0,1))

# Rank with ties receiving the same ranking
xL <- rank(xL, ties.method = "average")
# Scale ranks from 0 to 1
xL=xL-min(xL)
xL=xL/max(xL)
x[,i]=xL
}
}
# remove duplicated GN
idx=which(names(x) %in% "GN")
Rsum=rowSums(as.matrix(x[,-idx]))
idxL=order(Rsum,decreasing =TRUE)
x=x[idxL,]
x=x[!duplicated(x$GN),]


# merge with markers data
DF=EVqualityMS::ExoData[["Expression"]]
DF=dplyr::left_join(DF, x, by = "GN")

# remove NAs
DF[is.na(DF)]=0
mat=DF[,2:ncol(DF)]
rownames(mat)=DF$GN
#
if (IncludeClinical==FALSE){
iL=which(colnames(mat)=="BALexoNo" | colnames(mat)=="BALexoYes")
mat=mat[,-iL]
}
# IncludeAlbumine
if (IncludeAlbumine==FALSE){mat=mat[rownames(mat)!="ALB",]}
# plot
if (plotHeat==TRUE){
# set up column colors of labels
colL=rep("black",nrow(mat))
# Contaminants

for (i in 1:length(EVqualityMS::ExoData$Contaminants)){
id2=which(EVqualityMS::ExoData$Contaminants[i]==rownames(mat))
colL[id2]="red"
}
# TopMarker
for (i in 1:length(EVqualityMS::ExoData$TopMarker)){
id2=which(EVqualityMS::ExoData$TopMarker[i]==rownames(mat))
colL[id2]="green"
}
# 
for (i in 1:length(EVqualityMS::ExoData$lEV)){
id2=which(EVqualityMS::ExoData$lEV[i]==rownames(mat))
colL[id2]="blue"
}
# plot
x11(w=20,h=10)
p=ComplexHeatmap::Heatmap(t(mat), col = color, name = "ranked expression",heatmap_legend_param = list(at = c(0, 0.25, 0.5, 0.75, 1)),column_names_gp = ggfun::gpar(col = colL,fontsize = fontSizeCol),row_names_gp = grid::gpar(fontsize = fontSizeRow))

print(p)
}
# print Q index 
idxCL=which(x$GN %in% EVqualityMS::ExoData[["Contaminants"]])
idxEL=which(x$GN %in% EVqualityMS::ExoData[["TopMarker"]])
if (length(idxCL)>0){
qL=colSums(x[idxEL,-idx,drop=FALSE])/length(EVqualityMS::ExoData[["TopMarker"]])-colSums(x[idxCL,-idx,drop=FALSE])/length(EVqualityMS::ExoData[["Contaminants"]])
} else {
qL=colSums(x[idxEL,-idx,drop=FALSE])/length(EVqualityMS::ExoData[["TopMarker"]])
}
if (length(idxEL)==0){
qL=0
} 
cat("Overall quality metrics:\n")
print(qL)
# calculate contamination and exosome markers
if (length(idxCL)>0){
qC=colSums(x[idxCL,-idx,drop=FALSE])/length(EVqualityMS::ExoData[["Contaminants"]])
} else {qC=0}
if (length(idxEL)>0){
qE=colSums(x[idxEL,-idx,drop=FALSE])/length(EVqualityMS::ExoData[["TopMarker"]])
} else {qE=0}
# collect results
res=list()
res[["Expression"]]=mat
res[["Quality"]]=qL 
res[["Contaminants"]]=qC
res[["EV marker"]]=qE
if (plotHeat==TRUE){
res[["Heatmap"]]=p	
}
invisible(res)
}

 



#' Plot heatmap with markers from ISEV paper
#' @param DF (DF - data frame with a GN column)
#' @param Fac (Fac - should have length ncol(DF)-1)
#' @param color (color - vector with ten colors for type of marker)
#' @param CatColL (color - protein category)
#' @param IncludeAlbumine (bool )
#' @param plotCatDescription (bool )
#' @return list
#' @keywords Quality_control
#' @examples
#' ---
#' @export
heatmapEVqualityMISEV <- function(DF,Fac=NULL,color = NULL,CatColL=NULL,IncludeAlbumine = TRUE,plotCatDescription=TRUE){
#
#
if (is.null(color)){color=c("grey60","grey30","gray15","steelblue3","navy","seagreen1","seagreen4","palegreen","violetred4","red")}
if (is.null(CatColL)){CatColL=c("#F0A3FF", "#0075DC", "#993F00", "#4C005C", "#191919", "#005C31", 
"#2BCE48", "#FFCC99", "#808080", "#94FFB5", "#8F7C00", "#9DCC00", 
"#C20088", "#003380", "#FFA405")}

if (!is.data.frame(DF)){stop("x must be a data frame with a column named GN")}
#
x=DF
# median of sample groups
if (!is.null(Fac)){
idx=which(names(DF)=="GN")
mat=DF[,-idx]
FUN = function(x) {median(as.numeric(x), na.rm = T)}
mat=t(apply(mat, 1, function(x) {tapply(x, Fac, FUN = FUN)}))
x=as.data.frame(mat)
x$GN=DF[,idx]
#print(x)
#tmp1=x
}
# rank scale data

for (i in 1:ncol(x)){
if (names(x)[i]!="GN"){
xL=as.numeric(x[,i])
xL[is.na(xL)]=0
xL=scales::rescale(xL,to=c(0,1))
# Rank with ties receiving the same ranking
xL <- rank(xL, ties.method = "average")
# Scale ranks from 0 to 1
#xL <- (xL - 1) / (length(xL) - 1)
#print(xL)
#print(min(xL))
#print(max(xL))
#print("---------")
xL=xL-min(xL)
xL=xL/max(xL)
x[,i]=xL
}
}
#tmp2=x
# remove duplicated GN keep the one with most intensity
idx=which(names(x) %in% "GN")
Rsum=rowSums(as.matrix(x[,-idx]))
idxL=order(Rsum,decreasing =TRUE)
#Rsum[idxL]
x=x[idxL,]
x=x[!duplicated(x$GN),]
#str(as.matrix(x[,-idx]))

# merge with markers data
# DF=EVqualityMS::ExoData[["Expression"]]
#which(DFanno$GN=="ENO1")

DF=data.frame(GN=EVqualityMS::DFanno$GN)
DF= dplyr::left_join(DF, x, by = "GN")

# remove NAs
DF[is.na(DF)]=0
# move gene names to last column
mat=DF[,2:ncol(DF)]
rownames(mat)=DF$GN
# remove zero rows
idxL=which(rowSums(mat)==0)
if (length(idxL)>0){
mat=mat[-idxL,]	
}
# Include Albumine
if (IncludeAlbumine==FALSE){mat=mat[rownames(mat)!="ALB",]}
# plot
colL=rep("black",nrow(mat))
# Contaminants
conL=append(c("APOA1","APOA2","APOB","UMOD"),EVqualityMS::ExoData$Contaminants)
for (i in 1:length(conL)){
id2=which(conL[i]==rownames(mat))
colL[id2]="red"
}
# TopMarker
for (i in 1:length(EVqualityMS::ExoData$TopMarker)){
id2=which(EVqualityMS::ExoData$TopMarker[i]==rownames(mat))
colL[id2]="green"
}
# microvesicles
for (i in 1:length(EVqualityMS::ExoData$lEV)){
id2=which(EVqualityMS::ExoData$lEV[i]==rownames(mat))
colL[id2]="blue"
}
# cluster
clustL=rep("",nrow(mat))

for (i in 1:nrow(mat)){
#print(rownames(mat)[i])
id=which(rownames(mat)[i]==EVqualityMS::DFanno$GN)
#print(id)
if (id>=0){

if (plotCatDescription==TRUE){
clustL[i]=EVqualityMS::DFanno$CatDes[id]
} else {
clustL[i]=EVqualityMS::DFanno$Cat[id]
}
}
}
#
#print(clustL)
#

Ncol=length(unique(clustL))
CatColL=CatColL[1:Ncol]
names(CatColL)=unique(clustL)
#print(CatColL)

#
p1 = ComplexHeatmap::Heatmap(t(mat), name = "mat",col = color, column_names_gp = ggfun::gpar(col = colL),cluster_columns = ComplexHeatmap::cluster_within_group(t(mat),as.numeric(factor(clustL))),heatmap_legend_param = list(position = "top", title = "Regulation", legend_direction = "horizontal"),top_annotation = ComplexHeatmap::HeatmapAnnotation(Category = clustL, col = list(Category = CatColL)))
p1 = ComplexHeatmap::draw(p1, heatmap_legend_side = "bottom")

# print Q index 
idxCL=which(x$GN %in% EVqualityMS::ExoData[["Contaminants"]])
idxEL=which(x$GN %in% EVqualityMS::ExoData[["TopMarker"]])
if (length(idxCL)>0){
qL=colSums(x[idxEL,-idx,drop=FALSE])/length(EVqualityMS::ExoData[["TopMarker"]])-colSums(x[idxCL,-idx,drop=FALSE])/length(EVqualityMS::ExoData[["Contaminants"]])
} else {
qL=colSums(x[idxEL,-idx,drop=FALSE])/length(EVqualityMS::ExoData[["TopMarker"]])
}
if (length(idxEL)==0){
qL=0
} 
cat("Overall quality metrics:\n")
print(qL)
# calculate contamination and exosome markers
if (length(idxCL)>0){
qC=colSums(x[idxCL,-idx,drop=FALSE])/length(EVqualityMS::ExoData[["Contaminants"]])
} else {qC=0}
if (length(idxEL)>0){
qE=colSums(x[idxEL,-idx,drop=FALSE])/length(EVqualityMS::ExoData[["TopMarker"]])
} else {qE=0}
# collect results
res=list()
res[["Expression"]]=mat
res[["Quality"]]=qL 
res[["Contaminants"]]=qC
res[["EV marker"]]=qE
res[["Heatmap"]]=p1
#res[["tmp1"]]=tmp1
#res[["tmp2"]]=tmp2
invisible(res)
}


#' plot scatter of relative contamination and EV marker abundance
#' @param x (DF- quantitative values in columns and last column holds the gene names - GN)
#' @param Name (string - name of the new samples)
#' @param Title (string - plot title) 
#' @param src (string - which reference samples to include)
#' @param palette (colorL or palette)
#' @param jitter_b (boolean)
#' @param Psize (float)
#' @param scaleData (boolean)
#' @param NoDensity (boolean)
#' @param addRefLines (boolean)
#' @return ggplot
#' @keywords Quality_control
#' @examples
#' ---
#' @export
scatterPlot <- function(x=NULL,Name="new",Title="",src="all",palette="jco",jitter_b=TRUE,Psize=3,scaleData=TRUE,NoDensity=TRUE,addRefLines=FALSE){
#
#
#x = DFurine; Name = "Urine"; src = "plasma"
if (!is.null(x)){
res=EVqualityMS::HeatmapEVmarkers(x,plotHeat =FALSE)
if (Name[1]==""){
Name=names(res$Quality)
Name=sapply(strsplit(Name, "_"), function(x) x[1])
DF=data.frame(EVmarker=res$"EV marker",Contaminants=res$Contaminants,source=Name,sample="")	
}else {
DF=data.frame(EVmarker=res$"EV marker",Contaminants=res$Contaminants,source=Name,sample="")
}
print(DF)
}
DFref=EVqualityMS::DFqualityRef
# select reference samples
if (src!="all"){DFref=DFref[DFref$source==src,]}
#
if (!is.null(x)){
DFref=rbind(DFref,DF)	
}
#
# palette = c("red", "green", "blue")
## sp <- ggpubr::ggscatter(DFref, x =  "EVmarker", y ="Non-EV marker",color = "source", palette = palette,size = 3, alpha = 0.6)+ggpubr::border()
# Marginal density plot of x (top panel) and y (right panel)
## xplot <- ggpubr::ggdensity(DFref, x="EVmarker", fill = "source",palette = palette)
## xplot <- xplot + ggpubr::clean_theme() + ggpubr::rremove("legend")
# Arranging the plot using cowplot
## p <- cowplot::plot_grid(xplot, sp,ncol = 1, align = "hv", rel_heights = c(1, 1.5,1))
if (jitter_b==TRUE){
H=0.01
W=0.01	
} else {
H=0
W=0	
}
#data.frame(c(0,1,0,1),c(0,1,1,0))
if (scaleData==TRUE){
p=ggpubr::ggscatterhist(
  DFref, x = "EVmarker", y = "Contaminants",position=ggplot2::position_jitter(h=H,w=W),title = Title,
  color = "source", size = Psize,shape ="source", alpha = 0.6,font.tickslab = c(20, "plain"),font.legend = c(16, "plain"),font.x=c(16, "plain"),font.y=c(16, "plain"),
  palette = palette,
  margin.params = list(fill = "source", color = "black", size = 0.2)
  )
#p=p+  ggplot2::geom_point(data = DF, ggplot2::aes(x = "EVmarker", y ="Contaminants"), color = "red", size = 3)
} else {
# Main plot
p <-  ggplot2::ggplot(DFref, ggplot2::aes(x = EVmarker, y = Contaminants, color = source))+ggplot2::geom_point()+ggplot2::scale_color_manual(values=palette)+ ggplot2::xlim(0, 1)+ ggplot2::ylim(0, 1) +ggplot2::geom_point(size=3) # +ggplot2::xlab('EV marker') + ggplot2::ylab('Non-EV marker')
#
p = p + ggplot2::theme_bw()
p = p + ggplot2::theme(axis.text.x = ggplot2::element_text(family = "Arial", face = "plain", colour = "black", size = 20, hjust = 0.5, vjust = 0.5, angle = 0))
p = p + ggplot2::theme(axis.text.y = ggplot2::element_text(family = "Arial", face = "plain", colour = "black", size = 20, hjust = 0.5, vjust = 0.5, angle = 0))
p = p + ggplot2::theme(axis.title.x = ggplot2::element_text(family = "Arial", face = "plain", colour = "black", size = 22, hjust = 0.5, vjust = 0.5))
p = p + ggplot2::theme(axis.title.y = ggplot2::element_text(family = "Arial", face = "plain", colour = "black", size = 22, hjust = 0.5, vjust = 0.5))
p= p + ggplot2::labs(title=Title,x ="EV marker", y = "Non-EV marker")
#
p=p+ggplot2::theme(legend.position="bottom")+ ggplot2::theme(legend.text=ggplot2::element_text(size=16))
#
if (length(unique(DFref[,c("source")]))<=4 & NoDensity==FALSE){
# Marginal densities along x axis
xdens <-  cowplot::axis_canvas(p, axis = "x")+ggplot2::geom_density(data = DFref, ggplot2::aes(x = EVmarker, fill = source),alpha = 0.7, linewidth = 0.2)+ggplot2::scale_color_manual(values=palette)
#ggpubr::fill_palette("jco") #+ xlim(0, 1)
# Marginal densities along y axis

# Need to set coord_flip = TRUE, if you plan to use coord_flip()
ydens <-  cowplot::axis_canvas(p, axis = "y", coord_flip = TRUE)+ggplot2::geom_density(data = DFref, ggplot2::aes(x = Contaminants, fill = source),alpha = 0.7, size = 0.2)+ggplot2::coord_flip()+ggplot2::scale_color_manual(values=palette)
#ggpubr::fill_palette("jco") 

p <- cowplot::insert_xaxis_grob(p, xdens, grid::unit(.2, "null"), position = "top")
p <- cowplot::insert_yaxis_grob(p, ydens, grid::unit(.2, "null"), position = "right")
p=cowplot::ggdraw(p)

} 
}
if (addRefLines==TRUE){
p=p+ggplot2::geom_hline(yintercept=0.16, linetype="dashed", color = "red")	
p=p+ggplot2::geom_hline(yintercept=0.32, linetype="dashed", color = "red")	
p=p+ggplot2::geom_vline(xintercept=0, linetype="dashed", color = "red")	
p=p+ggplot2::geom_vline(xintercept=0.93, linetype="dashed", color = "red")	
}



# fix xlab and ylab
#p=p+ggplot2::labs(x = "EV marker", y = "Non-EV marker") # xlab("EV marker") + ylab("Non-EV marker")
# stat
print("EV marker:")
DF=data.frame(values=DFref$EVmarker,groups=DFref$source)
print(EVqualityMS::statFunc(DF))

print("Non-EV markers:")
DF=data.frame(values=DFref$Contaminants,groups=DFref$source)
print(EVqualityMS::statFunc(DF))

return(p) # invisible(res)
}



#' Plot heatmap with subpopulation  markers from  Kowal et al 2016 paper
#' @param DF (DF - data frame with a GN column)
#' @param Fac (Fac - should have length ncol(DF)-1)
#' @param color (color - vector with ten colors for type of marker)
#' @param CatColL (color - protein category)
#' @param IncludeAlbumine (bool )
#' @return list
#' @keywords Quality_control
#' @examples
#' ---
#' @export
heatmapEVsub <- function(DF,Fac=NULL,color = NULL,CatColL=NULL,IncludeAlbumine = TRUE){
#
#
if (is.null(color)){color=c("grey60","grey30","gray15","steelblue3","navy","seagreen1","seagreen4","palegreen","violetred4","red")}
if (is.null(CatColL)){CatColL=c("#F0A3FF", "#0075DC", "#993F00", "#4C005C", "#191919", "#005C31", 
"#2BCE48", "#FFCC99", "#808080", "#94FFB5", "#8F7C00", "#9DCC00", 
"#C20088", "#003380", "#FFA405")}

if (!is.data.frame(DF)){stop("x must be a data frame with a column named GN")}
#
x=DF
# median of sample groups
if (!is.null(Fac)){
idx=which(names(DF)=="GN")
mat=DF[,-idx]
FUN = function(x) {median(as.numeric(x), na.rm = T)}
mat=t(apply(mat, 1, function(x) {tapply(x, Fac, FUN = FUN)}))
x=as.data.frame(mat)
x$GN=DF[,idx]
print(x)
}
# rank scale data

for (i in 1:ncol(x)){
if (names(x)[i]!="GN"){
xL=as.numeric(x[,i])
xL[is.na(xL)]=0
xL=scales::rescale(xL,to=c(0,1))
# Rank with ties receiving the same ranking
xL <- rank(xL, ties.method = "average")
# Scale ranks from 0 to 1
xL=xL-min(xL)
xL=xL/max(xL)
x[,i]=xL
}
}
# remove duplicated GN
idx=which(names(x) %in% "GN")
Rsum=rowSums(as.matrix(x[,-idx]))
idxL=order(Rsum,decreasing =TRUE)
#Rsum[idxL]
x=x[idxL,]
x=x[!duplicated(x$GN),]

#str(as.matrix(x[,-idx]))

# merge with markers data
# DF=EVqualityMS::ExoData[["Expression"]]
#which(DFanno$GN=="ENO1")

DF=data.frame(GN=EVqualityMS::DFsubEV$GN)
DF= dplyr::left_join(DF, x, by = "GN")

# remove NAs
DF[is.na(DF)]=0
mat=DF[,2:ncol(DF)]
rownames(mat)=DF$GN
# remove zero rows
idxL=which(rowSums(mat)==0)
if (length(idxL)>0){
mat=mat[-idxL,]	
}
# IncludeAlbumine
if (IncludeAlbumine==FALSE){mat=mat[rownames(mat)!="ALB",]}
# plot
colL=rep("black",nrow(mat))
# Contaminants
#conL=append(c("APOA1","APOA2","APOB","UMOD"),EVqualityMS::ExoData$Contaminants)
#for (i in 1:length(conL)){
#id2=which(conL[i]==rownames(mat))
#colL[id2]="red"
#}
# TopMarker
#for (i in 1:length(EVqualityMS::ExoData$TopMarker)){
#id2=which(EVqualityMS::ExoData$TopMarker[i]==rownames(mat))
#colL[id2]="green"
#}
# microvesicles
#for (i in 1:length(EVqualityMS::ExoData$lEV)){
#id2=which(EVqualityMS::ExoData$lEV[i]==rownames(mat))
#colL[id2]="blue"
#}
# cluster
clustL=rep("",nrow(mat))

for (i in 1:nrow(mat)){
#print(rownames(mat)[i])
id=which(rownames(mat)[i]==EVqualityMS::DFsubEV$GN)
#print(id)
if (id>=0){
clustL[i]=EVqualityMS::DFsubEV$Group[id]
} # if id>0
} # i
#
#print(clustL)
percentage_coverages=calculate_coverage_percentages(clustL,EVqualityMS::DFsubEV$Group)
cat("Percentage marker coverages:\n")
for (term in names(percentage_coverages)) {
  cat(term, ": ", percentage_coverages[[term]], "%\n")
}
#

Ncol=length(unique(clustL))
CatColL=CatColL[1:Ncol]
names(CatColL)=unique(clustL)
#print(CatColL)
p1 = ComplexHeatmap::Heatmap(t(mat), name = "mat", column_names_gp = ggfun::gpar(col = colL),cluster_columns = ComplexHeatmap::cluster_within_group(t(mat),as.numeric(factor(clustL))),heatmap_legend_param = list(position = "top", title = "Scaled expression", legend_direction = "horizontal"),top_annotation = ComplexHeatmap::HeatmapAnnotation(Category = clustL, col = list(Category = CatColL)))
p1 = ComplexHeatmap::draw(p1, heatmap_legend_side = "bottom")

# print Q index 
idxCL=which(x$GN %in% EVqualityMS::ExoData[["Contaminants"]])
idxEL=which(x$GN %in% EVqualityMS::ExoData[["TopMarker"]])
if (length(idxCL)>0){
qL=colSums(x[idxEL,-idx,drop=FALSE])/length(EVqualityMS::ExoData[["TopMarker"]])-colSums(x[idxCL,-idx,drop=FALSE])/length(EVqualityMS::ExoData[["Contaminants"]])
} else {
qL=colSums(x[idxEL,-idx,drop=FALSE])/length(EVqualityMS::ExoData[["TopMarker"]])
}
if (length(idxEL)==0){
qL=0
} 
cat("Overall quality metrics:\n")
print(qL)
# calculate contamination and exosome markers
if (length(idxCL)>0){
qC=colSums(x[idxCL,-idx,drop=FALSE])/length(EVqualityMS::ExoData[["Contaminants"]])
} else {qC=0}
if (length(idxEL)>0){
qE=colSums(x[idxEL,-idx,drop=FALSE])/length(EVqualityMS::ExoData[["TopMarker"]])
} else {qE=0}
# collect results
res=list()
res[["Expression"]]=mat
res[["Quality"]]=qL 
res[["Contaminants"]]=qC
res[["EV marker"]]=qE
res[["Heatmap"]]=p1
invisible(res)
}



#' Plot heatmap with subpopulation  markers from cd9, cd63 and CD81 pulldown from Kowal et al 2016 paper
#' @param DF (DF - data frame with a GN column)
#' @param Fac (Fac - should have length ncol(DF)-1)
#' @param color (color - vector with ten colors for type of marker)
#' @param CatColL (color - protein category)
#' @param IncludeAlbumine (bool )
#' @return list
#' @keywords Quality_control
#' @examples
#' ---
#' @export
heatmapEVsubCD <- function(DF,Fac=NULL,color = NULL,CatColL=NULL,IncludeAlbumine = TRUE){
#
#
if (is.null(color)){color=c("grey60","grey30","gray15","steelblue3","navy","seagreen1","seagreen4","palegreen","violetred4","red")}
if (is.null(CatColL)){CatColL=c("#F0A3FF", "#0075DC", "#993F00", "#4C005C", "#191919", "#005C31", 
"#2BCE48", "#FFCC99", "#808080", "#94FFB5", "#8F7C00", "#9DCC00", 
"#C20088", "#003380", "#FFA405")}

if (!is.data.frame(DF)){stop("x must be a data frame with a column named GN")}
#
x=DF
# median of sample groups
if (!is.null(Fac)){
idx=which(names(DF)=="GN")
mat=DF[,-idx]
FUN = function(x) {median(as.numeric(x), na.rm = T)}
mat=t(apply(mat, 1, function(x) {tapply(x, Fac, FUN = FUN)}))
x=as.data.frame(mat)
x$GN=DF[,idx]
print(x)
}
# rank scale data

for (i in 1:ncol(x)){
if (names(x)[i]!="GN"){
xL=as.numeric(x[,i])
xL[is.na(xL)]=0
xL=scales::rescale(xL,to=c(0,1))
# Rank with ties receiving the same ranking
xL <- rank(xL, ties.method = "average")
# Scale ranks from 0 to 1
xL=xL-min(xL)
xL=xL/max(xL)
x[,i]=xL
}
}
# remove duplicated GN
idx=which(names(x) %in% "GN")
Rsum=rowSums(as.matrix(x[,-idx]))
idxL=order(Rsum,decreasing =TRUE)
#Rsum[idxL]
x=x[idxL,]
x=x[!duplicated(x$GN),]

#str(as.matrix(x[,-idx]))

# merge with markers data
# DF=EVqualityMS::ExoData[["Expression"]]
#which(DFanno$GN=="ENO1")

DF=data.frame(GN=EVqualityMS::DFsubCD$GN)
DF= dplyr::left_join(DF, x, by = "GN")

# remove NAs
DF[is.na(DF)]=0
mat=DF[,2:ncol(DF)]
rownames(mat)=DF$GN
# remove zero rows
idxL=which(rowSums(mat)==0)
if (length(idxL)>0){
mat=mat[-idxL,]	
}
# IncludeAlbumine
if (IncludeAlbumine==FALSE){mat=mat[rownames(mat)!="ALB",]}
# plot
colL=rep("black",nrow(mat))
# Contaminants
#conL=append(c("APOA1","APOA2","APOB","UMOD"),EVqualityMS::ExoData$Contaminants)
#for (i in 1:length(conL)){
#id2=which(conL[i]==rownames(mat))
#colL[id2]="red"
#}
# TopMarker
#for (i in 1:length(EVqualityMS::ExoData$TopMarker)){
#id2=which(EVqualityMS::ExoData$TopMarker[i]==rownames(mat))
#colL[id2]="green"
#}
# microvesicles
#for (i in 1:length(EVqualityMS::ExoData$lEV)){
#id2=which(EVqualityMS::ExoData$lEV[i]==rownames(mat))
#colL[id2]="blue"
#}
# cluster
clustL=rep("",nrow(mat))

for (i in 1:nrow(mat)){
#print(rownames(mat)[i])
id=which(rownames(mat)[i]==EVqualityMS::DFsubCD$GN)
#print(id)
if (id>=0){
clustL[i]=EVqualityMS::DFsubCD$Group[id]
} # if id>0
} # i
#
#print(clustL)
# print percentage of markers identified
#print(table(clustL))
#print(table(EVqualityMS::DFsubCD$Group))
# Print the coverage percentages for each term in vector1
percentage_coverages=calculate_coverage_percentages(clustL,EVqualityMS::DFsubCD$Group)
cat("Percentage marker coverages: \n")
for (term in names(percentage_coverages)) {
  cat(term, ": ", percentage_coverages[[term]], "%\n")
}

#
Ncol=length(unique(clustL))
CatColL=CatColL[1:Ncol]
names(CatColL)=unique(clustL)
#print(CatColL)
p1 = ComplexHeatmap::Heatmap(t(mat), name = "mat", column_names_gp = ggfun::gpar(col = colL),cluster_columns = ComplexHeatmap::cluster_within_group(t(mat),as.numeric(factor(clustL))),heatmap_legend_param = list(position = "top", title = "Scaled expression", legend_direction = "horizontal"),top_annotation = ComplexHeatmap::HeatmapAnnotation(Category = clustL, col = list(Category = CatColL)))
p1 = ComplexHeatmap::draw(p1, heatmap_legend_side = "bottom")

# print Q index 
idxCL=which(x$GN %in% EVqualityMS::ExoData[["Contaminants"]])
idxEL=which(x$GN %in% EVqualityMS::ExoData[["TopMarker"]])
if (length(idxCL)>0){
qL=colSums(x[idxEL,-idx,drop=FALSE])/length(EVqualityMS::ExoData[["TopMarker"]])-colSums(x[idxCL,-idx,drop=FALSE])/length(EVqualityMS::ExoData[["Contaminants"]])
} else {
qL=colSums(x[idxEL,-idx,drop=FALSE])/length(EVqualityMS::ExoData[["TopMarker"]])
}
if (length(idxEL)==0){
qL=0
} 
cat("Overall quality metrics:\n")
print(qL)
# calculate contamination and exosome markers
if (length(idxCL)>0){
qC=colSums(x[idxCL,-idx,drop=FALSE])/length(EVqualityMS::ExoData[["Contaminants"]])
} else {qC=0}
if (length(idxEL)>0){
qE=colSums(x[idxEL,-idx,drop=FALSE])/length(EVqualityMS::ExoData[["TopMarker"]])
} else {qE=0}
# collect results
res=list()
res[["Expression"]]=mat
res[["Quality"]]=qL 
res[["Contaminants"]]=qC
res[["EV marker"]]=qE
res[["Heatmap"]]=p1
invisible(res)
}


#' Import MaxQuant protein group file
#' @param filename_txt (string)
#' @param data_s (string)
#' @return DF
#' @keywords import
#' @examples
#' ---
#' @export
ImportMaxQuantProteinGroups <- function(filename_txt,data_s="iBAQ"){
#
#
DF=DFstatRM::loadFromFile(filename_txt, sep = "\t", fill = F)
# remove contaminants and revers
idxL=grep("^CON_",DF$Protein.IDs)
DF=DF[-idxL,]
idxL=grep("^REV_",DF$Protein.IDs)
DF=DF[-idxL,]
#import
if (any(grepl(data_s,names(DF)))){
print("Extracting iBAQ values")
idxL=grep(data_s,names(DF))	
DFt=DF[,idxL[3:length(idxL)]]
names(DFt)=sub("iBAQ."," ",names(DFt))
DFt$GN=ifelse(grepl(".* GN=([^ ]+).*", DF[,c("Fasta.headers")]), gsub(".* GN=([^ ]+).*", "\\1", DF[,c("Fasta.headers")]), "")
DF=DFt
} else {
# grep lfq
if (!any(grepl("LFQ.intensity.",names(DF)))){stop("No LFQ values in result file")}
print("Extracting LFQ.intensity.")
idxL=grep("LFQ.intensity.",names(DF))
#names(DF)[idxL]	
DFt=DF[,idxL]
names(DFt)=sub("LFQ.intensity."," ",names(DFt))
DFt$GN=ifelse(grepl(".* GN=([^ ]+).*", DF[,c("Fasta.headers")]), gsub(".* GN=([^ ]+).*", "\\1", DF[,c("Fasta.headers")]), "")
DF=DFt
}

return(DF) # invisible(DF)
}

