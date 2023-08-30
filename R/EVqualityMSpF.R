# _"/media/rune/Jochime2/RscripterCool/Rprojects/rmpackges/EVqualityMSp.R"



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
xL=scales::rescale(xL,to=c(0,1))
xL[is.na(xL)]=0
xL=rank(xL)/max(rank(xL))
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
DF=DF %>% dplyr::left_join( x, by = "GN")

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
for (i in it(EVqualityMS::ExoData$Contaminants)){
id2=which(EVqualityMS::ExoData$Contaminants[i]==rownames(mat))
colL[id2]="red"
}
# TopMarker
for (i in it(EVqualityMS::ExoData$TopMarker)){
id2=which(EVqualityMS::ExoData$TopMarker[i]==rownames(mat))
colL[id2]="green"
}
# 
for (i in it(EVqualityMS::ExoData$lEV)){
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
heatmapEVqualityISEV <- function(DF,Fac=NULL,color = NULL,CatColL=NULL,IncludeAlbumine = TRUE,plotCatDescription=TRUE){
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
FUN = function(x) {median(x, na.rm = T)}
mat=t(apply(mat, 1, function(x) {tapply(x, Fac, FUN = FUN)}))
x=as.data.frame(mat)
x$GN=DF[,idx]
}
# rank scale data

for (i in 1:ncol(x)){
if (names(x)[i]!="GN"){
xL=as.numeric(x[,i])
xL=scales::rescale(xL,to=c(0,1))
xL[is.na(xL)]=0
xL=rank(xL)/max(rank(xL))
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

DF=data.frame(GN=EVqualityMS::DFanno$GN)
DF=DF %>% dplyr::left_join( x, by = "GN")

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
conL=append(c("APOA1","APOA2","APOB","UMOD"),EVqualityMS::ExoData$Contaminants)
for (i in it(conL)){
id2=which(conL[i]==rownames(mat))
colL[id2]="red"
}
# TopMarker
for (i in it(EVqualityMS::ExoData$TopMarker)){
id2=which(EVqualityMS::ExoData$TopMarker[i]==rownames(mat))
colL[id2]="green"
}
# microvesicles
for (i in it(EVqualityMS::ExoData$lEV)){
id2=which(EVqualityMS::ExoData$lEV[i]==rownames(mat))
colL[id2]="blue"
}
# cluster
clustL=rep("",nrow(mat))

for (i in it(mat)){
print(rownames(mat)[i])
id=which(rownames(mat)[i]==EVqualityMS::DFanno$GN)
if (id>=0){

if (plotCatDescription==TRUE){
clustL[i]=EVqualityMS::DFanno$CatDes[id]
} else {
clustL[i]=EVqualityMS::DFanno$Cat[id]
}
}
}
Ncol=clustL %>% unique %>% length
CatColL=CatColL[1:Ncol]
names(CatColL)=unique(clustL)
print(CatColL)
p1 = ComplexHeatmap::Heatmap(t(mat), name = "mat", column_names_gp = ggfun::gpar(col = colL),cluster_columns = ComplexHeatmap::cluster_within_group(t(mat),as.numeric(factor(clustL))),heatmap_legend_param = list(position = "top", title = "Regulation", legend_direction = "horizontal"),top_annotation = ComplexHeatmap::HeatmapAnnotation(Category = clustL, col = list(Category = CatColL)))
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


#' plot scatter of relative contamination and EV marker abundance
#' @param x (DF- quantitative values in columns and last column holds the gene names - GN)
#' @param Name (string - name of the new samples)
#' @param src (string - which reference samples to include)
#' @param palette (colorL or palette)
#' @return ggplot
#' @keywords Quality_control
#' @examples
#' ---
#' @export
scatterPlot <- function(x=NULL,Name="new",src="all",palette="jco"){
#
#
#x = DFurine; Name = "Urine"; src = "plasma"
if (!is.null(x)){
res=EVqualityMS::HeatmapEVmarkers(x,plotHeat =FALSE)
DF=data.frame(EVmarker=res$"EV marker",Contaminants=res$Contaminants,source=Name,sample="")
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
## sp <- ggpubr::ggscatter(DFref, x =  "EVmarker", y ="Contaminants",color = "source", palette = palette,size = 3, alpha = 0.6)+ggpubr::border()
# Marginal density plot of x (top panel) and y (right panel)
## xplot <- ggpubr::ggdensity(DFref, x="EVmarker", fill = "source",palette = palette)
## xplot <- xplot + ggpubr::clean_theme() + ggpubr::rremove("legend")
# Arranging the plot using cowplot
## p <- cowplot::plot_grid(xplot, sp,ncol = 1, align = "hv", rel_heights = c(1, 1.5,1))
p=ggpubr::ggscatterhist(
  DFref, x = "EVmarker", y = "Contaminants",position=ggplot2::position_jitter(h=0.05,w=0.05),
  color = "source", size = "source",shape ="source", alpha = 0.6,font.tickslab = c(20, "plain"),font.legend = c(16, "plain"),font.x=c(16, "plain"),font.y=c(16, "plain"),
  palette = palette,
  margin.params = list(fill = "source", color = "black", size = 0.2)
  ) 

return(p) # invisible(res)
}

