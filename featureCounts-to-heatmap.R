setwd("~/Dropbox/Manuscripts/onyvax/R/") # OSX computer



# ######################################################################
#             LOAD sample counts (featureCounts)
# ######################################################################
matching <- read.delim("./samplelist.txt", header=F, comment.char = "#")

result = list()
for (i in 1:nrow(matching)) {
  
  temp <- read.delim(paste("./data/featureCounts/", matching[i,1], ".count", sep=""), comment.char = "#", header=T)

  rownames(temp) <- temp[,1]
  
  rownames(temp)
  
  temp <- temp[,c(1,7)]
  
  result[[as.character(matching[i,1])]] <- temp
}

# which gene names are common?
gene.common <- table(unlist(lapply(result, rownames)))
nosamples <- length(result)
gene.common <- names(gene.common)[gene.common==nosamples]  # change to # of samples!  

# only keep common genes in the list results (really a dfList)
result.common <- lapply(result, function(x) x[gene.common,])

# create the final table with counts 
count.table <- do.call(cbind.data.frame, lapply(result.common, function(x) x[,2])) # ,2 is count
rownames(count.table) <- gene.common

# ######################################################################
#             NORMALISE AND FILTER
# ######################################################################

# Normalization time!
library("edgeR")

tmm <- calcNormFactors(count.table, method = "TMM")
lib.size <- colSums(count.table)/10^6
eff.lib.size <- lib.size * tmm
count.table.scale <- t(t(count.table) / eff.lib.size)

min.count = apply(count.table.scale, 1, function(x) sum(x<3))

# get rid of rows with only zeros...
count.table.scale.no0 <- count.table.scale[ rowSums(count.table.scale)!=0, ] 
count.table.scale <- count.table.scale.no0 

# replace 0s with NA
count.table.scale.nona <- count.table.scale
count.table.scale.nona[count.table.scale.nona==0] = NA
count.table.scale <- count.table.scale.nona

length(count.table.scale)/nosamples # 17,430 genes 
counts <- count.table.scale

# 
# rename samples
counts <- as.data.frame(counts)
colnames(counts)
class(counts)
names(counts)[names(counts)=="11_2004"] <- "AQ0411"
names(counts)[names(counts)=="15_2004"] <- "AQ0415"
names(counts)[names(counts)=="20_2004"] <- "AQ0420"
names(counts)[names(counts)=="96_2003"] <- "AQ0396"
# save outout to an R object
save(counts,file = "./counts.Robj") 
# load("./counts.Robj")



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# HEAT MAP
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
library(gplots)
library(RColorBrewer)

# load normalised gene counts
load("./counts.Robj")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ONLY DRAW PARTICULAR GENES
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# tempinput <- c("NFKBIE","ATG16L2","RGS19","ADRBK1")

# custom gene list #1
tempinput <- readLines("KEGG_map05215_prostate_cancer.txt")
tempinput <- sort(tempinput) # sort if for convenience
# give your list a name for later
pdffn = "KEGG_map05215_prostate_cancer-heatmap.pdf"
# convert to df
countsDF <- as.data.frame(counts)
countsDF["AR",] # very high AR in LNCaP, very little in PC3 (but there), a small amount in BPH1 and 100X less than LNCaP in 96_2003

value.top = countsDF[tempinput,] # top values
head(value.top)
value.top["AR",]

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# HEAT MAP WITH ALL GENES
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
pdffn = "heatmap_onyvax-17430-genes-heatmap3.pdf"
load("./counts.Robj")
countsDF <- as.data.frame(counts)
# we want all genes!
tempinput <- as.character(row.names(countsDF))
value.top = countsDF[tempinput,] # top values for 17,300 genes



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# DRAW HEATMAP
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
value.topLOG = as.matrix(log2(value.top))
value.topLOG[is.na(value.topLOG)] <- 0 # replace NA with 0 or heatmap will choke
value.topLOG["AR",]
# -log10 if I want to draw the target on RHS...
class(value.topLOG) # matrix
library(gplots)
library(RColorBrewer)
cols <- rev(colorRampPalette(brewer.pal(10,"RdBu"))(256))
library(devtools)
source("heatmap.3.R")
#
any(is.na(value.topLOG)) # check for any pesky  NAs
value.topLOG <- na.omit(value.topLOG)
value.topLOG["AR",]
#
value.topLOG <- as.matrix(value.topLOG[!grepl('NA', rownames(value.topLOG)), ])
# remove NA rownames (sometimes introduced if you are using a different OS at home etc...)
#
#Define custom dist and hclust functions for use with heatmaps
mydist=function(c) {dist(c,method="euclidian")}
myclust=function(c) {hclust(c,method="average")}
pdf(file=pdffn)
#@ cols <- colorRampPalette(c("green", "black", "red"))(n = 256) # classical
#@ cols=rev(redgreen(256))
h <- heatmap.3(scale(value.topLOG), col=cols, scale="column", trace="none",density.info="none", dendrogram="column", key=TRUE,cexRow=0.5 , hclustfun=myclust, distfun=mydist, breaks=seq(-3,3,6/256) )
dev.off()
# , Colv=FALSE



