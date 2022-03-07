### SCRIPT: AUCell method enrichement analysis Tokio data in IBEX

## 03.03.22 Laura Sudupe , git @lsudupe

#https://bioconductor.org/packages/release/bioc/vignettes/AUCell/inst/doc/AUCell.html
#https://rdrr.io/bioc/AUCell/man/GeneSet-methods.html

###### libraries
library("Seurat")
library("GSEABase")
library("AUCell")
library("ggplot2")
library("escape")
library("dittoSeq")
library("dichromat")

##CHOOSE THE auxMaxRank
genes.porcentage <- 750 #porcentage depending number of genes in spatial


###### Read the gene sets
genes.velocity <- read.csv("../data/genes_velocity.csv", header = TRUE, sep = ",")
B <- read.table("../data/clustB_signature.txt", header = FALSE)
B <- as.vector(B$V1)
D1 <- as.vector(genes.velocity$Dynamics_1)
D2 <- as.vector(genes.velocity$Dynamics_2)
D2[D2 == ""] <- NA 
D2 <- D2[!is.na(D2)]

###### Create genesets
b.sig <- GeneSet(B, setName="geneSetB")
d1.sig <- GeneSet(D1, setName="geneSetD1")
d2.sig <- GeneSet(D2, setName="geneSetD2")
geneSets <- GeneSetCollection(b.sig, d1.sig, d2.sig)

###### Download  data
fibro <- readRDS("../data/integrated.sp.sct.fibro.rds")
no.fibro <- readRDS("../data/integrated.sp.sct.nofibro.rds")

objects <- c(fibro, no.fibro)
names(objects) <- c("fibro", "no.fibro")

###############################################333
for (i in 1:length(objects)){
          ###matrix
          a <- objects[[i]]
          a@assays[["integrated"]]@counts <- a@assays[["integrated"]]@data
          matrix <- a@assays[["integrated"]]@counts
          matrix <- as.matrix(matrix)
          ###### AUC score 
          cells_rankings <- AUCell_buildRankings(matrix, nCores=1)#, plotStats=TRUE)
          #rankings <- getRanking(cells_rankings)
          cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings, aucMaxRank=genes.porcentage)
          #extract AUC values
          auc_per_cell_all <- as.data.frame(t(getAUC(cells_AUC)))
          #calculate relation d1/d2
          d1 <- as.vector(auc_per_cell_all$geneSetD1)
          d2 <- as.vector(auc_per_cell_all$geneSetD2)
          d1.d2 <- log10(d1/d2)
          d1.d2 <-sapply(d1.d2, is.infinite)
          d1.d2[is.na(d1.d2)] = 0
          auc_per_cell_all$relation_log.d1.d2 <- d1.d2
          ##save meta
          a <- AddMetaData(a, auc_per_cell_all)
          ##save object
          saveRDS(a,file = paste0(names(objects[i]),".rds"))
          #feature_plot
          p1 <- SpatialFeaturePlot(a, features = "nCount_SCT", pt.size.factor = 10, combine = FALSE)
          fix.p1 <- scale_fill_continuous(limits = c(30000,40000), 
                                	      breaks = c(30000,40000),
                                          type = "viridis")
          p2 <- lapply(p1, function (x) x + fix.p1)

          pdf(paste(names(objects[i]),"count.dens.pdf",sep=""))
          print(CombinePlots(p2))
          dev.off()
}         

###################################################PLOTS

fibro <- readRDS("./fibro.rds")
no.fibro <- readRDS("./no.fibro.rds")

#Spatial plots
bl <- colorRampPalette(c("navy","royalblue","lightskyblue"))(200)                      
re <- colorRampPalette(c("mistyrose", "red2","darkred"))(200)

feature.list <- c("geneSetB", "geneSetD1","geneSetD2", "relation_log.d1.d2")

for (i in feature.list){
  p1 <- SpatialFeaturePlot(fibro, features = i, pt.size.factor = 10, combine = FALSE)
  fix.p1 <- scale_fill_gradientn(colours=c(bl,"white", re), na.value = "grey98",limits = c(-1,1))
  p2 <- lapply(p1, function (x) x + fix.p1)
  
  pdf(paste(i,"spatial.fibro.pdf",sep=""))
  print(CombinePlots(p2))
  dev.off()
}

#Density plots

dens <- c("geneSetB", "geneSetD1","geneSetD2", "relation_log.d1.d2", "nCount_SCT")

##fibro
for (i in dens){
  pdf(paste(i,"dens.fibro.pdf",sep=""))
  print(dittoRidgePlot(fibro, i, group.by = "orig.ident"))
  dev.off()
}

##no fibro
for (i in dens){
  pdf(paste(i,"dens.nofibro.pdf",sep=""))
  print(dittoRidgePlot(no.fibro, i, group.by = "orig.ident"))
  dev.off()
}

#Density plots in clusters

##fibro
for (i in dens){
  pdf(paste(i,"clusterdens.fibro.pdf",sep=""))
  print(ridgeEnrichment(fibro@meta.data, gene.set = i, group = "ident", facet = "orig.ident", add.rug = TRUE))
  dev.off()
}

##no fibro
for (i in dens){
  pdf(paste(i,"clusterdens.nofibro.pdf",sep=""))
  print(ridgeEnrichment(no.fibro@meta.data, gene.set = i, group = "ident", facet = "orig.ident", add.rug = TRUE))
  dev.off()
}






