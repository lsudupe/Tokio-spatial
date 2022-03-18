### SCRIPT: addModuleScore() function from Seurat for enrichement analysis Tokio data in IBEX

## 10.03.22 Laura Sudupe , git @lsudupe

## https://www.waltermuskovic.com/2021/04/15/seurat-s-addmodulescore-function/

###### libraries
library("Seurat")
library("ggplot2")
library(data.table)
library("escape")
library("GSEABase")
library("RColorBrewer")
library("dittoSeq")
library("dichromat")

###### Download  data
fibro <- readRDS("../data/integrated.sp.sct.fibro.new.rds")
no.fibro <- readRDS("../data/integrated.sp.sct.nofibro.new.rds")

objects <- c(fibro, no.fibro)
names(objects) <- c("fibro", "no.fibro")

###### Read the gene sets
genes.velocity <- read.csv("../data/genes_velocity.csv", header = TRUE, sep = ",")
B <- read.table("../data/clustB_signature.txt", header = FALSE)
B <- as.vector(B$V1)
D1 <- as.vector(genes.velocity$Dynamics_1)
D2 <- as.vector(genes.velocity$Dynamics_2)
D2[D2 == ""] <- NA 
D2 <- D2[!is.na(D2)]

##create list
sig.list <- list(B, D1, D2)
names(sig.list) <- c("B", "D1", "D2")
saveRDS(sig.list, "../data/sig.list.rds")

##our gene list
sig.list <- readRDS("../data/sig.list.rds")
B <- sig.list[1]
D1 <- sig.list[2]
D2 <- sig.list[3]

###############################################333
for (i in 1:length(objects)){
  ###matrix
  a <- objects[[i]]
  a@assays[["integrated"]]@counts <- a@assays[["integrated"]]@data
  ###### addModuleScore score 
  a <- AddModuleScore(a, features = B, assay = "integrated", name = "geneSetB_")
  a <- AddModuleScore(a, features = D1, assay = "integrated", name = "geneSetD1_")
  a <- AddModuleScore(a, features = D2, assay = "integrated", name = "geneSetD2_")
  #extract module Scores
  d1 <- as.vector(a@meta.data[["geneSetD1_1"]])
  d2 <- as.vector(a@meta.data[["geneSetD2_1"]])
  #calculate relation d1/d2
  d1.d2 <- log10(d1/d2)
  #is.na(d1.d2) <-sapply(d1.d2, is.infinite)
  d1.d2[is.na(d1.d2)] = 0
  a@meta.data$relation_log.d1.d2 <- d1.d2
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
###### no present genes
#b.no.present <- c("1500015O10Rik", "Wisp2", "Ctgf", "Cyr61", "Cpe", "Id1", "Palld", "Cdc42ep3", "Nr4a2", "Lbh", 
#                  "Gng8", "Fam180a", "Tcea3", "Cyp26a1", "Wisp1", "Nox4", "2310022B05Rik", "Myl12a", "Fam19a5", 
#                  "Bsg", "Cldn11", "Grb14", "Nuak1", "Cdkn2b", "Ttll5", "Maff", "Vamp5", "Bdnf", "Stc2", "Arid5b", 
#                  "Slc44a2", "Mpzl1", "Hcfc1r1", "Luzp1", "Smad7", "1810058I24Rik", "Tyms", "Klf2", "Tmem160", 
#                  "Eif3f", "Nenf", "Dnajb4", "Dut", "Rgs3", "Siva1", "Vapa", "Bri3", "Fhl2", "Sgk1", "Ier5", "Fam198b", 
#                  "Gm26532", "Pam", "Btg")
#d1.no.presnt <- c("Nox4", "Fhl2", "Smad7", "Nr4a2", "Grb14", "Nuak1", "Stc2", "Bsg", "Vamp5", "Slc44a2", "Tcea3", 
#                 "Cpe", "Gng8")
#d2.no.present <- c("Fam198b", "Ctgf", "Cdkn2b", "Cyr61", "Dut", "Bdnf", "Wisp1", "Tyms", "Luzp1", "Palld")

###################################################PLOTS

fibro <- readRDS("./fibro.rds")
no.fibro <- readRDS("./no.fibro.rds")

objects <- c(fibro, no.fibro)
names(objects) <- c("fibro", "no.fibro")

clusters <- c("geneSetB_1", "geneSetD1_1", "geneSetD2_1", "relation_log.d1.d2")

##### UMAP Plot scores
##fibro
for (i in clusters){
  a <- FeaturePlot(fibro, features = "geneSetB_1", label = TRUE, repel = TRUE) +
    scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
  b <- FeaturePlot(fibro, features = "geneSetD1_1", label = TRUE, repel = TRUE) +
    scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
  c <- FeaturePlot(fibro, features = "geneSetD2_1", label = TRUE, repel = TRUE) +
    scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
  d <- FeaturePlot(fibro, features = "relation_log.d1.d2", label = TRUE, repel = TRUE) +
    scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
  
  pdf(paste(i,"umap.scores.fibro.pdf",sep=""))
  print(cowplot::plot_grid(a,NULL,b,NULL,c,NULL,d, nrow = 5, ncol = 1, rel_heights = c(4,0.5,3,0.5,3)))
  dev.off()
}

##no-fibro
for (i in clusters){
  a <- FeaturePlot(no.fibro, features = "geneSetB_1", label = TRUE, repel = TRUE) +
    scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
  b <- FeaturePlot(no.fibro, features = "geneSetD1_1", label = TRUE, repel = TRUE) +
    scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
  c <- FeaturePlot(no.fibro, features = "geneSetD2_1", label = TRUE, repel = TRUE) +
    scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
  d <- FeaturePlot(no.fibro, features = "relation_log.d1.d2", label = TRUE, repel = TRUE) +
    scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
  
  pdf(paste(i,"umap.scores.no.fibro.pdf",sep=""))
  print(cowplot::plot_grid(a,NULL,b,NULL,c,NULL,d, nrow = 5, ncol = 1, rel_heights = c(4,0.5,3,0.5,3)))
  dev.off()
}


###### UMAP plot, enrichment, orig.ident, clusters name relation d1.d2
##fibro
for (i in clusters){
  u1 <- FeaturePlot(fibro, features = i, label = TRUE, repel = TRUE) +
    scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
  u2 <- DimPlot(fibro, group.by = "orig.ident")
  u3 <- DimPlot(fibro, reduction= "umap", label = T, pt.size = 1)
  pdf(paste(i,"umap.fibro.pdf",sep=""))
  print(cowplot::plot_grid(u1,NULL,u2,NULL,u3,NULL, nrow = 5, ncol = 1, rel_heights = c(4,0.5,3,0.5,3)))
  dev.off()
}

##no.fibro
for (i in clusters){
  u1 <- FeaturePlot(no.fibro, features = i, label = TRUE, repel = TRUE) +
    scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
  u2 <- DimPlot(no.fibro, group.by = "orig.ident")
  u3 <- DimPlot(no.fibro, reduction= "umap", label = T, pt.size = 1)
  pdf(paste(i,"umap.d1.d2.relations.fibro.pdf",sep=""))
  print(cowplot::plot_grid(u1,NULL,u2,NULL,u3,NULL, nrow = 5, ncol = 1, rel_heights = c(4,0.5,3,0.5,3)))
  dev.off()
}


###### Density plots
##fibro
for (i in clusters){
  pdf(paste(i,"dens.fibro.pdf",sep=""))
  print(dittoRidgePlot(fibro, i, group.by = "orig.ident"))
  dev.off()
}

##no fibro
for (i in clusters){
  pdf(paste(i,"dens.nofibro.pdf",sep=""))
  print(dittoRidgePlot(no.fibro, i, group.by = "orig.ident"))
  dev.off()
}

###### Density plots in clusters

##fibro
for (i in clusters){
  pdf(paste(i,"clusterdens.fibro.pdf",sep=""))
  print(ridgeEnrichment(fibro@meta.data, gene.set = i, group = "ident", facet = "orig.ident", add.rug = TRUE))
  dev.off()
}

##no fibro
for (i in clusters){
  pdf(paste(i,"clusterdens.nofibro.pdf",sep=""))
  print(ridgeEnrichment(no.fibro@meta.data, gene.set = i, group = "ident", facet = "orig.ident", add.rug = TRUE))
  dev.off()
}


###### Spatial plots

bl <- colorRampPalette(c("navy","royalblue","lightskyblue"))(200)                      
re <- colorRampPalette(c("mistyrose", "red2","darkred"))(200)

for (i in clusters){
  p1 <- SpatialFeaturePlot(fibro, features = i, pt.size.factor = 10, combine = FALSE)
  fix.p1 <- scale_fill_gradientn(colours=c(bl,"white", re), na.value = "grey98",limits = c(-3,6))
  p2 <- lapply(p1, function (x) x + fix.p1)
  
  pdf(paste(i,".pdf",sep=""))
  print(CombinePlots(p2))
  dev.off()
}

######d1.d2

p1 <- SpatialFeaturePlot(fibro, features = "relation_log.d1.d2", pt.size.factor = 10, alpha= 0.8,combine = FALSE)
fix.p1 <- scale_fill_gradientn(colours=c(bl,"white", re), na.value = "grey98",limits = c(-2,2))
p2 <- lapply(p1, function (x) x + fix.p1)

pdf("d1.d2.spatial.pdf")
CombinePlots(p2)
dev.off()
