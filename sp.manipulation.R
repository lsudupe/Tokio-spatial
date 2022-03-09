### SCRIPT: Check and manipulate spatial data

## 02.12.21 Laura Sudupe , git @lsudupe

library("Seurat")
library("SPOTlight")
library("base")
library("dplyr")
library("gt")
library("DT")
library("ggplot2")
library("data.table")

##Read the data
spatial <- readRDS("./ibex/deconvolution/combined.sp.sct.rds")
integrated <- readRDS("./ibex/integration/integrated.sct.rds")
saveRDS(integrated, "./ibex/integration/integrated.sct.rds")

DimPlot(spatial, group.by = c("seurat_clusters"), label = T) + ggtitle("UMAP_r0.4")

pdf(file.path("./",filename = "INTEGRATED.pdf"))
SpatialDimPlot(object = integrated ,group.by = c("seurat_clusters"), pt.size.factor = 30)
dev.off()


########SPLIT OBJECT
integrated.list <- SplitObject(integrated, split.by = "orig.ident")

sham <- integrated.list[["sham"]]
saveRDS(sham, "./ibex/integration/integrated_sham.sct.rds")
sham <- readRDS("./ibex/integration/integrated_sham.sct.rds")

day1 <- integrated.list[["day 1 post mi"]]
saveRDS(day1, "./ibex/integration/integrated_day1.sct.rds")
day1 <- readRDS("./ibex/integration/integrated_day1.sct.rds")

day7 <- integrated.list[["day7 post mi"]]
saveRDS(day7, "./ibex/integration/integrated_day7.sct.rds")
day7 <- readRDS("./ibex/integration/integrated_day7.sct.rds")

day14 <- integrated.list[["day14 post mi"]]
saveRDS(day14, "./ibex/integration/integrated_day14.sct.rds")
day14 <- readRDS("./ibex/integration/integrated_day14.sct.rds")

####SUBSET THE SPATIAL DATA

integrated <- readRDS("./ibex/integration/integrated.sct.rds")

integrated@meta.data[["ident"]] <- integrated@meta.data[["seurat_clusters"]]
Seurat::Idents(object = integrated) <- integrated@meta.data[["ident"]]
integrated <- subset(x = integrated, idents = c("3","4","5","7","8","9","10"))

integrated@meta.data[["ident"]] <- integrated@active.ident

saveRDS(integrated.fibro, file = "./ibex/integration/integrated.sp.sct.fibro.new.rds")
integrated.fibro <- readRDS("./ibex/integration/integrated.sp.sct.fibro.new.rds")
integrated.fibro$orig.ident[integrated.fibro$orig.ident == "day14 post mi"] <- "4.day 14"

########### SAVE FIBRO WITH OUT SHAM
meta <- integrated.fibro@meta.data
#take our sham data
meta <- meta[!grepl("sham", meta[,1]),]
integrated.fibro@meta.data <- meta
#Idents(integrated.fibro) <- "orig.ident"
#take out image
integrated.fibro@images[["sham.slice"]] <- NULL
#save
saveRDS(integrated.fibro, file = "./ibex/integration/integrated.sp.sct.fibro.nosham.rds")
#################

##
integrated@meta.data[["ident"]] <- integrated@meta.data[["seurat_clusters"]]
Seurat::Idents(object = integrated) <- integrated@meta.data[["ident"]]
integrated <- subset(x = integrated, idents = c("0","1","2","6"))

integrated@meta.data[["ident"]] <- integrated@active.ident

saveRDS(no.fibro, file = "./ibex/integration/integrated.sp.sct.nofibro.new.rds")
no.fibro <- readRDS("./ibex/integration/integrated.sp.sct.nofibro.new.rds")
no.fibro$orig.ident[no.fibro$orig.ident == "day14 post mi"] <- "1.day 14"
##
day1 <- readRDS("./ibex/integration/integrated_day1.sct.rds")
SpatialPlot(object = day1, 
            features = c("seurat_clusters"), 
            pt.size.factor = 15,
            min.cutoff = 0,
            max.cutoff = 100,
            repel = TRUE)
SpatialPlot(object = day1, group.by= c("seurat_clusters"), pt.size.factor = 15, ncol = 1)

day1@meta.data[["ident"]] <- day1@meta.data[["seurat_clusters"]]
Seurat::Idents(object = day1) <- day1@meta.data[["ident"]]
day1 <- subset(x = day1, idents = c("3","4","5","7"))

day1@meta.data[["ident"]] <- day1@active.ident

saveRDS(day1, file = "./ibex/integration/day1.sp.sct.fibro.rds")
day1 <- readRDS("./ibex/integration/day1.sp.sct.fibro.rds")

#############

SpatialPlot(object = day7.prueba, 
            features = c("seurat_clusters"), 
            pt.size.factor = 15,
            min.cutoff = 0,
            max.cutoff = 100,
            repel = TRUE)
SpatialPlot(object = day7, group.by= c("seurat_clusters"), pt.size.factor = 15, ncol = 1)

day7@meta.data[["ident"]] <- day7@meta.data[["seurat_clusters"]]
Seurat::Idents(object = day7) <- day7@meta.data[["ident"]]
day7 <- subset(x = day7, idents = c("3","4","5","7"))

day7@meta.data[["ident"]] <- day7@active.ident

saveRDS(day7, file = "./ibex/integration/day7.sp.sct.fibro.rds")
day7 <- readRDS("./ibex/integration/day7.sp.sct.fibro.rds")

#############

SpatialPlot(object = day14, group.by= c("seurat_clusters"), pt.size.factor = 15, ncol = 1)

day14@meta.data[["ident"]] <- day14@meta.data[["seurat_clusters"]]
Seurat::Idents(object = day14) <- day14@meta.data[["ident"]]
day14 <- subset(x = day14, idents = c("3","4","5","7"))

day14@meta.data[["ident"]] <- day14@active.ident

saveRDS(day14, file = "./ibex/integration/day14.sp.sct.fibro.rds")
day14 <- readRDS("./ibex/integration/day14.sp.sct.fibro.rds")

##############################BvsRest

BvsRestmarkers <- readRDS("./ibex/deconvolution/results/IRCF_BvsRest_markers.rds")
BvsRest<- subset(BvsRestmarkers, p_val_adj < 0.05 & 0.5 < avg_log2FC)
top5.BvsRest <- BvsRest %>%
  group_by(cluster) %>%
  top_n(n = 5,
        wt = avg_log2FC)

BvsD1vsD2markers <- readRDS("./ibex/deconvolution/results/IRCF_BvsD1vsD2_markers.rds")
BvsD1vsD2<- subset(BvsD1vsD2markers, p_val_adj < 0.05 & 0.5 < avg_log2FC)
top5.BvsD1vsD2 <- BvsD1vsD2 %>%
  group_by(cluster) %>%
  top_n(n = 5,
        wt = avg_log2FC)


############download feature, barcode, matrix for Jupyter notebook
library(DropletUtils)
library(rhdf5)

write10xCounts(path = "./filtered_feature_bc_matrix.h5", type = c("HDF5"), overwrite = T,version = "3",
               x = GetAssayData(day1,slot = "data",assay = "integrated"), genome = "mm10", 
               barcodes = colnames(GetAssayData(day1,slot = "data",assay = "integrated")), 
               gene.id = rownames(GetAssayData(day1,slot = "data",assay = "integrated")))

library(SeuratDisk)
library(SeuratData)
SaveH5Seurat(day1, overwrite = TRUE)
Convert("SeuratProject.h5Seurat", dest = "h5ad")

######################SAVE COORDINATES

day1.coor <- day1@images[["day1_slice.1"]]@coordinates
day1.coor <- cbind(feature = rownames(day1.coor), day1.coor)
rownames(day1.coor) <- 1:nrow(day1.coor)

a <- data.frame(do.call('rbind', strsplit(as.character(day1.coor$feature),'_',fixed=TRUE)))
day1.coor$feature <- a$X2

write.csv(day1.coor,"./day1.coordi.csv",row.names = FALSE, col.names= FALSE)
