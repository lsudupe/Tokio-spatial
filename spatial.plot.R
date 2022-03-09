p1 <- SpatialFeaturePlot(fibro, features = "geneSetD2", pt.size.factor = 10,alpha = 0.6, combine = FALSE)
fix.p1 <- scale_fill_continuous(limits = c(0,0.6),type = "viridis")
p2 <- lapply(p1, function (x) x + fix.p1)

pdf("./ssGSEA/prueba.pdf")
CombinePlots(p2)
dev.off()

p1 <- SpatialFeaturePlot(integrated, features = "Prelp", pt.size.factor = 10, combine = FALSE)
fix.p1 <- scale_fill_continuous(limits = c(-5.0,20.0), 
                                breaks = c(-5.0,20.0),
                                type = "viridis")
p2 <- lapply(p1, function (x) x + fix.p1)

pdf("Prelp.pdf")
CombinePlots(p2)
dev.off()

pdf("aa.pdf")
SpatialPlot(object = integrated, 
            features = c("Cthrc1"), 
            pt.size.factor = 10,
            min.cutoff = -10,
            max.cutoff = 20,
            combine = FALSE) 
dev.off()


p1 <- SpatialFeaturePlot(integrated, features = "Cthrc1", pt.size.factor = 10, alpha= 0.6,combine = FALSE)
fix.p1 <- scale_fill_gradientn(limits = c(-5.0,20.0),
                                breaks=c("Min","Max"),
                                labels=c("Min","Max"),
                               colours=topo.colors(7))
p2 <- lapply(p1, function (x) x + fix.p1)

pdf("f.pdf")
CombinePlots(p2)
dev.off()

pdf("integrated.ident.pdf")
SpatialPlot(integrated, group.by = c("seurat_clusters"), pt.size.factor = 10,label = TRUE, combine = FALSE)
dev.off()



bl <- colorRampPalette(c("navy","royalblue","lightskyblue"))(200)                      
re <- colorRampPalette(c("mistyrose", "red2","darkred"))(200)

p1 <- SpatialFeaturePlot(fibro, features = "d1.d2", pt.size.factor = 10, combine = FALSE)
fix.p1 <- scale_fill_gradientn(colours=c(bl,"white", re), na.value = "grey98",limits = c(-0.5,0.5))
p2 <- lapply(p1, function (x) x + fix.p1)

pdf("./ssGSEA/d1.d2.fibro.minu0.5.0.5.pdf")
CombinePlots(p2)
dev.off()

a <- c("Col14a1", "Nox4")


####select only the genes in comun, for the loop ploting
a <- integrated@assays[["integrated"]]@var.features
genesD1 <-  base::intersect(a,D1)

for (i in a){
p1 <- SpatialFeaturePlot(integrated, features = i, pt.size.factor = 10, alpha = 0.6, combine = FALSE)
fix.p1 <- scale_fill_continuous(type = "viridis")
p2 <- lapply(p1, function (x) x + fix.p1)

pdf(paste(i,".pdf",sep=""))
print(CombinePlots(p2))
dev.off()
}

feature.list <- c("geneSetB", "geneSetD1", "geneSetD2", "relation_log(d1.d2)")

for (i in feature.list){
  p1 <- SpatialFeaturePlot(fibro, features = i, pt.size.factor = 10, combine = FALSE)
  fix.p1 <- scale_fill_gradientn(colours=c(bl,"white", re), na.value = "grey98",limits = c(-1,1))
  p2 <- lapply(p1, function (x) x + fix.p1)
  
  pdf(paste(i,".pdf",sep=""))
  print(CombinePlots(p2))
  dev.off()
}

