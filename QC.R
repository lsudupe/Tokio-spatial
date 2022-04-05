## SCRIPT: QC of the seurat spatial objects TOKIO project

## 29.03.22 Laura Sudupe , git @lsudupe

#Libraries---------------------------------
library(Seurat)
library(dplyr)
library(ggplot2)

#Data--------------------------------------
control <- readRDS("./data/initial/sham.sp.rds")
infarto_1dpi <- readRDS("./data/initial/day1.sp.rds")
infarto_7dpi <- readRDS("./data/initial/day7.sp.rds")
infarto_14dpi <- readRDS("./data/initial/day14.sp.rds")

##Merge them
combined_Tokio <- merge(control, y = c(infarto_1dpi, infarto_7dpi, infarto_14dpi ), 
                     add.cell.ids = c("control","1_dpi","7_dpi", "14_dpi"), project = "Tokio")

##Visualization
# Visualize the number of spots counts per sample
png(file.path("./QC/results_initial",filename = "Number of spot per sample.png"))
combined_Tokio@meta.data%>% 
  ggplot(aes(x=sample, fill=sample)) + 
  geom_bar() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("NSpots")
dev.off()

# Visualize the distribution of genes detected per spot via boxplot
png(file.path("./QC/results_initial",filename = "genes detected per spot boxplot.png"))
combined_Tokio@meta.data %>% 
  ggplot(aes(x=sample, y=log10(nFeature_Spatial), fill=sample)) + 
  geom_boxplot() + 
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("Ngenes vs Npots")
dev.off()

# Visualize the distribution of mitochondrial gene expression detected per spot
png(file.path("./QC/results_initial",filename = "mito percentage per spot.png"))
combined_Tokio@meta.data %>% 
  ggplot(aes(color=sample, x=percent_mito, fill=sample)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  geom_vline(xintercept = 0.2)
dev.off()


##########SPATIAL plots

feature.list <- c("nCount_Spatial", "nFeature_Spatial","percent_mito", "log10GenesPerUMI")

for (i in feature.list){
  p1 <- SpatialFeaturePlot(combined_Tokio, features = i,pt.size.factor = 10,combine = FALSE)
  fix.p1 <- scale_fill_continuous(type = "viridis")
  p2 <- lapply(p1, function (x) x + fix.p1)
  
  pdf(paste("./QC/initial_results",i,"spatial_no_range.pdf",sep=""))
  print(CombinePlots(p2))
  dev.off()
}

##specific plots
p1 <- SpatialFeaturePlot(combined_Tokio, features = "nCount_Spatial",pt.size.factor = 10,alpha = 0.5,combine = FALSE)
fix.p1 <- scale_fill_continuous(limits = c(0,110000), 
                                breaks = c(0,110000),
                                type = "viridis")
p2 <- lapply(p1, function (x) x + fix.p1)

pdf(paste("./QC/results_initial/nCount_Spatial.pdf",sep=""))
CombinePlots(p2)
dev.off()

p1 <- SpatialFeaturePlot(combined_Tokio, features = "nFeature_Spatial",pt.size.factor = 10,alpha = 0.5,combine = FALSE)
fix.p1 <- scale_fill_continuous(limits = c(0,9100), 
                                breaks = c(0,9100),
                                type = "viridis")
p2 <- lapply(p1, function (x) x + fix.p1)

pdf(paste("./QC/results_initial/nFeature_Spatial.pdf",sep=""))
CombinePlots(p2)
dev.off()

p1 <- SpatialFeaturePlot(combined_Tokio, features = "percent_mito",pt.size.factor = 10,alpha = 0.5,combine = FALSE)
fix.p1 <- scale_fill_continuous(limits = c(0,100), 
                                breaks = c(0,50,100),
                                type = "viridis")
p2 <- lapply(p1, function (x) x + fix.p1)

pdf(paste("./QC/results_initial/percent_mito.pdf",sep=""))
CombinePlots(p2)
dev.off()


#####save combined
saveRDS(combined_Tokio,"./data/initial/combined_Tokio.rds")

combined.meta <- combined_Tokio@meta.data

#####Violin plots
a <- ggplot(combined.meta, aes(x=sample, y=nCount_Spatial, color=sample)) +
  geom_violin(trim=FALSE) +
  labs(title="nCount_Spatial Tokio",x="sample", y = "nCount_Spatial")

b <- ggplot(combined.meta, aes(x=sample, y=nFeature_Spatial, color=sample)) +
  geom_violin(trim=FALSE) +
  labs(title="nFeature_Spatial Tokio",x="sample", y = "nFeature_Spatial")

c <- ggplot(combined.meta, aes(x=sample, y=percent_mito, color=sample)) +
  geom_violin(trim=FALSE) +
  labs(title="percent_mito Tokio",x="sample", y = "percent_mito")

pdf(paste("./QC/results_initial/features_violin.pdf",sep=""))
cowplot::plot_grid(a,NULL,b,NULL,c, nrow = 5, ncol = 1, rel_heights = c(4,0.5,3,0.5,3))
dev.off()


###########count/feature correlation

# Visualize the correlation between genes detected and number of UMIs and determine whether strong presence of cells with low numbers of genes/UMIs
pdf(paste("./QC/results_initial/count_features_cor.pdf",sep=""))
combined.meta %>% 
  ggplot(aes(x=nCount_Spatial, y=nFeature_Spatial, color=percent_mito)) + 
  geom_point() + 
  scale_colour_gradient(low = "gray90", high = "black") +
  stat_smooth(method=lm) +
  scale_x_log10() + 
  scale_y_log10() + 
  theme_classic() +
  geom_vline(xintercept = 500) +
  geom_hline(yintercept = 250) +
  facet_wrap(~sample)
dev.off()


d <- FeatureScatter(control, feature1 = "nCount_Spatial", feature2 = "nFeature_Spatial")
e <- FeatureScatter(infarto_1dpi, feature1 = "nCount_Spatial", feature2 = "nFeature_Spatial")
f <- FeatureScatter(infarto_7dpi, feature1 = "nCount_Spatial", feature2 = "nFeature_Spatial")
g <- FeatureScatter(infarto_14dpi, feature1 = "nCount_Spatial", feature2 = "nFeature_Spatial")

d + e + f + g


