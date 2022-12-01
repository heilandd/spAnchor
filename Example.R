library(SPATA2)
library(tidyverse)
library(Seurat)
library(tidyverse)

path_visium=("~/Desktop/SpatialTranscriptomics/Visium/Visium")
setwd(path_visium)
meta.st <- read.csv( "feat.csv", sep=";")

setwd("~/Desktop/SpatialTranscriptomics/Visium/Visium/CancerCell_Revision/Modules_WT")
obj.wt <- meta.st %>% filter(Region=="T") %>% filter(Tumor=="IDH-WT") %>% pull(files) %>% unique()
obj.wt <-data.frame(files=obj.wt,
                    path=paste0("~/Desktop/SpatialTranscriptomics/Visium/Visium/All_SPATA_Revisions/", "Revision_",obj.wt,"_SPATA_CNV_Pred.RDS"),
                    imh.hrs=paste0("~/Desktop/SpatialTranscriptomics/Visium/Visium/",obj.wt,"/outs/spatial/tissue_hires_image.png"))
obj.wt$Seuratstep1=NA

nr <- 1
sample <- obj.wt$files[nr]

object <- readRDS(obj.wt$path[nr])

library(spAnchor)
object <- spAnchor::inferCellPositionsfromHE(object, "~/Desktop/SpatialTranscriptomics/Visium/Visium/275_T")
cell_counts <- spAnchor::runSegmentfromCoords(object)
object <- SPATA2::addFeatures(object, cell_counts$Feature_cells, overwrite = T)

colors <- readRDS("~/Desktop/ImmunoSpatial/Paper/colors_cell_deconv.RDS")
## Try mapping

cell_types <- getSingleCellDeconv(object, deconv_cell_types= colors$annotation_level_4, spot_extension=0.7)

colors %>% filter(annotation_level_4 %in% unique(cell_types$celltypes))
cell_types$celltypes <- factor(cell_types$celltypes, levels = colors$annotation_level_4)
ggplot(cell_types, aes(x,y,color=celltypes))+
  geom_point()+theme_classic()+
  scale_color_manual(values=colors$colors)+
  Seurat::NoLegend()+
  coord_fixed()


### Define optimal gene Space

#Load ref
ref <- readRDS('~/Desktop/TCGA_Data/GBmap/Azimuth/azimuth_core_GBmap.rds')
ref@meta.data$annotation_level_4 <- ref@meta.data$annotation_level_4 %>% str_replace_all(., " ", "_") %>% str_replace_all(., "-", "_") %>% str_replace_all(., "/", "_")
ref@meta.data$annotation_level_4 %>% table() %>% as.data.frame() %>% arrange(Freq)
ref.new <- subset(ref, cells=rownames(ref@meta.data)[ref@meta.data$annotation_level_4!="Neuron"])
rm(ref)

base::options(future.globals.maxSize = 6000 * 1024^2)

reference <-DownScaleSeurat(ref.new, "annotation_level_4", max=10000, min=50)
rm(ref.new)

base::options(future.globals.maxSize = 600 * 1024^2)
metaSpace <- spAnchor::defineMetaSpace(object, reference, cell_type_var="annotation_level_4", scDF=cell_types)
gc()


save <- readRDS("save.RDS")
cell_map <- spAnchor::runMappingGA(object=save$object,
                                             reference=save$reference,
                                             cell_type_var="annotation_level_4",
                                             scDF=save$cell_types,
                                             metaSpace=save$metaSpace, ram = 8, workers=5)


save <- list(object, reference, cell_types, metaSpace)
names(save) <- c("object", "reference", "cell_types", "metaSpace")
saveRDS(save, file="save.RDS")

































