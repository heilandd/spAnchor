
#Read in example data
setwd("~/Desktop/RecurrentGBM/Data/SPATA_files")

all <- readRDS("817_R_spAnchor_list.RDS")



library(scattermore)
cell_map_df <- as.data.frame(do.call(rbind, all$cell_map))
cell_map_df$match <- cell_map_df$match %>% unlist()

#Add UMAP coords
UMAP <- all$reference@reductions$umap@cell.embeddings %>% as.data.frame() %>% rownames_to_column("match")
meta <- all$reference@meta.data %>% as.data.frame() %>% rownames_to_column("match")
cell_map_df <- cell_map_df %>% left_join(., UMAP, by="match") %>% left_join(., meta, by="match")


#Plot UMAP

colors <- readRDS("~/Desktop/ImmunoSpatial/Paper/colors_cell_deconv.RDS")
col <- colors$colors
names(col) <- colors$annotation_level_4

ggplot(cell_map_df, mapping=aes(x=UMAP_1, y=UMAP_2, color=annotation_level_4))+
  scattermore::geom_scattermore(pointsize=2, interpolate=T, pixels = c(1024, 1024))+
  theme_void()+
  scale_color_manual(values=col)+
  coord_fixed()+
  Seurat::NoLegend()

ggplot(cell_map_df, mapping=aes(x=x, y=y, color=annotation_level_4))+
  scattermore::geom_scattermore(pointsize=2, interpolate=T, pixels = c(1024, 1024))+
  theme_void()+
  scale_color_manual(values=col)+
  coord_fixed()+
  Seurat::NoLegend()


DimPlot(all$reference, group.by = "annotation_level_4", label = T) + Seurat::NoLegend()

#Get a Fill SPATA single cell object

cell_map_df <- cell_map_df %>% filter(match !="Empty")
mat <- reference@assays$SCT@counts[, cell_map_df$match %>% unlist()]
colnames(mat) <- cell_map_df$cells

obj_scSPATA <- SPATA2::initiateSpataObject_CountMtr(coords_df = data.frame(barcodes=cell_map_df$cells,
                                                                           x=cell_map_df$x,
                                                                           y=cell_map_df$y),
                                                    count_mtr = mat, sample_name = "275_T_sc", image = getImage(object))

obj_scSPATA@used_genesets <- object@used_genesets


plotSurface(obj_scSPATA, color_by = "GFAP", pt_size = 0.8, smooth = F, normalize = T, display_image = F)
plotSurface(object, color_by = "GFAP", smooth = F, normalize = T, display_image = F)

plotSurface(obj_scSPATA, color_by = "CHI3L1", pt_size = 0.8, smooth = F, normalize = T, display_image = F, pt_clrsp="Reds")
plotSurface(object, color_by = "CHI3L1", smooth = F, normalize = T, display_image = F, pt_clrsp="Reds")

plotSurface(obj_scSPATA, color_by = "OLIG1", pt_size = 0.8, smooth = F, normalize = T, display_image = F, pt_clrsp="Reds")
plotSurface(object, color_by = "OLIG1", smooth = F, normalize = T, display_image = F, pt_clrsp="Reds")


p1=plotSurface(obj_scSPATA, color_by = "HM_HYPOXIA", pt_size = 0.8, smooth = F, 
               normalize = T, display_image = F, pt_clrsp="Reds", limit=c(0.4, 0.8), oob=scales::squish)+Seurat::NoLegend()
p1=p1+ggtitle("SPATADeconvolution")
p2=plotSurface(object, color_by = "HM_HYPOXIA", smooth = F, 
               normalize = T, display_image = F, pt_clrsp="Reds", limit=c(0.4, 0.8), oob=scales::squish)
p2=p2+ggtitle("Raw Data")

p1+p2

prom="GFAP"

p1=plotSurface(obj_scSPATA, color_by = prom, pt_size = 0.8, smooth = F, 
               normalize = T, display_image = F, pt_clrsp="Greens", limit=c(0.3, 1), oob=scales::squish)+Seurat::NoLegend()
p1=p1+ggtitle("SPATADeconvolution")
p2=plotSurface(object, color_by = prom, smooth = F, 
               normalize = F, display_image = F, pt_clrsp="Greens", limit=c(0, 3), oob=scales::squish)
p2=p2+ggtitle("Raw Data")

p1+p2

prom="CHI3L1"

p1=plotSurface(obj_scSPATA, color_by = prom, pt_size = 0.8, smooth = F, 
               normalize = T, display_image = F, pt_clrsp="Greens", limit=c(0.2, 0.8), oob=scales::squish)+Seurat::NoLegend()
p1=p1+ggtitle("SPATADeconvolution")
p2=plotSurface(object, color_by = prom, smooth = F, 
               normalize = T, display_image = F, pt_clrsp="Greens", limit=c(0.2, 1), oob=scales::squish)
p2=p2+ggtitle("Raw Data")

p1+p2

object %>% getGeneSets(index="Ne")

prom="Neftel_AClike"

p1=plotSurface(obj_scSPATA, color_by = prom, pt_size = 0.8, smooth = F, 
               normalize = T, display_image = F, pt_clrsp="Greens", limit=c(0.3, 0.9), oob=scales::squish)+Seurat::NoLegend()
p1=p1+ggtitle("SPATADeconvolution")
p2=plotSurface(object, color_by = prom, smooth = F, 
               normalize = T, display_image = F, pt_clrsp="Greens", limit=c(0.4, 1), oob=scales::squish)
p2=p2+ggtitle("Raw Data")

p1+p2

plotSurfaceInteractive(obj_scSPATA)



# Summarize again and compare

mat <- reference@assays$SCT@counts[, cell_map_df$match %>% unlist()]
colnames(mat) <- cell_map_df$cells

mat <- mat %>% as.matrix() %>% t() %>% as.data.frame() %>% mutate(barcodes=cell_map_df$barcodes)
mat <- mat %>% group_by(barcodes) %>% summarise_all(.funs = sum)
mat_out <- mat[,2:ncol(mat)] %>% t() %>% Matrix::Matrix(., sparse = T)
colnames(mat_out) <- mat$barcodes

obj_scSPATA_sum <- SPATA2::initiateSpataObject_CountMtr(coords_df = 
                                                          object %>% getCoordsDf() %>% 
                                                          dplyr::select(barcodes,x,y, row, col) %>% 
                                                          filter(barcodes %in% mat$barcodes),
                                                        count_mtr = mat_out, 
                                                        sample_name = "275_T_sum", 
                                                        image = getImage(object))


prom="CHI3L1"

p1=plotSurface(obj_scSPATA_sum, color_by = prom, smooth = T, 
               normalize = T, display_image = F, pt_clrsp="Greens", limit=c(0.2, 1), oob=scales::squish)+Seurat::NoLegend()
p1=p1+ggtitle("SPATADeconvolution")
p2=plotSurface(object, color_by = prom, smooth = T, 
               normalize = T, display_image = F, pt_clrsp="Greens", limit=c(0, 1), oob=scales::squish)
p2=p2+ggtitle("Raw Data")

p1+p2


# get matrix

ncol=object %>% getCoordsDf() %>% pull(col) %>% max()-object %>% getCoordsDf() %>% pull(col) %>% min()
nrow=object %>% getCoordsDf() %>% pull(row) %>% max()-object %>% getCoordsDf() %>% pull(row) %>% min()

gene="CHI3L1"

m1 <- joinWith(object, genes = gene) %>% dplyr::select(row, col, !!sym(gene)) %>% reshape2::acast(row~col, value.var = gene)
m1[is.na(m1)] <- 0

m2 <- joinWith(obj_scSPATA_sum, genes = gene) %>% dplyr::select(row, col, !!sym(gene)) %>% reshape2::acast(row~col, value.var = gene)
m2[is.na(m2)] <- 0
m1 <- m1[rownames(m2), colnames(m2)]


library(CCP)
CCP::p.perm(m1, m2)
cc1 <- CCA::matcor(m1, m2)

cc1$XYcor %>% mean()

plot(as.raster(cc1$XYcor %>% scales::rescale(., c(0,1))))

