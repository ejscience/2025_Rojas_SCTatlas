setwd("/Users/kriniogiannikou/Desktop/Ernesto_spatial")

install.packages("remotes")
remotes::install_github(
  "jbergenstrahle/STUtility"
)

install.packages("devtools")
library(devtools)

devtools::install_github(
  "jbergenstrahle/STUtility"
  , force = TRUE)


install.packages("fftwtools", repos="http://cran.r-project.org")
library(fftwtools)

install.packages("SeuratObject", type = "source")
library("SeuratObject")

# Load required libraries
library(ggplot2)
library(ggalluvial)
library(dplyr)
library(tidyr)
library(RColorBrewer)
library(sf)
library(ggnewscale)
library(semla)
library(tidyr)
library(tibble)
library(singlet)
library(SeuratData)
library(hdf5r)
library(ggplot2)
library(patchwork)
library(viridis)
library(RColorBrewer)
library(magick)
library(dplyr)
library(jsonlite)
library(grid)
library(plotly)
library(scico)
library(ggsci)
library(patchwork)
library(DT)
library(ggfittext)
library(TabulaMurisSenisData)
library(SingleCellExperiment)
library(Seurat)
library(SeuratWrappers)
library(purrr)
library(RcppML)
library(R.utils)
library(devtools)
library(SeuratWrappers)
library(NMF)
library(tidygraph)

# Ensure all paths are absolute and correctd
samples_SCT04manual <- "/Users/kriniogiannikou/Downloads/SCT04manual/filtered_feature_bc_matrix.h5"
imgs_SCT04manual <- "/Users/kriniogiannikou/Downloads/SCT04manual/spatial/tissue_lowres_image.png"
spotfiles_SCT04manual <- "/Users/kriniogiannikou/Downloads/SCT04manual/spatial/tissue_positions.csv"
json_SCT04manual <- "/Users/kriniogiannikou/Downloads/SCT04manual/spatial/scalefactors_json.json"

# Rebuild infoTable
infoTable_SCT04manual <- tibble(samples = samples_SCT04manual,
                                imgs = imgs_SCT04manual,
                                spotfiles = spotfiles_SCT04manual,
                                json = json_SCT04manual)
# Check the table
View(infoTable_SCT04manual)

# Load the 10x Visium data into a Seurat object
SCT04manual <- ReadVisiumData(infoTable_SCT04manual)
spatial_data_SCT04manual <- GetStaffli(SCT04manual)
SCT04manual <- LoadImages(SCT04manual)
ImagePlot(SCT04manual)

# Plot with semla
MapFeatures(SCT04manual, features = "nFeature_Spatial", 
            image_use = "raw", override_plot_dims = TRUE) & ThemeLegendRight()

MapFeaturesSummary(SCT04manual, features = "nFeature_Spatial", subplot_type = "histogram")

SCT04manual_filtered <- SubsetSTData(SCT04manual, expression = nFeature_Spatial > 30)
SCT04manual_filtered

#Map numeric features- The most basic usage is to map gene expression spatially
# Blue to Red color scale
cols <- brewer.pal(11, "RdBu")


q <- MapFeaturesSummary(SCT04manual, 
                        features = "nFeature_Spatial", 
                        subplot_type = "violin",
                        colors = cols)
q


q <- MapFeaturesSummary(SCT04manual, 
                        features = "CASK", 
                        subplot_type = "violin",
                        colors = cols)

p1 <- MapFeatures(SCT04manual, 
                  features = "DPPA5", 
                  max_cutoff = 1, # Full data
                  colors = cols)

p2 <- MapFeatures(SCT04manual, 
                  section_number = 1, # Show only the first sample
                  features = "DPPA5", 
                  max_cutoff = 0.95, # Max value set to the 95th percentile
                  min_cutoff = 0.05, # Min value set to the 5th percentile
                  colors = cols)

p1 | p2


SCT04manual <- NormalizeData(SCT04manual)
SCT04manual <- ScaleData(SCT04manual)
SCT04manual <- FindVariableFeatures(SCT04manual, nfeatures = 10000, verbose = FALSE)
SCT04manual <- RunPCA(SCT04manual, verbose = FALSE)
SCT04manual <- FindNeighbors(SCT04manual, reduction = "pca", dims = 1:30)
SCT04manual <- FindClusters(SCT04manual, resolution=0.8)
SCT04manual <- RunUMAP(SCT04manual, reduction = "pca", dims = 1:30, n.neighbors = 10)
#SCT04manual <- LoadImages(SCT04manual)

pdf("SCT04manual_plot-UMAP-0.8-resolution.pdf", width = 12, height = 8)
DimPlot(SCT04manual, reduction = "umap", pt.size = 1.1)
dev.off()

MapLabels(SCT04manual, column_name = "seurat_clusters", 
          image_use = "raw", pt_alpha = 0.6, pt_size = 2) +
  plot_layout(guides = "collect") &
  theme(legend.position = "right") &
  guides(fill = guide_legend(override.aes = list(size = 3), ncol = 2))


SCT04manual$cluster_9 <- ifelse(SCT04manual$seurat_clusters %in% "9", "9", NA)
MapLabels(SCT04manual, column_name = "cluster_9", override_plot_dims = TRUE, 
          image_use = "raw", drop_na = TRUE, pt_size = 2) +
  plot_layout(guides = "collect") &
  theme(legend.position = "right") &
  guides(fill = guide_legend(override.aes = list(size = 3), ncol = 2))

#Disconnect regions
SCT04manual <- DisconnectRegions(SCT04manual, column_name = "seurat_clusters", selected_groups = "9")

MapLabels(SCT04manual, column_name = "9_split", override_plot_dims = TRUE, 
          image_use = "raw", drop_na = TRUE, pt_size = 2) +
  plot_layout(guides = "collect") &
  theme(legend.position = "right") &
  guides(fill = guide_legend(override.aes = list(size = 3), ncol = 2))

#Radial distance
SCT04manual$cluster_3 <- ifelse(SCT04manual$seurat_clusters %in% "3", "3", NA)
MapLabels(SCT04manual, column_name = "cluster_3", override_plot_dims = TRUE, 
          image_use = "raw", drop_na = TRUE, pt_size = 2) +
  theme(legend.position = "right") &
  guides(fill = guide_legend(override.aes = list(size = 5), ncol = 2))

SCT04manual <- RadialDistance(SCT04manual, column_name = "seurat_clusters", selected_groups = "3")
MapFeatures(SCT04manual, features = "r_dist_3", center_zero = TRUE, pt_size = 1.5, 
            colors = RColorBrewer::brewer.pal(n = 11, name = "RdBu") |> rev(),
            override_plot_dims = TRUE)

SCT04manual$r_dist_3 <- (100/273)*SCT04manual$r_dist_3
SCT04manual <- RadialDistance(SCT04manual, column_name = "seurat_clusters", 
                              selected_groups = "3", convert_to_microns = TRUE)

# Convert radial distances 
SCT04manual$r_dist_3_sqrt <- sign(SCT04manual$r_dist_3)*sqrt(abs(SCT04manual$r_dist_3))
MapFeatures(SCT04manual, features = "r_dist_3_sqrt", center_zero = TRUE, pt_size = 1.5, 
            colors = RColorBrewer::brewer.pal(n = 11, name = "RdBu") |> rev(),
            override_plot_dims = TRUE)

#we can now explore the expression of certain genes as a function of distance from our ROI.
sel_genes <- c("TMEM178B", "PPP2R2B")

SCT04manual[[]] |> 
  bind_cols(FetchData(SCT04manual, vars = sel_genes)) |> 
  filter(r_dist_3 < 1e3) |> 
  pivot_longer(all_of(sel_genes), names_to = "variable", values_to = "value") |> 
  ggplot(aes(r_dist_3, value, color = variable)) +
  geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs")) +
  geom_vline(aes(xintercept = 0, color = "border"), linetype = "dashed")


MapFeatures(SubsetSTData(SCT04manual, expression = r_dist_3 < 1e3), 
            features = sel_genes, override_plot_dims = TRUE, scale_alpha = TRUE,
            image_use = "raw", pt_size = 2)

############
pdf("SCT04manual_plot-PECAM1-VWF.pdf", width = 12, height = 8)
MapFeatures(SCT04manual, features = c("PECAM1", "VWF"), image_use = "raw", pt_size = 2.5,
            override_plot_dims = TRUE, colors = RColorBrewer::brewer.pal(n = 11, name = "Spectral") |> rev())
dev.off()

pdf("SCT04manual_plot-stroma-COL1A2-PRRX1.pdf", width = 12, height = 8)
MapFeatures(SCT04manual, features = c("COL1A2", "PRRX1"), image_use = "raw", pt_size = 2.5,
            override_plot_dims = TRUE, colors = RColorBrewer::brewer.pal(n = 11, name = "Spectral") |> rev())
dev.off()

pdf("SCT04manual_plot-immune-IKZF1-PTPRC.pdf", width = 12, height = 8)
MapFeatures(SCT04manual, features = c("IKZF1", "PTPRC"), image_use = "raw", pt_size = 2.5,
            override_plot_dims = TRUE, colors = RColorBrewer::brewer.pal(n = 11, name = "Spectral") |> rev())
dev.off()

pdf("SCT04manual_plot-epithelia-BAIAP2L1-KRT19.pdf", width = 12, height = 8)
MapFeatures(SCT04manual, features = c("BAIAP2L1", "KRT19"), image_use = "raw", pt_size = 2.5,
            override_plot_dims = TRUE, colors = RColorBrewer::brewer.pal(n = 11, name = "Spectral") |> rev())
dev.off()

pdf("SCT04manual_plot-neuroendoderm-CTNNA2-NRXN1.pdf", width = 12, height = 8)
MapFeatures(SCT04manual, features = c("CTNNA2", "NRXN1"), image_use = "raw", pt_size = 2.5,
            override_plot_dims = TRUE, colors = RColorBrewer::brewer.pal(n = 11, name = "Spectral") |> rev())
dev.off()

pdf("SCT05_plot-Tcells.pdf", width = 12, height = 8)
MapFeatures(SCT05outs, features = c("IL2RG", "CD3E", "CD7", "NKG7"), image_use = "raw", pt_size = 1.0,
            override_plot_dims = TRUE, colors = RColorBrewer::brewer.pal(n = 11, name = "Spectral") |> rev())
dev.off()


############
pdf("SCT04manual_plot-endothelia-CTNNA2.pdf", width = 12, height = 8)
MapFeatures(SCT04manual, features = c("V", "PECAM1"), image_use = "raw", pt_size = 1.8,
            override_plot_dims = TRUE, colors = c("blue", "lightblue", "lightgray") |> rev())
dev.off()


MapFeatures(SCT04manual, arrange_features = "row", features = c("DOCK2", "PTPRC"), subplot_type = "violin", colors = RColorBrewer::brewer.pal(n = 11, name = "RdBu"), pt.size = 4 |> rev())

#Map multiple numeric features
#With the MapMultipleFeatures() function, it is possible to plot several numeric features within the same plot.
pdf("se_SCT04manual_plot-pluripotency-markers.pdf", width = 12, height = 8)

MapMultipleFeatures(se_SCT04manual,
                    features = c("ALPL", "SOX2", "KLF4", "LIN28A", "LIN28B", "TFAP2C",
                                 "PRDM1", "SOX17", "KIT", "DND1", "NANOS3", "DDX4", "DAZL"),
                    colors = c(
                      "#AA6599", "#88CCEE", "#332288", "#44AA99", "#FF6347",
                      "#6B5ACD", "#20B2AA", "#FFD700", "#8B0000", "#D2691E",
                      "#4DAF4A", "#984EB1", "#E41A1C", "#377EB8"
                    ),
                    ncol = 2.3,
                    pt_size = 1.6
)

dev.off()


MapMultipleFeatures(SCT04manual,
                    features = c("MKI67", "TOP2A", "CDK1"),
                    colors = c("#AA4499", "#88CCEE", "#31A354"),
                    ncol = 2,
                    pt_size = 1.8
)

dev.off()

MapMultipleFeatures(SCT04manual,
                    features = c("POU5F1", "NANOG", "ALPL", "SOX2", "KLF4", "LIN28A", "LIN28B", "TFAP2C",
                                 "PRDM1", "SOX17", "KIT", "DND1", "NANOS3", "DDX4", "DAZL"),
                    colors = c("#AA4499", "#88CCEE", "#332288", "#44AA99", 
                               "#FF6347", "#6A5ACD", "#20B2AA", "#FFD700", 
                               "#8B0000", "#D2691E"),
                    ncol = 2,
                    pt_size = 1.8
)

dev.off()

#errors
AnglePlot(SCT04manual, column_name = "cluster_3", selected_group = "5", pt_size = 2,
          crop_area = c(0.2, 0.3, 1, 1), image_use = "raw", radius = 0.3, drop_na = TRUE)

#errors
SCT04manual <- RadialDistance(SCT04manual, column_name = "seurat_clusters", selected_groups = "3", 
                              angles = c(200, 240))

MapFeatures(SCT04manual, features = "r_dist_3", pt_size = 2)

SCT04manual <- RadialDistance(SCT04manual, column_name = "seurat_clusters", selected_groups = "3", 
                              angles = c(120, 320), angles_nbreaks = 6)

MapLabels(SCT04manual, column_name = "intervals_3", pt_size = 2, drop_na = TRUE,
          colors = RColorBrewer::brewer.pal(n = 6, name = "Spectral")) &
  guides(fill = guide_legend(override.aes = list(size = 5), ncol = 3))

sel_genes <- c("CRISP3", "TSC2")

SCT04manual[[]] |> 
  bind_cols(FetchData(SCT04manual, vars = sel_genes)) |> 
  filter(r_dist_3 < 1e3) |> 
  pivot_longer(all_of(sel_genes), names_to = "variable", values_to = "value") |> 
  ggplot(aes(r_dist_3, value, color = variable)) +
  geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs")) +
  geom_vline(aes(xintercept = 0, color = "border"), linetype = "dashed") +
  facet_wrap(~ intervals_3, ncol = 3) +
  theme_minimal()

sel_genes <- c("PGM5-AS1", "CXCL14", "IGHG3")

SCT04manual[[]] |> 
  bind_cols(FetchData(SCT04manual, vars = sel_genes)) |> 
  filter(r_dist_3 < 1e3) |> 
  pivot_longer(all_of(sel_genes), names_to = "variable", values_to = "value") |> 
  group_by(variable) |> 
  mutate(value_centered = value - mean(value)) |> # Note that the values have been centered
  ggplot(aes(r_dist_3, value_centered, color = variable)) +
  geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs")) +
  geom_vline(aes(xintercept = 0, color = "border"), linetype = "dashed") +
  facet_wrap(~ intervals_3, ncol = 3) +
  theme_minimal()

#!!!!!
MapFeatures(SCT04manual, crop_area = c(0.45, 0.6, 0.8, 0.93),
            features = sel_genes, scale_alpha = TRUE, ncol = 3,
            image_use = "raw", pt_size = 1.5, colors = viridis::viridis(n = 9),
            max_cutoff = 0.99) &
  theme(legend.position = "right", legend.text = element_text(angle = 0))

#Region neightbors
SCT04manual <- RegionNeighbors(SCT04manual, 
                               column_name = "seurat_clusters", 
                               column_labels = "3")

MapLabels(SCT04manual, crop_area = c(0.45, 0.6, 0.8, 0.93),
          column_name = "nb_to_3", drop_na = TRUE,
          image_use = "raw", pt_size = 3)

SCT04manual <- RegionNeighbors(SCT04manual, 
                               column_name = "seurat_clusters", 
                               column_labels = "3", 
                               mode = "inner")

MapLabels(SCT04manual, 
          crop_area = c(0.45, 0.6, 0.8, 0.93),
          column_name = "inner_border_3",
          image_use = "raw", 
          pt_size = 3, 
          drop_na = TRUE)

SCT04manual <- RegionNeighbors(SCT04manual, column_name = "seurat_clusters", 
                               column_labels = "3", mode = "inner_outer")

MapLabels(SCT04manual, crop_area = c(0.45, 0.6, 0.8, 0.93),
          column_name = "nb_to_3",
          image_use = "raw", pt_size = 3, drop_na = TRUE)


SCT04manual <- RegionNeighbors(SCT04manual, column_name = "seurat_clusters", column_key = "(all)nb_to_",
                               column_labels = "3", mode = "all_inner_outer")

MapLabels(SCT04manual, crop_area = c(0.45, 0.6, 0.8, 0.93),
          column_name = "(all)nb_to_3",
          image_use = "raw", pt_size = 3, drop_na = TRUE)


border_markers <- FindMarkers(SCT04manual, ident.1 = "nb_to_3", 
                              ident.2 = "3", group.by = "nb_to_3")

# Filter results
border_markers_up <- border_markers |> 
  filter(p_val_adj < 0.01, avg_log2FC > 0) |> 
  arrange(-avg_log2FC) |> 
  slice_head(n = 10)

# Subset Seurat object to include border spots
SCT04manual_border_spots <- SubsetSTData(SCT04manual, expression = nb_to_3 %in% c("3", "nb_to_3"))

# Violin plot
VlnPlot(SCT04manual_border_spots, features = rownames(border_markers_up), group.by = "nb_to_3") &
  theme(axis.title = element_blank())


VlnPlot(SCT04manual, features = c("MKI67",  "TOP2A", "CDK1")) &
  theme(axis.title = element_blank())


# Example genes to plot
genes_to_plot <- c("TSC2", "SOX3", "MITF")

#*********************************
# Violin plot for the specified genes
VlnPlot(SCT04manual, features = genes_to_plot) +
  theme(axis.title = element_blank())

png("SCT04manual_FeaturePlot-normalized-SCTransform3.png", width=9, height=9, units="in", res=300)
VlnPlot(SCT04manual, features = c("SPARC", "TIMP1", "B2M", "COL1A2"), ncol = 2, group.by = "seurat_clusters")
dev.off()   


top_100_positive_corr_to_XAratio_genes = c("EDA", "CASK", "NHS", "OPHN1", "COL4A5", "KLHL13", "MID1", "PBX1", "RBMX", "PTPRD", "PLEKHA5", "GRIK2", "PTPN13", "LSAMP", "HPSE2", "COL4A6", 
                                           "ANKRD30BL", "RNF150", "SRGAP3", "NKAIN3", "HIF3A", "BMPR1A", "TMEM132C", "GPC3", "BOC", "GPC5", "MALRD1", "PLAG1", "CTPS2", "XPO1", "RUNX1T1", "KANTR", 
                                           "MAGED1", "SLC44A5", "NCAM1", "CSRNP3", "ARMCX4", "LRCH2", "SSBP2", "EPHA3", "PDE5A", "EPHA7", "HDAC8", "PAK3", "EYA1", "DCC", "GRM7", "GRIP1", "ZC4H2", "TET1",
                                           "NEO1", "ASXL3", "NLGN4X", "PTPRT", "EDNRA", "MEIS1", "IL1RAPL1", "IGSF1", "ROBO1", "IGF1R", "FGFR2", "GARNL3", "ADAMTS19", "MEIS2", "SESN3", "BBX", "NCOA1", "CHIC1", 
                                           "RARB", "ITGA8", "FREM1", "MAPK10", "ERBB4", "HMGA2", "TBL1X", "GRID2", "CDH11", "TBX5", "DDB2", "COL3A1", "AFF3", "TENM4", "AUTS2", "TRPC1", "CECR2", "NAALAD2", "POSTN", 
                                           "BMP5", "HS6ST2", "COL1A2", "NBEA", "HIBCH", "FBLN1", "STAG2", "ADGRV1", "PDZRN3", "SGCD", "SEMA3A", "SPINK5", "GPR173")

top_100_negative_corr_to_XAratio_genes = c("LDLR", "MCL1", "RFX2", "HIF1A", "PITPNC1", "FNDC3B", "LRRFIP1", "YBX3", "NAMPT", "NEDD9", "CPNE8", "TUBB6", "GRB10", "KLF9",
                                           "AKAP13", "LMNA", "EMP1", "MAP2K3", "CRY1", "PACSIN2", "FOSL1", "OSMR", "SNX9", "SLC7A1", "SIK3", "PFKFB3", "ITGB1", "MYO1E", "TACC1", "PFKP",
                                           "BACE2", "MGMT", "SDCBP", "TNFRSF10D", "DUSP1", "PDLIM1", "PER1", "B3GNT5", "FOSL2", "FOSB", "EPAS1", "ZSWIM6", "NFATC1", "CCNL1", "NECTIN2", "ELL", "KDM6B", "NFKBIA", 
                                           "ESYT2", "PDLIM5", "PLAUR", "MAFF", "IL6R", "FOS", "GBE1", "FKBP5", "DDX21", "ZFP36", "AHNAK", "KLF6", "AGO2", "COL18A1", "TIPARP", "TEX14", "ATP10A", "HIVEP2", "CFH", "MIDN",
                                           "YWHAG", "NFATC2", "NR4A1", "SPSB1", "MAPKAPK2", "RNF19A", "NPAS2", "JUNB", "MYADM", "TM4SF4", "ZFAND3", "ELF1", "DIAPH1", "UBC", "IL4R", "SBNO2", "YWHAZ", "ANXA5",
                                           "MACF1", "ITGA6", "ST3GAL1", "EVA1C", "RELB", "ITPKC", "BAIAP2", "SYNJ2", "TES", "KLF13", "LATS2", "MAGI1", "ELOVL5")


library(UCell)

ANGIOGENESIS = c("APOH", "APP", "CCND2", "COL3A1", "COL5A2", "CXCL6", "FGFR1", "FSTL1", "ITGAV", "JAG1", "JAG2", "KCNJ8", "LPL", "LRPAP1", "LUM", "MSX1", "NRP1", "OLR1", "PDGFA", "PF4", "PGLYRP1", "POSTN", "PRG2", "PTK2", "S100A4", "SERPINA5", "SLCO2A1", "SPP1", "STC1", "THBD", "TIMP1", "TNFRSF21", "VAV2", "VCAN", "VEGFA", "VTN")
APICAL_JUNCTION = c("ACTA1", "ACTB", "ACTC1", "ACTG1", "ACTG2", "ACTN1", "ACTN2", "ACTN3", "ACTN4", "ADAM15", "ADAM23", "ADAM9", "ADAMTS5", "ADRA1B", "AKT2", "AKT3", "ALOX15B", "AMH", "AMIGO1", "AMIGO2", "ARHGEF6", "ARPC2", "ATP1A3", "B4GALT1", "BAIAP2", "BMP1", "CADM2", "CADM3", "CALB2", "CAP1", "CD209", "CD274", "CD276", "CD34", "CD86", "CD99", "CDH1", "CDH11", "CDH15", "CDH3", "CDH4", "CDH6", "CDH8", "CDK8", "CDSN", "CERCAM", "CLDN11", "CLDN14", "CLDN15", "CLDN18", "CLDN19", "CLDN4", "CLDN5", "CLDN6", "CLDN7", "CLDN8", "CLDN9", "CNN2", "CNTN1", "COL16A1", "COL17A1", "COL9A1", "CRAT", "CRB3", "CTNNA1", "CTNND1", "CX3CL1", "DHX16", "DLG1", "DMP1", "DSC1", "DSC3", "EGFR", "EPB41L2", "EVL", "EXOC4", "FBN1", "FLNC", "FSCN1", "FYB1", "GAMT", "GNAI1", "GNAI2", "GRB7", "GTF2F1", "HADH", "HRAS", "ICAM1", "ICAM2", "ICAM4", "ICAM5", "IKBKG", "INPPL1", "INSIG1", "IRS1", "ITGA10", "ITGA2", "ITGA3", "ITGA9", "ITGB1", "ITGB4", "JAM3", "JUP", "KCNH2", "KRT31", "LAMA3", "LAMB3", "LAMC2", "LAYN", "LDLRAP1", "LIMA1", "MADCAM1", "MAP3K20", "MAP4K2", "MAPK11", "MAPK13", "MAPK14", "MDK", "MMP2", "MMP9", "MPZL1", "MPZL2", "MSN", "MVD", "MYH10", "MYH9", "MYL12B", "MYL9", "NECTIN1", "NECTIN2", "NECTIN3", "NECTIN4", "NEGR1", "NEXN", "NF1", "NF2", "NFASC", "NLGN2", "NLGN3", "NRAP", "NRTN", "NRXN2", "PARD6G", "PARVA", "PBX2", "PCDH1", "PDZD3", "PECAM1", "PFN1", "PIK3CB", "PIK3R3", "PKD1", "PLCG1", "PPP2R2C", "PTEN", "PTK2", "PTPRC", "RAC2", "RASA1", "RHOF", "RRAS", "RSU1", "SDC3", "SGCE", "SHC1", "SHROOM2", "SIRPA", "SKAP2", "SLC30A3", "SLIT2", "SORBS3", "SPEG", "SRC", "STX4", "SYK", "SYMPK", "TAOK2", "TGFBI", "THBS3", "THY1", "TIAL1", "TJP1", "TMEM8B", "TNFRSF11B", "TRAF1", "TRO", "TSC1", "TSPAN4", "TUBG1", "VASP", "VAV2", "VCAM1", "VCAN", "VCL", "VWF", "WASL", "WNK4", "YWHAH", "ZYX")
APICAL_SURFACE=c("ADAM10", "ADIPOR2", "AFAP1L2", "AKAP7", "APP", "ATP6V0A4", "ATP8B1", "B4GALT1", "BRCA1", "CD160", "CROCC", "CRYBG1", "CX3CL1", "DCBLD2", "EFNA5", "EPHB4", "FLOT2", "GAS1", "GATA3", "GHRL", "GSTM3", "HSPB1", "IL2RB", "IL2RG", "LYN", "LYPD3", "MAL", "MDGA1", "NCOA6", "NTNG1", "PCSK9", "PKHD1", "PLAUR", "RHCG", "RTN4RL1", "SCUBE1", "SHROOM2", "SLC22A12", "SLC2A4", "SLC34A3", "SRPX", "SULF2", "THY1", "TMEM8B")
APOPTOSIS=c("ADD1", "AIFM3", "ANKH", "ANXA1", "APP", "ATF3", "AVPR1A", "BAX", "BCAP31", "BCL10", "BCL2L1", "BCL2L10", "BCL2L11", "BCL2L2", "BGN", "BID", "BIK", "BIRC3", "BMF", "BMP2", "BNIP3L", "BRCA1", "BTG2", "BTG3", "CASP2", "CASP3", "CASP4", "CASP6", "CASP7", "CASP8", "CASP9", "CAV1", "CCNA1", "CCND1", "CCND2", "CD14", "CD2", "CD38", "CD44", "CD69", "CDC25B", "CDK2", "CDKN1A", "CDKN1B", "CFLAR", "CLU", "CREBBP", "CTH", "CTNNB1", "CYLD", "DAP", "DAP3", "DCN", "DDIT3", "DFFA", "DNAJA1", "DNAJC3", "DNM1L", "DPYD", "EBP", "EGR3", "EMP1", "ENO2", "ERBB2", "ERBB3", "EREG", "ETF1", "F2", "F2R", "FAS", "FASLG", "FDXR", "FEZ1", "GADD45A", "GADD45B", "GCH1", "GNA15", "GPX3", "GPX4", "GSN", "GSR", "GSTM1", "HGF", "HMGB2", "HMOX1", "HSPB1", "IER3", "IFITM3", "IFNB1", "IFNGR1", "IGF2R", "IGFBP6", "IL18", "IL1A", "IL1B", "IL6", "IRF1", "ISG20", "JUN", "LEF1", "LGALS3", "LMNA", "LUM", "MADD", "MCL1", "MGMT", "MMP2", "NEDD9", "NEFH", "PAK1", "PDCD4", "PDGFRB", "PEA15", "PLAT", "PLCB2", "PLPPR4", "PMAIP1", "PPP2R5B", "PPP3R1", "PPT1", "PRF1", "PSEN1", "PSEN2", "PTK2", "RARA", "RELA", "RETSAT", "RHOB", "RHOT2", "RNASEL", "ROCK1", "SAT1", "SATB1", "SC5D", "SLC20A1", "SMAD7", "SOD1", "SOD2", "SPTAN1", "SQSTM1", "TAP1", "TGFB2", "TGFBR3", "TIMP1", "TIMP2", "TIMP3", "TNF", "TNFRSF12A", "TNFSF10", "TOP2A", "TSPO", "TXNIP", "VDAC2", "WEE1", "XIAP")
COMPLEMENT=c("ACTN2", "ADAM9", "ADRA2B", "AKAP10", "ANG", "ANXA5", "APOA4", "APOBEC3F", "APOBEC3G", "APOC1", "ATOX1", "BRPF3", "C1QA", "C1QC", "C1R", "C1S", "C2", "C3", "C4BPB", "C9", "CA2", "CALM1", "CALM3", "CASP10", "CASP3", "CASP4", "CASP5", "CASP7", "CASP9", "CBLB", "CCL5", "CD36", "CD40LG", "CD46", "CD55", "CD59", "CDA", "CDH13", "CDK5R1", "CEBPB", "CFB", "CFH", "CLU", "COL4A2", "CP", "CPM", "CPQ", "CR1", "CR2", "CSRP1", "CTSB", "CTSC", "CTSD", "CTSH", "CTSL", "CTSO", "CTSS", "CTSV", "CXCL1", "DGKG", "DGKH", "DOCK10", "DOCK4", "DOCK9", "DPP4", "DUSP5", "DUSP6", "DYRK2", "EHD1", "ERAP2", "F10", "F2", "F3", "F5", "F7", "F8", "FCER1G", "FCN1", "FDX1", "FN1", "FYN", "GATA3", "GCA", "GMFB", "GNAI2", "GNAI3", "GNB2", "GNB4", "GNG2", "GNGT2", "GP9", "GPD2", "GRB2", "GZMA", "GZMB", "GZMK", "HNF4A", "HPCAL4", "HSPA1A", "HSPA5", "IL6", "IRF1", "IRF2", "IRF7", "ITGAM", "ITIH1", "JAK2", "KCNIP2", "KCNIP3", "KIF2A", "KLK1", "KLKB1", "KYNU", "L3MBTL4", "LAMP2", "LAP3", "LCK", "LCP2", "LGALS3", "LGMN", "LIPA", "LRP1", "LTA4H", "LTF", "LYN", "MAFF", "ME1", "MMP12", "MMP13", "MMP14", "MMP15", "MMP8", "MSRB1", "MT3", "NOTCH4", "OLR1", "PCLO", "PCSK9", "PDGFB", "PDP1", "PFN1", "PHEX", "PIK3CA", "PIK3CG", "PIK3R5", "PIM1", "PLA2G4A", "PLA2G7", "PLAT", "PLAUR", "PLEK", "PLG", "PLSCR1", "PPP2CB", "PPP4C", "PRCP", "PRDM4", "PREP", "PRKCD", "PRSS3", "PRSS36", "PSEN1", "PSMB9", "RABIF", "RAF1", "RASGRP1", "RBSN", "RCE1", "RHOG", "RNF4", "S100A12", "S100A13", "S100A9", "SCG3", "SERPINA1", "SERPINB2", "SERPINC1", "SERPINE1", "SERPING1", "SH2B3", "SIRT6", "SPOCK2", "SRC", "STX4", "TFPI2", "TIMP1", "TIMP2", "TMPRSS6", "TNFAIP3", "USP14", "USP15", "USP16", "USP8", "VCPIP1", "WAS", "XPNPEP1", "ZEB1", "ZFPM2")
E2F_TARGETS=c("ANP32E", "ASF1A", "ASF1B", "ATAD2", "AURKA", "AURKB", "BARD1", "BIRC5", "BRCA1", "BRCA2", "BRMS1L", "BUB1B", "CBX5", "CCNB2", "CCNE1", "CCP110", "CDC20", "CDC25A", "CDC25B", "CDCA3", "CDCA8", "CDK1", "CDK4", "CDKN1A", "CDKN1B", "CDKN2A", "CDKN2C", "CDKN3", "CENPE", "CENPM", "CHEK1", "CHEK2", "CIT", "CKS1B", "CKS2", "CNOT9", "CSE1L", "CTCF", "CTPS1", "DCK", "DCLRE1B", "DCTPP1", "DDX39A", "DEK", "DEPDC1", "DIAPH3", "DLGAP5", "DNMT1", "DONSON", "DSCC1", "DUT", "E2F8", "EED", "EIF2S1", "ESPL1", "EXOSC8", "EZH2", "GINS1", "GINS3", "GINS4", "GSPT1", "HELLS", "HMGA1", "HMGB2", "HMGB3", "HMMR", "HNRNPD", "HUS1", "ILF3", "ING3", "JPT1", "KIF18B", "KIF22", "KIF2C", "KIF4A", "KPNA2", "LBR", "LIG1", "LMNB1", "LUC7L3", "LYAR", "MAD2L1", "MCM2", "MCM3", "MCM4", "MCM5", "MCM6", "MCM7", "MELK", "MKI67", "MLH1", "MMS22L", "MRE11", "MSH2", "MTHFD2", "MXD3", "MYBL2", "MYC", "NAA38", "NASP", "NBN", "NCAPD2", "NME1", "NOLC1", "NOP56", "NUDT21", "NUP107", "NUP153", "NUP205", "ORC2", "ORC6", "PA2G4", "PAICS", "PAN2", "PCNA", "PDS5B", "PHF5A", "PLK1", "PLK4", "PMS2", "PNN", "POLA2", "POLD1", "POLD2", "POLD3", "POLE", "POLE4", "POP7", "PPM1D", "PPP1R8", "PRDX4", "PRIM2", "PRKDC", "PRPS1", "PSIP1", "PSMC3IP", "PTTG1", "RAD1", "RAD21", "RAD50", "RAD51AP1", "RAD51C", "RBBP7", "RFC1", "RFC2", "RFC3", "RNASEH2A", "RPA1", "RPA2", "RPA3", "RRM2", "SHMT1", "SLBP", "SMC1A", "SMC3", "SMC4", "SMC6", "SNRPB", "SPAG5", "SPC24", "SPC25", "SRSF1", "SRSF2", "SSRP1", "STAG1", "STMN1", "SUV39H1", "SYNCRIP", "TACC3", "TBRG4", "TCF19", "TFRC", "TIMELESS", "TIPIN", "TK1", "TMPO", "TOP2A", "TP53", "TRA2B", "TRIP13", "TUBB", "TUBG1", "UBE2T", "UBR7", "UNG", "USP1", "WDR90", "WEE1", "XPO1", "XRCC6", "ZW10")
GLYCOLYSIS=c("ABCB6", "ADORA2B", "AGL", "AGRN", "AK4", "AKR1A1", "ALDH7A1", "ALDH9A1", "ALDOB", "ALG1", "ANG", "ANGPTL4", "ANKZF1", "ARPP19", "ARTN", "AURKA", "B3GALT6", "B3GAT1", "B3GAT3", "B3GNT3", "B4GALT1", "B4GALT2", "B4GALT4", "B4GALT7", "BIK", "BPNT1", "CACNA1H", "CAPN5", "CASP6", "CD44", "CDK1", "CENPA", "CHPF", "CHPF2", "CHST1", "CHST12", "CHST2", "CHST4", "CHST6", "CITED2", "CLDN3", "CLDN9", "CLN6", "COG2", "COL5A1", "COPB2", "CTH", "CXCR4", "CYB5A", "DCN", "DDIT4", "DEPDC1", "DLD", "DPYSL4", "DSC2", "ECD", "EFNA3", "EGFR", "EGLN3", "ELF3", "ENO1", "ENO2", "ERO1A", "EXT1", "EXT2", "FAM162A", "FBP2", "FUT8", "G6PD", "GAL3ST1", "GALE", "GALK1", "GALK2", "GAPDHS", "GCLC", "GFPT1", "GLCE", "GLRX", "GMPPA", "GMPPB", "GNE", "GNPDA1", "GOT1", "GOT2", "GPC1", "GPC3", "GPC4", "GPR87", "GUSB", "GYS1", "GYS2", "HAX1", "HDLBP", "HK2", "HMMR", "HOMER1", "HS2ST1", "HS6ST2", "HSPA5", "IDH1", "IDUA", "IER3", "IGFBP3", "IL13RA1", "IRS2", "ISG20", "KDELR3", "KIF20A", "KIF2A", "LCT", "LDHC", "LHPP", "LHX9", "MDH1", "MDH2", "ME1", "ME2", "MED24", "MERTK", "MET", "MIF", "MIOX", "MPI", "MXI1", "NANP", "NASP", "NDST3", "NDUFV3", "NOL3", "NSDHL", "NT5E", "P4HA1", "P4HA2", "PAM", "PAXIP1", "PC", "PDK3", "PFKFB1", "PFKP", "PGAM1", "PGAM2", "PGK1", "PGLS", "PGM2", "PHKA2", "PKM", "PKP2", "PLOD1", "PLOD2", "PMM2", "POLR3K", "PPFIA4", "PPIA", "PPP2CB", "PRPS1", "PSMC4", "PYGB", "PYGL", "QSOX1", "RBCK1", "RPE", "RRAGD", "SAP30", "SDC1", "SDC2", "SDC3", "SLC16A3", "SLC25A10", "SLC25A13", "SLC35A3", "SLC37A4", "SOD1", "SOX9", "SPAG4", "SRD5A3", "STC1", "STC2", "STMN1", "TALDO1", "TFF3", "TGFA", "TGFBI", "TKTL1", "TPBG", "TPST1", "UGP2", "VCAN", "VEGFA", "VLDLR", "XYLT2", "ZNF292")
HEDGEHOG_SIGNALING=c("ACHE", "ADGRG1", "AMOT", "CDK5R1", "CDK6", "CELSR1", "CNTFR", "CRMP1", "DPYSL2", "ETS2", "GLI1", "HEY1", "HEY2", "L1CAM", "LDB1", "MYH9", "NF1", "NKX6-1", "NRCAM", "NRP1", "NRP2", "OPHN1", "PLG", "PML", "PTCH1", "RASA1", "RTN1", "SCG2", "SHH", "SLIT1", "THY1", "TLE1", "TLE3", "UNC5C", "VEGFA", "VLDLR")
IL2_STAT5_SIGNALING=c("ABCB1", "ADAM19", "AGER", "AHCY", "AHNAK", "AHR", "ALCAM", "AMACR", "ANXA4", "APLP1", "ARL4A", "BATF", "BATF3", "BCL2", "BCL2L1", "BHLHE40", "BMP2", "BMPR2", "CA2", "CAPG", "CAPN3", "CASP3", "CCND2", "CCND3", "CCNE1", "CCR4", "CD44", "CD48", "CD79B", "CD81", "CD83", "CD86", "CDC42SE2", "CDC6", "CDCP1", "CDKN1C", "CISH", "CKAP4", "COCH", "COL6A1", "CSF1", "CSF2", "CST7", "CTLA4", "CTSZ", "CXCL10", "CYFIP1", "DCPS", "DENND5A", "DHRS3", "DRC1", "ECM1", "EEF1AKMT1", "EMP1", "ENO3", "ENPP1", "EOMES", "ETFBKMT", "ETV4", "F2RL2", "FAH", "FAM126B", "FGL2", "FLT3LG", "FURIN", "GADD45B", "GALM", "GATA1", "GBP4", "GLIPR2", "GPR65", "GPR83", "GPX4", "GSTO1", "GUCY1B1", "HIPK2", "HK2", "HOPX", "HUWE1", "ICOS", "IFITM3", "IFNGR1", "IGF1R", "IGF2R", "IKZF2", "IKZF4", "IL10", "IL10RA", "IL13", "IL18R1", "IL1R2", "IL1RL1", "IL2RA", "IL2RB", "IL3RA", "IL4R", "IRF4", "IRF6", "IRF8", "ITGA6", "ITGAV", "ITIH5", "KLF6", "LCLAT1", "LIF", "LRIG1", "LRRC8C", "LTB", "MAFF", "MAP3K8", "MAP6", "MAPKAPK2", "MUC1", "MXD1", "MYC", "MYO1C", "MYO1E", "NCOA3", "NCS1", "NDRG1", "NFIL3", "NFKBIZ", "NOP2", "NRP1", "NT5E", "ODC1", "P2RX4", "P4HA1", "PDCD2L", "PENK", "PHLDA1", "PHTF2", "PIM1", "PLAGL1", "PLEC", "PLIN2", "PLPP1", "PLSCR1", "PNP", "POU2F1", "PRAF2", "PRKCH", "PRNP", "PTCH1", "PTGER2", "PTH1R", "PTRH2", "PUS1", "RABGAP1L", "RGS16", "RHOB", "RHOH", "RNH1", "RORA", "RRAGD", "S100A1", "SCN9A", "SELL", "SELP", "SERPINB6", "SERPINC1", "SH3BGRL2", "SHE", "SLC1A5", "SLC29A2", "SLC2A3", "SLC39A8", "SMPDL3A", "SNX14", "SNX9", "SOCS1", "SOCS2", "SPP1", "SPRED2", "SPRY4", "ST3GAL4", "SWAP70", "SYNGR2", "SYT11", "TGM2", "TIAM1", "TLR7", "TNFRSF18", "TNFRSF1B", "TNFRSF21", "TNFRSF4", "TNFRSF8", "TNFRSF9", "TNFSF10", "TNFSF11", "TRAF1", "TTC39B", "TWSG1", "UCK2", "UMPS", "WLS", "XBP1")
IL6_JAK_STAT3_SIGNALING=c("A2M", "ACVR1B", "ACVRL1", "BAK1", "CBL", "CCL7", "CCR1", "CD14", "CD36", "CD38", "CD44", "CD9", "CNTFR", "CRLF2", "CSF1", "CSF2", "CSF2RA", "CSF2RB", "CSF3R", "CXCL1", "CXCL10", "CXCL11", "CXCL13", "CXCL3", "CXCL9", "DNTT", "EBI3", "FAS", "GRB2", "HAX1", "HMOX1", "IFNAR1", "IFNGR1", "IFNGR2", "IL10RB", "IL12RB1", "IL13RA1", "IL15RA", "IL17RA", "IL17RB", "IL18R1", "IL1B", "IL1R1", "IL1R2", "IL2RA", "IL2RG", "IL3RA", "IL4R", "IL6", "IL7", "IL9R", "INHBE", "IRF1", "IRF9", "ITGA4", "ITGB3", "JUN", "LEPR", "LTB", "LTBR", "MAP3K8", "MYD88", "OSMR", "PDGFC", "PF4", "PIK3R5", "PIM1", "PLA2G2A", "PTPN1", "PTPN11", "PTPN2", "REG1A", "SOCS1", "SOCS3", "STAM2", "STAT1", "STAT2", "STAT3", "TGFB1", "TLR2", "TNF", "TNFRSF12A", "TNFRSF1A", "TNFRSF1B", "TNFRSF21", "TYK2")
INFLAMMATORY_RESPONSE=c("ABCA1", "ABI1", "ACVR1B", "ACVR2A", "ADGRE1", "ADM", "ADORA2B", "ADRM1", "AHR", "APLNR", "AQP9", "ATP2A2", "ATP2B1", "ATP2C1", "AXL", "BDKRB1", "BEST1", "BST2", "BTG2", "C3AR1", "C5AR1", "CALCRL", "CCL17", "CCL2", "CCL20", "CCL22", "CCL24", "CCL5", "CCL7", "CCR7", "CCRL2", "CD14", "CD40", "CD48", "CD55", "CD69", "CD70", "CD82", "CDKN1A", "CHST2", "CLEC5A", "CMKLR1", "CSF1", "CSF3", "CSF3R", "CX3CL1", "CXCL10", "CXCL11", "CXCL6", "CXCL8", "CXCL9", "CXCR6", "CYBB", "DCBLD2", "EBI3", "EDN1", "EIF2AK2", "EMP3", "EREG", "F3", "FFAR2", "FPR1", "FZD5", "GABBR1", "GCH1", "GNA15", "GNAI3", "GPC3", "GPR132", "GPR183", "HAS2", "HBEGF", "HIF1A", "HPN", "HRH1", "ICAM1", "ICAM4", "IFITM1", "IFNAR1", "IFNGR2", "IL10", "IL10RA", "IL12B", "IL15", "IL15RA", "IL18", "IL18R1", "IL18RAP", "IL1A", "IL1B", "IL1R1", "IL2RB", "IL4R", "IL6", "IL7R", "INHBA", "IRAK2", "IRF1", "IRF7", "ITGA5", "ITGB3", "ITGB8", "KCNA3", "KCNJ2", "KCNMB2", "KIF1B", "KLF6", "LAMP3", "LCK", "LCP2", "LDLR", "LIF", "LPAR1", "LTA", "LY6E", "LYN", "MARCO", "MEFV", "MEP1A", "MET", "MMP14", "MSR1", "MXD1", "MYC", "NAMPT", "NDP", "NFKB1", "NFKBIA", "NLRP3", "NMI", "NMUR1", "NOD2", "NPFFR2", "OLR1", "OPRK1", "OSM", "OSMR", "P2RX4", "P2RX7", "P2RY2", "PCDH7", "PDE4B", "PDPN", "PIK3R5", "PLAUR", "PROK2", "PSEN1", "PTAFR", "PTGER2", "PTGER4", "PTGIR", "PTPRE", "PVR", "RAF1", "RASGRP1", "RELA", "RGS1", "RGS16", "RHOG", "RIPK2", "RNF144B", "ROS1", "RTP4", "SCARF1", "SCN1B", "SELE", "SELENOS", "SELL", "SEMA4D", "SERPINE1", "SGMS2", "SLAMF1", "SLC11A2", "SLC1A2", "SLC28A2", "SLC31A1", "SLC31A2", "SLC4A4", "SLC7A1", "SLC7A2", "SPHK1", "SRI", "STAB1", "TACR1", "TACR3", "TAPBP", "TIMP1", "TLR1", "TLR2", "TLR3", "TNFAIP6", "TNFRSF1B", "TNFRSF9", "TNFSF10", "TNFSF15", "TNFSF9", "TPBG", "VIP")
INTERFERON_ALPHA_RESPONSE=c("ADAR", "B2M", "BATF2", "BST2", "C1S", "CASP8", "CCRL2", "CD47", "CD74", "CMPK2", "CMTR1", "CNP", "CSF1", "CXCL10", "CXCL11", "DDX60", "DHX58", "EIF2AK2", "ELF1", "EPSTI1", "GBP2", "GBP4", "GMPR", "HELZ2", "HERC6", "IFI27", "IFI30", "IFI35", "IFI44", "IFI44L", "IFIH1", "IFIT2", "IFIT3", "IFITM1", "IFITM2", "IFITM3", "IL15", "IL4R", "IL7", "IRF1", "IRF2", "IRF7", "IRF9", "ISG15", "ISG20", "LAMP3", "LAP3", "LGALS3BP", "LPAR6", "LY6E", "MOV10", "MVB12A", "MX1", "NCOA7", "NMI", "NUB1", "OAS1", "OASL", "OGFR", "PARP12", "PARP14", "PARP9", "PLSCR1", "PROCR", "PSMA3", "PSMB8", "PSMB9", "PSME1", "PSME2", "RIPK2", "RNF31", "RSAD2", "RTP4", "SAMD9", "SAMD9L", "SELL", "SLC25A28", "SP110", "STAT2", "TAP1", "TDRD7", "TENT5A", "TMEM140", "TRAFD1", "TRIM14", "TRIM21", "TRIM25", "TRIM26", "TRIM5", "TXNIP", "UBA7", "UBE2L6", "USP18")
INTERFERON_GAMMA_RESPONSE=c("ADAR", "APOL6", "ARID5B", "ARL4A", "AUTS2", "B2M", "BANK1", "BATF2", "BPGM", "BST2", "BTG1", "C1R", "C1S", "CASP3", "CASP4", "CASP7", "CASP8", "CCL2", "CCL5", "CCL7", "CD274", "CD38", "CD40", "CD69", "CD74", "CD86", "CDKN1A", "CFB", "CFH", "CIITA", "CMKLR1", "CMPK2", "CMTR1", "CSF2RB", "CXCL10", "CXCL11", "CXCL9", "DDX58", "DDX60", "DHX58", "EIF2AK2", "EIF4E3", "EPSTI1", "FAS", "FCGR1A", "FGL2", "FPR1", "GBP4", "GBP6", "GCH1", "GPR18", "GZMA", "HELZ2", "HERC6", "HIF1A", "ICAM1", "IDO1", "IFI27", "IFI30", "IFI35", "IFI44", "IFI44L", "IFIH1", "IFIT1", "IFIT2", "IFIT3", "IFITM2", "IFITM3", "IFNAR2", "IL10RA", "IL15", "IL15RA", "IL18BP", "IL2RB", "IL4R", "IL6", "IL7", "IRF1", "IRF2", "IRF4", "IRF5", "IRF7", "IRF8", "IRF9", "ISG15", "ISG20", "ISOC1", "ITGB7", "JAK2", "KLRK1", "LAP3", "LATS2", "LCP2", "LGALS3BP", "LY6E", "LYSMD2", "METTL7B", "MT2A", "MTHFD2", "MVP", "MX1", "MX2", "MYD88", "NAMPT", "NCOA3", "NFKB1", "NFKBIA", "NLRC5", "NMI", "NOD1", "NUP93", "OAS2", "OAS3", "OASL", "OGFR", "P2RY14", "PARP12", "PARP14", "PDE4B", "PELI1", "PFKP", "PIM1", "PLA2G4A", "PLSCR1", "PML", "PNP", "PSMA2", "PSMA3", "PSMB10", "PSMB2", "PSMB8", "PSMB9", "PSME1", "PSME2", "PTGS2", "PTPN1", "PTPN2", "PTPN6", "RAPGEF6", "RBCK1", "RIPK1", "RIPK2", "RNF213", "RNF31", "RSAD2", "RTP4", "SAMD9L", "SAMHD1", "SECTM1", "SELP", "SERPING1", "SLAMF7", "SLC25A28", "SOCS1", "SOCS3", "SOD2", "SP110", "SPPL2A", "SRI", "SSPN", "ST3GAL5", "ST8SIA4", "STAT1", "STAT2", "STAT3", "STAT4", "TAP1", "TAPBP", "TDRD7", "TNFAIP2", "TNFAIP3", "TNFAIP6", "TNFSF10", "TOR1B", "TRAFD1", "TRIM14", "TRIM21", "TRIM25", "TRIM26", "TXNIP", "UBE2L6", "UPP1", "USP18", "VAMP5", "VAMP8", "VCAM1", "XAF1", "XCL1", "ZBP1", "ZNFX1")
NOTCH_SIGNALING=c("APH1A", "ARRB1", "CCND1", "CUL1", "DLL1", "DTX1", "DTX2", "DTX4", "FBXW11", "FZD1", "FZD5", "FZD7", "HES1", "HEYL", "JAG1", "KAT2A", "LFNG", "MAML2", "NOTCH1", "NOTCH2", "NOTCH3", "PPARD", "PRKCA", "PSEN2", "PSENEN", "RBX1", "SAP30", "ST3GAL6", "TCF7L2", "WNT2", "WNT5A")
OXIDATIVE_PHOSPHORYLATION=c("ABCB7", "ACAA1", "ACAA2", "ACADM", "ACADSB", "ACADVL", "ACAT1", "ACO2", "AFG3L2", "AIFM1", "ALAS1", "ALDH6A1", "ATP1B1", "ATP5F1A", "ATP5F1B", "ATP5F1C", "ATP5F1D", "ATP5F1E", "ATP5MC1", "ATP5MC2", "ATP5MC3", "ATP5ME", "ATP5MF", "ATP5MG", "ATP5PB", "ATP5PD", "ATP5PF", "ATP5PO", "ATP6AP1", "ATP6V0B", "ATP6V0C", "ATP6V0E1", "ATP6V1C1", "ATP6V1D", "ATP6V1E1", "ATP6V1F", "ATP6V1G1", "ATP6V1H", "BAX", "BCKDHA", "BDH2", "CASP7", "COX10", "COX11", "COX15", "COX17", "COX4I1", "COX5B", "COX6A1", "COX6B1", "COX6C", "COX7A2", "COX7A2L", "COX7B", "COX7C", "COX8A", "CPT1A", "CS", "CYB5A", "CYB5R3", "CYC1", "CYCS", "DECR1", "DLAT", "DLD", "DLST", "ECH1", "ECHS1", "ECI1", "ETFA", "ETFB", "ETFDH", "FDX1", "FH", "FXN", "GLUD1", "GOT2", "GPI", "GPX4", "GRPEL1", "HADHA", "HADHB", "HCCS", "HSD17B10", "HSPA9", "HTRA2", "IDH1", "IDH2", "IDH3A", "IDH3B", "IDH3G", "IMMT", "ISCA1", "ISCU", "LDHB", "LRPPRC", "MAOB", "MDH1", "MDH2", "MFN2", "MGST3", "MPC1", "MTRF1", "MTRR", "MTX2", "NDUFA1", "NDUFA2", "NDUFA3", "NDUFA4", "NDUFA5", "NDUFA6", "NDUFA7", "NDUFA8", "NDUFA9", "NDUFAB1", "NDUFB1", "NDUFB2", "NDUFB3", "NDUFB4", "NDUFB5", "NDUFB6", "NDUFB7", "NDUFB8", "NDUFC1", "NDUFC2", "NDUFS1", "NDUFS2", "NDUFS3", "NDUFS4", "NDUFS6", "NDUFS7", "NDUFS8", "NDUFV1", "NDUFV2", "NNT", "NQO2", "OAT", "OGDH", "OPA1", "OXA1L", "PDHA1", "PDHB", "PDHX", "PDK4", "PDP1", "PHB2", "PHYH", "PMPCA", "POLR2F", "POR", "PRDX3", "RETSAT", "RHOT1", "RHOT2", "SDHA", "SDHB", "SLC25A11", "SLC25A12", "SLC25A20", "SLC25A4", "SLC25A5", "SLC25A6", "SUCLA2", "SUCLG1", "SUPV3L1", "SURF1", "TCIRG1", "TIMM10", "TIMM13", "TIMM17A", "TIMM50", "TIMM8B", "TIMM9", "TOMM22", "TOMM70", "UQCR10", "UQCR11", "UQCRB", "UQCRC1", "UQCRC2", "UQCRFS1", "UQCRH", "UQCRQ", "VDAC1", "VDAC2", "VDAC3")
P53_PATHWAY=c("ABAT", "ABCC5", "ABHD4", "ACVR1B", "ADA", "AEN", "AK1", "ALOX15B", "ANKRA2", "APAF1", "APP", "ATF3", "BAIAP2", "BAK1", "BAX", "BLCAP", "BMP2", "BTG1", "BTG2", "CCND2", "CCND3", "CCNG1", "CCNK", "CCP110", "CD81", "CD82", "CDH13", "CDK5R1", "CDKN1A", "CDKN2A", "CDKN2AIP", "CDKN2B", "CEBPA", "CGRRF1", "CLCA2", "COQ8A", "CSRNP2", "CTSD", "CTSF", "CYFIP2", "DCXR", "DDB2", "DDIT3", "DDIT4", "DEF6", "DGKA", "DNTTIP2", "DRAM1", "EI24", "ELP1", "EPHA2", "EPHX1", "EPS8L2", "ERCC5", "F2R", "FAM162A", "FAS", "FBXW7", "FDXR", "FGF13", "FOS", "FOXO3", "FUCA1", "GADD45A", "GLS2", "GM2A", "GPX2", "HBEGF", "HDAC3", "HEXIM1", "HINT1", "HMOX1", "HRAS", "HSPA4L", "IER3", "IER5", "IFI30", "IL1A", "INHBB", "IP6K2", "IRAK1", "ISCU", "ITGB4", "JAG2", "JUN", "KIF13B", "KLF4", "KLK8", "KRT17", "LDHB", "LIF", "MAPKAPK3", "MDM2", "MKNK2", "MXD1", "MXD4", "NDRG1", "NHLH2", "NINJ1", "NOL8", "NOTCH1", "NUDT15", "NUPR1", "OSGIN1", "PCNA", "PDGFA", "PERP", "PHLDA3", "PIDD1", "PITPNC1", "PLK2", "PLK3", "PLXNB2", "PMM1", "POLH", "POM121", "PPM1D", "PPP1R15A", "PRKAB1", "PRMT2", "PROCR", "PTPN14", "PTPRE", "RAB40C", "RACK1", "RAD51C", "RAD9A", "RALGDS", "RAP2B", "RB1", "RCHY1", "RETSAT", "RGS16", "RHBDF2", "RNF19B", "RRAD", "RRP8", "RXRA", "S100A10", "S100A4", "SAT1", "SDC1", "SEC61A1", "SERPINB5", "SERTAD3", "SESN1", "SFN", "SLC19A2", "SLC35D1", "SLC3A2", "SLC7A11", "SOCS1", "SP1", "SPHK1", "ST14", "STEAP3", "STOM", "TAP1", "TAX1BP3", "TCHH", "TCN2", "TGFA", "TGFB1", "TM4SF1", "TM7SF3", "TNFSF9", "TNNI1", "TOB1", "TP53", "TP63", "TPD52L1", "TPRKB", "TRAF4", "TRAFD1", "TRIAP1", "TRIB3", "TSC22D1", "TSPYL2", "TXNIP", "UPP1", "VAMP8", "VDR", "VWA5A", "WRAP73", "WWP1", "XPC", "ZBTB16", "ZFP36L1", "ZMAT3", "ZNF365")
PEROXISOME=c("ABCB1", "ABCB4", "ABCB9", "ABCC5", "ABCC8", "ABCD1", "ABCD2", "ABCD3", "ACAA1", "ACOT8", "ACOX1", "ACSL1", "ACSL4", "ACSL5", "ALB", "ALDH1A1", "ALDH9A1", "ATXN1", "BCL10", "CACNA1B", "CADM1", "CAT", "CDK7", "CEL", "CLN6", "CLN8", "CNBP", "CRABP1", "CRABP2", "CRAT", "CTPS1", "DHCR24", "DHRS3", "DIO1", "DLG4", "ECH1", "ECI2", "EHHADH", "ELOVL5", "EPHX2", "ERCC1", "ERCC3", "FABP6", "FADS1", "FDPS", "FIS1", "GNPAT", "GSTK1", "HAO2", "HMGCL", "HRAS", "HSD11B2", "HSD17B11", "HSD17B4", "HSD3B7", "IDE", "IDH1", "IDH2", "IDI1", "ISOC1", "ITGB1BP1", "LONP2", "MLYCD", "MSH2", "MVP", "NR1I2", "NUDT19", "PEX11A", "PEX11B", "PEX13", "PEX14", "PEX2", "PEX5", "PEX6", "PRDX1", "PRDX5", "RDH11", "RETSAT", "RXRG", "SCGB1A1", "SCP2", "SEMA3C", "SERPINA6", "SIAH1", "SLC23A2", "SLC25A17", "SLC25A19", "SLC25A4", "SLC27A2", "SLC35B2", "SMARCC1", "SOD1", "SOD2", "STS", "SULT2B1", "TOP2A", "TSPO", "TTR", "UGT2B17", "VPS4B", "YWHAH")
PI3K_AKT_MTOR_SIGNALING=c("ACACA", "ACTR2", "ACTR3", "ADCY2", "AKT1", "AKT1S1", "AP2M1", "ARF1", "ARHGDIA", "ARPC3", "ATF1", "CAB39", "CAB39L", "CALR", "CAMK4", "CDK1", "CDK2", "CDK4", "CDKN1A", "CDKN1B", "CFL1", "CLTC", "CSNK2B", "CXCR4", "DAPP1", "DDIT3", "DUSP3", "E2F1", "ECSIT", "EGFR", "EIF4E", "FASLG", "FGF17", "FGF22", "FGF6", "GNA14", "GNGT1", "GRB2", "GRK2", "GSK3B", "HRAS", "HSP90B1", "IL2RG", "IL4", "IRAK4", "ITPR2", "LCK", "MAP2K3", "MAP2K6", "MAP3K7", "MAPK1", "MAPK10", "MAPK8", "MAPK9", "MAPKAP1", "MKNK1", "MKNK2", "MYD88", "NCK1", "NFKBIB", "NGF", "NOD1", "PAK4", "PDK1", "PFN1", "PIK3R3", "PIKFYVE", "PIN1", "PITX2", "PLA2G12A", "PLCB1", "PLCG1", "PPP1CA", "PPP2R1B", "PRKAA2", "PRKAG1", "PRKAR2A", "PRKCB", "PTEN", "PTPN11", "RAC1", "RAF1", "RALB", "RIPK1", "RIT1", "RPTOR", "SFN", "SLA", "SLC2A1", "SMAD2", "SQSTM1", "STAT2", "TBK1", "THEM4", "TIAM1", "TNFRSF1A", "TRAF2", "TRIB3", "TSC2", "UBE2D3", "UBE2N", "VAV3", "YWHAB")
PROTEIN_SECRETION=c("ABCA1", "ADAM10", "ANP32E", "AP1G1", "AP2B1", "AP2M1", "AP2S1", "AP3B1", "ARCN1", "ARF1", "ARFGAP3", "ARFGEF1", "ARFGEF2", "ARFIP1", "ATP1A1", "ATP6V1B1", "ATP6V1H", "ATP7A", "BET1", "CAV2", "CD63", "CLCN3", "CLN5", "CLTA", "CLTC", "COG2", "COPB1", "COPB2", "COPE", "CTSC", "DNM1L", "DOP1A", "DST", "EGFR", "ERGIC3", "GALC", "GBF1", "GLA", "GOSR2", "ICA1", "IGF2R", "KIF1B", "LAMP2", "LMAN1", "M6PR", "MAPK1", "MON2", "NAPA", "NAPG", "OCRL", "PAM", "PPT1", "RAB14", "RAB22A", "RAB2A", "RAB5A", "RAB9A", "RER1", "SCAMP1", "SCAMP3", "SCRN1", "SEC22B", "SEC24D", "SEC31A", "SGMS1", "SH3GL2", "SNAP23", "SNX2", "SOD1", "SSPN", "STAM", "STX12", "STX16", "STX7", "TMED10", "TMED2", "TMX1", "TOM1L1", "TPD52", "TSG101", "TSPAN8", "USO1", "VAMP3", "VAMP4", "VAMP7", "VPS45", "VPS4B", "YIPF6", "YKT6", "ZW10")
REACTIVE_OXYGEN_SPECIES_PATHWAY=c("ABCC1", "ATOX1", "CAT", "CDKN2D", "EGLN2", "ERCC2", "FES", "FTL", "G6PD", "GCLC", "GCLM", "GLRX", "GLRX2", "GPX3", "GPX4", "GSR", "HHEX", "HMOX2", "IPCEF1", "JUNB", "LAMTOR5", "LSP1", "MBP", "MGST1", "MPO", "MSRA", "NDUFA6", "NDUFB4", "NDUFS2", "NQO1", "OXSR1", "PDLIM1", "PFKP", "PRDX1", "PRDX2", "PRDX4", "PRDX6", "PRNP", "PTPA", "SBNO2", "SCAF4", "SELENOS", "SOD1", "SOD2", "SRXN1", "STK25", "TXNRD1", "TXNRD2")
TGF_BETA_SIGNALING=c("ACVR1", "APC", "ARID4B", "BCAR3", "BMP2", "BMPR1A", "BMPR2", "CDH1", "CDK9", "CDKN1C", "CTNNB1", "ENG", "FNTA", "FURIN", "HDAC1", "HIPK2", "ID1", "ID2", "ID3", "IFNGR2", "JUNB", "KLF10", "LEFTY2", "LTBP2", "MAP3K7", "NCOR2", "NOG", "PMEPA1", "PPM1A", "PPP1CA", "PPP1R15A", "RAB31", "RHOA", "SERPINE1", "SKI", "SKIL", "SLC20A1", "SMAD1", "SMAD3", "SMAD6", "SMAD7", "SMURF1", "SMURF2", "SPTBN1", "TGFB1", "TGFBR1", "TGIF1", "THBS1", "TJP1", "TRIM33", "UBE2D3", "WWTR1", "XIAP")
TNFA_SIGNALING_VIA_NFKB=c("ABCA1", "ACKR3", "AREG", "ATF3", "ATP2B1", "B4GALT1", "B4GALT5", "BCL2A1", "BCL3", "BCL6", "BHLHE40", "BIRC2", "BIRC3", "BMP2", "BTG1", "BTG2", "BTG3", "CCL2", "CCL20", "CCL4", "CCL5", "CCN1", "CCND1", "CCNL1", "CCRL2", "CD44", "CD69", "CD80", "CD83", "CDKN1A", "CEBPB", "CEBPD", "CFLAR", "CLCF1", "CSF1", "CSF2", "CXCL1", "CXCL10", "CXCL11", "CXCL2", "CXCL3", "CXCL6", "DDX58", "DENND5A", "DNAJB4", "DRAM1", "DUSP1", "DUSP2", "DUSP4", "DUSP5", "EDN1", "EFNA1", "EGR1", "EGR2", "EGR3", "EHD1", "EIF1", "ETS2", "F2RL1", "F3", "FJX1", "FOS", "FOSB", "FOSL1", "FOSL2", "FUT4", "G0S2", "GADD45A", "GADD45B", "GCH1", "GEM", "GFPT2", "GPR183", "HBEGF", "HES1", "ICAM1", "ID2", "IER2", "IER3", "IER5", "IFIH1", "IFIT2", "IFNGR2", "IL12B", "IL15RA", "IL18", "IL1A", "IL1B", "IL23A", "IL6", "IL7R", "INHBA", "IRF1", "IRS2", "JAG1", "JUN", "JUNB", "KDM6B", "KLF10", "KLF2", "KLF4", "KLF6", "KLF9", "KYNU", "LAMB3", "LDLR", "LIF", "LITAF", "MAFF", "MAP2K3", "MAP3K8", "MARCKS", "MCL1", "MSC", "MXD1", "MYC", "NAMPT", "NFAT5", "NFE2L2", "NFIL3", "NFKB1", "NFKB2", "NFKBIA", "NFKBIE", "NINJ1", "NR4A1", "NR4A2", "NR4A3", "OLR1", "PANX1", "PDE4B", "PDLIM5", "PER1", "PFKFB3", "PHLDA1", "PHLDA2", "PLAU", "PLAUR", "PLEK", "PLK2", "PLPP3", "PMEPA1", "PNRC1", "PPP1R15A", "PTGER4", "PTGS2", "PTPRE", "PTX3", "RCAN1", "REL", "RELA", "RELB", "RHOB", "RIPK2", "RNF19B", "SAT1", "SDC4", "SERPINB2", "SERPINB8", "SERPINE1", "SGK1", "SLC16A6", "SLC2A3", "SLC2A6", "SMAD3", "SNN", "SOCS3", "SOD2", "SPHK1", "SPSB1", "SQSTM1", "STAT5A", "TANK", "TAP1", "TGIF1", "TIPARP", "TLR2", "TNC", "TNF", "TNFAIP2", "TNFAIP3", "TNFAIP6", "TNFAIP8", "TNFRSF9", "TNFSF9", "TNIP1", "TNIP2", "TRAF1", "TRIB1", "TRIP10", "TSC22D1", "TUBB2A", "VEGFA", "YRDC", "ZBTB10", "ZC3H12A", "ZFP36")
WNT_BETA_CATENIN_SIGNALING=c("ADAM17", "AXIN1", "AXIN2", "CCND2", "CSNK1E", "CTNNB1", "CUL1", "DKK1", "DKK4", "DLL1", "DVL2", "FRAT1", "FZD1", "FZD8", "GNAI1", "HDAC11", "HDAC2", "HDAC5", "HEY1", "HEY2", "JAG1", "JAG2", "KAT2A", "LEF1", "MAML1", "MYC", "NCOR2", "NCSTN", "NKD1", "NOTCH1", "NUMB", "PPARD", "PSEN2", "PTCH1", "RBPJ", "SKP2", "TCF7", "TP53", "WNT1", "WNT5B", "WNT6")

#setwd("/Users/kriniogiannikou/Downloads/Ernesto-spatial-10X-partA/Ernesto-spatial-10X-partB")

SCT04manual_module <- AddModuleScore_UCell(SCT04manual,
                                           features = list(WNT_BETA_CATENIN_SIGNALING),  # wrap in list!
                                           name = "WNT_BETA_CATENIN_SIGNALING"
)
head(  SCT04manual_module)
colnames(  SCT04manual_module@meta.data)

pdf("SCT04manual_UCell_plot-WNT_BETA_CATENIN_SIGNALING.pdf", width = 14, height = 9)

MapFeatures(
  object =   SCT04manual_module,
  features = "signature_1WNT_BETA_CATENIN_SIGNALING",
  pt_size = 2.7,
  pt_alpha = 0.7,
  colors = viridis(n = 9),
  scale = "shared",
  add_scalebar = TRUE,
  scalebar_position = c(0.8, 0.8),
  scalebar_height = 0.05
)

dev.off()

# Define your genes of interest
Cycling= c("MKI67",  "TOP2A", "CDK1")

# Check if these genes are present in the Seurat object
genes_present <- Cycling[Cycling %in% rownames(SCT04manual)]

# Create a DotPlot for the selected genes across all clusters
dot_plot <- DotPlot(SCT04manual, features = genes_present) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ggtitle("Expression of Cycling Genes Across Clusters") +
  scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0.5) # Adjust the color scale

# Print the plot
print(dot_plot)

Neuroectoderm=c("TMEM178B", "PPP2R2B", "DPP10", "CTNNA2", "DPP6", "NRXN1", "DOCK3", "FMN2", "CADPS", "DTNA")
Epithelia=c("LIPH", "TMC5", "ELF3", "GRHL2", "MACC1", "KRT19", "BAIAP2L1", "LYPD6B", "ABCC3", "CDH1")
Stroma=c("COL6A3", "COL5A2", "PRR16", "PRRX1", "COL1A2", "COL5A1", "SLIT3", "COL1A1", "CDH11", "EBF1")
Endothelia=c("VWF", "EGFL7", "DIPK2B", "CDH5", "ADGRL4", "FLT1", "PECAM1", "AFAP1L1", "SHANK3", "ERG")
Immune=c("CSF2RA", "IKZF1", "PIK3R5", "SAMSN1", "DOCK2", "PTPRC", "APBB1IP", "LCP2", "SRGN", "ARHGAP15")
Cycling= c("MKI67",  "TOP2A", "CDK1")


if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("UCell")
library(UCell)

# Create a gene set list with your positive and negative genes
gene_set <- list(Positive_Corr_XAratio = top_positive_genes,
                 Negative_Corr_XAratio = top_negative_genes)

# Add UCell module scores to your Seurat object
# Adjust the 'assay' and 'slot' parameters if your data is stored differently
SCT04manual_UCells <- AddModuleScore_UCell(se_SCT04manual,
                                           features = gene_set)

# The function will add new metadata columns, such as 'UCellScore_Positive' and 'UCellScore_Negative'

FeatureScatter_scCustom(seurat_object = nucSamples_female, feature1 = "X_to_Autosome_Ratio", feature2 = "Positive_Corr_XAratio_UCell", group.by = "BroadCategory", colors_use = pal_broad)
FeatureScatter_scCustom(seurat_object = nucSamples_female, feature1 = "X_to_Autosome_Ratio", feature2 = "Negative_Corr_XAratio_UCell", group.by = "BroadCategory", colors_use = pal_broad)

nucSamples_female <- AddModuleScore_UCell(
  obj      = sobj, # Be sure to change object here
  features = Hallmark_to_Test, # Just use the list that we made above to add the 
  name     = NULL         # using the list name directly (“XIST_Autosomal”) :contentReference[oaicite:0]{index=0}
)
#You can use the same Hallmark_to_Test from the R file attached (SpatialGeneSetsForPlotting_20250509.R)



pdf("SCT04manual_Neuroectoderm_FeaturePlot.pdf", width = 8, height = 6)
FeaturePlot(SCT04manual, dims = c(1, 2), features = Neuroectoderm, cols = c("lightgrey", "lightblue", "blue"), pt.size = 0.5)
dev.off()

DoHeatmap(SCT04manual, features = Neuroectoderm, group.colors = c("light green", "lightgray"), disp.min = -2, disp.max = 2)
dev.off()

DoHeatmap(SCT04manual, features = Cycling, group.colors = c("#AA4499", "#88CCEE", "#20B2AA"), disp.min = -2, disp.max = 2)


SCT04manual_module <- AddModuleScore(
  object = SCT04manual,
  features = top_100_positive_corr_to_XAratio_genes,
  rownames(x=SCT04manual),
  ctrl = 100,
  name = 'top_100_positive_corr_to_XAratio_genes_genes')
head(x = SCT04manual_module[])

FeaturePlot(SCT04manual, dims = c(1, 2), features = c("ABCC3", "CDH1"), cols = c("lightgrey", "lightblue", "blue"), pt.size = 0.5)


# Visualizing the module score (CAF_Features1) using MapFeatures
MapFeatures(
  object = SCT04manual_module,                  # Seurat object with module scores
  features = "top_100_positive_corr_to_XAratio_genes_genes1",                    # The name of the module score feature
  pt_size = 1.8,                                   # Point size (adjust based on your dataset)
  pt_alpha = 0.7,                                # Adjust opacity if needed
  colors = viridis(n = 9),                      # Color palette for visualization
  scale = "shared",                              # Scale colors for the entire dataset
  add_scalebar = TRUE,                           # Option to add a scale bar
  scalebar_position = c(0.8, 0.8),               # Position of scale bar
  scalebar_height = 0.05                         # Height of the scale bar relative to plot
)

FeaturePlot(SCT04manual_module, dims = c(1, 2), features = c("top_100_negative_corr_to_XAratio_genes_genes1"), cols = c("lightgrey", "lightblue", "blue"), pt.size = 1)

spatial_data_SCT04manual <- GetStaffli(SCT04manual)

png("SCT04manual_plot.png", width=7, height=8, units="in", res=300)
MapFeatures(SCT04manual, features = "top_100_positive_corr_to_XAratio_genes", cols = c("lightgrey", "lightblue", "blue"), pt.size = 1)
MapFeatures(SCT04manual, arrange_features = "row", features = c("TMEM178B"), subplot_type = "violin", colors = RColorBrewer::brewer.pal(n = 11, name = "RdBu"), pt.size = 3 |> rev())

pdf("SCT04manual_cycling markers-vs3.pdf", width = 14, height = 10)
#Map multiple numeric features
#With the MapMultipleFeatures() function, it is possible to plot several numeric features within the same plot.
MapMultipleFeatures(SCT04manual, 
                    features = c("MKI67", "CDK1", "TOP2A"),
                    colors = c("#FF6347", "#88CCEE", "#4DAF4A"),
                    ncol = 1.5, 
                    pt_size = 2)
dev.off()

SCT04manual_semla <- UpdateSeuratForSemla(SCT04manual)
SCT04manual <- LoadImages(SCT04manual)


# Plot the spatial data with ggplot2
ggplot(SCT04manual) +
  geom_sf(aes(fill = SCT04manual_module)) +  # Using module_score to fill color
  scale_fill_viridis_c() +  # Apply a color scale (viridis)
  theme_minimal() +         # Clean theme
  ggtitle("Spatial Heatmap of Module Scores")  # Title


# Set seed for reproducibility
set.seed(42)
SCT04manual <- RunNMF(SCT04manual)


###Spatial statistics for labelled spots
unique_clusters <- levels(SCT04manual$seurat_clusters)
cluster_colors <- setNames(ggsci::pal_d3("category20")(length(unique_clusters)), nm = unique_clusters)

pdf("SCT04manual_plot-spots-cluster.pdf", width = 16, height = 10)
MapLabelsSummary(SCT04manual, 
                 column_name = "seurat_clusters", 
                 bar_display = "count",
                 pt_size = 3.8,
                 override_plot_dims = F, 
                 bar_label_size = 2.5, 
                 colors = cluster_colors)
theme(plot.title = element_text(hjust=1, size = 14, face = "bold"))

dev.off()

spatial_network_SCT04manual <- GetSpatialNetwork(SCT04manual)[[1]]
SCT04manual$random_clusters <- SCT04manual$seurat_clusters |> as.character() |> sample(ncol(SCT04manual))
spatial_network_SCT04manual$cluster <- SCT04manual[[]][spatial_network_SCT04manual$from, "seurat_clusters", drop = TRUE] |> as.character()
spatial_network_SCT04manual$random_clusters <- SCT04manual[[]][spatial_network_SCT04manual$from, "random_clusters", drop = TRUE] |> as.character()

head(spatial_network, 4)

ggplot(spatial_network_SCT04manual, aes(x = x_start, xend = x_end, y = y_start, yend = y_end, color=cluster)) +
  geom_segment() +
  geom_point() +
  scale_color_manual(values = cluster_colors) +
  labs(title="Spatial network") +
  scale_y_reverse() + 
  theme_void() +
  theme(aspect.ratio = 1, plot.title = element_text(hjust=1))

res_lat_SCT04manual <- RunLabelAssortativityTest(SCT04manual, column_name = "seurat_clusters")

datatable(
  res_lat_SCT04manual |> arrange(desc(avg_k_scaled)),
  rownames = F,
  caption = "Label assortativity results"
)MapFeatures(SCT04manual, features = c("ABCC3", "CDH1"), image_use = "raw", pt_size = 1.8,override_plot_dims = TRUE, colors = RColorBrewer::brewer.pal(n = 11, name = "Spectral") |> rev())

ggplot(res_lat_SCT04manual, aes(x=reorder(label, avg_k_scaled), y=avg_k_scaled)) +
  geom_col(width = 0.8) +
  geom_hline(yintercept = 0, color="black", size=0.2) +
  labs(x="Label", y="Scaled average degree", title="Label assortativity score") +
  coord_flip() +
  theme_linedraw() +
  theme(plot.title = element_text(hjust=0.5, size=12, face = "bold"),
        axis.text = element_text(size=10),
        legend.position = "none",
        legend.title = element_text(size=10),
        legend.text = element_text(size=10))

p1 <- ggplot(spatial_network_SCT04manual, aes(x = x_start, xend = x_end, y = y_start, yend = y_end)) +
  geom_segment(color="grey90") +
  geom_point(data = subset(spatial_network_SCT04manual, cluster == "0"), color = cluster_colors["0"]) +
  labs(title="Cluster 0",
       caption = paste0("Observed avg k: ", round(subset(res_lat_SCT04manual, label=="1")[,"avg_k"], 2))
  )

p2 <- ggplot(spatial_network_SCT04manual, aes(x = x_start, xend = x_end, y = y_start, yend = y_end)) +
  geom_segment(color="grey90") +
  geom_point(data = subset(spatial_network_SCT04manual, random_clusters == "0"), color = cluster_colors["0"]) +
  labs(title="Cluster 0 randomly dispersed",
       caption = paste0("Random mean avg k: ", round(subset(res_lat_SCT04manual, label=="0")[,"min_avg_k_mean"], 2))
  )

(p1 | p2) & scale_y_reverse() & theme_void() & theme(aspect.ratio = 1, plot.title = element_text(hjust=0.5), plot.caption = element_text(hjust=0.5))


p1 <- ggplot(spatial_network_SCT04manual, aes(x = x_start, xend = x_end, y = y_start, yend = y_end)) +
  geom_segment(color="grey90") +
  geom_point(data = subset(spatial_network_SCT04manual, cluster == "0"), color = cluster_colors["0"]) +
  labs(title="Cluster 0",
       caption = paste0("Observed avg k: ", round(subset(res_lat_SCT04manual, label=="7")[,"avg_k"], 2))
  )

p2 <- ggplot(spatial_network_SCT04manual, aes(x = x_start, xend = x_end, y = y_start, yend = y_end)) +
  geom_segment(color="grey90") +
  geom_point(data = subset(spatial_network_SCT04manual, random_clusters == "0"), color = cluster_colors["0"]) +
  labs(title="Cluster 0 randomly dispersed",
       caption = paste0("Random mean avg k: ", round(subset(res_lat_SCT04manual, label=="0")[,"min_avg_k_mean"], 2))
  )

(p1 | p2) & scale_y_reverse() & theme_void() & theme(aspect.ratio = 1, plot.title = element_text(hjust=0.5), plot.caption = element_text(hjust=0.5))


perm_test_SCT04manual <- c(10, 20, 50, 100, seq(from=200, to=1000, by=200))
lat_perm_test_list_SCT04manual <- lapply(perm_test_SCT04manual, function(n_perm){
  f_time <- system.time(result <- RunLabelAssortativityTest(object = SCT04manual, 
                                                            column_name = "seurat_clusters", 
                                                            n_permutations = n_perm, 
                                                            verbose = F))
  result$time_s <- f_time[3]
  result$n_perm <- n_perm
  result
})
lat_perm_test_SCT04manual <- do.call(bind_rows, lat_perm_test_list_SCT04manual)


ggplot(lat_perm_test_SCT04manual, aes(x=n_perm, y = min_avg_k_mean)) +
  geom_point(size=1) +
  geom_errorbar(aes(ymin = min_avg_k_mean - min_avg_k_sd, ymax = min_avg_k_mean + min_avg_k_sd)) +
  geom_line(data = lat_perm_test_SCT04manual, mapping = aes(x=n_perm, y = avg_k), color = "orange", linewidth=0.5) +
  facet_wrap(~label, scales = "free_y") +
  labs(title="Effect of n_permutations on average degree in randomized network", 
       caption = "Black dots and error bars: average degree mean and s.d from permutations. Orange line: actual observed average degree.",
       y="Average degree", x="n permutations") +
  theme_bw()


mean(lat_perm_test_SCT04manual$min_avg_k_sd)

t_plot_SCT04manual <- lat_perm_test_SCT04manual
t_plot_SCT04manual <- t_plot_SCT04manual |> 
  group_by(n_perm) |> 
  summarise(.groups = "keep", time_s= max(time_s))

ggplot(t_plot_SCT04manual, aes(x=n_perm, y = time_s)) +
  geom_line(size=.5) +
  geom_point(size=2) +
  labs(title="Function run time with increasing n permutations",
       y="time (s)", x="n permutations") +
  theme_bw()

#Neighborhood Enrichment Analysis
res_net_SCT04manual <- RunNeighborhoodEnrichmentTest(SCT04manual, column_name = "seurat_clusters", n_permutations = 1000) 

datatable(
  res_net_SCT04manual |> arrange(desc(abs(z_score))), 
  rownames = F, 
  caption = "Neighborhood enrichment analysis output"
)

cluster_labels_SCT04manual <- paste0("Label_", sort(unique(SCT04manual$seurat_clusters)))
hm_plot_data <- res_net_SCT04manual

hm_plot_data$label_1 <- factor(hm_plot_data$label_1, levels = cluster_labels_SCT04manual)
hm_plot_data$label_2 <- factor(hm_plot_data$label_2, levels = cluster_labels_SCT04manual)

ggplot(hm_plot_data, aes(label_1, label_2, fill= z_score)) + 
  geom_tile(color = "white",
            lwd = 1.5,
            linetype = 1) +
  scale_fill_gradient2(low = "#0474BA",
                       mid = "grey90",
                       high = "#F17720") +
  labs(x="", y="", title="Neighborhood enrichment", fill="Z-score") +
  theme_bw() +
  coord_fixed() +
  theme(panel.grid.major = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(hjust=0.5, size=12, face = "bold"),
        axis.text = element_text(size=10),
        legend.title = element_text(size=10), 
        legend.text = element_text(size=10))


top_n_plot_SCT04manual <- 10
plot_data_SCT04manual <- res_net_SCT04manual |> 
  mutate(across(where(is.factor), as.character)) |> 
  group_by(grp = paste0(pmin(label_1, label_2), "-", pmax(label_1, label_2))) |> 
  slice(1)  |> 
  ungroup()  |> 
  select(-grp)

plot_data_SCT04manual <- rbind(
  plot_data_SCT04manual |> arrange(z_score) |> head(top_n_plot_SCT04manual) |> filter(z_score < 0),
  plot_data_SCT04manual |> arrange(desc(z_score)) |> head(top_n_plot_SCT04manual) |> filter(z_score > 0)
)

plot_data_SCT04manual$direction <- ifelse(plot_data_SCT04manual$z_score>0, "over-represented", "under-represented")
colors_direction_fill <- setNames(c("#F17720", "#0474BA"), nm = c("over-represented", "under-represented"))

ggplot(plot_data_SCT04manual, aes(x=reorder(label_label, z_score), y=z_score, fill = direction)) +
  geom_col(width = 0.6) +
  labs(x="", y="Z-score", title="Top enriched label pairs", fill="Enrichment") +
  scale_fill_manual(values = colors_direction_fill) +
  geom_hline(yintercept = 0, color="black", size=0.5) +
  coord_flip() +
  theme_linedraw() +
  theme(plot.title = element_text(hjust=0.5, size=12, face = "bold"),
        axis.text = element_text(size=10),
        panel.grid = element_blank(), 
        legend.position = "top", 
        legend.title = element_text(size=10), 
        legend.text = element_text(size=10))


p1 <- ggplot(spatial_network_SCT04manual, aes(x = x_start, xend = x_end, y = y_start, yend = y_end)) +
  geom_segment(color="grey90") +
  geom_point(data = subset(spatial_network_SCT04manual, cluster == "1"), color = cluster_colors["1"]) +
  geom_point(data = subset(spatial_network_SCT04manual, cluster == "12"), color = cluster_colors["12"]) +
  labs(title="Clusters 1 and 12")

p2 <- ggplot(spatial_network_SCT04manual, aes(x = x_start, xend = x_end, y = y_start, yend = y_end)) +
  geom_segment(color="grey90") +
  geom_point(data = subset(spatial_network_SCT04manual, cluster == "3"), color = cluster_colors["3"]) +
  geom_point(data = subset(spatial_network_SCT04manual, cluster == "4"), color = cluster_colors["4"]) +
  labs(title="Clusters 3 and 4")

(p1 | p2) & scale_y_reverse() & theme_void() & theme(aspect.ratio = 1, plot.title = element_text(hjust=0.5))

p1 <- ggplot(spatial_network_SCT04manual, aes(x = x_start, xend = x_end, y = y_start, yend = y_end)) +
  geom_segment(color="grey90") +
  geom_point(data = subset(spatial_network_SCT04manual, cluster == "0"), color = cluster_colors["0"]) +
  geom_point(data = subset(spatial_network_SCT04manual, cluster == "1"), color = cluster_colors["1"]) +
  labs(title="Clusters 0 and 1")

p2 <- ggplot(spatial_network_SCT04manual, aes(x = x_start, xend = x_end, y = y_start, yend = y_end)) +
  geom_segment(color="grey90") +
  geom_point(data = subset(spatial_network_SCT04manual, cluster == "1"), color = cluster_colors["1"]) +
  geom_point(data = subset(spatial_network_SCT04manual, cluster == "2"), color = cluster_colors["2"]) +
  labs(title="Clusters 1 and 2")

(p1 | p2) & scale_y_reverse() & theme_void() & theme(aspect.ratio = 1, plot.title = element_text(hjust=0.5))


#Explore parameter settings: n_permutations
perm_test_SCT04manual <- c(10, 20, 50, 100, seq(from=200, to=1000, by=200))
nea_perm_test_list_SCT04manual <- lapply(perm_test_SCT04manual, function(n_perm){
  f_time <- system.time(result <- RunNeighborhoodEnrichmentTest(object = SCT04manual, 
                                                                column_name = "seurat_clusters", 
                                                                n_permutations = n_perm, 
                                                                verbose = F))
  result$time_s <- f_time[3]
  result$n_perm <- n_perm
  result
})
nea_perm_test_SCT04manual <- do.call(bind_rows, nea_perm_test_list_SCT04manual)

ggplot(nea_perm_test_SCT04manual, aes(x=n_perm, y = z_score, color=label_label)) +
  geom_point(size=1) +
  geom_line() +
  scale_x_continuous(breaks = perm_test_SCT04manual) +
  labs(title="Effect of n_permutations on z-score", 
       caption = "Each line corresponds to the z-score of a label pair",
       y="Z-score", x="n permutations") +
  theme_bw() +
  theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1))

######################
t_plot_SCT04manual <- nea_perm_test_SCT04manual
t_plot_SCT04manual <- t_plot_SCT04manual |> 
  group_by(n_perm) |> 
  summarise(.groups = "keep", time_s= max(time_s))

ggplot(t_plot_SCT04manual, aes(x=n_perm, y = time_s)) +
  geom_line(size=.5) +
  geom_point(size=2) +
  scale_x_continuous(breaks = perm_test_SCT04manual) +
  labs(title="Function run time with increasing n permutations",
       y="time (s)", x="n permutations") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

#Explore parameter settings: column_labels
cluster_selection_SCT04manual <- c("1", "12", "3", "4")
res_net_medulla_SCT04manual <- RunNeighborhoodEnrichmentTest(SCT04manual, 
                                                             column_name = "seurat_clusters", 
                                                             column_labels = cluster_selection_SCT04manual, 
                                                             n_permutations = 1000)

datatable(
  res_net_medulla_SCT04manual |> arrange(desc(abs(z_score))), 
  rownames = F, 
  caption = "Neighborhood enrichment analysis output for label subset"
)


cluster_labels_SCT04manual <- paste0("Label_", cluster_selection_SCT04manual)
hm_plot_data_SCT04manual <- res_net_medulla_SCT04manual

hm_plot_data$label_1 <- factor(hm_plot_data$label_1, levels = cluster_labels_SCT04manual)
hm_plot_data$label_2 <- factor(hm_plot_data$label_2, levels = cluster_labels_SCT04manual)

ggplot(hm_plot_data_SCT04manual, aes(label_1, label_2, fill= z_score)) + 
  geom_tile(color = "white",
            lwd = 1.5,
            linetype = 1) +
  scale_fill_gradient2(low = "#0474BA",
                       mid = "grey90",
                       high = "#F17720") +
  labs(x="", y="", title="Neighborhood enrichment", fill="Z-score") +
  theme_bw() +
  coord_fixed() +
  theme(panel.grid.major = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(hjust=0.5, size=12, face = "bold"),
        axis.text = element_text(size=10),
        legend.title = element_text(size=10), 
        legend.text = element_text(size=10))


ggplot(subset(spatial_network_SCT04manual, cluster %in% cluster_selection_SCT04manual), aes(x = x_start, xend = x_end, y = y_start, yend = y_end)) +
  geom_segment(color="grey90") +
  geom_point(mapping = aes(x = x_start, y = y_start, color = cluster)) +
  scale_color_manual(values = cluster_colors[cluster_selection_SCT04manual]) +
  labs(title="SCT associated clusters") & 
  scale_y_reverse() & 
  theme_void() & 
  theme(aspect.ratio = 1, plot.title = element_text(hjust=0.5))



#.....
PlotFeatureLoadings(SCT04manual, 
                    dims = selected_factors, 
                    reduction = "nmf", 
                    nfeatures = 30, 
                    mode = "heatmap", 
                    gradient_colors = viridis::magma(n = 11, direction = -1))

MapLabels(SCT04manual, column_name = "seurat_clusters", ncol = 1) &
  theme(legend.position = "right")

MapFeatures(SCT04manual, 
            features = "PC_1", 
            center_zero = TRUE, 
            section_number = 1, 
            pt_size = 2,
            colors = RColorBrewer::brewer.pal(n = 11, name = "RdBu") |> rev())

#summarizing the percentage of spots for each cluster in the section
MapLabelsSummary(SCT04manual, 
                 column_name = "seurat_clusters", 
                 ncol = 1, 
                 bar_display = "count", #remove
                 section_number = 1) &
  theme(legend.position = "none")

# Plot with Seurat
SpatialFeaturePlot(SCT04manual, features = "nFeature_Spatial")


SCT04manual <- LoadImages(SCT04manual, verbose = FALSE)
cols <- RColorBrewer::brewer.pal(11, "Spectral") |> rev()

p <- MapFeatures(SCT04manual, 
                 features = c("MITF", "SOX2"), 
                 image_use = "raw", 
                 colors = cols)
p

p <- MapFeatures(SCT04manual, 
                 features = c("MITF", "SOX2"), 
                 image_use = "raw", 
                 colors = cols,
                 scale_alpha = TRUE)
p


plot_ly(adjusted_coords, x = ~x, y = ~y, z = ~z, type = "scatter3d", mode = "markers", color = se_merged$seurat_clusters, size = 1)

# Ensure all paths are absolute and correct
samples_SCT04manual <- "/Users/kriniogiannikou/Downloads/SCT04manual/filtered_feature_bc_matrix.h5"
imgs_SCT04manual<- "/Users/kriniogiannikou/Downloads/SCT04manual/spatial/tissue_lowres_image.png"
spotfiles_SCT04manual <- "/Users/kriniogiannikou/Downloads/SCT04manual/spatial/tissue_positions.csv"
json_SCT04manual <- "/Users/kriniogiannikou/Downloads/SCT04manual/spatial/scalefactors_json.json"

# Rebuild infoTable
infoTable_SCT04manual <- tibble(samples = samples_SCT04manual,
                                imgs = imgs_SCT04manual,
                                spotfiles = spotfiles_SCT04manual,
                                json = json_SCT04manual)
# Check the table
View(infoTable_SCT04manual)

# Load the 10x Visium data into a Seurat object
se_SCT04manual <- ReadVisiumData(infoTable_SCT04manual)
spatial_data_se_SCT04manual <- GetStaffli(se_SCT04manual)
SpatialFeaturePlot(spatial_data_SCT04manual, features = "nFeature_Spatial")

# Ensure all paths are absolute and correct
samples_SCT04manual <- "/Users/kriniogiannikou/Downloads/SCT04manual/filtered_feature_bc_matrix.h5"
imgs_SCT04manual<- "/Users/kriniogiannikou/Downloads/SCT04manual/spatial/tissue_lowres_image.png"
spotfiles_SCT04manual <- "/Users/kriniogiannikou/Downloads/SCT04manual/spatial/tissue_positions.csv"
json_SCT04manual <- "/Users/kriniogiannikou/Downloads/SCT04manual/spatial/scalefactors_json.json"

# Rebuild infoTable
infoTable_SCT04manual <- tibble(samples = samples_SCT04manual,
                                imgs = imgs_SCT04manual,
                                spotfiles = spotfiles_SCT04manual,
                                json = json_SCT04manual)
# Check the table
View(infoTable_SCT04manual)

# Load the 10x Visium data into a Seurat object
se_SCT04manual <- ReadVisiumData(infoTable_SCT04manual)
spatial_data_SCT04manual <- GetStaffli(se_SCT04manual)

se_SCT04manual <- LoadImages(se_SCT04manual)
ImagePlot(se_SCT04manual)


############
# Step 1: Load the RNA-seq data (assuming a tab-delimited text file)
# Replace 'your_rnaseq_data.txt' with your actual data file path
data <- read.table("/Users/kriniogiannikou/Downloads/Mehrdad-RNA-seq/DEseq2_NIPU_all/NIPU-all-PCA-input.txt", header = TRUE, sep = "\t")

# Step 2: Identify the most variable genes (top 3000 by variance)
variances <- apply(data, 1, var)  # Calculate variance for each gene
top_variable_genes <- names(sort(variances, decreasing = TRUE))[1:3000]

# Step 3: Identify the most highly expressed genes (top 300 by mean expression)
means <- rowMeans(data)  # Calculate the mean expression for each gene
top_highly_expressed_genes <- names(sort(means, decreasing = TRUE))[1:300]

# Step 4: Extract expression data for the selected genes

# For top 3000 most variable genes
data_variable_genes <- data[top_variable_genes, ]

# For top 300 most highly expressed genes
data_highly_expressed_genes <- data[top_highly_expressed_genes, ]

# Step 5: Standardize the data (important for UMAP)
# Scaling the data, samples as columns
data_variable_genes_scaled <- scale(t(data_variable_genes))
data_highly_expressed_genes_scaled <- scale(t(data_highly_expressed_genes))

# Step 6: Apply UMAP for most variable genes
umap_variable <- umap(data_variable_genes_scaled)  # UMAP dimensionality reduction

# Apply UMAP for most highly expressed genes
umap_highly_expressed <- umap(data_highly_expressed_genes_scaled)  # UMAP dimensionality reduction

# Step 7: Visualize the UMAP result for the most variable genes
umap_variable_df <- data.frame(UMAP1 = umap_variable$layout[,1], UMAP2 = umap_variable$layout[,2])
ggplot(umap_variable_df, aes(x = UMAP1, y = UMAP2)) +
  geom_point(aes(color = "black"), size = 2) +
  labs(title = "UMAP of Most Variable Genes (Top 3000)", x = "UMAP 1", y = "UMAP 2") +
  theme_minimal()

# Step 8: Visualize the UMAP result for the most highly expressed genes
umap_highly_expressed_df <- data.frame(UMAP1 = umap_highly_expressed$layout[,1], UMAP2 = umap_highly_expressed$layout[,2])
ggplot(umap_highly_expressed_df, aes(x = UMAP1, y = UMAP2)) +
  geom_point(aes(color = "black"), size = 2) +
  labs(title = "UMAP of Most Highly Expressed Genes (Top 300)", x = "UMAP 1", y = "UMAP 2") +
  theme_minimal()

str(umap_variable)
str(umap_highly_expressed)

sum(is.na(data_variable_genes_scaled))
sum(is.infinite(data_variable_genes_scaled))

plot(umap_variable$layout[, 1], umap_variable$layout[, 2])

str(umap_variable)

dim(data_variable_genes_scaled)

library(devtools)

devtools::install_github("YingMa0107/CARD", upgrade = "never", build_vignettes = FALSE)

usethis::create_github_token()

gitcreds::gitcreds_set()

library(CARD)


#########

# CARD analysis
library(devtools)
library(CARD)
library(Seurat)
library(scCustomize)
library(ggprism)
library(colorRamp2)
library(tidyverse)
library(ComplexHeatmap)
library(ggpattern)
library(viridis)
library(scales)
library(dplyr)
library(pheatmap)
library(MuSiC)

get_legend <- function(my_plot) {
  tmp <- ggplot_gtable(ggplot_build(my_plot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

# Immune (Blue) shades
Immune_shades <- c("#85C1E9", "#5DADE2", "#3498DB", "#2980B9", "#1F618D")

# Neuroectoderm (Red) shades
Neuroectoderm_shades <- c("#E57373", "#D64A4A", "#B53A3A", "#A83232")

# Epithelia (Green) shades
Epithelia_shades <- c("#A8E6A3", "#7CC47F", "#61B14B", "#88C992", "#4E9835", "#009A44", "#3CB371")

# Stroma (Purple) shades
Stroma_shades <- c("#D1A1DA", "#C882D5", "#A259C0", "#8E38A7", "#754C96", "#A67BB4", "#B46BB2")

# Endothelia (Orange) shades
Endothelia_shades <- c("#F9A14B", "#F47D2D", "#D36E20", "#FFB14C")

# Combine all shades into one palette
pal1 <- c(
  Epithelia_shades,
  Stroma_shades,
  Endothelia_shades,
  Immune_shades,
  Neuroectoderm_shades
)

# Assign names to the colors
names(pal1) <- c(
  "BPIFB1+_epi",
  "OCT4+_epi",
  "LGR5+_epi",
  "Epithelia_spare",
  "cycl_epi",
  "cilia_epi",  
  "enterochromaffin", # Extra in case you need it
  
  "fibroblast",
  "COL12A1+_myofibro",
  "sm_musc",
  "chondrocyte",
  "cycl_mesench",
  "PAX7+_musc_MSC",
  "sk_musc",
  
  "vasc_endo",
  "vein_endo",
  "lymph_endo",
  "cycl_endo",
  
  "mac",
  "infl_mac",
  "t_cell",
  "mast",
  "nkt",
  
  "GFAP+_astro_radGlia",
  "neuro",
  "oligo",
  "CUX2+_cilia_astro_radGlia"
)

# Check the final palette
pal1


get_legend <- function(my_plot) {
  tmp <- ggplot_gtable(ggplot_build(my_plot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

# Immune (Blue) shades
Immune_shades <- c("#85C1E9", "#5DADE2", "#3498DB", "#2980B9", "#1F618D")

# Neuroectoderm (Red) shades
Neuroectoderm_shades <- c("#E57373", "#D64A4A", "#B53A3A", "#A83232")

# Epithelia (Green) shades
Epithelia_shades <- c("#A8E6A3", "#7CC47F", "#61B14B", "#88C992", "#4E9835", "#009A44", "#3CB371")

# Stroma (Purple) shades
Stroma_shades <- c("#D1A1DA", "#C882D5", "#A259C0", "#8E38A7", "#754C96", "#A67BB4", "#B46BB2")

# Endothelia (Orange) shades
Endothelia_shades <- c("#F9A14B", "#F47D2D", "#D36E20", "#FFB14C")

# Combine all shades into one palette
pal1 <- c(
  Epithelia_shades,
  Stroma_shades,
  Endothelia_shades,
  Immune_shades,
  Neuroectoderm_shades
)

# Assign names to the colors
names(pal1) <- c(
  "BPIFB1+_epi",
  "OCT4+_epi",
  "LGR5+_epi",
  "Epithelia_spare",
  "cycl_epi",
  "cilia_epi",  
  "enterochromaffin", # Extra in case you need it
  
  "fibroblast",
  "COL12A1+_myofibro",
  "sm_musc",
  "chondrocyte",
  "cycl_mesench",
  "PAX7+_musc_MSC",
  "sk_musc",
  
  "vasc_endo",
  "vein_endo",
  "lymph_endo",
  "cycl_endo",
  
  "mac",
  "infl_mac",
  "t_cell",
  "mast",
  "nkt",
  
  "GFAP+_astro_radGlia",
  "neuro",
  "oligo",
  "CUX2+_cilia_astro_radGlia"
)

# Check the final palette
pal1


pal_broad <- c("Stroma" =  "#9467bd", "Epithelia" = "#2ca02c", "Endothelia" = "#ff7f0e", "Neuroectoderm" = "#d62728", "Immune" = "#1f77b4")
pal_solidcystic <- c("solid" = "#696969", "cystic" = "#D3D3D3")
pal_epithrichpoor <- c("EpitheliaRich" = "orange", "EpitheliaPoor" = "blue", "epithelia_rich" = "orange", "epithelia_poor"= "blue")
pal_samples <- c(
  "SCT02" = "#66C2A5",
  "SCT06" = "#FC8D62",
  "SCT09" = "#8DA0CB",
  "SCT01" = "#E78AC3",
  "SCT03" = "#A6D854",
  "SCT10" = "#FFD92F"
)
set.seed(12345)

nucSamples <- readRDS("/Users/kriniogiannikou/Downloads/nucSamplesForKrinio_deconvolution_20250123-EJ.rds")

SCT04edge <- Load10X_Spatial(data.dir = "/Users/kriniogiannikou/Downloads/SCT04redo")
SCT04middle <- Load10X_Spatial(data.dir = "/Users/kriniogiannikou/Downloads/SCT04manual")
SCT05 <- Load10X_Spatial(data.dir = "/Users/kriniogiannikou/Downloads/SCT05outs")

# 1) extract the raw counts matrix
sc_count <- GetAssayData(
  object = nucSamples,
  assay  = "RNA"
)
stopifnot(inherits(sc_count, "dgCMatrix"))
sc_count[1:4, 1:4]

# 2) build your metadata data.frame
sc_meta <- nucSamples@meta.data

# add a cellID column (so that it shows up when you print the table)
sc_meta$cellID     <- rownames(sc_meta)

# rename the cell‐type & sample columns to match CARD’s tutorial
sc_meta$cellType   <- sc_meta$HiResNamedOrdered    # assuming your meta has “CellType”
sc_meta$sampleInfo <- sc_meta$orig.ident      # assuming your meta has “Sample”

# keep only the three required columns (CELL ID, cellType, sampleInfo)
sc_meta <- sc_meta[, c("cellID", "cellType", "sampleInfo")]

# sanity check: rownames(sc_meta) must line up with colnames(sc_count)
stopifnot(identical(rownames(sc_meta), colnames(sc_count)))

# view the first few rows
head(sc_meta, 4)

spatial_count_SCT04edge <- GetAssayData(SCT04edge, assay = "Spatial", slot = "counts") # grabbing the counts here 
# 2) extract the per‐spot XY coordinates from the centroids slot
cent <- SCT04edge@images[["slice1"]]@boundaries$centroids
coords_mat <- cent@coords                 # numeric matrix n_spots × 2
barcodes   <- cent@cells                  # character vector of spot IDs

# Option A. negate y at construction time
spatial_location_SCT04edge <- data.frame(
  x         = coords_mat[,2], # The X axis seems to be in column 2 
  y         = - coords_mat[,1], # KEY FINDING - here we multiply y by -1 to flip the Y and we found the y-axis is in column 1
  row.names = barcodes,
  stringsAsFactors = FALSE
)

# sanity check: rownames(spatial_location) must match your spot barcodes
if (!all(rownames(spatial_location_SCT04edge) == colnames(spatial_count_SCT04edge))) {
  stop("Spot IDs in spatial_location and spatial_count do not line up.")
}

CARD_obj_SCT04edge = createCARDObject(
  sc_count = sc_count,
  sc_meta = sc_meta,
  spatial_count = spatial_count_SCT04edge,
  spatial_location = spatial_location_SCT04edge,
  ct.varname = "cellType",
  ct.select = unique(sc_meta$cellType),
  sample.varname = "sampleInfo",
  minCountGene = 100,
  minCountSpot = 5) 

CARD_obj_SCT04edge = CARD_deconvolution(CARD_object = CARD_obj_SCT04edge)
print(CARD_obj_SCT04edge@Proportion_CARD[1:2,])

colors = pal1
p1 <- CARD.visualize.pie(
  proportion = CARD_obj_SCT04edge@Proportion_CARD,
  spatial_location = CARD_obj_SCT04edge@spatial_location, 
  colors = colors) + 
  
  coord_fixed() # VERY IMPORTANT: need to add coord_fixed so that the shape is maintained after we flipped the y-axis

print(p1)

# pull CARD proportions into long form
prop_df  <- as.data.frame(CARD_obj_SCT04edge@Proportion_CARD)
prop_df$spot <- rownames(prop_df)

prop_long <- prop_df %>%
  pivot_longer(-spot, names_to="celltype", values_to="prop")

# 1) sum up each celltype across all spots
prop_sum <- prop_long %>%
  group_by(celltype) %>%
  summarise(total_prop = sum(prop), .groups="drop")

# 2) if you want relative % of the whole tissue, divide by grand‐total
prop_sum <- prop_sum %>%
  mutate(frac = total_prop / sum(total_prop))

# 3a) absolute “sum of proportions” stacked bar
ggplot(prop_sum, aes(x = 1, y = total_prop, fill = celltype)) +
  geom_col() +
  scale_y_continuous(expand = c(0,0)) +
  scale_fill_manual(values = pal1) +
  labs(x=NULL, y="Sum of spot‐proportions", fill="Cell Type") +
  theme_classic() +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

# 3b) relative % stacked bar
ggplot(prop_sum, aes(x = 1, y = frac, fill = celltype)) +
  geom_col() +
  scale_y_continuous(labels = percent_format(), expand = c(0,0)) +
  scale_fill_manual(values = pal1) +
  labs(x=NULL, y="Fraction of tissue", fill="Cell Type") +
  theme_classic() +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank())


#####
spatial_count_SCT04middle <- GetAssayData(SCT04middle, assay = "Spatial", slot = "counts") # grabbing the counts here 
# 2) extract the per‐spot XY coordinates from the centroids slot
cent <- SCT04middle@images[["slice1"]]@boundaries$centroids
coords_mat <- cent@coords                 # numeric matrix n_spots × 2
barcodes   <- cent@cells                  # character vector of spot IDs

# Option A. negate y at construction time
spatial_location_SCT04middle <- data.frame(
  x         = coords_mat[,2], # The X axis seems to be in column 2 
  y         = -coords_mat[,1], # KEY FINDING - here we multiply y by -1 to flip the Y and we found the y-axis is in column 1
  row.names = barcodes,
  stringsAsFactors = FALSE
)


# sanity check: rownames(spatial_location) must match your spot barcodes
if (!all(rownames(spatial_location_SCT04middle) == colnames(spatial_count_SCT04middle))) {
  stop("Spot IDs in spatial_location and spatial_count do not line up.")
}


CARD_obj_SCT04middle = createCARDObject(
  sc_count = sc_count,
  sc_meta = sc_meta,
  spatial_count = spatial_count_SCT04middle,
  spatial_location = spatial_location_SCT04middle,
  ct.varname = "cellType",
  ct.select = unique(sc_meta$cellType),
  sample.varname = "sampleInfo",
  minCountGene = 100,
  minCountSpot = 5) 

CARD_obj_SCT04middle = CARD_deconvolution(CARD_object = CARD_obj_SCT04middle)

print(CARD_obj_SCT04middle@Proportion_CARD[1:2,])

colors = pal1
p1 <- CARD.visualize.pie(
  proportion = CARD_obj_SCT04middle@Proportion_CARD,
  spatial_location = CARD_obj_SCT04middle@spatial_location, 
  colors = colors) + 
  coord_fixed() # VERY IMPORTANT: need to add coord_fixed so that the shape is maintained after we flipped the y-axis

print(p1)

# pull CARD proportions into long form
prop_df  <- as.data.frame(CARD_obj_SCT04middle@Proportion_CARD)
prop_df$spot <- rownames(prop_df)

prop_long <- prop_df %>%
  pivot_longer(-spot, names_to="celltype", values_to="prop")

# 1) sum up each celltype across all spots
prop_sum <- prop_long %>%
  group_by(celltype) %>%
  summarise(total_prop = sum(prop), .groups="drop")

# 2) if you want relative % of the whole tissue, divide by grand‐total
prop_sum <- prop_sum %>%
  mutate(frac = total_prop / sum(total_prop))

# 3a) absolute “sum of proportions” stacked bar
ggplot(prop_sum, aes(x = 1, y = total_prop, fill = celltype)) +
  geom_col() +
  scale_y_continuous(expand = c(0,0)) +
  scale_fill_manual(values = pal1) +
  labs(x=NULL, y="Sum of spot‐proportions", fill="Cell Type") +
  theme_classic() +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

# 3b) relative % stacked bar
ggplot(prop_sum, aes(x = 1, y = frac, fill = celltype)) +
  geom_col() +
  scale_y_continuous(labels = percent_format(), expand = c(0,0)) +
  scale_fill_manual(values = pal1) +
  labs(x=NULL, y="Fraction of tissue", fill="Cell Type") +
  theme_classic() +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank())


spatial_count_SCT05 <- GetAssayData(SCT05, assay = "Spatial", slot = "counts") # grabbing the counts here 
# 2) extract the per‐spot XY coordinates from the centroids slot
cent <- SCT05@images[["slice1"]]@boundaries$centroids
coords_mat <- cent@coords                 # numeric matrix n_spots × 2
barcodes   <- cent@cells                  # character vector of spot IDs

# Option A. negate y at construction time
spatial_location_SCT05 <- data.frame(
  x         = coords_mat[,2], # The X axis seems to be in column 2 
  y         = -coords_mat[,1], # KEY FINDING - here we multiply y by -1 to flip the Y and we found the y-axis is in column 1
  row.names = barcodes,
  stringsAsFactors = FALSE
)


# sanity check: rownames(spatial_location) must match your spot barcodes
if (!all(rownames(spatial_location_SCT05) == colnames(spatial_count_SCT05))) {
  stop("Spot IDs in spatial_location and spatial_count do not line up.")
}

CARD_obj_SCT05 = createCARDObject(
  sc_count = sc_count,
  sc_meta = sc_meta,
  spatial_count = spatial_count_SCT05,
  spatial_location = spatial_location_SCT05,
  ct.varname = "cellType",
  ct.select = unique(sc_meta$cellType),
  sample.varname = "sampleInfo",
  minCountGene = 100,
  minCountSpot = 5) 

CARD_obj_SCT05 = CARD_deconvolution(CARD_object = CARD_obj_SCT05)

print(CARD_obj_SCT05@Proportion_CARD[1:2,])

colors = pal1
p1 <- CARD.visualize.pie(
  proportion = CARD_obj_SCT05@Proportion_CARD,
  spatial_location = CARD_obj_SCT05@spatial_location, 
  colors = colors) + 
  coord_fixed() # VERY IMPORTANT: need to add coord_fixed so that the shape is maintained after we flipped the y-axis

print(p1)

# pull CARD proportions into long form
prop_df  <- as.data.frame(CARD_obj_SCT05@Proportion_CARD)
prop_df$spot <- rownames(prop_df)

prop_long <- prop_df %>%
  pivot_longer(-spot, names_to="celltype", values_to="prop")

# 1) sum up each celltype across all spots
prop_sum <- prop_long %>%
  group_by(celltype) %>%
  summarise(total_prop = sum(prop), .groups="drop")

# 2) if you want relative % of the whole tissue, divide by grand‐total
prop_sum <- prop_sum %>%
  mutate(frac = total_prop / sum(total_prop))

# 3a) absolute “sum of proportions” stacked bar
ggplot(prop_sum, aes(x = 1, y = total_prop, fill = celltype)) +
  geom_col() +
  scale_y_continuous(expand = c(0,0)) +
  scale_fill_manual(values = pal1) +
  labs(x=NULL, y="Sum of spot‐proportions", fill="Cell Type") +
  theme_classic() +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank())


# 3b) relative % stacked bar
ggplot(prop_sum, aes(x = 1, y = frac, fill = celltype)) +
  geom_col() +
  scale_y_continuous(labels = percent_format(), expand = c(0,0)) +
  scale_fill_manual(values = pal1) +
  labs(x=NULL, y="Fraction of tissue", fill="Cell Type") +
  theme_classic() +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

p8 <- CARD.visualize.Cor(CARD_obj_SCT04middle@Proportion_CARD,colors = NULL) # if not provide, we will use the default colors
print(p8)


## select the cell type that we are interested
ct.visualize = c( "vasc_endo",
                  "artery_endon",
                  "vein_endo",
                  "lymph_endo",
                  "cycl_endo")

## visualize the spatial distribution of the cell type proportion
pdf("SCT04edge_endo_CARD.pdf", width = 7, height = 6)
CARD.visualize.prop(
  proportion = CARD_obj_SCT04edge@Proportion_CARD,        
  spatial_location = CARD_obj_SCT04edge@spatial_location, 
  ct.visualize = ct.visualize,                 ### selected cell types to visualize
  colors = c("lightblue","lightyellow","red"), ### if not provide, we will use the default colors
  NumCols = 3,                                 ### number of columns in the figure panel
  pointSize = 0.5)                             ### point size in ggplot2 scatterplot  

dev.off()

## select the cell type that we are interested
ct.visualize = c( "vasc_endo",
                  "artery_endon",
                  "vein_endo",
                  "lymph_endo",
                  "cycl_endo")

## visualize the spatial distribution of the cell type proportion
pdf("test_CARD.pdf", width = 7, height = 6)
CARD.visualize.prop(
  proportion = CARD_obj_SCT05@Proportion_CARD,        
  spatial_location = CARD_obj_SCT05@spatial_location, 
  ct.visualize = ct.visualize,                 ### selected cell types to visualize
  colors = c("lightblue","lightyellow","red"), ### if not provide, we will use the default colors
  NumCols = 3,                                 ### number of columns in the figure panel
  pointSize = 0.5)                             ### point size in ggplot2 scatterplot  

dev.off()

#######

# Define palettes
Immune_shades <- c("#85C1E9", "#5DADE2", "#3498DB", "#2980B9", "#1F618D")
Neuroectoderm_shades <- c("#E57373", "#D64A4A", "#B53A3A", "#A83232")
Epithelia_shades <- c("#A8E6A3", "#7CC47F", "#61B14B", "#88C992", "#4E9835", "#009A44", "#3CB371")
Stroma_shades <- c("#D1A1DA", "#C882D5", "#A259C0", "#8E38A7", "#754C96", "#A67BB4", "#B46BB2")
Endothelia_shades <- c("#F9A14B", "#F47D2D", "#D36E20", "#FFB14C")

pal1 <- c(
  Epithelia_shades,
  Stroma_shades,
  Endothelia_shades,
  Immune_shades,
  Neuroectoderm_shades
)

library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)

# Helper function to reshape and summarize CARD proportions
get_prop_summary <- function(card_obj, sample_name) {
  prop_df <- as.data.frame(card_obj@Proportion_CARD) %>%
    mutate(spot = rownames(.)) %>%
    pivot_longer(-spot, names_to = "celltype", values_to = "prop") %>%
    group_by(celltype) %>%
    summarise(total_prop = sum(prop), .groups = "drop") %>%
    mutate(frac = total_prop / sum(total_prop),
           sample = sample_name)
  return(prop_df)
}

# Extract and summarize from each CARD object
prop_SCT05       <- get_prop_summary(CARD_obj_SCT05, "SCT05")
prop_SCT04edge   <- get_prop_summary(CARD_obj_SCT04edge, "SCT04_edge")
prop_SCT04middle <- get_prop_summary(CARD_obj_SCT04middle, "SCT04_middle")

# Combine all data
prop_all <- bind_rows(prop_SCT05, prop_SCT04edge, prop_SCT04middle)

# (Optional) reorder factor levels for consistent x-axis
prop_all$sample <- factor(prop_all$sample, levels = c("SCT05", "SCT04_edge", "SCT04_middle"))

# Relative % stacked bar plot
ggplot(prop_all, aes(x = sample, y = frac, fill = celltype)) +
  geom_col() +
  scale_y_continuous(labels = percent_format(), expand = c(0,0)) +
  scale_fill_manual(values = pal1) +  # use your predefined palette
  labs(x = NULL, y = "Fraction of tissue", fill = "Cell Type") +
  theme_classic() +
  theme(axis.text.x = element_text(size = 10),
        axis.ticks.x = element_blank())


###
# Load required package
library(RColorBrewer)

# 1. Extract number of cell types from the CARD object
n_celltypes <- ncol(CARD_obj_SCT04edge@Proportion_CARD)

# 2. Generate a color palette with enough distinct colors
colors <- colorRampPalette(brewer.pal(12, "Set3"))(n_celltypes)

# 3. Plot CARD pie chart with specified radius
p1 <- CARD.visualize.pie(
  proportion = CARD_obj_SCT04edge@Proportion_CARD,
  spatial_location = CARD_obj_SCT04edge@spatial_location,
  colors = colors,
  radius = 0.52  # adjust or set to NULL as needed
)

# 4. Display the plot
print(p1)


##### Broad categories####

nucSamples <- readRDS("/Users/kriniogiannikou/Downloads/nucSamplesForKrinio_deconvolution_20250123-EJ.rds")

SCT04edge <- Load10X_Spatial(data.dir = "/Users/kriniogiannikou/Downloads/SCT04redo")
SCT05 <- Load10X_Spatial(data.dir = "/Users/kriniogiannikou/Downloads/SCT05outs")
SCT04middle <- Load10X_Spatial(data.dir = "/Users/kriniogiannikou/Downloads/SCT04manual")

# 1) extract the raw counts matrix
sc_count <- GetAssayData(
  object = nucSamples,
  assay  = "RNA"
)
stopifnot(inherits(sc_count, "dgCMatrix"))
sc_count[1:4, 1:4]
#> 4 x 4 sparse Matrix of class "dgCMatrix"
#>                        SCT02_… SCT02_… SCT02_… SCT02_…
#> AL627309.1                    .       .       .       .
#> AL627309.5                    .       .       .       .
#> AL627309.4                    .       .       .       .
#> LINC01409                     .       .       .       .

# 2) build your metadata data.frame
sc_meta <- nucSamples@meta.data

# add a cellID column (so that it shows up when you print the table)
sc_meta$cellID     <- rownames(sc_meta)

# rename the cell‐type & sample columns to match CARD’s tutorial
sc_meta$cellType   <- sc_meta$BroadCategory    # assuming your meta has “CellType”
sc_meta$sampleInfo <- sc_meta$orig.ident      # assuming your meta has “Sample”

# keep only the three required columns (CELL ID, cellType, sampleInfo)
sc_meta <- sc_meta[, c("cellID", "cellType", "sampleInfo")]

# sanity check: rownames(sc_meta) must line up with colnames(sc_count)
stopifnot(identical(rownames(sc_meta), colnames(sc_count)))

# view the first few rows
head(sc_meta, 4)


spatial_count_SCT04middle <- GetAssayData(SCT04middle, assay = "Spatial", slot = "counts") # grabbing the counts here 
# 2) extract the per‐spot XY coordinates from the centroids slot
cent <- SCT04middle@images[["slice1"]]@boundaries$centroids
coords_mat <- cent@coords                 # numeric matrix n_spots × 2
barcodes   <- cent@cells                  # character vector of spot IDs

# Option A. negate y at construction time
spatial_location_SCT04middle <- data.frame(
  x         = coords_mat[,2], # The X axis seems to be in column 2 
  y         = - coords_mat[,1], # KEY FINDING - here we multiply y by -1 to flip the Y and we found the y-axis is in column 1
  row.names = barcodes,
  stringsAsFactors = FALSE
)


# sanity check: rownames(spatial_location) must match your spot barcodes
if (!all(rownames(spatial_location_SCT04middle) == colnames(spatial_count_SCT04middle))) {
  stop("Spot IDs in spatial_location and spatial_count do not line up.")
}

# 1) extract the raw counts matrix
sc_count <- GetAssayData(
  object = nucSamples,
  assay  = "RNA"
)
stopifnot(inherits(sc_count, "dgCMatrix"))
sc_count[1:4, 1:4]


# 2) build your metadata data.frame
sc_meta <- nucSamples@meta.data

# add a cellID column (so that it shows up when you print the table)
sc_meta$cellID     <- rownames(sc_meta)

# rename the cell‐type & sample columns to match CARD’s tutorial
sc_meta$cellType   <- sc_meta$BroadCategory    # assuming your meta has “CellType”
sc_meta$sampleInfo <- sc_meta$orig.ident      # assuming your meta has “Sample”

# keep only the three required columns (CELL ID, cellType, sampleInfo)
sc_meta <- sc_meta[, c("cellID", "cellType", "sampleInfo")]

# sanity check: rownames(sc_meta) must line up with colnames(sc_count)
stopifnot(identical(rownames(sc_meta), colnames(sc_count)))

# view the first few rows
head(sc_meta, 4)


CARD_obj_SCT04middle = createCARDObject(
  sc_count = sc_count,
  sc_meta = sc_meta,
  spatial_count = spatial_count_SCT04middle,
  spatial_location = spatial_location_SCT04middle,
  ct.varname = "cellType",
  ct.select = unique(sc_meta$cellType),
  sample.varname = "sampleInfo",
  minCountGene = 10,
  minCountSpot = 3) 

CARD_obj_SCT04middle = CARD_deconvolution(CARD_object = CARD_obj_SCT04middle)

print(CARD_obj_SCT04middle@Proportion_CARD[1:2,])

colors = pal_broad
p1 <- CARD.visualize.pie(
  proportion = CARD_obj_SCT04middle@Proportion_CARD,
  spatial_location = CARD_obj_SCT04middle@spatial_location, 
  colors = colors) + 
  
  coord_fixed() # VERY IMPORTANT: need to add coord_fixed so that the shape is maintained after we flipped the y-axis

pdf("CARD_SCT04middle_piechart.pdf", width = 15, height = 9)
print(p1)
dev.off()

# pull CARD proportions into long form
prop_df  <- as.data.frame(CARD_obj_SCT04middle@Proportion_CARD)
prop_df$spot <- rownames(prop_df)

prop_long <- prop_df %>%
  pivot_longer(-spot, names_to="celltype", values_to="prop")

# 1) sum up each celltype across all spots
prop_sum <- prop_long %>%
  group_by(celltype) %>%
  summarise(total_prop = sum(prop), .groups="drop")

# 2) if you want relative % of the whole tissue, divide by grand‐total
prop_sum <- prop_sum %>%
  mutate(frac = total_prop / sum(total_prop))

pdf("SCT04middle_stacked_bar_plot_CARD_proportions.pdf", width = 6, height = 5)

# 3a) absolute “sum of proportions” stacked bar
ggplot(prop_sum, aes(x = 1, y = total_prop, fill = celltype)) +
  geom_col() +
  scale_y_continuous(expand = c(0,0)) +
  scale_fill_manual(values = pal_broad) +
  labs(x=NULL, y="Sum of spot‐proportions", fill="Cell Type") +
  theme_classic() +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

dev.off() 

pdf("SCT04middle_stacked_bar_plot_CARD_percentage.pdf", width = 6, height = 5)

# 3b) relative % stacked bar
ggplot(prop_sum, aes(x = 1, y = frac, fill = celltype)) +
  geom_col() +
  scale_y_continuous(labels = percent_format(), expand = c(0,0)) +
  scale_fill_manual(values = pal_broad) +
  labs(x=NULL, y="Fraction of tissue", fill="Cell Type") +
  theme_classic() +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

dev.off() 

######
spatial_count_SCT05 <- GetAssayData(SCT05, assay = "Spatial", slot = "counts") # grabbing the counts here 
# 2) extract the per‐spot XY coordinates from the centroids slot
cent <- SCT05@images[["slice1"]]@boundaries$centroids
coords_mat <- cent@coords                 # numeric matrix n_spots × 2
barcodes   <- cent@cells                  # character vector of spot IDs

# Option A. negate y at construction time
spatial_location_SCT05 <- data.frame(
  x         = coords_mat[,2], # The X axis seems to be in column 2 
  y         = -coords_mat[,1], # KEY FINDING - here we multiply y by -1 to flip the Y and we found the y-axis is in column 1
  row.names = barcodes,
  stringsAsFactors = FALSE
)


# sanity check: rownames(spatial_location) must match your spot barcodes
if (!all(rownames(spatial_location_SCT05) == colnames(spatial_count_SCT05))) {
  stop("Spot IDs in spatial_location and spatial_count do not line up.")
}

CARD_obj_SCT05 = createCARDObject(
  sc_count = sc_count,
  sc_meta = sc_meta,
  spatial_count = spatial_count_SCT05,
  spatial_location = spatial_location_SCT05,
  ct.varname = "cellType",
  ct.select = unique(sc_meta$cellType),
  sample.varname = "sampleInfo",
  minCountGene = 100,
  minCountSpot = 5) 


CARD_obj_SCT05 = CARD_deconvolution(CARD_object = CARD_obj_SCT05)
print(CARD_obj_SCT05@Proportion_CARD[1:2,])


colors = pal_broad
p1 <- CARD.visualize.pie(
  proportion = CARD_obj_SCT05@Proportion_CARD,
  spatial_location = CARD_obj_SCT05@spatial_location, 
  colors = colors) + 
  
  coord_fixed() # VERY IMPORTANT: need to add coord_fixed so that the shape is maintained after we flipped the y-axis

pdf("CARD_SCT05_piechart.pdf", width = 10, height = 6)
print(p1)
dev.off()

# pull CARD proportions into long form
prop_df  <- as.data.frame(CARD_obj_SCT05@Proportion_CARD)
prop_df$spot <- rownames(prop_df)

prop_long <- prop_df %>%
  pivot_longer(-spot, names_to="celltype", values_to="prop")

# 1) sum up each celltype across all spots
prop_sum <- prop_long %>%
  group_by(celltype) %>%
  summarise(total_prop = sum(prop), .groups="drop")

# 2) if you want relative % of the whole tissue, divide by grand‐total
prop_sum <- prop_sum %>%
  mutate(frac = total_prop / sum(total_prop))

pdf("SCT05_stacked_bar_plot_CARD_proportions.pdf", width = 6, height = 5)

# 3a) absolute “sum of proportions” stacked bar
ggplot(prop_sum, aes(x = 1, y = total_prop, fill = celltype)) +
  geom_col() +
  scale_y_continuous(expand = c(0,0)) +
  scale_fill_manual(values = pal_broad) +
  labs(x=NULL, y="Sum of spot‐proportions", fill="Cell Type") +
  theme_classic() +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

dev.off()


pdf("SCT05_stacked_bar_plot_CARD_percentage.pdf", width = 6, height = 5)

# 3b) relative % stacked bar
ggplot(prop_sum, aes(x = 1, y = frac, fill = celltype)) +
  geom_col() +
  scale_y_continuous(labels = percent_format(), expand = c(0,0)) +
  scale_fill_manual(values = pal_broad) +
  labs(x=NULL, y="Fraction of tissue", fill="Cell Type") +
  theme_classic() +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

dev.off()

######
spatial_count_SCT04edge <- GetAssayData(SCT04edge, assay = "Spatial", slot = "counts") # grabbing the counts here 
# 2) extract the per‐spot XY coordinates from the centroids slot
cent <- SCT04edge@images[["slice1"]]@boundaries$centroids
coords_mat <- cent@coords                 # numeric matrix n_spots × 2
barcodes   <- cent@cells                  # character vector of spot IDs

# Option A. negate y at construction time
spatial_location_SCT04edge <- data.frame(
  x         = coords_mat[,2], # The X axis seems to be in column 2 
  y         = -coords_mat[,1], # KEY FINDING - here we multiply y by -1 to flip the Y and we found the y-axis is in column 1
  row.names = barcodes,
  stringsAsFactors = FALSE
)

# sanity check: rownames(spatial_location) must match your spot barcodes
if (!all(rownames(spatial_location_SCT04edge) == colnames(spatial_count_SCT04edge))) {
  stop("Spot IDs in spatial_location and spatial_count do not line up.")
}


CARD_obj_SCT04edge = createCARDObject(
  sc_count = sc_count,
  sc_meta = sc_meta,
  spatial_count = spatial_count_SCT04edge,
  spatial_location = spatial_location_SCT04edge,
  ct.varname = "cellType",
  ct.select = unique(sc_meta$cellType),
  sample.varname = "sampleInfo",
  minCountGene = 100,
  minCountSpot = 5) 


CARD_obj_SCT04edge = CARD_deconvolution(CARD_object = CARD_obj_SCT04edge)
print(CARD_obj_SCT04edge@Proportion_CARD[1:2,])


colors = pal_broad
p1 <- CARD.visualize.pie(
  proportion = CARD_obj_SCT04edge@Proportion_CARD,
  spatial_location = CARD_obj_SCT04edge@spatial_location, 
  colors = colors) + 
  
  coord_fixed() # VERY IMPORTANT: need to add coord_fixed so that the shape is maintained after we flipped the y-axis

pdf("CARD_SCT04edge_piechart.pdf", width = 10, height = 6)
print(p1)
dev.off()

# pull CARD proportions into long form
prop_df  <- as.data.frame(CARD_obj_SCT04edge@Proportion_CARD)
prop_df$spot <- rownames(prop_df)

prop_long <- prop_df %>%
  pivot_longer(-spot, names_to="celltype", values_to="prop")

# 1) sum up each celltype across all spots
prop_sum <- prop_long %>%
  group_by(celltype) %>%
  summarise(total_prop = sum(prop), .groups="drop")

# 2) if you want relative % of the whole tissue, divide by grand‐total
prop_sum <- prop_sum %>%
  mutate(frac = total_prop / sum(total_prop))

pdf("SCT04edge_stacked_bar_plot_CARD_proportions.pdf", width = 6, height = 5)

# 3a) absolute “sum of proportions” stacked bar
ggplot(prop_sum, aes(x = 1, y = total_prop, fill = celltype)) +
  geom_col() +
  scale_y_continuous(expand = c(0,0)) +
  scale_fill_manual(values = pal_broad) +
  labs(x=NULL, y="Sum of spot‐proportions", fill="Cell Type") +
  theme_classic() +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

dev.off()


pdf("SCT04edge_stacked_bar_plot_CARD_percentage.pdf", width = 6, height = 5)

# 3b) relative % stacked bar
ggplot(prop_sum, aes(x = 1, y = frac, fill = celltype)) +
  geom_col() +
  scale_y_continuous(labels = percent_format(), expand = c(0,0)) +
  scale_fill_manual(values = pal_broad) +
  labs(x=NULL, y="Fraction of tissue", fill="Cell Type") +
  theme_classic() +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

dev.off()

pdf("SCT04edge_correlation_plot_CARD.pdf", width = 6, height = 5)
p8 <- CARD.visualize.Cor(CARD_obj_SCT04edge@Proportion_CARD,colors = NULL) # if not provide, we will use the default colors
print(p8)
dev.off()

# Filter the ct.visualize vector to include only valid cell types
ct.visualize <- ct.visualize[ct.visualize %in% available_ct]

## select the cell type that we are interested
ct.visualize = c( "Stroma",
                  "Epithelia",
                  "Endothelia",
                  "Neuroectoderm",
                  "Immune")

## visualize the spatial distribution of the cell type proportion
pdf("SCT04edge_broad_categories_CARD.pdf", width = 7, height = 6)
CARD.visualize.prop(
  proportion = CARD_obj_SCT04edge@Proportion_CARD,        
  spatial_location = CARD_obj_SCT04edge@spatial_location, 
  ct.visualize = ct.visualize,
  colors = c("lightblue", "lightyellow", "red")[1:length(ct.visualize)],
  NumCols = 3,
  pointSize = 0.5
)
dev.off()


ct.visualize = c("Stroma", "Epithelia", "Endothelia", "Neuroectoderm", "Immune")

pdf("SCT05_broad_categories_CARD-vs2.pdf", width = 7, height = 6)
CARD.visualize.prop(
  proportion = CARD_obj_SCT05@Proportion_CARD,
  spatial_location = CARD_obj_SCT05@spatial_location,
  ct.visualize = ct.visualize,
  colors = colorRampPalette(c("lightblue", "lightyellow", "red"))(100),  # continuous gradient
  NumCols = 3,
  pointSize = 0.4
)
dev.off()

# Define color palettes
Immune_shades <- c("#85C1E9", "#5DADE2", "#3498DB", "#2980B9", "#1F618D")
Neuroectoderm_shades <- c("#E57373", "#D64A4A", "#B53A3A", "#A83232")
Epithelia_shades <- c("#A8E6A3", "#7CC47F", "#61B14B", "#88C992", "#4E9835", "#009A44", "#3CB371")
Stroma_shades <- c("#D1A1DA", "#C882D5", "#A259C0", "#8E38A7", "#754C96", "#A67BB4", "#B46BB2")
Endothelia_shades <- c("#F9A14B", "#F47D2D", "#D36E20", "#FFB14C")

# Combine into one palette
pal1 <- c(
  Epithelia_shades,
  Stroma_shades,
  Endothelia_shades,
  Immune_shades,
  Neuroectoderm_shades
)

markers_for_subtypes <- list(
  GFAP_plus_astro_radGlia = c("GFAP"),
  CUX2_plus_cilia_astro_radGlia = c("CUX2", "PCDH11X"),
  neuro = c("GAD2", "SRRM3", "MYT1L", "MIAT", "CACNA1B"),
  oligo = c("OLIG1", "OLIG2"),
  enterochromaffin = c("CHGA"),
  
  ALL_EPI = c("EPCAM", "ELF3"),
  BPIFB1_plus_epi = c("BPIFB1"),
  OCT4_plus_epi = c("POU5F1"),
  LGR5_plus_epi = c("LGR5"),
  cilia_epi = c("DNAAF1", "DNAH3"),
  cycl_epi = c(),
  
  ALL_STROMA = c("COL14A1", "POSTN"),
  fibroblast = c(),
  COL12A1_plus_myofibro = c("COL12A1"),
  PAX7_plus_musc_MSC = c("PAX7"),
  cycl_mesench = c(),
  chondrocyte = c("CHI3L1", "COMP"),
  sm_musc = c("NOTCH3", "ACTA2"),
  sk_musc = c("ASB5", "DES", "NCAM1"),
  
  vasc_endo = c("PECAM1", "VWF", "CDH5"),
  artery_endon = c(),
  vein_endo = c("ACKR1", "SELP"),
  lymph_endo = c("TBX1", "PROX1"),
  cycl_endo = c(),
  
  mac = c("PTPRC", "MSR1",  "CD163"), #removed "CD86",
  t_cell = c("IL2RG", "CD3E", "CD7", "NKG7"),
  nkt = c("CD247", "CD96"),
  infl_mac = c("AIF1", "ITGB2"), # removed"MNDA", "FGL2", "S100A9"
  mast = c("CHAT", "GCSAML" ), #removed "MS4A2", "CPA3"
  
  ALL_Cycling = c("MKI67")
)

all_markers_for_subtypes <- unique(unlist(markers_for_subtypes))

#STACKBAR plot

# Load libraries
library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)

# Define the broad lineage palette
pal_broad <- c(
  "Stroma" = "#9467bd",
  "Epithelia" = "#2ca02c",
  "Endothelia" = "#ff7f0e",
  "Neuroectoderm" = "#d62728",
  "Immune" = "#1f77b4"
)

# --------------------------------------------
# Helper function to summarize CARD proportions
# --------------------------------------------
get_prop_summary <- function(card_obj, sample_name) {
  as.data.frame(card_obj@Proportion_CARD) %>%
    mutate(spot = rownames(.)) %>%
    pivot_longer(-spot, names_to = "celltype", values_to = "prop") %>%
    group_by(celltype) %>%
    summarise(total_prop = sum(prop), .groups = "drop") %>%
    mutate(frac = total_prop / sum(total_prop),
           sample = sample_name)
}

# --------------------------------------------
# Summarize from each CARD object
# --------------------------------------------
prop_SCT05       <- get_prop_summary(CARD_obj_SCT05, "SCT05")
prop_SCT04edge   <- get_prop_summary(CARD_obj_SCT04edge, "SCT04_edge")
prop_SCT04middle <- get_prop_summary(CARD_obj_SCT04middle, "SCT04_middle")

# Combine all summaries
prop_all <- bind_rows(prop_SCT05, prop_SCT04edge, prop_SCT04middle)

# Ensure sample order
prop_all$sample <- factor(prop_all$sample, levels = c("SCT05", "SCT04_edge", "SCT04_middle"))

# --------------------------------------------
# Plot and Save as PDF using broad palette
# --------------------------------------------
pdf("stacked_barplot_CARD.pdf", width = 6, height = 4)

ggplot(prop_all, aes(x = sample, y = frac, fill = celltype)) +
  geom_col(width = 0.7) +
  scale_y_continuous(labels = percent_format(), expand = c(0, 0)) +
  scale_fill_manual(values = pal_broad) +
  labs(x = NULL, y = "Fraction of tissue", fill = "Cell Type") +
  theme_classic() +
  theme(
    axis.text.x = element_text(size = 10),
    axis.ticks.x = element_blank()
  )

dev.off()
