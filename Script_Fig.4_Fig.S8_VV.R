
######
#### "Type 1 and type 2 dendritic cell subsets are sufficient to maintain intestinal immune tolerance via integrin αvβ8-mediated TGF-β activation" by This and Brichart-Vernos et al. 
#### Script related to analysis of Figure 4 and Supplemental Figure S8
#### Author: Venla Väänänen (venlan@sund.ku.dk)

library(Seurat)
library(ggplot2)
library(tidyverse)
library(tidyr)
library(dplyr)
library(ggpubr)
library(scales)
library(BPCells)


mycols <- rev(c('#a50026','#d73027','#f46d43','#fdae61','#fee090','#e0f3f8','#abd9e9','#74add1','#4575b4','#313695'))
mycols_b <- c("#bdbdbd","#d9d9d9","#FDDBC7","#F4A582","#D6604D","#B2182B","#67001F")

## Load the preprocessed object
MNP_MLN_obj8 <- readRDS("MNP_MLN_obj8.rds")

# Fig. 4A
DimPlot(MNP_MLN_obj8, group.by = 'Cell_annotations', reduction = 'umap.harmony', raster = TRUE, raster.dpi = c(700,700), pt.size = 5)+ coord_fixed(ratio = 1.25)

# Fig. 4B

DC1score <- c("XCR1","CADM1","CLEC9A")
DC1.score <- list(DC1score)
MNP_MLN_obj8 <- AddModuleScore(MNP_MLN_obj8, features = DC1.score, ctrl = 5, name = 'DC1_smallscore')

CCR7DCscore <- c('CCR7', 'LAMP3', 'CD40')
CCR7DC.score <- list(CCR7DCscore)
MNP_MLN_obj8 <- AddModuleScore(MNP_MLN_obj8, features = CCR7DC.score, ctrl = 5, name = 'CCR7DC_score')

CLEC10A_pos <- c("CEBPD","IL1R2","FCER1A","FCN1","IL1B","CD1C","CEBPB","CLEC10A")
CLEC10A_pos.score <- list(CLEC10A_pos)
MNP_MLN_obj8 <- AddModuleScore(MNP_MLN_obj8, features = CLEC10A_pos.score, ctrl = 5, name = 'CLEC10A_pos.score')

CLEC10A_neg <- c("CLEC4A","LTB","IL22RA2","CD3E","ZFP36L2","CD300A","AREG","PPP1R14A","LILRA4")
CLEC10A_neg.score <- list(CLEC10A_neg)
MNP_MLN_obj8 <- AddModuleScore(MNP_MLN_obj8, features = CLEC10A_neg.score, ctrl = 5, name = 'CLEC10A_neg.score')

RDCscore <- c("RORC", "LTB", "CLEC4A", "SOX4", "SPI1", "CCR6", "PRDM16")
RDC.score <- list(RDCscore)
MNP_MLN_obj8 <- AddModuleScore(MNP_MLN_obj8, features = RDC.score, ctrl = 5, name = 'RORC_DC_score')


FeaturePlot(MNP_MLN_obj8, features = 'DC1_smallscore1', reduction = 'umap.harmony', raster = TRUE, raster.dpi = c(700,700), pt.size = 5)+ coord_fixed(ratio = 1.25) +scale_colour_gradientn(colours=mycols)
FeaturePlot(MNP_MLN_obj8, features = 'CCR7DC_score1', reduction = 'umap.harmony', raster = TRUE, raster.dpi = c(700,700), pt.size = 5)+ coord_fixed(ratio = 1.25) +scale_colour_gradientn(colours=mycols)
FeaturePlot(MNP_MLN_obj8, features = 'CLEC10A_pos.score1', reduction = 'umap.harmony', raster = TRUE, raster.dpi = c(700,700), pt.size = 5)+ coord_fixed(ratio = 1.25) +scale_colour_gradientn(colours=mycols)
FeaturePlot(MNP_MLN_obj8, features = 'CLEC10A_neg.score1', reduction = 'umap.harmony', raster = TRUE, raster.dpi = c(700,700), pt.size = 5)+ coord_fixed(ratio = 1.25) +scale_colour_gradientn(colours=mycols)
FeaturePlot(MNP_MLN_obj8, features = 'RORC_DC_score1', reduction = 'umap.harmony', raster = TRUE, raster.dpi = c(700,700), pt.size = 5)+ coord_fixed(ratio = 1.25) +scale_colour_gradientn(colours=mycols)


# Fig. 4C
# Make the DotPlot
gene_of_interest <- "ITGB8"

# Compute the basic DotPlot data
p <- DotPlot(
  MNP_MLN_obj8,
  features = gene_of_interest,
  group.by = "CellType_Patient",
  assay = "RNA"
)

# Split composite grouping back into cell type and patient columns
p$data <- p$data %>%
  tidyr::separate(id, into = c("CellType", "Patient"), sep = "_")

# Define the desired cell type order (top -> bottom)
desired_order <- c("CCR7+ DC", "CLEC10A- cDC2", "CLEC10A+ cDC2", "cDC1", "RORgt+ DC")

# Set factor levels — reversed so top cell type appears at top of plot
p$data$CellType <- factor(p$data$CellType, levels = rev(desired_order))

# Scale expression within each patient (column-wise normalization)
p$data <- p$data %>%
  group_by(Patient) %>%
  mutate(avg.exp.patient.scaled = scales::rescale(avg.exp, to = c(0, 1))) %>%
  ungroup()

# Plot using the new per-patient scaled expression
ggplot(p$data, aes(x = Patient, y = CellType)) +
  geom_point(aes(size = pct.exp, color = avg.exp.patient.scaled)) +
  scale_color_gradientn(
    colors = c("#d9d9d9", "#D6604D", "#67001F"),
    values = scales::rescale(c(0, 0.5, 1))
  ) +
  theme_bw() +
  labs(
    title = paste("Expression of", gene_of_interest),
    x = "Patient",
    y = "Cell Type",
    color = "Scaled\n(per patient)",
    size = "% Expressing"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank()
  )


# Fig. S8A

MNP_MLN_obj8@meta.data
Idents(MNP_MLN_obj8) <- MNP_MLN_obj8@meta.data$Cell_annotations
MNP_MLN_obj8_ave <- AverageExpression(MNP_MLN_obj8, assays = 'RNA', return.seurat = TRUE)

table(MNP_MLN_obj8_ave$Cell_annotations)
my_levels <- c('CCR7+ DC','CLEC10A- cDC2','CLEC10A+ cDC2','cDC1','RORgt+ DC')
levels(MNP_MLN_obj8_ave) <- my_levels

genes3 <- c('CD40','CD80','CD86','RELB','CD83','CD274','PDCD1LG2','CD200','FAS','SOCS1','SOCS2','CCR7','MYO1G','FSCN1','MARCKS','MARCKSL1') # maturation and regulatory
DoHeatmap(MNP_MLN_obj8_ave, features = c(genes3), draw.lines = FALSE)+scale_fill_gradientn(colours=mycols)


# Fig. S8B

# To order highes expressing cells on top (same as order = T but keeping rastering)
df <- FetchData(MNP_MLN_obj8, vars = c("ITGB8", "umapharmony_1", "umapharmony_2"))
df <- df[order(df$ITGB8, decreasing = FALSE), ]   # low first, high last (so high is on top)

ggplot(df, aes(x = umapharmony_1, y = umapharmony_2, color = ITGB8)) +
  geom_point(size = 1, alpha = 1, raster = TRUE) +
  coord_fixed(1.25) +
  scale_color_gradientn(colours = mycols_b) +
  theme_classic()













