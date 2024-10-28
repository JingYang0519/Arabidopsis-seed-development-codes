### hdWGCNA --------------------------------------------------------------------

library(Seurat)
library(cowplot)
library(patchwork)
library(WGCNA)
library(hdWGCNA)
library(monocle3)
library(SeuratWrappers)
library(ggplot2)
library(viridis)
library(dplyr)

# set random seed for reproducibility
set.seed(12345)

# optionally enable multithreading
enableWGCNAThreads(nThreads = 8)

seurat_obj <- readRDS("f03.subset.AnnoLabel.merge.fix.rds")
seurat_obj <- SetupForWGCNA(
  seurat_obj,
  gene_select = "fraction", # the gene selection approach
  fraction = 0.05, # fraction of cells that a gene needs to be expressed in order to be included
  wgcna_name = 'trajectory' # the name of the hdWGCNA experiment
)

seurat_obj$DetailAnno <- gsub("Seed.coat", "Seed.Coat", seurat_obj$DetailAnno)
print(sort(unique(seurat_obj$DetailAnno)))

seurat_obj <- MetacellsByGroups(
  seurat_obj = seurat_obj,
  group.by = c("Batch","DetailAnno"), # specify the columns in seurat_obj@meta.data to group by
  reduction = 'umap', # select the dimensionality reduction to perform KNN on
  k = 25, # nearest-neighbors parameter
  max_shared = 10, # maximum number of shared cells between two metacells
  ident.group = "DetailAnno" # set the Idents of the metacell seurat object
)

# normalize metacell expression matrix:
seurat_obj <- NormalizeMetacells(seurat_obj)

print(paste0("Current time1: ", Sys.time()))

#### Analysis ------------------------------------------------------------------

seurat_obj <- SetDatExpr(
  seurat_obj,
  group_name = c("Embryo","Embryo_Cotyledon","Embryo_Provasculature,.QC","Embryo_RAM",
                 "Endosperm","Endosperm_CZE","Endosperm_MCE","Endosperm_PEN",
                 "Seed.Coat","Seed.Coat_CZSC","Seed.Coat_Endothelium","Seed.Coat_Outer.integument"),
  group.by='DetailAnno', assay = "RNA")

# Test different soft powers:
seurat_obj <- TestSoftPowers(
  seurat_obj,
  networkType = 'signed' # you can also use "unsigned" or "signed hybrid"
)

# plot the results:
plot_list <- PlotSoftPowers(seurat_obj)

# assemble with patchwork
pdf("p03.soft-power.threshold.pdf",10,10)
print(wrap_plots(plot_list, ncol=2))
dev.off()

# Using soft_power = 6
seurat_obj <- ConstructNetwork(
  seurat_obj,
  tom_name = "trajectory",
  overwrite_tom=TRUE
)

saveRDS(seurat_obj, "f04.AnnoLabel.ConstructNetwork.rds")

seurat_obj <- readRDS("f04.AnnoLabel.ConstructNetwork.rds")
print(paste0("Current time2: ", Sys.time()))

pdf("p04.PlotDendrogram.pdf",10,10)
print(PlotDendrogram(seurat_obj, main='Embryo Endosperm Seedcoat hdWGCNA Dendrogram'))
dev.off()

seurat_obj <- ModuleEigengenes(seurat_obj)
seurat_obj <- ModuleConnectivity(seurat_obj)

saveRDS(seurat_obj, "f05.AnnoLabel.ModuleEigengenesConnectivity.rds")

# Rename before you run: 
# modules <- GetModules(seurat_obj)
# mods <- levels(modules$module)
# mods <- mods[mods!='grey']
# 
# meta <- seurat_obj@meta.data
# seurat_obj@meta.data <- cbind(meta, MEs)


# # make dotplot
# p <- DotPlot(
#   seurat_obj,
#   group.by='DetailAnno',
#   features = rev(mods)
# ) + RotatedAxis() +
#   scale_color_gradient2(high='red', mid='grey95', low='blue') + xlab('') + ylab('') +
#   theme(
#     plot.title = element_text(hjust = 0.5),
#     axis.line.x = element_blank(),
#     axis.line.y = element_blank(),
#     panel.border = element_rect(colour = "black", fill=NA, size=1)
#   )
# 
# pdf("p05.dotplot_MEs.pdf", 8, 4)
# print(p)
# dev.off()

###### -------------------------------------------------------------------------
library(Seurat)
library(cowplot)
library(patchwork)
library(WGCNA)
library(hdWGCNA)
library(monocle3)
library(SeuratWrappers)
library(ggplot2)
library(viridis)
library(dplyr)

# set random seed for reproducibility
set.seed(12345)

# optionally enable multithreading
enableWGCNAThreads(nThreads = 8)

seurat_obj <- readRDS("f05.AnnoLabel.ModuleEigengenesConnectivity.rds")

p  <- PlotModuleTrajectory(
  seurat_obj, ncol = 3,
  pseudotime_col = c('embryo_pseudotime', 'endosperm_pseudotime', 'seedcoat_pseudotime'),
  group_colors = c("#357B6E","#63a4cf","#F4A370"))

pdf("p06.PlotModuleTrajectory.pdf", 8, 8)
print(p)
dev.off()

p  <- PlotModuleTrajectory(
  seurat_obj, ncol = 3,
  pseudotime_col = c('pseudotime'))

pdf("p06.PlotModuleTrajectory.pseudobulk.pdf", 8, 8)
print(p)
dev.off()

plot_list <- ModuleFeaturePlot(
  seurat_obj,
  features='hMEs', # plot the hMEs
  order=TRUE # order so the points with highest hMEs are on top
)
p <- wrap_plots(plot_list, ncol=5)
pdf("p07.hME.projection.pdf", 12, 6)
print(p)
dev.off()

## See ModuleNetworks dir
ModuleNetworkPlot(seurat_obj)

seurat_obj <- RunModuleUMAP(
  seurat_obj,
  n_hubs = 10, # number of hub genes to include for the UMAP embedding
  n_neighbors=15, # neighbors parameter for UMAP
  min_dist=0.1 # min distance between points in UMAP space
)

ModuleUMAPPlot.cry <- function(
    seurat_obj,
    sample_edges = TRUE, # TRUE if we sample edges randomly, FALSE if we take the top edges
    edge_prop = 0.2,
    label_hubs = 5, # how many hub genes to label?
    edge.alpha=0.25,
    vertex.label.cex=0.5,
    label_genes = NULL,
    return_graph = FALSE, # this returns the igraph object instead of plotting
    keep_grey_edges = TRUE,
    wgcna_name=NULL,
    ...
){
  
  if(is.null(wgcna_name)){wgcna_name <- seurat_obj@misc$active_wgcna}
  
  # get the TOM
  TOM <- GetTOM(seurat_obj, wgcna_name)
  
  # get modules,
  modules <- GetModules(seurat_obj, wgcna_name)
  
  # get the UMAP df:
  umap_df <- GetModuleUMAP(seurat_obj, wgcna_name)
  mods <- levels(umap_df$module)
  mods <- mods[mods != 'grey']
  
  # subset the TOM:
  subset_TOM <- TOM[umap_df$gene, umap_df$gene[umap_df$hub == 'hub']]
  
  # genes to label:
  # hub_labels <- selected_modules %>% group_by(module) %>% top_n(label_hubs, wt=kME) %>% .$gene_name
  hub_list <- lapply(mods, function(cur_mod){
    cur <- subset(modules, module == cur_mod)
    cur[,c('gene_name', paste0('kME_', cur_mod))] %>%
      top_n(label_hubs) %>% .$gene_name
  })
  names(hub_list) <- mods
  hub_labels <- as.character(unlist(hub_list))
  print('hub labels')
  print(hub_labels)
  print(label_genes)
  if(is.null(label_genes)){
    label_genes <- hub_labels
  } else{
    if(!any(label_genes %in% umap_df$gene)){
      stop("Some genes in label_genes not found in the UMAP.")
    }
    label_genes <- unique(c(label_genes, hub_labels))
  }
  print(label_genes)
  
  # subset module df by genes in the UMAP df:
  selected_modules <- modules[umap_df$gene,]
  selected_modules <- cbind(selected_modules, umap_df[,c('UMAP1', 'UMAP2', 'hub', 'kME')])
  
  selected_modules$label <- ifelse(selected_modules$gene_name %in% label_genes, selected_modules$gene_name, '')
  selected_modules$fontcolor <- ifelse(selected_modules$color == 'black', 'gray50', 'black')
  
  # set frome color
  # same color as module for all genes, black outline for the selected hub genes
  selected_modules$framecolor <- ifelse(selected_modules$gene_name %in% label_genes, ifelse(selected_modules$color == 'black', 'white', 'black'), selected_modules$color)
  
  # melt TOM into long format
  edge_df <- subset_TOM %>% reshape2::melt()
  print(dim(edge_df))
  
  # set color of each edge based on value:
  edge_df$color <- future.apply::future_sapply(1:nrow(edge_df), function(i){
    gene1 = as.character(edge_df[i,'Var1'])
    gene2 = as.character(edge_df[i,'Var2'])
    
    col1 <- selected_modules[selected_modules$gene_name == gene1, 'color']
    col2 <- selected_modules[selected_modules$gene_name == gene2, 'color']
    
    if(col1 == col2){
      col = col1
    } else{
      col = 'grey90'
    }
    col
  })
  
  # keep grey edges?
  if(!keep_grey_edges){
    edge_df <- edge_df %>% subset(color != 'grey90')
  }
  
  # subset edges:
  groups <- unique(edge_df$color)
  if(sample_edges){
    # randomly sample
    temp <- do.call(rbind, lapply(groups, function(cur_group){
      cur_df <- edge_df %>% subset(color == cur_group)
      n_edges <- nrow(cur_df)
      cur_sample <- sample(1:n_edges, round(n_edges * edge_prop))
      cur_df[cur_sample,]
    }))
  } else{
    
    # get top strongest edges
    temp <- do.call(rbind, lapply(groups, function(cur_group){
      cur_df <- edge_df %>% subset(color == cur_group)
      n_edges <- nrow(cur_df)
      cur_df %>% dplyr::top_n(round(n_edges * edge_prop), wt=value)
    }))
  }
  
  edge_df <- temp
  print(dim(edge_df))
  
  # scale edge values between 0 and 1 for each module
  edge_df <- edge_df %>% group_by(color) %>% mutate(value=scale01(value))
  
  # edges & vertices are plotted in igraph starting with the first row, so re-order s.t. strong edges are on bottom, all gray on the top of the table:
  edge_df <- edge_df %>% arrange(value)
  edge_df <- rbind(
    subset(edge_df, color == 'grey90'),
    subset(edge_df, color != 'grey90')
  )
  
  # set alpha of edges based on kME
  edge_df$color_alpha <- ifelse(
    edge_df$color == 'grey90',
    alpha(edge_df$color, alpha=edge_df$value/2),
    alpha(edge_df$color, alpha=edge_df$value)
  )
  
  # re-order vertices so hubs are plotted on top
  selected_modules <- rbind(
    subset(selected_modules , hub == 'other'),
    subset(selected_modules , hub != 'other')
  )
  
  # re-order vertices so labeled genes are on top
  selected_modules <- rbind(
    subset(selected_modules , label == ''),
    subset(selected_modules , label != '')
  )
  
  # setup igraph:
  g <- igraph::graph_from_data_frame(
    edge_df,
    directed=FALSE,
    vertices=selected_modules
  )
  
  print('making net')
  print(head(edge_df))
  print(head(selected_modules))
  
  if(return_graph){return(g)}
  
  plot(
    g,
    layout=  as.matrix(selected_modules[,c('UMAP1', 'UMAP2')]),
    # edge.color=adjustcolor(igraph::E(g)$color, alpha.f=edge.alpha),
    edge.color=adjustcolor(igraph::E(g)$color_alpha, alpha.f=edge.alpha),
    vertex.size=igraph::V(g)$kME * 3,
    edge.curved=0,
    edge.width=0.5,
    vertex.color=igraph::V(g)$color,
    vertex.label=igraph::V(g)$label,
    vertex.label.dist=1.1,
    vertex.label.degree=-pi/4,
    vertex.label.family='Helvetica', #vertex.label.font=vertex_df$font,
    vertex.label.font = 3,
    vertex.label.color = igraph::V(g)$fontcolor,
    vertex.label.cex=0,
    vertex.frame.color=igraph::V(g)$framecolor,
    margin=0
  )
  
}


pdf("p08.ModuleUMAPPlot.v5.pdf", 8, 8)
ModuleUMAPPlot.cry(
  seurat_obj,
  label_genes = c("AT1G61720",
                  "AT5G48880",
                  "AT1G07720",
                  "AT2G40170",
                  "AT3G51810",
                  "AT1G48130",
                  "AT1G24360",
                  "AT5G47670",
                  "AT5G46290",
                  "AT2G38530",
                  "AT5G49190",
                  "AT5G64080",
                  "AT1G08560",
                  "AT3G08900",
                  "AT1G71250",
                  "AT2G42840",
                  "AT1G27950",
                  "AT4G20140",
                  "AT5G49360",
                  "AT1G23200",
                  "AT3G20210"),
  edge.alpha=0.25,
  sample_edges=TRUE,
  edge_prop=0.1, # proportion of edges to sample (20% here)
  label_hubs = 0,# how many hub genes to plot per module?
  keep_grey_edges=FALSE
)
dev.off()

### Pseudotime dimplot ---------------------------------------------------------
meta <- seurat_obj@meta.data
meta <- arrange(meta, desc(CellType))

names(meta) <- gsub("pseudotime","ppseudotime",names(meta))
names(meta) <- gsub("embryo_ppseudotime","epseudotime",names(meta))
names(meta) <- gsub("endosperm_ppseudotime","npseudotime",names(meta))
names(meta) <- gsub("seedcoat_ppseudotime","spseudotime",names(meta))

p1 <- meta %>%
  ggplot(aes(x=UMAP1, y=UMAP2, color=epseudotime)) +
  scale_color_gradientn(colors=rev(magma(300))[20:250], na.value='grey') +
  geom_point(size = 0.1) + theme_classic()

meta1 <- meta[which(meta$CellType != "Endosperm"),]
meta2 <- meta[which(meta$CellType == "Endosperm"),]
meta3 <- rbind(meta1, meta2)

p2 <- meta3 %>%
  ggplot(aes(x=UMAP1, y=UMAP2, color=npseudotime)) +
  scale_color_gradientn(colors=rev(magma(300))[20:250], na.value='grey') +
  geom_point(size = 0.1) + theme_classic()

meta1 <- meta[which(meta$CellType != "Seed_coat"),]
meta2 <- meta[which(meta$CellType == "Seed_coat"),]
meta3 <- rbind(meta1, meta2)

p3 <- meta3 %>%
  ggplot(aes(x=UMAP1, y=UMAP2, color=spseudotime)) +
  scale_color_gradientn(colors=rev(magma(300))[20:250], na.value='grey') +
  geom_point(size = 0.1) + theme_classic()

p4 <- meta %>%
  ggplot(aes(x=UMAP1, y=UMAP2, color=ppseudotime)) +
  # scale_color_gradientn(colors=rev(viridis(10))[4:9], na.value='grey') +
  scale_color_gradientn(colors=rev(magma(300))[20:250], na.value='grey') +
  geom_point(size = 0.1) + theme_classic()

library(cowplot)
pdf('p02.umap_pseudotime.v2.pdf', 16, 3)
plot_grid(p4, p1, p2, p3, ncol =4)
dev.off()

### Module dotplot -------------------------------------------------------------

MEs <- GetMEs(seurat_obj)
modules <- GetModules(seurat_obj)
mods <- levels(modules$module)
mods <- mods[mods!='grey']

meta <- seurat_obj@meta.data
seurat_obj@meta.data <- cbind(meta, MEs)
seurat_obj <- subset(seurat_obj, DetailAnno != "Unknown")

library(MetBrewer)
p <- DotPlot(
  seurat_obj,
  group.by='DetailAnno',
  features = rev(mods)
) + RotatedAxis() +
  scale_color_gradientn(colours = colorRampPalette(c("#3dacdb","#DBDCE1","#ffa55c","#ff745c","#ba061c"))(100)) + 
  # scale_color_gradient2(high='red', mid='grey95', low='blue') + 
  xlab('') + ylab('') +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.line.x = element_blank(),
    axis.line.y = element_blank(),
    panel.border = element_rect(colour = "black", fill=NA, size=1)
  )

pdf("p05.dotplot_MEs.v2.pdf", 7, 5)
print(p)
dev.off()

