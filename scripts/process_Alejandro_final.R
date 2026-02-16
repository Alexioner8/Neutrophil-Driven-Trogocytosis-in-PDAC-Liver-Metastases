source("/mnt/raidexttmp/Alejandro/functions.R")
setwd("/mnt/raidexttmp/Alejandro/Larissa_10X/Seurat")

#
# this searches for all input matrices
#
files=list()
for (x in c("LB1", "LB2", "LB3"))
{
file<- Sys.glob(paste(paste("/mnt/raidexttmp/Alejandro/Larissa_10X/raw_feature_bc_matrix", x, sep=""), "/matrix.mtx.gz", sep=""))

files[x]= file
}




#
# here we iterate over all input matrices and load them into allfiles (gene expression) and allABs (antibody capture)
#
allfiles.raw = list()
allABs.raw = list()
for (file in files)
{
  samplename = sub(".*matrix","",str_split(dirname(file), "/")[[1]][6]) #estos numeros indican el nombre de la carpeta
  foldername = dirname(file)
  
  print(paste(samplename, foldername))
  
  h5file = Read10X(foldername,unique.features = TRUE)

  if (is.null(names(h5file)))
  {
      print(paste("WITHOUT AB", samplename))
    allfiles.raw[[samplename]] = h5file
  } else {
      print(paste("WITH AB", samplename))
    allfiles.raw[[samplename]] = h5file$`Gene Expression`
    allABs.raw[[samplename]] = h5file$`Antibody Capture`
  }

  print(paste(samplename, nrow(allfiles.raw[[samplename]]), "x", ncol(allfiles.raw[[samplename]]), "genes x cells"))
}

length(allfiles.raw)
names(allfiles.raw)
length(allABs.raw)
names(allABs.raw)


#
# here we create a list of seurat object. each entry corresponds to an input matrix from above
#

objlist = list()
for (x in names(allfiles.raw))
{

    matrix = allfiles.raw[[x]]
    
    # this creates a Seurat object from the count matrix. it sets the object's project to x and prepends the sample name to all cells
    # the patternlist.mouse contains patterns for mt and RP-genes
    filteredObj = makeSeuratObj(matrix, x, patternList.mouse)
    
    # this creates log-normalized count matrices in RNA assay
    filteredObj <- NormalizeData(filteredObj, verbose = FALSE)
    # this calculates the most (2000) variable features per data set. variable features are features which show a high variance between all cells of a sample
    filteredObj <- FindVariableFeatures(filteredObj, verbose = FALSE)
    
    objlist[[x]] = filteredObj

    print(x)
    print(filteredObj)
    
    
}

names(objlist)

objlist.raw = objlist

objlist <- lapply(X = objlist.raw, FUN = function(obj) {
  # mt content: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6072887/
  print(paste("Seurat obj project", obj@project.name))
  print(obj)
  obj <- subset(obj, subset = nFeature_RNA > 100 & nFeature_RNA < 6000 & nCount_RNA > 300)
  obj <- subset(obj, subset = percent.mt < 15)
  print(obj)
  
  return(obj)
})


for (name in names(objlist))
{
  p=VlnPlot(objlist[[name]], features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0)
  save_plot(p, paste("/mnt/raidexttmp/Alejandro/Larissa_10X/Seurat/QC", paste(name, "filtered_violins_qc", sep="_"), sep="/"), fig.width=10, fig.height=6)

  p=VlnPlot(objlist[[name]], features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0, combine=F)
  p[[1]] = p[[1]] + scale_y_continuous(limits = c(0, 1000), breaks = seq(0,1000,100))
  p[[2]] = p[[2]] + scale_y_continuous(limits = c(0, 1000), breaks = seq(0,1000,100))
  p[[3]] = p[[3]] + scale_y_continuous(limits = c(0, 100), breaks = seq(0,100,5))
  p = combine_plot_grid_list(plotlist=p, ncol=3)
  save_plot(p, paste("/mnt/raidexttmp/Alejandro/Larissa_10X/Seurat/QC", paste(name, "filtered_violins_detail_qc", sep="_"), sep="/"), fig.width=18, fig.height=6)
  
  
  plot1 <- FeatureScatter(objlist[[name]], feature1 = "nCount_RNA", feature2 = "percent.mt")
  plot2 <- FeatureScatter(objlist[[name]], feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
  save_plot(plot1 + plot2, paste("/mnt/raidexttmp/Alejandro/Larissa_10X/Seurat/QC", paste(name, "filtered_scatter_ncount_mt", sep="_"), sep="/"), fig.width=10, fig.height=6)
  
  plot1 <- FeatureScatter(objlist[[name]], feature1 = "nCount_RNA", feature2 = "percent.rp")
  plot2 <- FeatureScatter(objlist[[name]], feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
  save_plot(plot1 + plot2, paste("/mnt/raidexttmp/Alejandro/Larissa_10X/Seurat/QC", paste(name, "filtered_scatter_ncount_rp", sep="_"), sep="/"), fig.width=10, fig.height=6)
}




# Define your sample names
samples <- c("LB1", "LB2", "LB3")

# Loop over each sample
for (sample_name in samples) {
  
  # Get Seurat object and HTO matrix
  obj <- objlist[[sample_name]]
  hto_matrix <- allABs.raw[[sample_name]]

  # ⚠️ Remove rows 9 and 10 (HTOs not used)
  hto_matrix <- hto_matrix[1:8, ]
  
  # Prefix HTO matrix column names with sample name
  colnames(hto_matrix) <- paste0(sample_name, "_", colnames(hto_matrix))
  
  # Subset HTO matrix to match cells in the Seurat object
  obj.htos <- hto_matrix[, colnames(obj)]
  
  # Add HTO assay to Seurat object
  obj[["HTO"]] <- CreateAssayObject(counts = obj.htos)
  
  
  # Normalize and demultiplex
  obj <- NormalizeData(obj, assay = "HTO", normalization.method = "CLR")
  obj <- HTODemux(obj, assay = "HTO", positive.quantile = 0.99)
  
  # Print classification table
  cat("\n", sample_name, "HTO classification:\n")
  print(table(obj$HTO_classification.global))
  
  # Plot ridge plot
  p <- RidgePlot(obj, assay = "HTO", features = rownames(obj[["HTO"]]), ncol = 4)
  
  # Save plot
  ggsave(
    filename = paste0("/mnt/raidexttmp/Alejandro/Larissa_10X/Seurat/HTOridge_", sample_name, ".png"),
    plot = p,
    width = 40, height = 20, units = "cm", dpi = 300
  )
  
  # Store updated Seurat object back into objlist
  objlist[[sample_name]] <- obj
}




# Define HTO label mapping per sample
hto_label_map <- list(
  LB1 = c("#1" = "C1", "#2" = "C2", "#3" = "C3", "#4" = "S1",
          "#5" = "S2", "#6" = "S3", "#7" = "L1", "#8" = "L2"),
  LB2 = c("#1" = "S4", "#2" = "S5", "#3" = "S6", "#4" = "L3",
          "#5" = "L4", "#6" = "L5", "#7" = "C4", "#8" = "C5"),
  LB3 = c("#1" = "L6", "#2" = "L7", "#3" = "L8", "#4" = "C6",
          "#5" = "C7", "#6" = "C8", "#7" = "S8", "#8" = "S9")
)

# Loop through samples
singletlist=list()
for (sample_name in names(hto_label_map)) {
  
  obj <- objlist[[sample_name]]
  
  # Keep only singlets
  obj_singlets <- obj[, obj$HTO_classification.global != "Doublet"]
  
  # Get current mapping
  label_map <- hto_label_map[[sample_name]]
  
  # Initialize classification vector
  featVec <- obj_singlets$HTO_classification
  
  # Replace HTO labels with biological names
  featVec <- unname(label_map[featVec])  # `unname()` avoids NAs from unmatched keys
  
  # Fill in "Negative" manually
  featVec[is.na(featVec) & obj_singlets$HTO_classification == "Negative"] <- "Negative"
  
  # Add new metadata column
  obj_singlets$CSclassification <- featVec
  
  # Store updated Seurat object back into objlist
  singletlist[[sample_name]] <- obj_singlets
}


# Clear large unused objects to free up memory
rm(allfiles.raw)
rm(allABs.raw)
rm(objlist.raw)
# Run garbage collection to reclaim memory
gc()




####Integration

prepareFinalList = function(finalList)
{



print("cells per experiment")
print(mapply(sum, lapply(finalList, function(x) {dim(x)[2]})))
print("total cells")
print(sum(mapply(sum, lapply(finalList, function(x) {dim(x)[2]}))))


objlist = list()
for (objname in names(finalList))
{

    x = finalList[[objname]]  
    print(objname)
    print(paste("Seurat obj project", x@project.name))
    
    Project(x) = objname
    print(paste("Seurat obj project", x@project.name))

    DefaultAssay(x) = "RNA"
    x <- NormalizeData(x, verbose = FALSE)
    x <- FindVariableFeatures(x, verbose = FALSE)

    x$library = objname

    objlist[[objname]] = x
}

features <- SelectIntegrationFeatures(object.list = objlist, nfeatures = 2000)

objlist <- lapply(X = objlist, FUN = function(x) {

    print(paste("Seurat obj project", x@project.name))
    print(x)



    x <- ScaleData(x, features = features, verbose = FALSE, assay="RNA", vars.to.regress = c('percent.rp', 'percent.mt', "nCount_RNA"))
    x <- RunPCA(x, verbose = FALSE, reduction.name="pca", assay="RNA")
    #x <- suppressWarnings(SCTransform(x, verbose = FALSE,vars.to.regress = c('percent.rp', 'percent.mt', "nCount_RNA","S.Score", "G2M.Score")))
    
    plot1 <- ElbowPlot(x)
    save_plot(plot1, paste("/mnt/raidexttmp/Alejandro/Larissa_10X/Seurat/QC", paste(x@project.name, "elbowplot_dimensionality", sep="_"), sep="/"), fig.width=10, fig.height=6)


    x$project = x@project.name

    return(x)
})

return(objlist)

}



analyseFinalList = function(objlist, intname) 
{

dir.create(intname, recursive = TRUE)


#
# integrate based on RNA/GEX assay
#
objSamples = objlist
print(objSamples)

objSamples = lapply(objSamples, function(x) {
  DefaultAssay(x) <- 'RNA'
  x <- RunPCA(x, verbose = FALSE, reduction.name="pca",  assay="RNA")
#  DefaultAssay(x) <- 'SCT'
#  print(x@reductions$pca)
  return(x)
})
print("GEX integration features")
print(objSamples)

features_gex <- SelectIntegrationFeatures(object.list = objSamples, nfeatures = 2000)
objlist.anchors <- FindIntegrationAnchors(object.list = objSamples,  reduction = "rpca", dims = 1:30, anchor.features = features_gex) 
obj.list.integrated <- IntegrateData(new.assay.name = "integratedgex", anchorset = objlist.anchors, dims = 1:30, verbose=T) 
print("GEX integration done")

#
# integrated GEX viz
#
obj.list.integrated = ScaleData(obj.list.integrated, assay="integratedgex")
obj.list.integrated <- RunPCA(obj.list.integrated, reduction.name="igpca", assay="integratedgex")
obj.list.integrated <- RunUMAP(obj.list.integrated, reduction = "igpca", dims = 1:30, reduction.key = "UMAP_",)
p=DimPlot(obj.list.integrated, reduction="umap", shuffle = T, seed = 1, group.by= "orig.ident")
save_plot(p, paste(intname, "wnn_ig_dimplot", sep="/"), 12, 8)

p=DimPlot(obj.list.integrated, reduction="igpca", group.by= "orig.ident")
save_plot(p, paste(intname, "wnn_pca_ig_dimplot", sep="/"), 12, 8)

p=DimPlot(obj.list.integrated, shuffle = T, seed = 1, group.by= "CSclassification")
save_plot(p, "/mnt/raidexttmp/Alejandro/Larissa_10X/Seurat/CSclassification_dimplot", 12, 8)
p=DimPlot(obj.list.integrated, shuffle = T, seed = 1, group.by= "HTO_classification", split.by= "HTO_classification", ncol=4)
save_plot(p, "/mnt/raidexttmp/Alejandro/Larissa_10X/Seurat/wnn_HTO_dimplot", 48, 24)




return(obj.list.integrated)

}


finalList_sample = prepareFinalList(singletlist)
integratedList_LB = analyseFinalList(finalList_sample, "LB_int")



source("/mnt/raidexttmp/Alejandro/functions.R")
setwd("/mnt/raidexttmp/Alejandro/Larissa_10X/Seurat")
#saveRDS(integratedList_LB, file = "/mnt/raidexttmp/Alejandro/Larissa_10X/Seurat/LB10X.rds")
#integratedList_LB=readRDS("/mnt/raidexttmp/Alejandro/Larissa_10X/Seurat/LB10X.rds")



##############################
cellList = colnames(integratedList_LB)

featVec <- vector(mode="character", length=length(cellList))
featVec = as.character(integratedList_LB$CSclassification)

featVec[featVec == "C1"] = "Control"
featVec[featVec == "C2"] = "Control"
featVec[featVec == "C3"] = "Control"
featVec[featVec == "C4"] = "Control"
featVec[featVec == "C5"] = "Control"
featVec[featVec == "C6"] = "Control"
featVec[featVec == "C7"] = "Control"
featVec[featVec == "C8"] = "Control"
featVec[featVec == "S1"] = "5day_tumor"
featVec[featVec == "S2"] = "5day_tumor"
featVec[featVec == "S3"] = "5day_tumor"
featVec[featVec == "S4"] = "5day_tumor"
featVec[featVec == "S5"] = "5day_tumor"
featVec[featVec == "S6"] = "5day_tumor"
featVec[featVec == "S8"] = "5day_tumor"
featVec[featVec == "S9"] = "5day_tumor"
featVec[featVec == "L1"] = "10day_tumor"
featVec[featVec == "L2"] = "10day_tumor"
featVec[featVec == "L3"] = "10day_tumor"
featVec[featVec == "L4"] = "10day_tumor"
featVec[featVec == "L5"] = "10day_tumor"
featVec[featVec == "L6"] = "10day_tumor"
featVec[featVec == "L7"] = "10day_tumor"
featVec[featVec == "L8"] = "10day_tumor"


integratedList_LB$Pheno=featVec





##############

#cluster identification
DefaultAssay(integratedList_LB) = "integratedgex"
obj.in <- FindNeighbors(integratedList_LB, reduction="igpca", dims = 1:30)
DefaultAssay(obj.in) <- "integratedgex"
obj.in <- FindClusters(obj.in, resolution = 0.2, algorithm = 4)

#new.cluster.ids <- c("SMCs",)
#names(new.cluster.ids) <- levels(obj.in)
#obj.in <- RenameIdents(obj.in, new.cluster.ids)

p=DimPlot(obj.in, pt.size = 0.5, label=T, reduction = "umap")
save_plot(p, paste("/mnt/raidexttmp/Alejandro/Larissa_10X/Seurat", "Dimplot_umap_ident", sep="/"), 15, 10)

p=DimPlot(obj.in, shuffle = T, seed = 1, group.by= "Pheno", split.by= "Pheno", ncol=2)
save_plot(p, "/mnt/raidexttmp/Alejandro/Larissa_10X/Seurat/Group_dimplot", 15*2, 10*2)

obj.pos= obj.in[,obj.in$Pheno!="Negative"]

p=DimPlot(obj.pos, shuffle = T, seed = 1, group.by= "Pheno", split.by= "Pheno", ncol=3)
save_plot(p, "/mnt/raidexttmp/Alejandro/Larissa_10X/Seurat/Grouppos_dimplot", 15*3, 10)


#######
# Extract metadata
metadata <- obj.pos@meta.data

# Count and compute percentage per cluster
metadata_summary <- metadata %>%
  group_by(seurat_clusters, Pheno) %>%
  summarise(cell_count = n(), .groups = 'drop') %>%
  group_by(seurat_clusters) %>%
  mutate(
    percent = cell_count / sum(cell_count) * 100,
    percent_label = paste0(round(percent, 1), "%"),
    ymax = cumsum(percent),
    ymin = ymax - percent,
    y_label = 100-(ymin + (percent/2))  # middle of each segment
  )

# Plot
p <- ggplot(metadata_summary, aes(x = seurat_clusters, y = percent, fill = Pheno)) +
  geom_bar(stat = "identity", position = "stack") +
  geom_text(aes(y = y_label, label = percent_label), color = "black", size = 3.5) +
  theme_minimal() +
  scale_fill_brewer(palette = "Set2") +
  labs(
    title = "Cluster Composition by Phenotype (%)",
    x = "Seurat Cluster",
    y = "Percentage of Cells",
    fill = "Pheno"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    text = element_text(size = 12),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

# Save the plot
save_plot(
  p,
  "/mnt/raidexttmp/Alejandro/Larissa_10X/Seurat/Barplot_annotation_percent",
  fig.width = 15,
  fig.height = 10
)

# Extract metadata
metadata <- obj.pos@meta.data

# Get Seurat color palette (same as DimPlot)
# Get colors from your existing UMAP plot
umap_plot <- DimPlot(obj.pos, group.by = "seurat_clusters", shuffle = TRUE)
umap_data <- ggplot_build(umap_plot)$data[[1]]

# Extract mapping between cluster identity and color
umap_color_map <- data.frame(
  cluster = umap_data$group,
  color = umap_data$colour
) %>%
  dplyr::distinct(cluster, color)

# Convert cluster to same type as in metadata (to ensure consistent sorting)
umap_color_map$cluster <- as.character(umap_color_map$cluster)

# Get ordered cluster levels from Seurat object
cluster_levels <- levels(obj.pos$seurat_clusters)

# Reorder colors to match cluster levels
ordered_colors <- umap_color_map$color[match(cluster_levels, umap_color_map$cluster)]
names(ordered_colors) <- cluster_levels


# Count and compute percentage per Pheno
metadata_summary <- metadata %>%
  group_by(Pheno, seurat_clusters) %>%
  summarise(cell_count = n(), .groups = 'drop') %>%
  group_by(Pheno) %>%
  mutate(
    percent = cell_count / sum(cell_count) * 100,
    percent_label = paste0(round(percent, 1), "%")
  )

# Plot using the Seurat palette
p <- ggplot(metadata_summary, aes(x = Pheno, y = percent, fill = seurat_clusters)) +
  geom_bar(stat = "identity", position = "stack", color = "white") +
  geom_text(aes(label = percent_label),
            color = "black", size = 3.2, position = position_stack(vjust = 0.5)) +
  scale_fill_manual(values = ordered_colors) +
  labs(
    title = "Cellular Composition per Phenotype",
    x = "Phenotype (Condition)",
    y = "Percentage of Cells per Cluster",
    fill = "Seurat Cluster"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
    axis.title = element_text(face = "bold"),
    legend.position = "right",
    panel.grid = element_blank()
  )

# Save the plot
save_plot(
  p,
  "/mnt/raidexttmp/Alejandro/Larissa_10X/Seurat/Barplot_cluster_percent_per_pheno_SeuratColors",
  fig.width = 8,
  fig.height = 10
)
 
########
DefaultAssay(obj.in) <- "RNA"
DefaultAssay(obj.pos) <- "RNA"
obj.in = JoinLayers(obj.in)
obj.pos = JoinLayers(obj.pos)


#~/R.sh
.libPaths("~/R/x86_64-suse-linux-gnu-library/4.3")
source("/mnt/raidexttmp/Alejandro/functions.R")
setwd("/mnt/raidexttmp/Alejandro/Larissa_10X/Seurat")
#saveRDS(obj.in, file = "/mnt/raidexttmp/Alejandro/Larissa_10X/Seurat/LB10X.rds")
#obj.in=readRDS("/mnt/raidexttmp/Alejandro/Larissa_10X/Seurat/LB10X.rds")

obj.pos= obj.in[,obj.in$Pheno!="Negative"]

cellList = colnames(obj.pos)

featVec <- vector(mode="character", length=length(cellList))
featVec = as.character(obj.pos$seurat_clusters)

featVec[featVec == "1"] = "Macrophages"
featVec[featVec == "2"] = "Neutrophils"
featVec[featVec == "3"] = "B-cells 1"
featVec[featVec == "4"] = "T-cells 1"
featVec[featVec == "5"] = "B-cells 2"
featVec[featVec == "6"] = "Kupffer cells"
featVec[featVec == "7"] = "Endothelial cells"
featVec[featVec == "8"] = "T-cells 2"
featVec[featVec == "9"] = "Regulatory Kupffer cells"
featVec[featVec == "10"] = "Fibroblasts"
featVec[featVec == "11"] = "Dendritic cells"
featVec[featVec == "12"] = "NK cells"
featVec[featVec == "13"] = "Proliferating Macrophages"
featVec[featVec == "14"] = "Cxcr2 B-cells"


obj.pos$label=featVec


###########
deResTT = makeDEResults(obj.pos, group.by="seurat_clusters", assay="RNA", test="wilcox")
write.table(deResTT,"/mnt/raidexttmp/Alejandro/Larissa_10X/Seurat/expr_test_clusters_t.tsv", sep="\t", row.names=F, quote = F)
write_xlsx(deResTT, "/mnt/raidexttmp/Alejandro/Larissa_10X/Seurat/expr_test_clusters_t.xlsx")


#deResTT<-read_tsv("/mnt/raidexttmp/Alejandro/Larissa_10X/Seurat/expr_test_clusters_t.tsv")
DefaultAssay(obj.pos) <- "RNA"
markers.use.tt= subset(deResTT , avg_log2FC>0&p_val_adj<0.05&!startsWith(gene, "mt-")&!startsWith(gene, "rp"))

# Ensure clusters are ordered numerically
# Convert clusterID to numeric-ordered factor
markers.use.tt$clusterID <- factor(
  markers.use.tt$clusterID,
  levels = sort(as.numeric(unique(markers.use.tt$clusterID)))
)

finalMarkers.use.tt = markers.use.tt %>% arrange(p_val_adj, desc(abs(pct.1)*abs(avg_log2FC))) %>% group_by(clusterID) %>% dplyr::slice(1:15) 
finalMarkers.use.tt


data_dupli= finalMarkers.use.tt[!duplicated(finalMarkers.use.tt[ , "gene"]), ]


events= data_dupli %>% count(clusterID)
inser=cumsum(events$n)+0.5-events$n
insert=replace(inser, inser==0.5, 0)

xmi<- insert
xmin<- xmi[c(FALSE, TRUE)]
xma<- insert+events$n
xmax<- xma[c(FALSE, TRUE)]
ymi<- 0*events$n
ymin<- ymi[c(FALSE, TRUE)]
yma<- rep(length(events$n)+0.5, each=length(events$n))
ymax<- yma[c(FALSE, TRUE)]

 


p_dp_genes_idents = DotPlot(obj.pos, features = data_dupli$gene, assay="RNA", dot.scale = 5, group.by="seurat_clusters")+
    theme(axis.text.x = element_text(angle = 90, hjust = 1.0, vjust = 0.5))+
    annotate("rect", xmin=xmin, xmax=xmax, ymin=ymin , ymax=ymax, alpha=0.2, fill="blue") #rep(c("blue", "grey"), times= length(events$n)/2)
save_plot(p_dp_genes_idents, "/mnt/raidexttmp/Alejandro/Larissa_10X/Seurat/dotplot_cluster_genes_colored", 40, 8)








#####
DefaultAssay(obj.pos) = "RNA"
Feat= FeaturePlot(obj.pos, features = ("Ly6g"), pt.size = 0.5, max.cutoff=2)
save_plot(Feat, "/mnt/raidexttmp/Alejandro/Larissa_10X/Seurat/Featureplot_Ly6g", fig.width=15, fig.height=10)

DefaultAssay(obj.pos) = "RNA"
Feat= FeaturePlot(obj.pos, features = ("Ly6g"), pt.size = 0.5, split.by = "Pheno",  max.cutoff=2)
save_plot(Feat, "/mnt/raidexttmp/Alejandro/Larissa_10X/Seurat/Featureplot_Ly6g_split", fig.width=15*3, fig.height=10)




###### DEanalysis and volcanos

####################

deGroup = "Pheno"

retList_master = list()
retListW_master = list()

condlist <- list(
  c("Control", "5day_tumor"),
  c("Control", "10day_tumor"),
  c("5day_tumor", "10day_tumor")
)

for (x in condlist)
{
  allConds = x  
  print("all conds")
  print(allConds)

  all.cells = list()
  for (cond in allConds) {
    all.cells[[cond]] = cellIDForClusters(obj.pos, deGroup, c(cond))
  }

  for (i in 1:(length(allConds)-1))
  {
    for (j in (i+1):length(allConds))
    {
      print(paste(i,"<->",j))

      condI = allConds[i]
      condJ = allConds[j]

      condNameI = tolower(str_replace_all(condI, " ", "_"))
      condNameJ = tolower(str_replace_all(condJ, " ", "_"))
      compName = paste(condNameJ, condNameI, sep="_")

      message("Running DE for: ", compName)

      deMarkers = compareCellsByCluster(obj.pos, all.cells[[condJ]], all.cells[[condI]],
                                        condNameJ, condNameI, group.by="label",
                                        outfolder="/mnt/raidexttmp/Alejandro/Larissa_10X/Seurat/de", fcCutoff=0.1)

      makeVolcanos(deMarkers, paste("DE", condJ, "vs", condI),
                   paste("/mnt/raidexttmp/Alejandro/Larissa_10X/Seurat/de_volcano", compName, sep="/"),
                   turnExpression=F, FCcutoff=0.1, pCutoff=0.05)

      deMarkersW = compareCellsByCluster(obj.pos, all.cells[[condJ]], all.cells[[condI]],
                                         condNameJ, condNameI, group.by="label",
                                         outfolder="/mnt/raidexttmp/Alejandro/Larissa_10X/Seurat/dewilcox", test="wilcox", fcCutoff=0.1)

      makeVolcanos(deMarkersW, paste("DE", condJ, "vs", condI),
                   paste("/mnt/raidexttmp/Alejandro/Larissa_10X/Seurat/dewilcox_volcano", compName, sep="/"),
                   turnExpression=F, FCcutoff=0.1, pCutoff=0.05)

      # Store results in master list
      retList_master[[compName]] = deMarkers
      retListW_master[[compName]] = deMarkersW
    }
  }
}

# Save combined results
saveRDS(retList_master, file="/mnt/raidexttmp/Alejandro/Larissa_10X/Seurat/de/der_Clusters.rds")
saveRDS(retListW_master, file="/mnt/raidexttmp/Alejandro/Larissa_10X/Seurat/dewilcox/der_Clusters.rds")


#Violins Marie Kupffer cells Adrenoceptors
KC=obj.pos[,obj.pos$Pheno=="Control"]
KC=KC[,KC$seurat_clusters=="6"]

# Extract data for the desired features
features <- c("Adra1a", "Adra1b", "Adra1d", "Adra2a", "Adra2b", "Adra2c", "Adrb1", "Adrb2", "Adrb3")
p <- VlnPlot(
  object = KC,
  features = features,
  group.by = "label",  # or another metadata column (e.g. "Pheno", "Condition")
  pt.size = 1,                   # hides individual points for clarity
  combine = TRUE,                # one panel with all features
  assay = "RNA"
) 
# Save the figure
save_plot(
  p, "/mnt/raidexttmp/Alejandro/Marie_kupffer_cells/KC_Adreno_Violin",
  fig.width = 18,
  fig.height = 10
)


################################################
library(clusterProfiler)
library(org.Mm.eg.db)
library(org.Hs.eg.db)
library(ReactomePA)
library(DOSE)
library(enrichplot)

library(clusterProfiler)
library(org.Mm.eg.db)
library(ggplot2)
library(forcats)

# Function to save plots
save_plot <- function(plot, filename, width, height) {
  ggsave(filename = paste0(filename, ".png"), plot = plot, width = width, height = height)
}

# Load your DE results
inFile = "/mnt/raidexttmp/Alejandro/Larissa_10X/Seurat/dewilcox/der_Clusters.rds"
deList = readRDS(inFile)

baseDir = dirname(inFile)

# Loop over comparisons
for (comparison in names(deList)) {
  
  # Loop over clusters
  for (cluster_id in names(deList[[comparison]])) {
    resDF = deList[[comparison]][[cluster_id]]
    
    # Skip if no data
    if (is.null(resDF) || nrow(resDF) == 0) {
      message("Skipping ", comparison, " cluster ", cluster_id, " (empty)")
      next
    }
    
    # Ensure gene column exists
    if (!"gene" %in% colnames(resDF)) {
      resDF$gene = rownames(resDF)
    }
    
    # Filter significant genes
    siggenes = resDF[resDF$p_val_adj < 0.05, ]
    if (nrow(siggenes) == 0) {
      message("Skipping ", comparison, " cluster ", cluster_id, " (no significant genes)")
      next
    }
    
    # Convert SYMBOL to ENTREZID
    geneNames = bitr(siggenes$gene, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)
    siggenes = siggenes[siggenes$gene %in% geneNames$SYMBOL, ]
    
    # Create gene vector for GSEA
    geneVector = siggenes$avg_log2FC
    names(geneVector) = as.numeric(geneNames$ENTREZID[match(siggenes$gene, geneNames$SYMBOL)])
    
    # Remove NA and sort
    geneVector = geneVector[!is.na(geneVector)]
    geneVector = sort(geneVector, decreasing = TRUE)
    
    if (length(geneVector) < 10) {
      message("Skipping ", comparison, " cluster ", cluster_id, " (too few genes)")
      next
    }
    
    # Run GSEA
    gseGOres = tryCatch(
      gseGO(geneList = geneVector,
            OrgDb = org.Mm.eg.db,
            ont = "BP",
            pvalueCutoff = 0.05,
            pAdjustMethod = "BH",
            maxGSSize = 2000),
      error = function(e) {
        message("Error in GSEA for ", comparison, " cluster ", cluster_id, ": ", e$message)
        return(NULL)
      }
    )
    
    if (is.null(gseGOres) || nrow(gseGOres) == 0) {
      message("No GSEA results for ", comparison, " cluster ", cluster_id)
      next
    }
    
    # Prepare dataframe for plotting
    yfinal = data.frame(
      NES = gseGOres$NES,
      qvalues = gseGOres$qvalue,
      Description = gseGOres$Description
    )
    yfinal = yfinal[order(yfinal$NES), ]
    
    # Plot
    p = ggplot(yfinal, aes(NES, fct_reorder(Description, NES))) +
      geom_segment(aes(xend = 0, yend = Description)) +
      geom_point(aes(color = qvalues, size = NES)) +
      scale_color_viridis_c(guide = guide_colorbar(reverse = TRUE)) +
      scale_size_continuous(range = c(1, 5)) +
      theme_minimal() +
      xlab("NES") +
      ylab(NULL) +
      ggtitle(paste("GSEA enrichment for", comparison, "Cluster", cluster_id))
    
    # Save plot
    outFile = file.path(baseDir, paste0("GSEA_", comparison, "_Cl", cluster_id))
    save_plot(p, outFile, 8, 6)
    
    message("Saved plot: ", outFile)
  }
}





######################################################
#########Neutrophil subclustering
######################################################

Neut=obj.pos[,obj.pos$label=="Neutrophils"]
Neut <- FindNeighbors(Neut, dims = 1:30, reduction = "igpca")
DefaultAssay(Neut) <- "integratedgex"
Neut <- FindClusters(Neut, resolution = 0.2, algorithm = 4)
Neut <- RunUMAP(Neut, dims = 1:30, reduction = "igpca")

p=DimPlot(Neut, pt.size = 1, label=T, reduction = "umap")
save_plot(p, paste("/mnt/raidexttmp/Alejandro/Larissa_10X/Seurat", "Dimplot_Neutrophil_umap_ident", sep="/"), 12, 10)


###########
deResTT = makeDEResults(Neut, group.by="seurat_clusters", assay="RNA", test="wilcox")
write.table(deResTT,"/mnt/raidexttmp/Alejandro/Larissa_10X/Seurat/expr_Neut_test_clusters_t.tsv", sep="\t", row.names=F, quote = F)
write_xlsx(deResTT, "/mnt/raidexttmp/Alejandro/Larissa_10X/Seurat/expr_Neut_test_clusters_t.xlsx")


#deResTT<-read_tsv("/mnt/raidexttmp/Alejandro/Larissa_10X/Seurat/expr_Neut_test_clusters_t.tsv")
DefaultAssay(Neut) <- "RNA"
markers.use.tt= subset(deResTT , avg_log2FC>0&p_val_adj<0.05&!startsWith(gene, "mt-")&!startsWith(gene, "rp"))

# Ensure clusters are ordered numerically
# Convert clusterID to numeric-ordered factor
markers.use.tt$clusterID <- factor(
  markers.use.tt$clusterID,
  levels = sort(as.numeric(unique(markers.use.tt$clusterID)))
)

finalMarkers.use.tt = markers.use.tt %>% arrange(p_val_adj, desc(abs(pct.1)*abs(avg_log2FC))) %>% group_by(clusterID) %>% dplyr::slice(1:15) 
finalMarkers.use.tt


data_dupli= finalMarkers.use.tt[!duplicated(finalMarkers.use.tt[ , "gene"]), ]


events= data_dupli %>% count(clusterID)
inser=cumsum(events$n)+0.5-events$n
insert=replace(inser, inser==0.5, 0)

xmi<- insert
xmin<- xmi[c(FALSE, TRUE)]
xma<- insert+events$n
xmax<- xma[c(FALSE, TRUE)]
ymi<- 0*events$n
ymin<- ymi[c(FALSE, TRUE)]
yma<- rep(length(events$n)+0.5, each=length(events$n))
ymax<- yma[c(FALSE, TRUE)]

 


p_dp_genes_idents = DotPlot(Neut, features = data_dupli$gene, assay="RNA", dot.scale = 5, group.by="seurat_clusters")+
    theme(axis.text.x = element_text(angle = 90, hjust = 1.0, vjust = 0.5))+
    annotate("rect", xmin=xmin, xmax=xmax, ymin=ymin , ymax=ymax, alpha=0.2, fill="blue") #rep(c("blue", "grey"), times= length(events$n)/2)
save_plot(p_dp_genes_idents, "/mnt/raidexttmp/Alejandro/Larissa_10X/Seurat/dotplot_cluster_Neutro_genes_colored", 16, 5)



p=DimPlot(Neut, pt.size = 1, label=F, reduction = "umap", group.by="Pheno", split.by="Pheno")
save_plot(p, paste("/mnt/raidexttmp/Alejandro/Larissa_10X/Seurat", "Neutrophil_umap_split", sep="/"), 10*3, 10)


# Extract metadata
metadata <- Neut@meta.data

# Get Seurat color palette (same as DimPlot)
# Get colors from your existing UMAP plot
umap_plot <- DimPlot(Neut, group.by = "seurat_clusters", shuffle = TRUE)
umap_data <- ggplot_build(umap_plot)$data[[1]]

# Extract mapping between cluster identity and color
umap_color_map <- data.frame(
  cluster = umap_data$group,
  color = umap_data$colour
) %>%
  dplyr::distinct(cluster, color)

# Convert cluster to same type as in metadata (to ensure consistent sorting)
umap_color_map$cluster <- as.character(umap_color_map$cluster)

# Get ordered cluster levels from Seurat object
cluster_levels <- levels(Neut$seurat_clusters)

# Reorder colors to match cluster levels
ordered_colors <- umap_color_map$color[match(cluster_levels, umap_color_map$cluster)]
names(ordered_colors) <- cluster_levels


# Count and compute percentage per Pheno
metadata_summary <- metadata %>%
  group_by(Pheno, seurat_clusters) %>%
  summarise(cell_count = n(), .groups = 'drop') %>%
  group_by(Pheno) %>%
  mutate(
    percent = cell_count / sum(cell_count) * 100,
    percent_label = paste0(round(percent, 1), "%")
  )

# Plot using the Seurat palette
p <- ggplot(metadata_summary, aes(x = Pheno, y = percent, fill = seurat_clusters)) +
  geom_bar(stat = "identity", position = "stack", color = "white") +
  geom_text(aes(label = percent_label),
            color = "black", size = 3.2, position = position_stack(vjust = 0.5)) +
  scale_fill_manual(values = ordered_colors) +
  labs(
    title = "Cellular Composition per Phenotype",
    x = "Phenotype (Condition)",
    y = "Percentage of Cells per Cluster",
    fill = "Seurat Cluster"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
    axis.title = element_text(face = "bold"),
    legend.position = "right",
    panel.grid = element_blank()
  )

# Save the plot
save_plot(
  p,
  "/mnt/raidexttmp/Alejandro/Larissa_10X/Seurat/Barplot_cluster_Neutrophils",
  fig.width = 8,
  fig.height = 10
)
 


###### DEanalysis and volcanos


####################

deGroup = "Pheno"

retList_master = list()
retListW_master = list()

condlist <- list(
  c("Control", "5day_tumor"),
  c("Control", "10day_tumor"),
  c("5day_tumor", "10day_tumor")
)

for (x in condlist)
{
  allConds = x  
  print("all conds")
  print(allConds)

  all.cells = list()
  for (cond in allConds) {
    all.cells[[cond]] = cellIDForClusters(Neut, deGroup, c(cond))
  }

  for (i in 1:(length(allConds)-1))
  {
    for (j in (i+1):length(allConds))
    {
      print(paste(i,"<->",j))

      condI = allConds[i]
      condJ = allConds[j]

      condNameI = tolower(str_replace_all(condI, " ", "_"))
      condNameJ = tolower(str_replace_all(condJ, " ", "_"))
      compName = paste(condNameJ, condNameI, sep="_")

      message("Running DE for: ", compName)

      deMarkers = compareCellsByCluster(Neut, all.cells[[condJ]], all.cells[[condI]],
                                        condNameJ, condNameI, group.by="seurat_clusters",
                                        outfolder="/mnt/raidexttmp/Alejandro/Larissa_10X/Seurat/de/Neutrophils", fcCutoff=0.1)

      makeVolcanos(deMarkers, paste("DE", condJ, "vs", condI),
                   paste("/mnt/raidexttmp/Alejandro/Larissa_10X/Seurat/de_volcano/Neutrophils", compName, sep="/"),
                   turnExpression=F, FCcutoff=0.1, pCutoff=0.05)

      deMarkersW = compareCellsByCluster(Neut, all.cells[[condJ]], all.cells[[condI]],
                                         condNameJ, condNameI, group.by="seurat_clusters",
                                         outfolder="/mnt/raidexttmp/Alejandro/Larissa_10X/Seurat/dewilcox/Neutrophils", test="wilcox", fcCutoff=0.1)

      makeVolcanos(deMarkersW, paste("DE", condJ, "vs", condI),
                   paste("/mnt/raidexttmp/Alejandro/Larissa_10X/Seurat/dewilcox_volcano/Neutrophils", compName, sep="/"),
                   turnExpression=F, FCcutoff=0.1, pCutoff=0.05)

      # Store results in master list
      retList_master[[compName]] = deMarkers
      retListW_master[[compName]] = deMarkersW
    }
  }
}

# Save combined results
saveRDS(retList_master, file="/mnt/raidexttmp/Alejandro/Larissa_10X/Seurat/de/der_Neut_Clusters.rds")
saveRDS(retListW_master, file="/mnt/raidexttmp/Alejandro/Larissa_10X/Seurat/dewilcox/der_Neut_Clusters.rds")


############Scores for adhesion, cytotoxicity, etc

cytotoxicity <- c(
  "Cybb","Ncf1","Ncf2","Ncf4","Rac2","Mpo","Elane","Prtn3","Ctsg",
  "Ltf","Camp","Lcn2","Mmp8","Mmp9","Pglyrp1","Padi4","S100a8","S100a9"
)

adhesion <- c(
  "Itgal","Itgam","Itgb2","Selplg","Sell","Tln1","Fermt3","Rap1b","Vav1",
  "Syk","Fgr","Hck","Fcer1g","Marcksl1","Arhgap15","Arhgdib",
  "Cxcr2","Ccrl2","Cxcr4"
)

suppression <- c(
  "Cd274","Il1rn","Isg15","Ifit1","Ifit3","Oas1a","Chil3","Il1r2",
  "Wfdc21","Arg1","Plin2","Bcl2l1","G0s2","Morrbid","Ccrl2","Cd9"
)

DefaultAssay(Neut) <- "RNA"

cytotoxicity <- intersect(cytotoxicity, rownames(Neut))
adhesion     <- intersect(adhesion, rownames(Neut))
suppression  <- intersect(suppression, rownames(Neut))


Neut <- AddModuleScore(
  Neut,
  features = list(cytotoxicity),
  name = "Cytotoxicity"
)

Neut <- AddModuleScore(
  Neut,
  features = list(adhesion),
  name = "Adhesion"
)

Neut <- AddModuleScore(
  Neut,
  features = list(suppression),
  name = "Suppression"
)

p = VlnPlot(
  Neut,
  features = c("Cytotoxicity1","Adhesion1","Suppression1"),
  group.by = "Pheno"
)

save_plot(
  p,
  "/mnt/raidexttmp/Alejandro/Larissa_10X/Seurat/Neutrophil_scores_Cytotoxicity_Adhesion_Suppression",
  fig.width = 18,
  fig.height = 6
)


cellList = colnames(Neut)

featVec <- vector(mode="character", length=length(cellList))
featVec = as.character(Neut$CSclassification)

featVec[featVec == "C1"] = "Control 1"
featVec[featVec == "C2"] = "Control 2"
featVec[featVec == "C3"] = "Control 3"
featVec[featVec == "C4"] = "Control 4"
featVec[featVec == "C5"] = "Control 5"
featVec[featVec == "C6"] = "Control 6"
featVec[featVec == "C7"] = "Control 7"
featVec[featVec == "C8"] = "Control 8"
featVec[featVec == "S1"] = "5day_tumor 1"
featVec[featVec == "S2"] = "5day_tumor 2"
featVec[featVec == "S3"] = "5day_tumor 3"
featVec[featVec == "S4"] = "5day_tumor 4"
featVec[featVec == "S5"] = "5day_tumor 5"
featVec[featVec == "S6"] = "5day_tumor 6"
featVec[featVec == "S8"] = "5day_tumor 7"
featVec[featVec == "S9"] = "5day_tumor 9"
featVec[featVec == "L1"] = "10day_tumor 1"
featVec[featVec == "L2"] = "10day_tumor 2"
featVec[featVec == "L3"] = "10day_tumor 3"
featVec[featVec == "L4"] = "10day_tumor 4"
featVec[featVec == "L5"] = "10day_tumor 5"
featVec[featVec == "L6"] = "10day_tumor 6"
featVec[featVec == "L7"] = "10day_tumor 7"
featVec[featVec == "L8"] = "10day_tumor 8"


Neut$PhenoMouse=featVec


###
cellList = colnames(Neut)

featVec <- vector(mode="character", length=length(cellList))
featVec = as.character(Neut$integratedgex_snn_res.0.2)

featVec[featVec == "1"] = "Neut_subcluster_1"
featVec[featVec == "2"] = "Neut_subcluster_2"
featVec[featVec == "3"] = "Neut_subcluster_3"

Neut$Neutsubcluster=featVec


library(dplyr)

df_mouse <- data.frame(
  Mouse = Neut$PhenoMouse,
  Pheno = Neut$Pheno,
  Cytotoxicity = Neut$Cytotoxicity1,
  Adhesion     = Neut$Adhesion1,
  Suppression  = Neut$Suppression1
) %>%
  group_by(Mouse, Pheno) %>%
  summarise(
    Cytotoxicity_mean = mean(Cytotoxicity, na.rm = TRUE),
    Adhesion_mean     = mean(Adhesion, na.rm = TRUE),
    Suppression_mean  = mean(Suppression, na.rm = TRUE),
    .groups = "drop"
  )

df_mouse$Pheno <- factor(
  df_mouse$Pheno,
  levels = c("Control", "5day_tumor", "10day_tumor")
)


cols <- c(
  Control      = "#2166AC",
  `5day_tumor`  = "#B2182B",
  `10day_tumor` = "#1B9E77"
)

wilcox_all <- function(df, variable) {

  p1 <- wilcox.test(
    as.formula(paste(variable, "~ Pheno")),
    data = subset(df, Pheno %in% c("Control", "5day_tumor"))
  )$p.value

  p2 <- wilcox.test(
    as.formula(paste(variable, "~ Pheno")),
    data = subset(df, Pheno %in% c("Control", "10day_tumor"))
  )$p.value

  p3 <- wilcox.test(
    as.formula(paste(variable, "~ Pheno")),
    data = subset(df, Pheno %in% c("5day_tumor", "10day_tumor"))
  )$p.value

  c(
    ctrl_5d  = p1,
    ctrl_10d = p2,
    d5_d10   = p3
  )
}


p_cyto <- wilcox_all(df_mouse, "Cytotoxicity_mean")
p_adh  <- wilcox_all(df_mouse, "Adhesion_mean")
p_sup  <- wilcox_all(df_mouse, "Suppression_mean")

make_labs <- function(pvals) {
  paste0(
    c("Ctrl vs 5d: ",
      "Ctrl vs 10d: ",
      "5d vs 10d: "),
    format.pval(pvals, digits = 3)
  )
}

library(ggplot2)

add_sig_bars <- function(pvals, ymax) {

  data.frame(
    x1 = c(1, 1, 2),
    x2 = c(2, 3, 3),

    # vertical placement of the bars
    y = ymax * c(1.20, 1.45, 1.70),

    label = paste0(
      "Wilcoxon p = ",
      format.pval(pvals, digits = 3)
    )
  )
}

plot_module <- function(df, variable, title, pvals) {

  ymax <- max(df[[variable]], na.rm = TRUE)

  sig_df <- add_sig_bars(pvals, ymax)

  ggplot(df, aes(x = Pheno, y = .data[[variable]], fill = Pheno)) +

    geom_violin(trim = FALSE, alpha = 0.4) +

    geom_boxplot(
      width = 0.15,
      outlier.shape = NA,
      alpha = 0.7
    ) +

    geom_jitter(
      aes(color = Pheno),
      width = 0.1,
      size = 3
    ) +

    ## horizontal significance lines
    geom_segment(
      data = sig_df,
      aes(x = x1, xend = x2, y = y, yend = y),
      inherit.aes = FALSE,
      linewidth = 0.8
    ) +

    ## p-value labels
    geom_text(
      data = sig_df,
      aes(x = (x1 + x2) / 2, y = y*1.05, label = label),
      inherit.aes = FALSE,
      size = 5,
      fontface = "bold"
    ) +

    scale_fill_manual(values = cols) +
    scale_color_manual(values = cols) +

    ggtitle(title) +
    ylab("Mean module score per mouse") +

    theme_classic(base_size = 16) +
    theme(
      legend.position = "none",
      plot.title = element_text(face = "bold")
    ) 
}



p1 <- plot_module(
  df_mouse,
  "Cytotoxicity_mean",
  "Neutrophil cytotoxicity score",
  p_cyto
)

save_plot(
  p1,
  "/mnt/raidexttmp/Alejandro/Larissa_10X/Seurat/Neutrophil_cytotoxicity_scores",
  fig.width = 10,
  fig.height = 6
)

p2 <- plot_module(
  df_mouse,
  "Adhesion_mean",
  "Neutrophil adhesion score",
  p_adh
)
save_plot(
  p2,
  "/mnt/raidexttmp/Alejandro/Larissa_10X/Seurat/Neutrophil_Adhesion_scores",
  fig.width = 10,
  fig.height = 6
)

p3 <- plot_module(
  df_mouse,
  "Suppression_mean",
  "Neutrophil suppressive score",
  p_sup
)
save_plot(
  p3,
  "/mnt/raidexttmp/Alejandro/Larissa_10X/Seurat/Neutrophil_Suppression_scores",
  fig.width = 10,
  fig.height = 6
)

##################
#Neutrophil trajectories (pseudotime)
##################
library(SeuratWrappers)
library(monocle3)

cds <- as.cell_data_set(Neut)

colData(cds)$Neutsubcluster <- Neut$Neutsubcluster
colData(cds)$Pheno <- Neut$Pheno

reducedDims(cds)$UMAP <- Neut@reductions$umap@cell.embeddings
cds <- cluster_cells(
  cds,
  reduction_method = "UMAP",
  cluster_method = "leiden",
  resolution = 0.0005
)


cds <- learn_graph(
  cds,
  use_partition = F
)

root_cells <- colnames(Neut)[Neut$Neutsubcluster == "Neut_subcluster_3"]

cds <- order_cells(
  cds,
  reduction_method = "UMAP",
  root_cells = root_cells
)

p=plot_cells(
  cds,
  color_cells_by = "Neutsubcluster",
  label_groups_by_cluster = F,
  label_cell_groups = T,
  label_leaves = T,
  label_branch_points = T,
  label_roots = FALSE,
  group_label_size = 3
)
save_plot(p, "/mnt/raidexttmp/Alejandro/Larissa_10X/Seurat/Neutrophil_pseudotime_trajectory", 7, 6)






#################################
#Macrophage subclustering
##################################

Macro=obj.pos[,obj.pos$label %in% c("Macrophages")]
Macro <- FindNeighbors(Macro, dims = 1:30, reduction = "igpca")
DefaultAssay(Macro) <- "integratedgex"
Macro <- FindClusters(Macro, resolution = 0.2, algorithm = 4)
Macro <- RunUMAP(Macro, dims = 1:30, reduction = "igpca")

p=DimPlot(Macro, pt.size = 1, label=T, reduction = "umap")
save_plot(p, paste("/mnt/raidexttmp/Alejandro/Larissa_10X/Seurat", "Dimplot_Macro_umap_ident", sep="/"), 12, 10)


###########
deResTT = makeDEResults(Macro, group.by="seurat_clusters", assay="RNA", test="wilcox")
write.table(deResTT,"/mnt/raidexttmp/Alejandro/Larissa_10X/Seurat/expr_Macro_test_clusters_t.tsv", sep="\t", row.names=F, quote = F)
write_xlsx(deResTT, "/mnt/raidexttmp/Alejandro/Larissa_10X/Seurat/expr_Macro_test_clusters_t.xlsx")


#deResTT<-read_tsv("/mnt/raidexttmp/Alejandro/Larissa_10X/Seurat/expr_Macro_test_clusters_t.tsv")
DefaultAssay(Macro) <- "RNA"
markers.use.tt= subset(deResTT , avg_log2FC>0&p_val_adj<0.05&!startsWith(gene, "mt-")&!startsWith(gene, "rp"))

# Ensure clusters are ordered numerically
# Convert clusterID to numeric-ordered factor
markers.use.tt$clusterID <- factor(
  markers.use.tt$clusterID,
  levels = sort(as.numeric(unique(markers.use.tt$clusterID)))
)

finalMarkers.use.tt = markers.use.tt %>% arrange(p_val_adj, desc(abs(pct.1)*abs(avg_log2FC))) %>% group_by(clusterID) %>% dplyr::slice(1:15) 
finalMarkers.use.tt


data_dupli= finalMarkers.use.tt[!duplicated(finalMarkers.use.tt[ , "gene"]), ]


events= data_dupli %>% count(clusterID)
inser=cumsum(events$n)+0.5-events$n
insert=replace(inser, inser==0.5, 0)

xmi<- insert
xmin<- xmi[c(FALSE, TRUE)]
xma<- insert+events$n
xmax<- xma[c(FALSE, TRUE)]
ymi<- 0*events$n
ymin<- ymi[c(FALSE, TRUE)]
yma<- rep(length(events$n)+0.5, each=length(events$n))
ymax<- yma[c(FALSE, TRUE)]

 


p_dp_genes_idents = DotPlot(Macro, features = data_dupli$gene, assay="RNA", dot.scale = 5, group.by="seurat_clusters")+
    theme(axis.text.x = element_text(angle = 90, hjust = 1.0, vjust = 0.5))+
    annotate("rect", xmin=xmin, xmax=xmax, ymin=ymin , ymax=ymax, alpha=0.2, fill="blue") #rep(c("blue", "grey"), times= length(events$n)/2)
save_plot(p_dp_genes_idents, "/mnt/raidexttmp/Alejandro/Larissa_10X/Seurat/dotplot_cluster_Macro_genes_colored", 16, 5)



p=DimPlot(Macro, pt.size = 1, label=F, reduction = "umap", group.by="Pheno", split.by="Pheno")
save_plot(p, paste("/mnt/raidexttmp/Alejandro/Larissa_10X/Seurat", "Macro_umap_split", sep="/"), 10*3, 10)

DefaultAssay(Macro) = "RNA"
Feat= FeaturePlot(Macro, features = ("Cd68"), pt.size = 1, max.cutoff=2)
save_plot(Feat, "/mnt/raidexttmp/Alejandro/Larissa_10X/Seurat/Featureplot_CD68_Macro", fig.width=15, fig.height=10)

DefaultAssay(Macro) = "RNA"
Feat= FeaturePlot(Macro, features = ("Adgre1"), pt.size = 1, max.cutoff=2)
save_plot(Feat, "/mnt/raidexttmp/Alejandro/Larissa_10X/Seurat/Featureplot_F4-80_Macro", fig.width=15, fig.height=10)


# Extract metadata
metadata <- Macro@meta.data

# Get Seurat color palette (same as DimPlot)
# Get colors from your existing UMAP plot
umap_plot <- DimPlot(Macro, group.by = "seurat_clusters", shuffle = TRUE)
umap_data <- ggplot_build(umap_plot)$data[[1]]

# Extract mapping between cluster identity and color
umap_color_map <- data.frame(
  cluster = umap_data$group,
  color = umap_data$colour
) %>%
  dplyr::distinct(cluster, color)

# Convert cluster to same type as in metadata (to ensure consistent sorting)
umap_color_map$cluster <- as.character(umap_color_map$cluster)

# Get ordered cluster levels from Seurat object
cluster_levels <- levels(Macro$seurat_clusters)

# Reorder colors to match cluster levels
ordered_colors <- umap_color_map$color[match(cluster_levels, umap_color_map$cluster)]
names(ordered_colors) <- cluster_levels


# Count and compute percentage per Pheno
metadata_summary <- metadata %>%
  group_by(Pheno, seurat_clusters) %>%
  summarise(cell_count = n(), .groups = 'drop') %>%
  group_by(Pheno) %>%
  mutate(
    percent = cell_count / sum(cell_count) * 100,
    percent_label = paste0(round(percent, 1), "%")
  )

# Plot using the Seurat palette
p <- ggplot(metadata_summary, aes(x = Pheno, y = percent, fill = seurat_clusters)) +
  geom_bar(stat = "identity", position = "stack", color = "white") +
  geom_text(aes(label = percent_label),
            color = "black", size = 3.2, position = position_stack(vjust = 0.5)) +
  scale_fill_manual(values = ordered_colors) +
  labs(
    title = "Cellular Composition per Phenotype",
    x = "Phenotype (Condition)",
    y = "Percentage of Cells per Cluster",
    fill = "Seurat Cluster"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
    axis.title = element_text(face = "bold"),
    legend.position = "right",
    panel.grid = element_blank()
  )

# Save the plot
save_plot(
  p,
  "/mnt/raidexttmp/Alejandro/Larissa_10X/Seurat/Barplot_cluster_Macrophages",
  fig.width = 8,
  fig.height = 10
)
 


###### DEanalysis and volcanos


####################

deGroup = "Pheno"

retList_master = list()
retListW_master = list()

condlist <- list(
  c("Control", "5day_tumor"),
  c("Control", "10day_tumor"),
  c("5day_tumor", "10day_tumor")
)

for (x in condlist)
{
  allConds = x  
  print("all conds")
  print(allConds)

  all.cells = list()
  for (cond in allConds) {
    all.cells[[cond]] = cellIDForClusters(Macro, deGroup, c(cond))
  }

  for (i in 1:(length(allConds)-1))
  {
    for (j in (i+1):length(allConds))
    {
      print(paste(i,"<->",j))

      condI = allConds[i]
      condJ = allConds[j]

      condNameI = tolower(str_replace_all(condI, " ", "_"))
      condNameJ = tolower(str_replace_all(condJ, " ", "_"))
      compName = paste(condNameJ, condNameI, sep="_")

      message("Running DE for: ", compName)

      deMarkers = compareCellsByCluster(Macro, all.cells[[condJ]], all.cells[[condI]],
                                        condNameJ, condNameI, group.by="seurat_clusters",
                                        outfolder="/mnt/raidexttmp/Alejandro/Larissa_10X/Seurat/de/Macrophages", fcCutoff=0.1)

      makeVolcanos(deMarkers, paste("DE", condJ, "vs", condI),
                   paste("/mnt/raidexttmp/Alejandro/Larissa_10X/Seurat/de_volcano/Macrophages", compName, sep="/"),
                   turnExpression=F, FCcutoff=0.1, pCutoff=0.05)

      deMarkersW = compareCellsByCluster(Macro, all.cells[[condJ]], all.cells[[condI]],
                                         condNameJ, condNameI, group.by="seurat_clusters",
                                         outfolder="/mnt/raidexttmp/Alejandro/Larissa_10X/Seurat/dewilcox/Macrophages", test="wilcox", fcCutoff=0.1)

      makeVolcanos(deMarkersW, paste("DE", condJ, "vs", condI),
                   paste("/mnt/raidexttmp/Alejandro/Larissa_10X/Seurat/dewilcox_volcano/Macrophages", compName, sep="/"),
                   turnExpression=F, FCcutoff=0.1, pCutoff=0.05)

      # Store results in master list
      retList_master[[compName]] = deMarkers
      retListW_master[[compName]] = deMarkersW
    }
  }
}

# Save combined results
saveRDS(retList_master, file="/mnt/raidexttmp/Alejandro/Larissa_10X/Seurat/de/der_Macro_Clusters.rds")
saveRDS(retListW_master, file="/mnt/raidexttmp/Alejandro/Larissa_10X/Seurat/dewilcox/der_Macro_Clusters.rds")






######################################################
#########Cell Chat v2
######################################################
#.libPaths("~/R/x86_64-suse-linux-gnu-library/4.3")
#devtools::install_github("jinworks/CellChat")
library(CellChat)

obj = obj.pos

cellchatlist = list()
for (x in c("Control", "5day_tumor", "10day_tumor"))
{
obj.i = obj[,obj$Pheno == x]
data.input <- obj.i[["RNA"]]$data # normalized data matrix
Idents(obj.i) <- "label"
labels <- Idents(obj.i) # cell labels
meta <- data.frame(labels = labels, row.names = names(labels)) # create a dataframe of the cell labels

cellchat <- createCellChat(object = obj.i, group.by = "ident", assay = "RNA")

CellChatDB <- CellChatDB.mouse # use CellChatDB.mouse if running on mouse data
cellchat@DB <- CellChatDB

# subset the expression data of signaling genes for saving computation cost
cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database
future::plan("sequential") 
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)

cellchat <- computeCommunProb(cellchat, type = "triMean")
cellchat <- filterCommunication(cellchat, min.cells = 10)

cellchat <- computeCommunProbPathway(cellchat)

# Aggregate
cellchatlist[[x]] <- aggregateNet(cellchat)
}



# Extract group names and group sizes (same for all)
groups <- rownames(cellchatlist[[1]]@net$weight)
ng <- length(groups)
groupSize <- as.numeric(table(cellchatlist[[1]]@idents))

# Compute a global max weight across all three objects for consistent scaling
max.weight <- max(sapply(cellchatlist, function(x) max(x@net$weight)))

# --- OUTPUT PDF ---
pdf("/mnt/raidexttmp/Alejandro/Larissa_10X/Seurat/Cellchat_outgoing_3timepoints.pdf",
    width = 15, height = 5 * ceiling(ng))

# Layout:
# Each row = one cell group
# Each column = one timepoint
par(mfrow = c(ng, length(cellchatlist)), xpd = TRUE)

# Loop over each sender group
for (i in seq_len(ng)) {
  sender <- groups[i]
  
  # Loop across each CellChat object (3 timepoints)
  for (tp in names(cellchatlist)) {
    mat <- cellchatlist[[tp]]@net$weight
    
    # Mat with outgoing edges only from this sender
    mat2 <- matrix(0, nrow = ng, ncol = ng, dimnames = dimnames(mat))
    mat2[sender, ] <- mat[sender, ]
    
    # Plot
    netVisual_circle(
      mat2,
      vertex.weight = groupSize,
      weight.scale = TRUE,
      edge.weight.max = max.weight,   # same scaling across all panels
      label.edge = FALSE,
      title.name = paste(sender, "\n", tp)
    )
  }
}

dev.off()


##############

# Prefixed cluster labels
macro_labels <- paste0("Macro_", Macro$seurat_clusters)
neut_labels  <- paste0("Neut_", Neut$seurat_clusters)

# Name the vectors with cell barcodes
names(macro_labels) <- colnames(Macro)
names(neut_labels)  <- colnames(Neut)

combined_labels <- c(macro_labels, neut_labels)

cells_in_obj <- intersect(names(combined_labels), colnames(obj.pos))

filtered_labels <- combined_labels[cells_in_obj]

obj.pos$newlabel <- NA
obj.pos$newlabel[cells_in_obj] <- filtered_labels

obj = obj.pos[,!is.na(obj.pos$newlabel)]


obj <- FindVariableFeatures(obj)
obj <- ScaleData(obj)
obj <- RunPCA(obj, npcs = 20)

obj <- FindNeighbors(obj, dims = 1:20)
obj <- FindClusters(obj, resolution = 0.2, algorithm = 4)
obj <- RunUMAP(obj, dims = 1:20)

p=DimPlot(obj, pt.size = 1, label=T, reduction = "umap", group.by="newlabel")
save_plot(p, paste("/mnt/raidexttmp/Alejandro/Larissa_10X/Seurat", "Dimplot_APCs&Neut_umap_ident", sep="/"), 12, 10)

cellchatlist = list()
for (x in c("Control", "5day_tumor", "10day_tumor"))
{
obj.i = obj[,obj$Pheno == x]
data.input <- obj.i[["RNA"]]$data # normalized data matrix
Idents(obj.i) <- "newlabel"
labels <- Idents(obj.i) # cell labels
meta <- data.frame(labels = labels, row.names = names(labels)) # create a dataframe of the cell labels

cellchat <- createCellChat(object = obj.i, group.by = "ident", assay = "RNA")

CellChatDB <- CellChatDB.mouse # use CellChatDB.mouse if running on mouse data
cellchat@DB <- CellChatDB

# subset the expression data of signaling genes for saving computation cost
cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database
future::plan("sequential") 
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)

cellchat <- computeCommunProb(cellchat, type = "triMean")
cellchat <- filterCommunication(cellchat, min.cells = 10)

cellchat <- computeCommunProbPathway(cellchat)

# Aggregate
cellchatlist[[x]] <- aggregateNet(cellchat)
}



# Extract group names and group sizes (same for all)
groups <- rownames(cellchatlist[[1]]@net$weight)
ng <- length(groups)
groupSize <- as.numeric(table(cellchatlist[[1]]@idents))

# Compute a global max weight across all three objects for consistent scaling
max.weight <- max(sapply(cellchatlist, function(x) max(x@net$weight)))

# --- OUTPUT PDF ---
pdf("/mnt/raidexttmp/Alejandro/Larissa_10X/Seurat/Cellchat_subcluster_3timepoints.pdf",
    width = 15, height = 5 * ceiling(ng))

# Layout:
# Each row = one cell group
# Each column = one timepoint
par(mfrow = c(ng, length(cellchatlist)), xpd = TRUE)

# Loop over each sender group
for (i in seq_len(ng)) {
  sender <- groups[i]
  
  # Loop across each CellChat object (3 timepoints)
  for (tp in names(cellchatlist)) {
    mat <- cellchatlist[[tp]]@net$weight
    
    # Mat with outgoing edges only from this sender
    mat2 <- matrix(0, nrow = ng, ncol = ng, dimnames = dimnames(mat))
    mat2[sender, ] <- mat[sender, ]
    
    # Plot
    netVisual_circle(
      mat2,
      vertex.weight = groupSize,
      weight.scale = TRUE,
      edge.weight.max = max.weight,   # same scaling across all panels
      label.edge = FALSE,
      title.name = paste(sender, "\n", tp)
    )
  }
}

dev.off()


###Heatmap


outdir <- "/mnt/raidexttmp/Alejandro/Larissa_10X/Seurat/Cellchat_heatmaps_by_condition/"
dir.create(outdir, showWarnings = FALSE)

for (tp in names(cellchatlist)) {

  cat("Plotting heatmap for:", tp, "\n")

  cc <- cellchatlist[[tp]]

  outfile <- file.path(outdir, paste0("CellChat_heatmap_", tp, ".pdf"))
  pdf(outfile, width = 12, height = 8)

  p <- netVisual_heatmap(
    cc,
    measure = "count",
    color.heatmap = "Reds",
    title.name = paste("Communication Heatmap -", tp),
    cluster.rows = TRUE,
    cluster.cols = TRUE,
  )

  print(p)
  dev.off()
}


###Bubbleplot interactions

outdir <- "/mnt/raidexttmp/Alejandro/Larissa_10X/Seurat/Cellchat_interactions_by_condition/"
dir.create(outdir, showWarnings = FALSE)

for (tp in names(cellchatlist)) {

  cat("Plotting Bubble for:", tp, "\n")

  cc <- cellchatlist[[tp]]

  outfile <- file.path(outdir, paste0("CellChat_Bubble_", tp, ".pdf"))
  pdf(outfile, width = 6, height = 6)

  p <- netVisual_bubble(cc, sources.use = "Neut_2", targets.use = c("Macro_1", "Macro_2", "Macro_3"), remove.isolate = FALSE)

  print(p)
  dev.off()
}


# remove datasets with no Neut_2 → Macro_1 interactions

extract_interactions <- function(cc, source, target, name) {
  df <- subsetCommunication(cc)
  df <- df[df$source == source & df$target == target, ]
  df$condition <- name
  df
}

df.all <- map2_df(
  cellchatlist,
  names(cellchatlist),
  ~ extract_interactions(.x, "Neut_2", "Macro_2", .y)
)


pdf("/mnt/raidexttmp/Alejandro/Larissa_10X/Seurat/Cellchat_interactions_by_condition/Bubble_Comparison_Neut2_to_Macro2.pdf",
    width = 6, height = 7)

ggplot(df.all, aes(
  x = factor(condition, levels = c("Control", "5day_tumor", "10day_tumor")),
  y = interaction_name_2,
  size = prob,
  color = pval
)) +
  geom_point(alpha = 0.9) +
  scale_size(range = c(3, 10)) +
  theme_bw(base_size = 12) +
  labs(
    title = "Neut_2 → Macro_2 communication",
    x = NULL,                      # removes x-axis title
    y = "Ligand–Receptor Pair",
    size = "Probability",
    color = "P-value"
  ) +
  theme(
    axis.title.x = element_blank() # ensures x-axis title is off
  )

dev.off()


