#------------------------- palettes ----------------------

ggplotColours <- function(n = 6, h = c(0, 360) + 15){
  if ((diff(h) %% 360) < 1) h[2] <- h[2] - 360/n
  hcl(h = (seq(h[1], h[2], length = n)), c = 100, l = 65)
}

palette <- c(
  '#F0A3FF', #Amethyst
  '#0075DC', #Blue
  '#993F00', #Caramel
  '#4C005C', #Damson
  '#191919', #Ebony
  '#005C31', #Forest
  '#2BCE48', #Green
  '#FFCC99', #Honeydew
  '#808080', #Iron
  '#94FFB5', #Jade
  '#8F7C00', #Khaki
  '#9DCC00', #Lime
  '#C20088', #Mallow
  '#003380', #Navy
  '#FFA405', #Orpiment
  '#FFA8BB', #Pink
  '#426600', #Quagmire
  '#FF0010', #Red
  '#5EF1F2', #Sky
  '#00998F', #Turquoise
  '#E0FF66', #Uranium
  '#740AFF', #Violet
  '#990000', #Wine
  '#FFFF80', #Xanthin
  '#FFFF00', #yellow
  '#FF5005' #Zinnia
)

#--------------------------------- CelliD -------------------------------

mapHGT <- function(obj, HGT, name, conf=2) {
  HGT_prediction <- rownames(HGT)[apply(HGT, 2, which.max)]
  HGT_conf <- apply(HGT, 2, max)
  HGT_signif <- ifelse(HGT_conf>conf, yes = HGT_prediction, "unassigned")
  obj[[name]] <- HGT_signif
  obj[[paste0(name,'_conf')]] <- HGT_conf
  return(obj)
}
HGTtable <- function(obj, name, lev='orig.ident', norm=TRUE) {
  bp_data <- table(tissues_merged_new[[name]][,1], tissues_merged_new[[lev]][,1])
  if(norm) bp_data <- t(t(bp_data)/colSums(bp_data))
  return(bp_data)
}

filterLowFreq <- function(x, n, lab='unassigned') {
  cfreq <- table(x)
  toremove <- names(cfreq[cfreq<n])
  x[x %in% toremove] <- "unassigned"
  return(x)  
}

#---------------------------- Seurat ------------------------------------


DoAverageHeatmap <- function(object, features, range_min=-2, range_max=2, 
                             cluster_rows = TRUE,
                             cluster_palette = NULL, 
                             ann_palette=c("snow2", "royalblue4"),
                             data_palette=c("magenta", "black", "yellow"), ...
) {
  stopifnot(require(pheatmap))
  if(is.null(cluster_palette)) {
    cluster_palette <- seq(levels(Idents(object)))
    names(cluster_palette) <- levels(Idents(object))
  }
  row_data <- data.frame(expr=rowMeans(object@assays$RNA@data[features,]),
                         row.names = make.unique(features))
  avg_expr <- AverageExpression(object, 
                    features = features, assays = 'RNA', 
                    slot = 'scale.data')$RNA[features,]
  # if(!cluster_rows)
  #   avg_expr <- avg_expr[order(row_data$expr, decreasing = T),]
  col_data <- data.frame(cluster=levels(Idents(object)), row.names = colnames(avg_expr))
  out <- pheatmap(avg_expr,
                  cluster_rows = cluster_rows, 
                  annotation_row = row_data, 
                  annotation_col = col_data,
                  annotation_colors = list(
                    expr=colorRampPalette(ann_palette)(100),
                    cluster=cluster_palette
                  ),
                  color = colorRampPalette(data_palette)(100),
                  breaks = seq(range_min,range_max,length.out=101),
                  treeheight_col = 10, treeheight_row = 10, ...)
  out$avg_expr <- avg_expr
  out <- out
}

DoAverageHeatmapByGroups <- 
  function(object, features, range_min=-2, range_max=2, 
           cluster_rows = TRUE, 
           cluster_palette = NULL,
           cluster_groups = NULL,
           group_palette = NULL,
           cluster_by_group = NULL,
           ann_palette=c("snow2", "royalblue4"),
           data_palette=c("magenta", "black", "yellow"), ...
  ) {
    stopifnot(require(pheatmap))
    if(is.null(cluster_palette)) {
      cluster_palette <- seq(levels(Idents(object)))
      names(cluster_palette) <- levels(Idents(object))
    }
    row_data <- data.frame(expr=rowMeans(object@assays$RNA@data[features,]),
                           row.names = features)
    avg_expr <- AverageExpression(object, 
                                  features = features, assays = 'RNA', 
                                  slot = 'scale.data')$RNA
    if(!cluster_rows)
      avg_expr <- avg_expr[order(row_data$expr, decreasing = T),]
    col_data <- data.frame(cluster=levels(Idents(object)), row.names = colnames(avg_expr))
    ann_colors <- list(
      expr=colorRampPalette(ann_palette)(100),
      cluster=cluster_palette)
    cluster_cols <- T
    if(!is.null(cluster_groups)) {
      col_data$groups <- NA
      for(group_name in names(cluster_groups))
        col_data[cluster_groups[[group_name]],'groups'] <- group_name
      if(!is.null(group_palette))
        ann_colors$groups <- group_palette
      if(!is.null(cluster_by_group)) {
        cluster_cols <- F
        group_as_factor <- factor(col_data$groups, levels = names(cluster_groups))
        ordered_clusters <- col_data$cluster[order(group_as_factor)]
        avg_expr <- avg_expr[,ordered_clusters]
      }
    }
    out <- pheatmap(avg_expr, cluster_cols = cluster_cols,
                    cluster_rows = cluster_rows, 
                    annotation_row = row_data, 
                    annotation_col = col_data,
                    annotation_colors = ann_colors,
                    color = colorRampPalette(data_palette)(100),
                    breaks = seq(range_min,range_max,length.out=101),
                    treeheight_col = 10, treeheight_row = 10, ...)
    out$avg_expr <- avg_expr
    out <- out
  }

viewGeneList <- function(gene_list, groups_order=NULL, h=1.3, z=.8, ncol=2, label=FALSE) {
  gene_list <- intersect(rownames(glioma_integrated_RNA), gene_list)
  fp <- FeaturePlot(glioma_integrated_RNA, features = gene_list, 
                    max.cutoff = 10, label = label, ncol=ncol)
  avg_heat <- DoAverageHeatmap(glioma_integrated_RNA, gene_list, 
                               cluster_palette = glioma_integrated@misc$cols_seurat, 
                               silent=TRUE)
  if(!is.null(groups_order))
    cluster_groups <- cluster_groups[groups_order]
  avg_heat_by_group <- 
    DoAverageHeatmapByGroups(glioma_integrated_RNA, gene_list, 
                             cluster_palette = glioma_integrated@misc$cols_seurat, 
                             cluster_groups = cluster_groups,
                             group_palette = glioma_integrated@misc$cols_groups,
                             cluster_by_group = T,
                             silent=TRUE)
  cluster_split <- with(avg_heat, split(tree_col$labels, cutree(tree_col, h = h)))
  most_represented_clusters <- names(which(apply(avg_heat$avg_expr > z, 2, any)))
  selected_clusters <- unlist(cluster_split[sapply(cluster_split, function(x) any(x %in% most_represented_clusters))])
  return(list(fp=fp, 
              avg_heat=avg_heat, 
              avg_heat_by_group=avg_heat_by_group,
              selected_clusters=selected_clusters))
}
