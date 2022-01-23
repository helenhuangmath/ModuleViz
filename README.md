# Heatmap Plot
Functions of heatmap plot 

3 different ways to plot heatmap 

1. pheatmap::pheatmap() 

pdf("heatmap.pdf")
p=pheatmap(scaled_mat, scale="none", kmeans_k = 5,show_rownames = F)
dev.off()


2. ComplexHeatmap::Heatmap()

pdf(paste0("ComplexHeatmap_",comp,"_fdr05_FC2.pdf"), width=5, height=12)
Heatmap(scaled_mat,show_row_names = F)
dev.off()


3. heatmaply

https://cran.r-project.org/web/packages/heatmaply/vignettes/heatmaply.html

heatmaply( as.matrix(dt[,c(1:12)]),scale="row", 
	scale_fill_gradient_fun = ggplot2::scale_fill_gradient2(low = "blue", high = "red",), # default is viridis
	file = paste0("heatmap_",comp,".pdf"), 
  height = 900, showticklabels = c(TRUE, FALSE) 
  )



4. ggheatmap

ggheatmap(
  dt,
  scale = "column",
  row_side_colors = dt[, c("cyl", "gear")]
)

