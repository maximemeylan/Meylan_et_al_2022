### Meylan et al 2021
### produce figure 1 and figure S1
### max.meylan@gmail.com

load("~/projects/visium/results/2021-04-12_seurat_processed.RData")
library(Seurat)
library(dplyr)
library(patchwork)
library(ggplot2)
library(gridExtra)
#ADD annotation 
for(slide in slide_list){
  #correct ident
  spatial_list[[slide]]$orig.ident <- slide
  #adding all annotation available
  annot <- list.files(paste0("/data/visium_ccRCC/annotations/",slide),full.names = T)
  for(x in annot){
    annot_table <- read.csv(x)
    rownames(annot_table) <- annot_table[,1]
    annot_table$Barcode <- NULL
    spatial_list[[slide]] <- AddMetaData(object= spatial_list[[slide]],
                                         metadata = annot_table,
                                         col.name = paste0(colnames(annot_table),"_annot"))
  }
}

#hide spatial axes for seurat plots 
hide_axis <- theme(axis.title.x=element_blank(),
                   axis.text.x=element_blank(),
                   axis.ticks.x=element_blank(),
                   axis.title.y=element_blank(),
                   axis.text.y=element_blank(),
                   axis.ticks.y=element_blank())

#figure 1 main 
b_1 <- spatial_list$b_1
Idents(b_1)  <-as.numeric(Idents(b_1)) 
#H&E
he_b_1 <- SpatialPlot(b_1,repel = F,label = F,image.alpha=1,alpha = c(0,0), pt.size.factor =0.000001) + geom_point(alpha=0)+
  NoLegend() +
  ggtitle("TLS positive")+  
  theme(text=element_text(size=14))+ 
  theme(text=element_text(face = "bold"))  
# MCP scores
plot_list_mcp <- lapply(cell_types[c(5,1,2,4,6,8,9,10)],function(x){
  plot_list_mcp <- SpatialPlot(b_1,
                               features = x,
                               image.alpha=0,
                               pt.size.factor = 1.8) +
    scale_fill_viridis_c(option="C") + 
    DarkTheme() +
    hide_axis +    
    theme(text=element_text(size=14))+ 
    theme(text=element_text(face = "bold"))+
    theme(legend.text=element_text(size=7))
  
})
#CA9
ca9 <- SpatialPlot(b_1,
            features = "CA9",
            image.alpha=0,
            pt.size.factor = 1.8) +
  scale_fill_viridis_c(option="C") + 
  DarkTheme() +
  hide_axis +    
  theme(text=element_text(size=14))+ 
  theme(text=element_text(face = "bold"))+
  theme(legend.text=element_text(size=7))
plot_list_mcp[[9]] <- ca9
plot_list_mcp[[10]] <- he_b_1
lay <- rbind(c(10,9,1,2,3,4),
             c(10,9,5,6,7,8))
png(paste0(res_folder,Sys.Date(),"_","FIGURE1_.png"),width = 450,height=250,res = 300,units = "mm",bg="transparent")
grid.arrange(grobs = plot_list_mcp, layout_matrix = lay)
dev.off()

#figure 1 main 
tls_neg_b_6 <- spatial_list$b_6
Idents(tls_neg_b_6)  <-as.numeric(Idents(tls_neg_b_6)) 
#H&E
he_b_1 <- SpatialPlot(tls_neg_b_6,repel = F,label = F,image.alpha=1,alpha = c(0,0), pt.size.factor =0.000001) +
  NoLegend() +
  ggtitle("TLS negative")+  
  theme(text=element_text(size=14))+ 
  theme(text=element_text(face = "bold"))
# MCP scores
plot_list_mcp <- lapply(cell_types[c(5,1,2,4,6,8,9,10)],function(x){
  plot_list_mcp <- SpatialPlot(tls_neg_b_6,
                               features = x,
                               image.alpha=0,
                               pt.size.factor = 1.8) +
    scale_fill_viridis_c(option="C") + 
    DarkTheme() +
    hide_axis +    
    theme(text=element_text(size=14))+ 
    theme(text=element_text(face = "bold"))+ 
    theme(legend.text=element_text(size=8))

})
#CA9
ca9 <- SpatialPlot(tls_neg_b_6,
                   features = "CA9",
                   image.alpha=0,
                   pt.size.factor = 1.8) +
  scale_fill_viridis_c(option="C") + 
  DarkTheme() +
  hide_axis +    
  theme(text=element_text(size=14))+ 
  theme(text=element_text(face = "bold"))+
  theme(legend.text=element_text(size=8))
plot_list_mcp[[9]] <- ca9
plot_list_mcp[[10]] <- he_b_1

lay <- rbind(c(10,9,1,2,3,4),
             c(10,9,5,6,7,8))
png(paste0(res_folder,Sys.Date(),"_","FIGURE1__neg.png"),width = 500,height=250,res = 300,units = "mm",bg="transparent")
grid.arrange(grobs = plot_list_mcp, layout_matrix = lay)
dev.off()

pdf(paste0(res_folder,Sys.Date(),"_","FIGURE1_supp.pdf"),width = 16,height=6,bg="transparent")
sapply(names(spatial_list),function(id){
  print(id)
  #figure 1 main 
  spatial_obj <- spatial_list[[id]]
  Idents(spatial_obj)  <-as.numeric(Idents(spatial_obj)) 
  #H&E
  he <- SpatialPlot(spatial_obj,repel = F,label = F,image.alpha=1,alpha = c(0,0), pt.size.factor =0.00000) +
    NoLegend() +
    ggtitle(id)+  
    theme(text=element_text(size=10))+ 
    theme(text=element_text(face = "bold"))
  #MCP scores
  mcp_keep <- intersect(cell_types[c(5,1,2,4,6,8,9,10)],colnames(spatial_obj@meta.data))
  plot_list_mcp <- lapply(mcp_keep,function(x){
    plot_list_mcp <- SpatialPlot(spatial_obj,
                                 features = x,
                                 image.alpha=0,
                                 pt.size.factor = 1.8) +
      scale_fill_viridis_c(option="C") + 
      DarkTheme() +
      hide_axis +    
      theme(text=element_text(size=10))+ 
      theme(text=element_text(face = "bold"))+ 
      theme(legend.text=element_text(size=8))
  })
  #CA9
  if("CA9" %in% rownames(spatial_obj)){
    ca9 <- SpatialPlot(spatial_obj,
                       features = "CA9",
                       image.alpha=0,
                       pt.size.factor = 1.8) +
      scale_fill_viridis_c(option="C") + 
      DarkTheme() +
      hide_axis +    
      theme(text=element_text(size=10))+ 
      theme(text=element_text(face = "bold"))+
      theme(legend.text=element_text(size=8))
  }else{
    ca9 <- SpatialPlot(spatial_obj,repel = F,label = F,image.alpha=0,alpha = c(0,0), pt.size.factor =0.00000) + NoLegend() + ggtitle("CA9 not expressed")
  }
  plot_list_mcp[[9]] <- ca9
  plot_list_mcp[[10]] <- he
  
  test_null <- sapply(plot_list_mcp,is.null)
  if(any(test_null)){
    check_null <- which(test_null)
    plot_list_mcp[[check_null]] <- SpatialPlot(spatial_obj,repel = F,label = F,image.alpha=0,alpha = c(0,0), pt.size.factor =0.00000) + NoLegend() 
  }
  lay <- rbind(c(10,9,1,2,3,4),
               c(10,9,5,6,7,8))
  grid.arrange(grobs = plot_list_mcp, layout_matrix = lay)
})
dev.off()

