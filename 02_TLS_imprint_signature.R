### Meylan et al 2021
### produce TLS imprint signature, figure 2, figure S3 and figure S4
### max.meylan@gmail.com

### TLS signature visium
load("~/projects/visium/results/2021-04-12_seurat_processed.RData")
source("~/script/maxime.utils.r")
library(Seurat)
library(dplyr)
library(patchwork)
library(ggplot2)
library(gridExtra)
library(stringi)
library(ROCit)
library(png)
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
# Differential expression in each samples
TLS_pos_ids <- c("b_1","a_15","a_4","b_17")
TLS_list_all_samples <- lapply(TLS_pos_ids,function(y){
  x <- spatial_list[[y]]
  to_use <- "TLS"
  if(length(to_use)!=0){
    tls_MAST <- FindMarkers(x,group.by = "TLS_2_cat_annot",ident.1 = to_use,test.use = "MAST")
  }else{
    return(NA)
  }
})

#perform discovery on 3 tumors and validate on 1
hide_ig <- F
cross_validation <- sapply(TLS_pos_ids[c(3,2,1,4)],function(x){
  print(x)
  sub <- TLS_pos_ids[TLS_pos_ids!=x]
  gene_list <- lapply(TLS_list_all_samples[TLS_pos_ids!=x],function(y){
    rownames(y)[which(y$avg_log2FC > 1 & y$p_val_adj < 0.05)]
  })
  tls_sign <- unique(unlist(gene_list))
  intersect_tls_sign <- intersect(tls_sign,rownames(spatial_list[[x]]))
  if(hide_ig){
    intersect_tls_sign <- grep(x =intersect_tls_sign, pattern = "IG",value = T,invert = T)
    intersect_tls_sign <- grep(x =intersect_tls_sign, pattern = "JCHAIN",value = T,invert = T)
  }
  #ROC
  rocs <- sapply(intersect_tls_sign,function(gene){
    roc_empirical <- rocit(score = spatial_list[[x]][["SCT"]]@data[gene,],
                           class =  spatial_list[[x]]$TLS_2_cat_annot=="TLS",
                           negref = F)
    roc_empirical$AUC
  })
  return(list(
    "tls_sign"=intersect_tls_sign,
    "auc"=sort(rocs)))
})

#get list with threshold
tls_sign <- sort(table(unlist(sapply(cross_validation["auc",],function(auc) names(which(auc > 0.7))))))
tls_sign <- names(which(tls_sign >=3))

#saveRDS(tls_sign,"~/projects/visium/results/tls_signature_29.rds")
tls_sign <- readRDS("~/projects/visium/results/tls_signature_29.rds")

# venn diagramm 
tls_sign_over_07 <- sapply(cross_validation["auc",],function(auc) auc[which(auc > 0.7)])
library(VennDiagram)
toto <- sapply(tls_sign_over_07,function(x) names(x))
pl <- venn.diagram(
  x = toto,
  category.names = names(toto),
  filename = NULL,cex=4,
  fill=c(viridis(4)),
  output=T
)
png(paste0(res_folder,Sys.Date(),"_FIGURE_TLS_signature_supp_venn.png"),width = 7,height=5,units = "in",res = 500,bg = 'transparent')
grid.draw(pl)
dev.off()
#visualise on all tumours 
opt<-"C"
val <-lapply(names(spatial_list)[1:12],function(x){
  print(x)
  intersect_tls_sign <- intersect(tls_sign,rownames(spatial_list[[x]])) 
  spatial_list[[x]] <- AddMetaData(spatial_list[[x]],apply(as.matrix(spatial_list[[x]][["SCT"]]@data[intersect_tls_sign,]),2,mean),col.name = "TLS_signature") 
  #ROC 
  if(any(spatial_list[[x]]$TLS_2_cat_annot=="TLS")){
    roc_empirical <- rocit(score = spatial_list[[x]]$TLS_signature,
                           class =  spatial_list[[x]]$TLS_2_cat_annot=="TLS",
                           negref = F)
    auc <- roc_empirical$AUC
  }else{
    auc <- NA
  }
  SpatialPlot(spatial_list[[x]],image.alpha = 0, features = "TLS_signature") +
    scale_fill_viridis_c(option=opt,limits = c(0,4.1)) + 
    DarkTheme() +
    hide_axis +
    ggtitle(label =paste0("AUC = ", round(auc,2)))
})
grid.arrange(grobs=val,ncol=6)

png(paste0(res_folder,Sys.Date(),"_FIGURE_TLS_signature_supp_heatmap.png"),width = 24,height=10,units = "in",res = 500,bg = 'transparent')
grid.arrange(grobs=val,ncol=6)
dev.off()
#figure supp TLS signature 
#plot volcano
names(TLS_list_all_samples) <- TLS_pos_ids 
all_volcanoes <- lapply(TLS_pos_ids[c(3,2,1,4)], function(ids){
  diff_exp <- TLS_list_all_samples[[ids]]
  diff_exp <- diff_exp[,c("avg_log2FC","p_val_adj")]
  colnames(diff_exp) <- c("log2FC","pval")
  plot_volcano(dataset = diff_exp,
               pval_FC=c(0.05,2),
               #label_to_id = tls_sign,
               hide_label = F,
               ylabel = "-log10 (Adjusted P value)") + xlim(lims=c(-5,5)) +
    ggtitle(label=ids) +
    theme(text=element_text(size=15))+ 
    theme(text=element_text(face = "bold"))
  })
png(paste0(res_folder,Sys.Date(),"_FIGURE_TLS_signature_supp.png"),width = 10,height=10,units = "in",res = 500)
grid.arrange(grobs=all_volcanoes,n=2)
dev.off()

###  BOXPLOT 
library(stringr)
toto <- data.frame(unlist(cross_validation["auc",]))
toto$id <- str_split(rownames(toto),"\\.",simplify = T)[,1]
colnames(toto)[1] <- "values"
png(paste0(res_folder,Sys.Date(),"_FIGURE_TLS_BOXPLOT_supp.png"),width = 5,height=3,units = "in",res = 500)
ggplot(data = toto, aes(x=id, y=values,fill=id)) +
  geom_boxplot(alpha=0.7)+
  geom_point(color = "black", size = 1) +
  scale_fill_manual(values =viridis(4))+
  ylab("AUC")+
  xlab("TLS+ tumours")+ 
  theme_linedraw()+
  geom_hline(yintercept = 0.7,colour="red") +
  theme(text=element_text(size=15))+ 
  theme(text=element_text(face = "bold"))
dev.off()

#main IMM15P7118
spat <- "b_1"
spatial_list[[spat]] <- AddMetaData(spatial_list[[spat]],apply(as.matrix(spatial_list[[spat]][["SCT"]]@data[c("IGKC","JCHAIN"),]),2,mean),col.name = "Pan_IG") 
intersect_tls_sign <- intersect(tls_sign,rownames(spatial_list[[spat]])) 
spatial_list[[spat]] <- AddMetaData(spatial_list[[spat]],apply(as.matrix(spatial_list[[spat]][["SCT"]]@data[intersect_tls_sign,]),2,mean),col.name = "TLS_Imprint") 
#ROC 
if(any(spatial_list[[spat]]$TLS_2_cat_annot=="TLS")){
  roc_empirical <- rocit(score = spatial_list[[spat]]$TLS_Imprint,
                         class =  spatial_list[[spat]]$TLS_2_cat_annot=="TLS",
                         negref = F)
  auc <- roc_empirical$AUC
}else{
  auc <- NA
}

TLS_Imprint <- as.vector(scale(spatial_list[[spat]]$TLS_Imprint))
Pan_IG <- as.vector(scale(spatial_list[[spat]]$Pan_IG))

spatial_list[[spat]]@assays$SCT@scale.data <- rbind(spatial_list[[spat]]@assays$SCT@scale.data,TLS_Imprint)
spatial_list[[spat]]@assays$SCT@scale.data <- rbind(spatial_list[[spat]]@assays$SCT@scale.data,Pan_IG)
sign_list <- c("TLS_Imprint",
               "MZB1",
               "Pan_IG",
               "IGHG1",
               "IGHA1",
               "IGHM",
               "COL1A1",
               "CXCL12")
spat_cor_int <- intersect(sign_list,rownames(spatial_list[[spat]]@assays$SCT@scale.data))
spat_cor_ext <- setdiff(sign_list,rownames(spatial_list[[spat]]@assays$SCT@scale.data))

all_plots <-lapply(sign_list,function(sign){
  p <- SpatialPlot(spatial_list[[spat]],image.alpha = 0, features = sign)
  p <- p+ scale_fill_viridis_c(option="C")
  p <- p+ DarkTheme()
  if(sign=="TLS_Imprint"){
    p <- p+ annotate("text", 
                     x=max(GetTissueCoordinates(object = spatial_list[[spat]])[,2]) - 10, 
                     y=max(GetTissueCoordinates(object = spatial_list[[spat]])[,1]) -20,
                     size=4,
                     vjust=1,
                     hjust=1,
                     col="white",
                     label = paste0("AUC = ", round(auc,2)))
  }
  p <- p + hide_axis +theme(legend.title=element_text(size=15))
})
png(paste0(res_folder,Sys.Date(),x,"_FIGURE2_b_1.png"),width = 15,height=12,units = "in",res = 500)
grid.arrange(grobs=all_plots[c(3,6,4,5,7,8,1,2)],ncol=6)
dev.off()

#figure supp 
opt <- "C"
plot_list <- list()
for(spat in names(spatial_list)[c(1:2,4:12)]){
  print(spat)
  library("patchwork")
  pat <- list.files(path = "~/projects/visium/data/ihc",pattern =spat,full.names =T)
  ihc <- readPNG(pat,native=T)
  
  spatial_list[[spat]] <- AddMetaData(spatial_list[[spat]],apply(as.matrix(spatial_list[[spat]][["SCT"]]@data[c("IGKC","JCHAIN"),]),2,mean),col.name = "Pan_IG") 
  intersect_tls_sign <- intersect(tls_sign,rownames(spatial_list[[spat]])) 
  spatial_list[[spat]] <- AddMetaData(spatial_list[[spat]],apply(as.matrix(spatial_list[[spat]][["SCT"]]@data[intersect_tls_sign,]),2,mean),col.name = "TLS_Imprint") 
  #ROC 
  if(any(spatial_list[[spat]]$TLS_2_cat_annot=="TLS")){
    print("compute_roc")
    roc_empirical <- rocit(score = spatial_list[[spat]]$TLS_Imprint,
                           class =  spatial_list[[spat]]$TLS_2_cat_annot=="TLS",
                           negref = F)
    auc <- roc_empirical$AUC
  }else{
    auc <- NA
  }
  #compute spatial cor 
  TLS_Imprint <- as.vector(scale(spatial_list[[spat]]$TLS_Imprint))
  Pan_IG <- as.vector(scale(spatial_list[[spat]]$Pan_IG))
  
  spatial_list[[spat]]@assays$SCT@scale.data <- rbind(spatial_list[[spat]]@assays$SCT@scale.data,TLS_Imprint)
  spatial_list[[spat]]@assays$SCT@scale.data <- rbind(spatial_list[[spat]]@assays$SCT@scale.data,Pan_IG)
  sign_list <- c("TLS_Imprint",
                 "MZB1",
                 "Pan_IG",
                 "IGHM",
                 "IGHG1",
                 "IGHA1",
                 "COL1A1",
                 "CXCL12")
  spat_cor_int <- intersect(sign_list,rownames(spatial_list[[spat]]@assays$SCT@scale.data))
  spat_cor_ext <- setdiff(sign_list,rownames(spatial_list[[spat]]@assays$SCT@scale.data))
  
  # spatial_cor <- RunMoransI(data = spatial_list[[spat]]@assays$SCT@scale.data[spat_cor_int,],
  #                               pos =  GetTissueCoordinates(object = spatial_list[[spat]]))
  # spatial_cor[spat_cor_ext,] <-c(NA,NA)
  he <- SpatialPlot(spatial_list[[spat]],repel = F,label = F,image.alpha=1,alpha = c(0,0), pt.size.factor =0) + NoLegend()
  blank <- SpatialPlot(spatial_list[[spat]],repel = F,label = F,image.alpha=0,alpha = c(0,0), pt.size.factor =0) + NoLegend()+
    inset_element(p = ihc,
                  left = 0,
                  bottom = 0,
                  right = 1,
                  top = 1)
  
  all_plots <-lapply(sign_list,function(sign){
      p <- SpatialPlot(spatial_list[[spat]],image.alpha = 0, features = sign)
      p <- p+ scale_fill_viridis_c(option="C")
      p <- p+ DarkTheme()
      # p <- p+ annotate("text", 
      #                  x=max(GetTissueCoordinates(object = spatial_list[[spat]])[,2]) - 10, 
      #                  y = max(GetTissueCoordinates(object = spatial_list[[spat]])[,1]) +30,
      #                  size=4,
      #                  vjust=1,
      #                  hjust=1,
      #                  col="white"),
      #                  label = paste0("Spatial cor: ", round(spatial_cor[sign,"observed"],digits = 2)))
      if(sign=="TLS_Imprint" & spat %in% names(spatial_list)[c(1,2,4,6:9)]){
        p <- p+ annotate("text", 
                   x=min(GetTissueCoordinates(object = spatial_list[[spat]])[,2]) + 100, 
                   y=min(GetTissueCoordinates(object = spatial_list[[spat]])[,1]) +30,
                   size=4,
                   vjust=1,
                   hjust=1,
                   col="white",
                   label = paste0("AUC = ", round(auc,2)))
      }
    p <- p + hide_axis +theme(legend.title=element_text(size=15))
  })
  plot_list[[spat]]  <- (he | blank | all_plots[[1]] | all_plots[[2]] |all_plots[[3]] ) / ( all_plots[[4]] | all_plots[[5]] | all_plots[[6]] |all_plots[[7]]| all_plots[[8]])
}
pdf(paste0(res_folder,Sys.Date(),"_FIGURE_4_supp.pdf"),width = 15,height=15)
plot_list
dev.off()
