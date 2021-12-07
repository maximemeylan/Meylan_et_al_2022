### Meylan et al 2021
### processing of clonotypes called by MiXCR from visium data
### max.meylan@gmail.com

load("~/projects/visium/results/2021-04-12_seurat_processed.RData")
source("~/script/maxime.utils.r")

library(Seurat)
library(dplyr)
library(patchwork)
library(ggplot2)
library(gridExtra)
library(stringi)
library(viridis)
library(reshape2)
library(scales)
library(parallel)
library(stringr)
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
    if(colnames(annot_table)=="TLS_2_cat"){
      annot_table[annot_table[,"TLS_2_cat"]!="TLS","TLS_2_cat"]<-"NO_TLS"
    }
    spatial_list[[slide]] <- AddMetaData(object= spatial_list[[slide]],
                                         metadata = annot_table,
                                         col.name = paste0(colnames(annot_table),"_annot"))
  }
}


#bulk clonotype analysis
clonotype_list <- NULL
to_use <- slide_list[1:12]
#import 10X spatial barcodes
barcodes <- read.csv(paste0("/data/visium_ccRCC/processed_data/",slide_list[1],"/outs/raw_feature_bc_matrix/barcodes.tsv.gz"),header = F)
barcodes <- gsub(pattern = "-1",replacement = "",barcodes$V1)
for(x in to_use){
  print(x)
  dir_fastq <- paste0("/home/maxmey/projects/visium/mixcr/results/",x,"/")
  temp <- read.csv(paste0(dir_fastq,x,"repertoire.clonotypes.ALL.txt"),sep = '\t')
  temp$sampleid <- x
  temp$spatial_barcodes <- NA
  temp$nspots <- NA
  temp$type <- substr(temp$allVHitsWithScore, start = 1, stop = 3)
  clust <- makeCluster(12)
  clusterExport(clust, c("dir_fastq","barcodes","temp"))
  res <- parLapply(cl=clust,temp$cloneId,function(clone){
    #for(clone in temp$cloneId){
    fastq_clone_path <- paste(dir_fastq,"clones_fastq/",unique(temp$sampleid),"_reads_cln",clone,"_R1.fastq.gz",sep = "")
    if(file.exists(fastq_clone_path)){
      r1 <- read.csv(fastq_clone_path,sep = "\n")[c(TRUE, FALSE, FALSE, FALSE), ]
      res <-NULL
      for(barcode in barcodes){
        if(any(grepl( barcode, r1, fixed = TRUE))){
          res <- c(res,paste0(barcode,"-1"))
          temp[temp$cloneId==clone,"spatial_barcodes"] <- paste0(unlist(res),collapse = ',')
          temp[temp$cloneId==clone,"nspots"] <- length(unique(unlist(res)))
        }
      }
    }
    temp[temp$cloneId==clone,]
  })
  temp <- data.frame(do.call(rbind, res))
  clonotype_list[[x]] <- temp
  #annotate seurat object
  IGH_barcodes <- unique(unlist(strsplit(temp[temp$type=="IGH","spatial_barcodes"],split = ",")))
  spot_info <- sapply(IGH_barcodes, function(bc){
    bc_cdr3 <- temp[grep(bc, temp$spatial_barcodes),"aaSeqImputedCDR3"]
    counts <- length(unique(bc_cdr3))
    max_length  <- max(nchar(bc_cdr3),na.rm = T)
    average_length  <- mean(nchar(bc_cdr3),na.rm = T)
    temp_igh <- temp[temp$type=="IGH",]
    major_isotype <- unique(substr(temp_igh[grep(bc, temp_igh$spatial_barcodes),"allCHitsWithScore"],start = 1,stop =4))
    major_isotype_names <- paste0(names(table(major_isotype)),collapse = "")
    return(list("unique_counts_IGH"= counts,
                "average_cdr3_length_IGH"= average_length,
                "max_cdr3_length_IGH"= max_length,
                "isotypes"= major_isotype_names))
  })
  for(info in rownames(spot_info)){
    spatial_list[[x]] <- AddMetaData(object= spatial_list[[x]],
                                     metadata = unlist(spot_info[info,]),
                                     col.name = info)
  }
  #IGL 
  IGL_barcodes <- unique(unlist(strsplit(temp[temp$type=="IGL","spatial_barcodes"],split = ",")))
  spot_info_igl <- sapply(IGL_barcodes, function(bc){
    bc_cdr3 <- temp[grep(bc, temp$spatial_barcodes),"aaSeqImputedCDR3"]
    counts <- length(unique(bc_cdr3))
    max_length  <- max(nchar(bc_cdr3),na.rm = T)
    average_length  <- mean(nchar(bc_cdr3),na.rm = T)
    temp_igl <- temp[temp$type=="IGL",]
    major_isotype <- unique(substr(temp_igl[grep(bc, temp_igl$spatial_barcodes),"allCHitsWithScore"],start = 1,stop =4))
    major_isotype_names <- paste0(names(table(major_isotype)),collapse = "")
    return(list("unique_counts_IGL"= counts,
                "average_cdr3_length_IGL"= average_length,
                "max_cdr3_length_IGL"= max_length,
                "isotypes_IGL"= major_isotype_names))
  })
  for(info in rownames(spot_info_igl)){
    spatial_list[[x]] <- AddMetaData(object= spatial_list[[x]],
                                     metadata = unlist(spot_info_igl[info,]),
                                     col.name = info)
  }
}
clonotype_merged <- Reduce(rbind,clonotype_list)

## generate table 
unique_CDR3 <- sapply(unique(clonotype_merged$type),function(x) nrow(clonotype_merged[clonotype_merged$type==x,]))
counts_CDR3 <- sapply(unique(clonotype_merged$type),function(x) sum(clonotype_merged[clonotype_merged$type==x,"cloneCount"]))


#mutation analysis 
# targetFrom | targetTo | targetLength | queryFrom | queryTo | mutations | alignmentScore
clusterExport(clust, "stri_split")
clusterExport(clust, "str_count")
mutations <- parApply(cl = clust,clonotype_merged,1, function(x){
  #compute J alignments
  J_alignment <- x["allJAlignments"]
  #J_target_length <- as.numeric(stri_split(str = J_alignment,fixed='|')[[1]][3])
  #may be biaised by length of RNA fragment mapping
  J_target_length <- as.numeric(stri_split(str = J_alignment,fixed='|')[[1]][5]) - as.numeric(stri_split(str = J_alignment,fixed='|')[[1]][4])
  J_subs <- str_count(J_alignment,pattern = "S")
  J_ins <- str_count(J_alignment,pattern = "I")
  J_del <- str_count(J_alignment,pattern = "D")
  J_all <- (J_subs+J_ins+J_del) / J_target_length
  #compute V alignments
  V_alignment <- x["allVAlignments"]
  V_alignment <- gsub(";",",",V_alignment)
  #take first alignment
  V_alignment <- unlist(stri_split(str = V_alignment,fixed=','))
  #remove ""
  V_alignment <- V_alignment[!V_alignment == ""]
  V_qual <- sapply(stri_split(str = V_alignment,fixed="|"),function(splitst) as.numeric(splitst[length(splitst)]))
  
  vquerylengths <- sapply(V_alignment,function(align){
    align_split <- stri_split(str = align,fixed='|')[[1]]
    if(any(align_split!="")){
      V_query_length <- as.numeric(align_split[5]) - as.numeric(align_split[4])
    }else{
      0
    }
  })
  V_query_length <- sum(vquerylengths)
  
  #quality and length filter
  V_alignment <- V_alignment[V_qual > 0 ] 
  #quality to return 
  V_qual <- sapply(stri_split(str = V_alignment,fixed="|"),function(splitst) as.numeric(splitst[length(splitst)]))
  
  #if no aligment match quality
  V_target_length <- as.numeric(nchar(x["targetSequences"]))
  if(identical(character(0),V_alignment)){
    return(c('J_target_length' =J_target_length,
             "J_mutation_freq"=NA,
             "J_mutation_perc"=NA,
             'V_target_length' =V_target_length,
             "V_query_length"=NA,
             "V_mutation_freq"=NA,
             "V_qual"=NA,
             "V_mutation_perc"=NA))
  }else{
    V_subs <- str_count(x["allVAlignments"],pattern = "S")
    V_ins <- str_count(x["allVAlignments"],pattern = "I")
    V_del <- str_count(x["allVAlignments"],pattern = "D")
    V_all <- (V_subs+V_ins+V_del) / V_query_length
    return(c('J_target_length' =J_target_length,
             "J_mutation_freq"=(J_subs+J_ins+J_del),
             "J_mutation_perc"=J_all,
             'V_target_length' =V_target_length,
             "V_query_length"=V_query_length,
             "V_mutation_freq"=(V_subs+V_ins+V_del),
             "V_qual"=median(V_qual),
             "V_mutation_perc"=V_all))
  }
})
stopCluster(clust)
mutations <- t(mutations)
clonotype_merged_mut <- cbind(clonotype_merged,mutations)


# table per sample 
toto <- sapply(slide_list[1:12],function(id){
  print(id)
  ## generate table 
  unique_CDR3 <- sapply(c("IGH","IGK","IGL"),function(x){
    nrow(clonotype_merged[clonotype_merged$type==x & clonotype_merged$sampleid==id,])
  }) 
  counts_CDR3 <- sapply(c("IGH","IGK","IGL"),function(x){
    sum(clonotype_merged[clonotype_merged$type==x & clonotype_merged$sampleid==id ,"cloneCount"])
  }) 
  
  if(any(spatial_list[[id]]$TLS_2_cat_annot =="TLS")){
    tls_status <- "TLS positive"
  }else{
    tls_status <- "TLS negative"
  }
  return(c("clonotypes"=unique_CDR3,
           'counts'=counts_CDR3,
           "TLS_status"=tls_status))
})

#supplementary
get_order <- order(as.numeric(toto["counts.IGL",]),decreasing = T)
grid.table(t(toto[,get_order]))
write.csv(x = t(toto[,get_order]),file = "~/projects/visium/results/supplementary_table_igs_visium.csv")


for(clon in names(clonotype_list)){
  clonotype_list[[clon]] <- clonotype_merged_mut[clonotype_merged_mut$sampleid==clon,]
}
# add clonotype information to spatial data
for(sub in names(spatial_list)[1:10]){
  print(sub)
  temp <- clonotype_merged_mut[clonotype_merged_mut$sampleid==sub,]
  IGH_barcodes <- na.omit(unique(unlist(strsplit(temp[temp$type=="IGL","spatial_barcodes"],split = ","))))
  spot_info <- sapply(IGH_barcodes, function(bc){
    bc_V_perc <- temp[grep(bc, temp$spatial_barcodes),"V_mutation_perc"]
    bc_V_freq <- temp[grep(bc, temp$spatial_barcodes),"V_mutation_freq"]
    
    bc_J_perc <- temp[grep(bc, temp$spatial_barcodes),"J_mutation_perc"]
    bc_J_freq  <- temp[grep(bc, temp$spatial_barcodes),"J_mutation_freq"]
    
    average_V_perc  <- mean(bc_V_perc,na.rm = T)*100
    average_V_freq  <- mean(bc_V_freq,na.rm = T)
    
    median_V_perc  <- median(bc_V_perc,na.rm = T)*100
    median_V_freq  <- median(bc_V_freq,na.rm = T)
    
    
    cv_V_perc  <- sd(bc_V_perc,na.rm = T)/mean(bc_V_perc, na.rm=TRUE)*100
    cv_V_freq  <- sd(bc_V_freq,na.rm = T)/mean(bc_V_freq, na.rm=TRUE)*100
    
    average_J_perc  <- mean(bc_J_perc,na.rm = T)*100
    average_J_freq  <- mean(bc_J_freq,na.rm = T)
    
    max_V_perc <- max(bc_V_perc,na.rm = T)*100
    min_V_perc <- min(bc_V_perc,na.rm = T)*100
    
    max_J_perc <- max(bc_J_perc,na.rm = T)*100
    
    return(list("average_V_perc"= average_V_perc,
                "average_V_freq"= average_V_freq,
                "average_J_perc"= average_J_perc,
                "average_J_freq"= average_J_freq,
                "max_V_perc"=max_V_perc,
                "min_V_perc"=min_V_perc,
                "max_J_perc"=max_J_perc,
                "median_V_perc"=median_V_perc,
                "median_V_freq"=median_V_freq,
                "cv_V_perc"=cv_V_perc,
                "cv_V_freq"=cv_V_freq))
  })
  for(info in rownames(spot_info)){
    spatial_list[[sub]] <- AddMetaData(object= spatial_list[[sub]],
                                       metadata = unlist(spot_info[info,]),
                                       col.name = info)
  }
}

# get best aligments score
clonotype_merged_mut$best_V <- sapply(clonotype_merged_mut$allVHitsWithScore,function(x)strsplit(x,split = "\\*")[[1]][1])
clonotype_merged_mut$best_V <- factor(clonotype_merged_mut$best_V,levels=rev(names(sort(table(clonotype_merged_mut$best_V),decreasing = T))))
clonotype_merged_mut$best_D <- sapply(clonotype_merged_mut$allDHitsWithScore,function(x)strsplit(x,split = "\\*")[[1]][1])
clonotype_merged_mut$best_J <- sapply(clonotype_merged_mut$allJHitsWithScore,function(x)strsplit(x,split = "\\*")[[1]][1])
clonotype_merged_mut$best_C <- sapply(clonotype_merged_mut$allCHitsWithScore,function(x)strsplit(x,split = "\\*")[[1]][1])
clonotype_merged_mut$gene_used <- paste0(clonotype_merged_mut[,"best_V"],"_",clonotype_merged_mut[,"best_D"],"_",clonotype_merged_mut[,"best_J"])
clonotype_merged_mut$gene_used <- gsub(pattern = "_NA",replacement = "",clonotype_merged_mut$gene_used )
clonotype_merged_mut$full_length_id_nseq<-  apply(clonotype_merged_mut[,c("nSeqImputedFR1","nSeqImputedCDR1",
                                                                          "nSeqImputedFR2","nSeqImputedCDR2",
                                                                          "nSeqImputedFR3","nSeqImputedCDR3",
                                                                          "nSeqImputedFR4")],1,function(x) toupper(paste0(x,collapse = "")))

clonotype_merged_mut$putative_nseq_CDR3_FR4<-  apply(clonotype_merged_mut[,c("nSeqImputedCDR3",
                                                                    "nSeqImputedFR4")],1,function(x) paste0(x,collapse = ""))

clonotype_merged_mut$putative_nseq_CDR3_FR4 <- sapply(clonotype_merged_mut$putative_nseq_CDR3_FR4,function(x){
  paste0(unlist(str_extract_all(x, "([A-Z]+(?=[^a-z]))")),collapse = "")
} )

clonotype_merged_mut$full_length_id_aaSeq<-  apply(clonotype_merged_mut[,c("aaSeqImputedFR1","aaSeqImputedCDR1",
                                                                          "aaSeqImputedFR2","aaSeqImputedCDR2",
                                                                          "aaSeqImputedFR3","aaSeqImputedCDR3",
                                                                          "aaSeqImputedFR4")],1,function(x) toupper(paste0(x,collapse = "")))
                
spot_info_sequences <-NULL     

for(sub in names(spatial_list)[1:10]){
  print(sub)
  temp <- clonotype_merged_mut[clonotype_merged_mut$sampleid==sub,]
  temp_IGH <- clonotype_merged_mut[clonotype_merged_mut$sampleid==sub & clonotype_merged_mut$type=="IGH",]
  temp_IGL <- clonotype_merged_mut[clonotype_merged_mut$sampleid==sub & clonotype_merged_mut$type=="IGL",]
  
  IGH_barcodes <- na.omit(unique(unlist(strsplit(temp[temp$type=="IGH","spatial_barcodes"],split = ","))))
  IGL_barcodes <- na.omit(unique(unlist(strsplit(temp[temp$type=="IGL","spatial_barcodes"],split = ","))))
  spot_info_IGH <- lapply(intersect(IGL_barcodes,IGH_barcodes), function(bc){
    seq_IGH <- temp_IGH[grep(bc, temp_IGH$spatial_barcodes),"full_length_id_nseq"]
    seq_IGL <- temp_IGL[grep(bc, temp_IGL$spatial_barcodes),"full_length_id_nseq"]
    res <- list("seq_IGH"=seq_IGH,
             "seq_IGL"=seq_IGL)
  })
}

res <- sapply(spot_info_IGH,function(x){
  if(temp_IGH$full_length_id_nseq[5] %in% x$seq_IGH){
    x$seq_IGL
  }else{
    NA
  }
})

#mutation rates for each sample
var_to_check <- c("cloneCount","V_mutation_perc","V_mutation_freq","J_mutation_perc","J_mutation_freq")
melted <- melt_df(clonotype_merged_mut[,c("sampleid",var_to_check)],var_to_group = "sampleid")
plots <- sapply(var_to_check,function(x){
  plot_group_boxplot(data.m = melted,
                     variable = x,
                     log_scale = F,
                     violin=T,
                     compare_groups = F,
                     specify_col = viridis(12),
                     plot_outlier =T,
                     title_size = 8,
                     alpha=0.001,
                     add_jitter=T,
                     labs=c("","",x))
})


#get all mutation rates measured in each zone (may be multiple values per spot)
plot_list <- list()
for(sub in slide_list[1:12]){
  print(sub)
  zones <- rev(unique(spatial_list[[sub]]$TLS_2_cat_annot))
  mutation_per_zone <- sapply(zones,function(zone){
    temp <- clonotype_merged_mut[clonotype_merged_mut$sampleid==sub & clonotype_merged_mut$type=="IGL",]
    zone_bc <- names(which(spatial_list[[sub]]$TLS_2_cat_annot==zone))
    mutation_zone <- sapply(zone_bc,function(x){
      temp[grepl(pattern = x,x =  temp$spatial_barcodes),"V_mutation_freq"]
    })
    unlist(mutation_zone)
  })
  names(mutation_per_zone) <- zones
  n <- max(length(mutation_per_zone[["NO_TLS"]]), length(mutation_per_zone[["TLS"]]))
  length(mutation_per_zone[["TLS"]]) <- n                      
  length(mutation_per_zone[["NO_TLS"]]) <- n
  
  perc_TLS_df <- cbind(mutation_per_zone[["NO_TLS"]],mutation_per_zone[["TLS"]])
  perc_TLS_df <- t(perc_TLS_df)
  perc_TLS_df <- data.frame(perc_TLS_df,check.names = F)
  perc_TLS_df$zone <- c("Tumor area","TLS")
  merged_melted <- na.omit(melt_df(perc_TLS_df,var_to_group = "zone"))
  merged_melted$variable <-"zone"
  merged_melted$value <- merged_melted$value+1
  merged_melted$groups <- factor(merged_melted$groups,levels=rev(c("Tumor area","TLS")))
  p <- plot_group_boxplot(data.m = merged_melted,
                          variable = "zone",
                          log_scale =T,
                          violin=T,
                          compare_groups = T,
                          specify_col = rev(viridis(6,option = "B")[c(3,5)]),
                          plot_outlier = T,
                          trim = T,
                          title_size = 8,
                          alpha=0.0000005,
                          add_jitter=T,
                          labs=c(sub,"","Mutation count + 1"))
  plot_list[[sub]]<- p$plot
}

#FIGURE 4 alternative VIRIDIS 
p1 <- SpatialPlot(spatial_list[[x]], image.alpha = 0, features = c("unique_counts_IGL")) + hide_axis + scale_fill_viridis_c(option=opt) + DarkTheme() 
p2 <- SpatialPlot(spatial_list[[x]], image.alpha = 0,features =c("median_V_freq"))+ hide_axis + scale_fill_viridis_c(option=opt) + DarkTheme()
pimm <- (p1 | p2 )

p4 <- SpatialPlot(spatial_list[[y]], image.alpha = 0, features = c("unique_counts_IGL")) + hide_axis + scale_fill_viridis_c(option=opt) + DarkTheme() + hide_axis
p6 <- SpatialPlot(spatial_list[[y]], image.alpha = 0,features =c("median_V_freq")) + hide_axis + scale_fill_viridis_c(option=opt) + DarkTheme() + hide_axis
pexhau <- (p4 | p6)

p7 <- SpatialPlot(spatial_list[[z]], image.alpha = 0, features = c("unique_counts_IGL"))+ hide_axis + scale_fill_viridis_c(option=opt) + DarkTheme() + hide_axis
p9 <- SpatialPlot(spatial_list[[z]], image.alpha = 0,features =c("median_V_freq"))+ hide_axis + scale_fill_viridis_c(option=opt) + DarkTheme() + hide_axis
pbionikk <- (p7 | p9 )

p <- (p1 | p2 | p3 ) / (p4 | p5 | p6) / (p7 | p8 | p9 )
boxplots <- grid.arrange(grobs=plot_list[c("b_1",
                                          "a_3",
                                          "c_57")],ncol=1) 

ggsave(paste0(res_folder,Sys.Date(),"_","FIGURE4_plots.png"),
       boxplots,
       width = 2.4,height=7,dpi = 300,scale=1)


#clustering with levenstein distance
library(stringdist)
library(dendextend)
library(pheatmap)
sub <-"b_1"
clon_sub <- clonotype_merged_mut[clonotype_merged_mut$sampleid==sub & clonotype_merged_mut$type=="IGH",]
rownames(clon_sub) <- paste0("X",clon_sub$cloneId)
all_dist_mat <- stringsimmatrix(clon_sub[,"nSeqImputedCDR3"],method = "lv")
all_dist <- stringdistmatrix(clon_sub[,"nSeqImputedCDR3"],method="lv")

colnames(all_dist_mat) <- paste0("IGH-",colnames(all_dist_mat))
rownames(all_dist_mat) <- paste0("IGH-",rownames(all_dist_mat))

distance_threshold <- 0.15
#define clonality by identity of V and J + similar length CDR3 
all_hits <- apply(clon_sub,1,function(ig1){
  V_query <- strsplit(x = ig1["allVHitsWithScore"],split = ",")
  V_query <- strsplit(x = V_query[[1]],split = "-")
  V_query <- unique(sapply(V_query, "[[", 1))

  
  J_query <- strsplit(x = ig1["allJHitsWithScore"],split = ",")
  J_query <- strsplit(x = J_query[[1]],split = "\\*")
  J_query <- unique(sapply(J_query, "[[", 1))

  hist_ig1 <- apply(clon_sub,1,function(ig2){
    if(ig1["cloneId"]==ig2["cloneId"]){
      return(NULL)
    }
    V_query_2 <- strsplit(x = ig2["allVHitsWithScore"],split = ",")
    V_query_2 <- strsplit(x = V_query_2[[1]],split = "-")
    V_query_2 <- unique(sapply(V_query_2, "[[", 1))

    
    J_query_2 <- strsplit(x = ig2["allJHitsWithScore"],split = ",")
    J_query_2 <- strsplit(x = J_query_2[[1]],split = "\\*")
    J_query_2 <- unique(sapply(J_query_2, "[[", 1))

    
    V_intersect <- intersect(V_query,V_query_2)
    J_intersect <- intersect(J_query,J_query_2)
    
    #test for same VJ and same length CDR3s
    identity_test <- all(!identical(V_intersect,character(0)) & !identical(J_intersect,character(0))) & nchar(ig1["nSeqImputedCDR3"]) == nchar(ig2["nSeqImputedCDR3"])
    #if test is true, measure CDR3 distances
    if(identity_test){
      cdr3_dist <- stringdist(a=ig1["nSeqImputedCDR3"],
                              b=ig2["nSeqImputedCDR3"],
                              method = "lv")
      if(cdr3_dist <= nchar(ig1["nSeqImputedCDR3"])*distance_threshold){
        return(list("V"=V_intersect,"J"=J_intersect,"cdr3_dist"=cdr3_dist))
      }
    }
  })
  
})
all_hits <- Filter(Negate(is.null), lapply(all_hits, function(x) Filter(Negate(is.null), x)))
all_hits <- unique(sapply(names(all_hits),function(x) sort(c(x,names(all_hits[[x]])))))
all_hits <- all_hits[order(sapply(all_hits,length),decreasing = T)]

#remove completely overlapping families
to_keep <- sapply(all_hits,function(x){
 res <- sapply(all_hits,function(y){
   all(x %in% y)
 })
 sum(res) >= 2 
})
all_hits <- all_hits[!to_keep]

#merge clone 1 and 2 with 14/15 overlap 
all_hits[[1]]  <- c(all_hits[[1]],setdiff(all_hits[[2]],all_hits[[1]]))
all_hits[[2]] <- NULL
#add clonal family to the df  
#order matters !
for(i in length(all_hits):1){
  clon_sub[all_hits[[i]],"clonal_family"] <- i
}
# create the ighclust variables for the dendogram
igh_clust <- setNames(clon_sub$clonal_family,1:length(clon_sub$clonal_family))

#show histograms to evaluate best leuvensthein threhsold 
toto <- as.dendrogram(hclust(all_dist))

clust_to_highlight <- which(table(igh_clust) > 1)
ordered_ids <- order.dendrogram(toto)
names(igh_clust) <-clon_sub$aaSeqImputedCDR3

#add cluster info + make unique 
unique_ids <- paste0( names(igh_clust),"_",igh_clust)
baba <- ifelse(duplicated(unique_ids),
       paste0(unique_ids,"'"),unique_ids)
bubu <- ifelse(duplicated(baba), paste0(unique_ids,"''"),baba)
cucu <- ifelse(duplicated(bubu), paste0(unique_ids,"'''"),bubu)
unique_ids <- cucu        

#add to the df
clon_sub[,"clonal_ids"] <- unique_ids
clon_sub[,"clonal_relationship"] <- igh_clust

#scale for label cex
cex_value <-clon_sub[,"cloneCount"][ordered_ids]
cex_size <- rescale(log(cex_value), to=c(0.2,2))

#chose colors for leave spatialization 
map2color <- function(x,n=20){
  viridis(n,option="C")[cut(x,n)]
}
cex_col_value <- clon_sub[,"nspots"][ordered_ids]
cex_col <- map2color(log(cex_col_value))
#cex_col <- map2color(cex_col_value)
cex_col[is.na(cex_col)] <- "#808080"
cex_col <- cex_col

d1 <- color_branches(as.dendrogram(toto)
                     ,clusters = igh_clust[ordered_ids],
                     groupLabels = F)
d1 <- dendextend::set(d1,"labels_cex", 0.5)
d1 <- dendextend::set(d1,"leaves_pch", 19)
d1 <- dendextend::set(d1,"leaves_cex", cex_size)
d1 <- dendextend::set(d1,"leaves_col", cex_col)

d1 <- dendextend::set(d1,"labels", unique_ids[ordered_ids])
d1 <- dendextend::set(d1,"branches_lwd",1.7)
label_cols <-sapply(get_leaves_edgePar(d1),function(x){
  ifelse(is.null(x$col),'#000000',x$col)
})
d1 <- dendextend::set(d1,"labels_col", label_cols)
par(mar = c(5,5,5,5),xpd = NA)

circlize_dendrogram(d1,
                    facing = "outside",
                    labels = T,
                    labels_track_height = 0.05,
                    dend_track_height = 0.89)
# histogramms 
fafa <- clon_sub[,]
neq_var <- c("nSeqImputedCDR3")
to_id <- 1
thresh <- nchar(fafa[which(igh_clust==to_id),neq_var])[1]*0.15
hist(as.matrix(all_dist)[which(igh_clust==to_id),],
     breaks=100,
     main = paste0("Clonal Igs: ",to_id),
     xlab = "levenstheim distance")
abline(v=thresh,col="red")
text(paste0("threshold = ",thresh),x = thresh + 2,y=15,col="red")


pdf(paste0(res_folder,Sys.Date(),"clonal_relationship_2.O.pdf"),width =10,height = 10)
par(mar = c(8,8,8,8),xpd = NA)
circlize_dendrogram(d1,
                    facing = "outside",
                    labels = T,
                    labels_track_height = 0.05,
                    dend_track_height = 0.89)
plot.new()
legend(      
  legend = c(80,25,10,1),
  title = "Absolute count",
  x = 0,
  box.col="white",
  y=0.5,
  cex =1,
  pch = 19,
  pt.cex = cex_size[order(cex_value,decreasing = T)][c(1,3,9,115)]) 
legend(
       legend = c(20,10,5,1,"NA"),
       title = "Unique spatial spot count",
       x = 0,
       box.col="white",
       y=0.3,
       pch = 19,
       pt.cex = 2,
       cex = 1,
       col = c(map2color(log(c(20,10,5,1))),"#808080")
)

dev.off()

# Create spatial plot identifying clonal families
df_clone_clust <- clon_sub[,c("spatial_barcodes",
                              "clonal_ids",
                              "clonal_family",
                              "V_mutation_perc")]

plots <- lapply(as.character(c(1:10)),function(clone_num){
  print(clone_num)
  to_select <- df_clone_clust[df_clone_clust$clonal_family==clone_num,"spatial_barcodes"]
  if(!identical(to_select,character(0))){
    to_select <- strsplit(x = to_select,split = ",")
    names(to_select) <-df_clone_clust[df_clone_clust$clonal_family==clone_num,"clonal_ids"]
    
    to_select <- sapply(to_select,function(x) colnames(spatial_list[[sub]])[match(x= x,table = colnames(spatial_list[[sub]]))])
    to_select <- to_select[names(which(!is.na(to_select)))]  
    
    #count number of matching spatial barcode 
    match_bc <- sum(sapply(to_select,function(x) any(x %in% colnames(spatial_list[[sub]]))))
    if(match_bc >1){
      #scatterpie df 
      library(scatterpie)
      get_used_bc <- na.omit(unique(unlist(to_select)))
      scatter_df <- GetTissueCoordinates(object = spatial_list[[sub]])[,]
      scatter_df[,names(to_select)] <- 0
      for(cln_id in 1:length(to_select)){
        clone <- to_select[cln_id]
        clean_bc <- as.character(na.omit(unlist(unique(clone))))
        scatter_df[clean_bc,names(clone)] <- scatter_df[clean_bc,names(clone)] + 1
      }
      scatter_df$region <- factor(1:nrow(scatter_df))
      #transform 
      y <- scatter_df$imagerow
      y_mid <- min(y)+0.5*(max(y)-min(y))
      coords_y <- y-y_mid # centre sur le milieu
      coords_y<- -coords_y # on inverse tout
      coords_y <- coords_y+y_mid #on repasse sur une echelle 
      scatter_df$imagerow <-coords_y 
      if(clone_num=="8"){
        plus_smth <- 1
      }else{
        plus_smth <- 0
      }
      p <- SpatialPlot(spatial_list[[sub]],
                       image.alpha = 0.7,
                       pt.size.factor =0,
                       cells.highlight = to_select,
                       label.size = 3,
                       cols.highlight = c(rainbow(match_bc-plus_smth),transparent_col))
      p <- p  + geom_scatterpie(data=scatter_df, aes(y=imagerow, x=imagecol,group=region),
                                pie_scale = 0.66, 
                                cols=colnames(scatter_df[3:(ncol(scatter_df)-1)])) 
      if(clone_num==0){
        p <- SpatialDimPlot(spatial_list[[sub]],
                            pt.size.factor = ifelse(points,4,0),
                            cells.highlight = points,
                            cols.highlight = c("grey","black"),
                            stroke = ifelse(points,2,0)
        )
        p <- p  + geom_scatterpie(data=scatter_df, aes(y=imagerow, x=imagecol,group=region),
                                  pie_scale = 0.66, 
                                  cols=colnames(scatter_df[3:(ncol(scatter_df)-1)])) 
      }
      p <- p + ggtitle(paste0("Clonal immunoglobulin family: ",clone_num))
    }
  }
  
})
plots <- plots[!sapply(plots,is.null)]
p1 <- grid.arrange(grobs=plots[1:5] ,ncol=3)
par(mar = c(8,8,8,8),xpd = NA)
ggsave(paste0(res_folder,Sys.Date(),"_","FIGURE4_clonal_igs_position.png"),
       p1,
       width = 24,height=11.67,dpi = 300)


                                      
