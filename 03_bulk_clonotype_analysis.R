### Meylan et al 2021
### processing of clonotypes called by MiXCR from bulk RNAseq data
### max.meylan@gmail.com
library(Seurat)
library(dplyr)
library(patchwork)
library(ggplot2)
library(gridExtra)
library(stringi)
library(viridis)
library(pastecs)
library(reshape2)
library(RColorBrewer)
source("~/script/maxime.utils.r")
res_folder <- "~/projects/visium/results/"
#bulk clonotype analysis
clonotype_list <- NULL
to_use <- c(readLines("/home/maxmey/projects/visium/mixcr/exhauIMM_id_list.txt"),
            readLines("/home/maxmey/projects/visium/mixcr/exhauHEGP_id_list.txt"),
            readLines("/home/maxmey/projects/visium/mixcr/visium_complete_id_list.txt"),
            readLines("/home/maxmey/projects/visium/mixcr/visium_complete_id_list_2.txt")
)
file_exhau <- list.files(path =c("/home/maxmey/projects/visium/mixcr/results/bulk_exhau/",
                                 "/home/maxmey/projects/visium/mixcr/results/bulk_visium_complete/",
                                 "/home/maxmey/projects/visium/mixcr/results/bulk_visium_complete_2/"),
                         pattern = "repertoire.clonotypes.ALL.txt",
                         recursive = T,
                         full.names = T)
for(x in file_exhau){
  temp <- read.csv(x,sep = '\t')
  temp$type <- substr(temp$allVHitsWithScore, start = 1, stop = 3)
  temp$isotype <- substr(temp[,"allCHitsWithScore"],start = 1,stop =4)
  
  #edit names 
  new_name <- strsplit(x,split = "/")[[1]][length(strsplit(x,split = "/")[[1]])]
  new_name <- gsub(pattern = "repertoire.clonotypes.ALL.txt",replacement = "",new_name)
  new_name <- gsub(pattern = "ExCRF_",replacement = "",new_name)
  new_name <- gsub(pattern = "exhauCRF_",replacement = "",new_name)
  new_name <- gsub(pattern = "exhauCRF-",replacement = "",new_name)
  new_name <- gsub(pattern = "IMM_",replacement = "",new_name)
  new_name <- gsub(pattern = "_RNA",replacement = "",new_name)
  
  temp$sampleid <- new_name
  clonotype_list[[new_name]] <- temp
}
clonotype_merged <- Reduce(rbind,clonotype_list)

#filter out tonsil
clonotype_merged <- clonotype_merged[clonotype_merged$sampleid!="Tonsil",]

## generate table
unique_CDR3 <- sapply(unique(clonotype_merged$type),function(x) length(unique(clonotype_merged[clonotype_merged$type==x,"aaSeqImputedCDR3"])))
counts_CDR3 <- sapply(unique(clonotype_merged$type),function(x) sum(clonotype_merged[clonotype_merged$type==x,"cloneCount"]))
unique_CDR3_per_sample <- sapply(unique(clonotype_merged$sampleid),function(x) length(unique(clonotype_merged[clonotype_merged$sampleid==x & clonotype_merged$type=="IGH","aaSeqImputedCDR3"])))
total_IGH_per_sample <- sapply(unique(clonotype_merged$sampleid),function(x) sum(clonotype_merged[clonotype_merged$sampleid==x & clonotype_merged$type=="IGH","cloneCount"]))

#mutation analysis
# targetFrom | targetTo | targetLength | queryFrom | queryTo | mutations | alignmentScore
library(stringr)
library(parallel)
time1 <- Sys.time()
clust <- makeCluster(12)

#clusterExport(clust, "clonotype_merged")
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
  #
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
  
  #quality filter
  V_alignment <- V_alignment[V_qual > 750 ] 
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
time2 <- Sys.time()
time2 - time1
stopCluster(clust)
mutations <- t(mutations)
clonotype_bulk_m <- cbind(clonotype_merged,mutations)

#take only best aligments score
time1 <- Sys.time()
clust <- makeCluster(32)
clusterExport(clust, "stri_split")
clusterExport(clust, "str_count")
res <- parApply(cl = clust,clonotype_bulk_m,1, function(x){
  best_V <- strsplit(x["allVHitsWithScore"],split = "\\*")[[1]][1]
  best_V <- factor(best_V,levels=rev(names(sort(table(best_V),decreasing = T))))
  best_D <- strsplit(x["allDHitsWithScore"],split = "\\*")[[1]][1]
  best_J <- strsplit(x["allJHitsWithScore"],split = "\\*")[[1]][1]
  best_C <- strsplit(x["allCHitsWithScore"],split = "\\*")[[1]][1]
  gene_used <- paste0(best_V,"_",best_D,"_",best_J)
  gene_used <- gsub(pattern = "_NA",replacement = "",gene_used)
  
  full_length_id_nseq<-  toupper(paste0(x[c("nSeqImputedFR1","nSeqImputedCDR1",
                                            "nSeqImputedFR2","nSeqImputedCDR2",
                                            "nSeqImputedFR3","nSeqImputedCDR3",
                                            "nSeqImputedFR4")]
                                        ,collapse = ""))
  
  full_length_id_aaseq<- toupper(paste0(x[c("aaSeqImputedFR1","aaSeqImputedCDR1",
                                            "aaSeqImputedFR2","aaSeqImputedCDR2",
                                            "aaSeqImputedFR3","aaSeqImputedCDR3",
                                            "aaSeqImputedFR4")],
                                        collapse = ""))
  return(c("best_V"=best_V,
         "best_D"=best_D,
         "best_J"=best_J,
         "best_C"=best_C,
         "gene_used"=gene_used,
         "full_length_id_nseq"=full_length_id_nseq,
         "full_length_id_aaseq"=full_length_id_aaseq))
})
time2 <- Sys.time()
time2 - time1
stopCluster(clust)
res <- t(res)
clonotype_bulk_m <- cbind(clonotype_bulk_m,res)

var_to_check <- c("cloneCount","V_mutation_perc","V_mutation_freq","V_target_length")
rna_imm_hegp <- readRDS("~/projects/visium/data/batch_corrected_RNA_seq_TPM.rds")

TLS_signature <- readRDS("~/projects/visium/results/tls_signature_29.rds") 
TLS_signature <- intersect(TLS_signature,rownames(rna_imm_hegp))
#define threshold
all_tum_ids <- colnames(rna_imm_hegp)[colnames(rna_imm_hegp) != "Tonsil"]
TLS_scores <- apply(rna_imm_hegp[TLS_signature,all_tum_ids],2,mean)
thresh_tlspos <- quantile(TLS_scores,probs = 0.5)

TLS_positivity <- ifelse(TLS_scores > thresh_tlspos,"high","low")
clonotype_bulk_m$TLS_positivity <- factor(ifelse(colMeans(rna_imm_hegp[TLS_signature,clonotype_bulk_m$sampleid]) > thresh_tlspos,"high","low"))

# Average mutation rates for each patients then compute pvalues boxplot
var_to_check <- c("cloneCount","V_mutation_freq","V_mutation_perc")
res_stats <- sapply(unique(clonotype_bulk_m$sampleid),function(pid){
  temp <- clonotype_bulk_m[clonotype_bulk_m$sampleid==pid & clonotype_bulk_m$type=="IGL" ,var_to_check]
  res <- stat.desc(temp)
  res <- as.vector(res["mean",])
})
res_stats <- data.frame(apply(res_stats,1,function(x) unlist(x)),check.rows = F)
res_stats[,"TLS_positivity"] <- sapply(rownames(res_stats), function(pid) unique(clonotype_bulk_m[clonotype_bulk_m$sampleid==pid  ,"TLS_positivity"]))
res_stats[,"richness"] <- sapply(rownames(res_stats), function(pid) nrow(clonotype_bulk_m[clonotype_bulk_m$sampleid==pid & clonotype_bulk_m$type=="IGL" ,]))

colnames(res_stats) <- c("Average IGL Clone Count","Average VL gene mutation count","Average VL gene mutation rate","TLS_positivity","IGL Richness")
#mutation rates according to TLS positivity 
melted <- melt_df(res_stats,var_to_group = "TLS_positivity")
melted$value <- melted$value
melted <- na.omit(melted)
plots <-sapply(colnames(res_stats)[c(1:3,5)],function(x){
  plot_group_boxplot(data.m = melted,
                     variable = x,
                     log_scale = T,
                     violin=F,
                     trim = T,
                     compare_groups = F,
                     specify_col = brewer.pal(3,"Pastel2")[c(2,1)],
                     plot_outlier =T,
                     title_size = 8,
                    # ylim = c(0,0.13),
                     alpha=0.3,
                     add_jitter=T,
                     labs=c("","TLS signature",x))
})
p1 <- plots[1,colnames(res_stats)[c(1:3,5)]]
p1$`Average VL gene mutation count`
res_stats <- sapply(unique(clonotype_bulk_m$sampleid),function(pid){
  temp <- clonotype_bulk_m[clonotype_bulk_m$sampleid==pid & clonotype_bulk_m$type=="IGH" ,var_to_check]
  res <- stat.desc(temp)
  res <- as.vector(res["mean",])
})
res_stats <- data.frame(apply(res_stats,1,function(x) unlist(x)),check.rows = F)
res_stats[,"TLS_positivity"] <- sapply(rownames(res_stats), function(pid) unique(clonotype_bulk_m[clonotype_bulk_m$sampleid==pid  ,"TLS_positivity"]))
res_stats[,"richness"] <- sapply(rownames(res_stats), function(pid) nrow(clonotype_bulk_m[clonotype_bulk_m$sampleid==pid & clonotype_bulk_m$type=="IGH" ,]))
colnames(res_stats) <- c("Average IGH Clone Count","Average VH gene mutation count","Average VH gene mutation rate","TLS_positivity","Richness")
#mutation rates according to TLS positivity 
library(RColorBrewer)
melted <- melt_df(res_stats,var_to_group = "TLS_positivity")
melted$value <- melted$value
melted <- na.omit(melted)
plots <-sapply(colnames(res_stats)[c(1:3,5)],function(x){
  plot_group_boxplot(data.m = melted,
                     variable = x,
                     log_scale = T,
                     violin=F,
                     trim = T,
                     compare_groups = F,
                     specify_col = brewer.pal(3,"Pastel2")[c(2,1)],
                     plot_outlier =T,
                     title_size = 8,
                     #ylim = c(0,0.13),
                     alpha=0.3,
                     add_jitter=T,
                     labs=c("","TLS signature",x))
})
p2 <- plots[1,colnames(res_stats)[c(1:3,5)]]
p2$`Average VH gene mutation count`

count_richness_IGH <- p2$`Average IGH Clone Count` + p2$Richness
mutations_IGL_IGH <- p1$`Average VL gene mutation count` + p2$`Average VH gene mutation count`

ggsave(paste0(res_folder,Sys.Date(),"_","FIGURE5_clone_count.png"),
       plot=count_richness_IGH,  
       width = 130,height=100,units = "mm")

ggsave(paste0(res_folder,Sys.Date(),"_","FIGURE5_clone_mutation.png"),
       plot=mutations_IGL_IGH,
       width = 120,height=120,units = "mm")

library(immunarch)
rep_dir <-list.files(c("/home/maxmey/projects/visium/mixcr/results/bulk_exhau/",
                       "/home/maxmey/projects/visium/mixcr/results/bulk_visium_complete/",
                       "/home/maxmey/projects/visium/mixcr/results/bulk_visium_complete_2/"),
                     pattern = "IGH.txt",
                     full.names = T,
                     recursive = T)
#remove tonsil 
rep_dir <- grep(pattern = "Tonsil",x = rep_dir,invert = T,value = T)
immdata_mixcr <- immunarch::repLoad(rep_dir)

#saveRDS(immdata_mixcr,file = "/home/maxmey/projects/visium/mixcr/results/bulk_exhau/immunoarch_object_all.rds")
#immdata_mixcr<- readRDS("/home/maxmey/projects/visium/mixcr/results/bulk_exhau/immunoarch_object_all.rds")
names(immdata_mixcr$data) <- gsub(pattern = "-IGH","",names(immdata_mixcr$data))
names(immdata_mixcr$data) <- gsub(pattern = "repertoire.clonotypes.IGH","",names(immdata_mixcr$data))
names(immdata_mixcr$data) <- gsub(pattern = "ExCRF_","",names(immdata_mixcr$data))
names(immdata_mixcr$data) <- gsub(pattern = "exhauCRF-","",names(immdata_mixcr$data))
names(immdata_mixcr$data) <- gsub(pattern = "_RNA","",names(immdata_mixcr$data))
names(immdata_mixcr$data) <- gsub(pattern = "exhauCRF_","",names(immdata_mixcr$data))
names(immdata_mixcr$data) <- gsub(pattern = "IMM_","",names(immdata_mixcr$data))

immdata_mixcr$meta$Sample <- names(immdata_mixcr$data)
immdata_mixcr$meta$TLS <- factor(TLS_positivity[immdata_mixcr$meta$Sample])

#immdata_mixcr$data <- immdata_mixcr$data[names(factor(TLS_positivity[immdata_mixcr$meta$Sample]))]
#immdata_mixcr$meta <- immdata_mixcr$meta[match(immdata_mixcr$meta$Sample,names(sort(TLS_scores))),]

exp_len <- repExplore(immdata_mixcr$data, .method = "len", .col = "aa")
exp_vol <- repExplore(immdata_mixcr$data, .method = "volume")
to_keep <- order(exp_vol$Volume,decreasing = T)[1:10]
exp_cnt <- repExplore(immdata_mixcr$data[to_keep], .method = "count")

p1 <- vis(exp_vol)
exp_vol$Sample <- factor(exp_vol$Sample,levels=c(names(sort(TLS_scores))))
vis(exp_vol)
p1 
p1 <- vis(exp_vol, .by = "TLS", .meta = immdata_mixcr$meta)
p1

p2 <- vis(exp_cnt)
p2 <- vis(exp_cnt, .by = "TLS", .meta = immdata_mixcr$meta)
p2
p4 <- vis(exp_len, .by = "TLS", .meta = immdata_mixcr$meta)
exp_cnt <- repExplore(immdata_mixcr$data, .method = "count")
p2 <- vis(exp_cnt)
p2

#clonality prop
tls_neg_ids <- names(grep(x = TLS_positivity,pattern =  "low",value = T))
tls_pos_ids <- names(grep(x = TLS_positivity,pattern =  "high",value = T))
tls_pos_ids <- grep(pattern = "Tonsil",x = tls_pos_ids,invert = T,value = T)
tls_order<- names(sort(TLS_scores))


percs <- repClonality(immdata_mixcr$data, .method = "rare")
abso <- sapply(rownames(percs),function(x){
  percs[x,] * sum(clonotype_bulk_m[clonotype_bulk_m$sampleid==x,"cloneCount"])
})
abso <- t(abso)
class(abso) <- class(percs)
vis(abso) + scale_x_discrete(limits=tls_order)# + scale_y_continuous(trans = "log10",breaks = base_breaks(), labels = prettyNum)

#recompute clones per clone clount
clon_sum <- sapply(unique(clonotype_bulk_m$sampleid),function(id){
  sub <- clonotype_bulk_m[clonotype_bulk_m$sampleid==id & clonotype_bulk_m$type=="IGH",]
  one_count <- (which(sub$cloneCount == 1))
  two_three <- (which(sub$cloneCount == 2 | sub$cloneCount == 3))
  four_ten <- (which(sub$cloneCount>=4 & sub$cloneCount <= 10))
  eleven_thirty <- (which(sub$cloneCount>=11 & sub$cloneCount <= 30))
  thirtyone_hundred <- (which(sub$cloneCount>=31 & sub$cloneCount <= 100))
  over_hundred <- (which(sub$cloneCount >= 101))
  total <- sum(sub$cloneCount)
  
  a <- sum(sub[one_count,"cloneCount"])
  b <- sum(sub[two_three,"cloneCount"])
  c <- sum(sub[four_ten,"cloneCount"])
  d <- sum(sub[eleven_thirty,"cloneCount"])
  e <- sum(sub[thirtyone_hundred,"cloneCount"])
  f <- sum(sub[over_hundred,"cloneCount"])

  c("1"=a,
    "2-3"=b,
    "4-10"=c,
    "11-30"=d,
    "31-100"=e,
    "101-MAX"=f)
})
clon_sum <- data.frame(t(clon_sum),check.names = F)
clon_sum$id <- rownames(clon_sum)
#clon_sum[,1:6] <- log10(clon_sum[,1:6]+10)

melted_clon_sum <- melt(clon_sum,"id")
melted_clon_sum$TLS_signature <- factor(TLS_positivity[melted_clon_sum$id],levels=c("high","low"))
melted_clon_sum$id <- factor(melted_clon_sum$id,levels=rev(tls_order))
melted_clon_sum$variable <- factor(melted_clon_sum$variable,levels=rev(levels(melted_clon_sum$variable)))
# absolute values no logs scale
p1 <- ggplot(melted_clon_sum, aes(x = id, y = value,fill=variable)) +
  geom_bar(stat = 'identity', position = 'stack',color="black")+
  theme_bw()+
  ylab(label = "Absolute IGH clonotpe count")+
  theme(axis.text.x = element_text(angle = 90))+
  scale_fill_manual(values = rev(viridis(6,option = "C")))+ 
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
  theme(legend.position = "none")+
p1

# percentages values no logs scale
clon_sum <- sapply(unique(clonotype_bulk_m$sampleid),function(id){
    sub <- clonotype_bulk_m[clonotype_bulk_m$sampleid==id & clonotype_bulk_m$type=="IGH",]
    one_count <- (which(sub$cloneCount == 1))
    two_three <- (which(sub$cloneCount == 2 | sub$cloneCount == 3))
    four_ten <- (which(sub$cloneCount>=4 & sub$cloneCount <= 10))
    eleven_thirty <- (which(sub$cloneCount>=11 & sub$cloneCount <= 30))
    thirtyone_hundred <- (which(sub$cloneCount>=31 & sub$cloneCount <= 100))
    over_hundred <- (which(sub$cloneCount >= 101))
    total <- sum(sub$cloneCount)
    
    a <- sum(sub[one_count,"cloneCount"])
    b <- sum(sub[two_three,"cloneCount"])
    c <- sum(sub[four_ten,"cloneCount"])
    d <- sum(sub[eleven_thirty,"cloneCount"])
    e <- sum(sub[thirtyone_hundred,"cloneCount"])
    f <- sum(sub[over_hundred,"cloneCount"])
    
    c("1"=a,
      "2-3"=b,
      "4-10"=c,
      "11-30"=d,
      "31-100"=e,
      "101-MAX"=f) / total
  })
clon_sum <- data.frame(t(clon_sum),check.names = F)
clon_sum$id <- rownames(clon_sum)
melted_clon_sum <- melt(clon_sum,"id")
melted_clon_sum$TLS_signature <- factor(TLS_positivity[melted_clon_sum$id],levels=c("high","low"))

melted_clon_sum$id <- factor(melted_clon_sum$id,levels=rev(tls_order))
melted_clon_sum$variable <- factor(melted_clon_sum$variable,levels=rev(levels(melted_clon_sum$variable)))

# absolute values 
p2 <- ggplot(melted_clon_sum, aes(x = id, y = value,fill=variable)) +
  geom_bar(stat = 'identity', position = 'stack',color="black")+
  theme_bw()+
  ylab(label = "Percentage of occupied IGH repertoire" )+
  theme(axis.text.x = element_text(angle = 90))+
  scale_fill_manual(values = rev(viridis(6,option = "C")))+
  guides(fill=guide_legend(title="Clonotype counts"))

ggsave(paste0(res_folder,Sys.Date(),"_","FIGURE3_repertoireBCR.png"),
       plot=p1/p2,  
       width = 160,height=150,units = "mm")


