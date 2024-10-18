#load output
rm(list = ls())
library(reshape2);
library(diffcyt);
library(ggpubr);
library(readxl);
library(randomcoloR);
library(pals);
library(ggplot2);
library(colorspace);
library(dplyr);
library(scales);
library(tibble);
library(lme4);
library(multcomp);
library(mvtnorm);
library(flowCore);
library(umap);
library(matrixStats); 
library(premessa); 
library(limma);
library(ComplexHeatmap);
library(circlize); 
library(stringr); 
library(edgeR);
library(corrplot);
library(ggridges);
library(Hmisc);
library(pheatmap);
library(RColorBrewer);
library(tidyr);
library(openxlsx)
#set levels
setwd("E:/WJHLab/cytofAtlas/updated_data/J13108/finalrun_revised")
output <- readRDS("output113108_wPubidFcs.rds")
md <- read.xlsx("E:/WJHLab/cytofAtlas/updated_data/J13108/Config/metadata_wPubID.xlsx")
current_names <- names(output$fcs@frames)
new_names <- md$file_name[match(current_names,md$file_name_old)]

for (i in seq_along(current_names)) {
  # Get the current object from the environment
  obj <- get(current_names[i], envir = output$fcs@frames)
  
  # Assign the object back with the correctly ordered new name
  assign(new_names[i], obj, envir = output$fcs@frames)
  
  # Remove the object with the old name
  rm(list = current_names[i], envir = output$fcs@frames)
}
md <- md[,!names(md)=="file_name_old"]
md <- md[,!names(md)=="sample_id_old"]

md$sample_id_old <- sub("\\.fcs$", "", md$file_name_old)

indices <- match(output$sample_ids,md$sample_id_old)
output$sample_ids <- md$sample_id[indices]
output$meta_data <- md

saveRDS(output,"output_forCellxGene.rds")

oldoutput <- readRDS("backup_output.rds")
write.xlsx(md,"E:/WJHLab/cytofAtlas/updated_data/J13108/Config/metadata_wPubID.xlsx")

saveRDS(output,)
clusterlevels<-c("Th",
                 "ThN",
                 "Th1",
                 "ThCM",
                 "ThEM_Acv",
                 "ThEM",
                 "Tc",
                 "Tc1",
                 "TcN",
                 "Tc_Acv",
                 "TcCM",
                 "TcEFF",
                 "TcEFF_Acv",
                 "Treg",
                 "DNT",
                 "NK",
                 "B",
                 "UA")
clusterMergeFile = 'J13108_merging_revisedforintegration.xlsx'

timepointlevels=c('0wk','7wk')
ptlevels<-unique(md$Pub_ID)

clusterlevels2<-c("Th","Th_Acv","ThN","Th1","ThCM_Acv","ThEM_Acv","ThEM","Tc","Tc1","TcN","Tc_Acv","TcEFF_Acv","Treg","DNT","NK","B","UA")
clustercolors <- as.character(c(cols25(n=25),alphabet(n=19)))
#For clustering heatmap merged
plot_clustering_heatmap_wrapper2 <- function(fcs, cell_clustering, nclusters=40,
                                             color_clusters=clustercolors,
                                             subtype_markers,
                                             fileName = 'clusteringheatmap.pdf'){
  require(matrixStats);require(dplyr);require(RColorBrewer);require(pheatmap);require(readxl);require(flowCore);require(scales)
  ## Will output the heatmap object and print it 
  color_clusters=color_clusters
  #get expression
  expr <- fsApply(fcs, exprs);expr <-expr[,subtype_markers]
  ## Scale expression of all markers to values between 0 and 1
  rng <- colQuantiles(expr, probs = c(0.01, 0.99))
  expr01 <- t((t(expr) - rng[, 1]) / (rng[, 2] - rng[, 1]))
  expr01[expr01 < 0] <- 0; expr01[expr01 > 1] <- 1;expr01 <-expr01[,subtype_markers]
  ## Calculate the mean expression##################################################
  pdf(fileName, width=8, height=11) 
  expr_mean <- data.frame(expr, cell_clustering = cell_clustering, check.names = FALSE) %>%
    group_by(cell_clustering) %>% summarize_all(funs(mean))
  expr01_mean <- data.frame(expr01, cell_clustering = cell_clustering, check.names = FALSE) %>%
    group_by(cell_clustering) %>% summarize_all(funs(mean))
  
  ## Calculate cluster frequencies
  
  clustering_table <- as.numeric(table(cell_clustering))
  clustering_prop <- round(clustering_table / sum(clustering_table) * 100, 2)
  
  ## Sort the cell clusters with hierarchical clustering
  
  d <- dist(expr_mean[, colnames(expr)], method = "euclidean")
  cluster_rows <- hclust(d, method = "average")
  
  expr_heat <- as.matrix(expr01_mean[, colnames(expr01)])
  rownames(expr_heat) <- expr01_mean$cell_clustering
  
  ## Colors for the heatmap
  
  #legend_breaks = seq(from = 0, to = 1, by = 0.2)
  #labels_row <- expr01_mean$cell_clustering
  
  labels_row <- paste0(expr01_mean$cell_clustering, " ")
  
  ## Annotation for the original clusters
  
  annotation_row <- data.frame(Cluster = factor(expr01_mean$cell_clustering))
  rownames(annotation_row) <- rownames(expr_heat)
  color_clusters1 <- color_clusters[1:nlevels(annotation_row$Cluster)]
  names(color_clusters1) <- levels(annotation_row$Cluster)
  annotation_colors <- list(Cluster = color_clusters1)
  
  p <- pheatmap(expr_heat, 
                color = rev(brewer.rdylbu(100)), 
                cluster_cols = T,
                cluster_rows = F, 
                labels_row = labels_row,
                #scale="column",
                display_numbers = F, 
                number_color = "black",
                fontsize = 9, fontsize_number = 6,  
                #legend_breaks = legend_breaks,
                annotation_row = annotation_row, 
                annotation_colors = annotation_colors,
                cellwidth = 8,
                cellheight = 8,
                border_color = NA,
                annotation_legend = F
  )
  dev.off() 
  print('Colors:')
  print(color_clusters)
  print(p);return(p)
}


plot_clustering_heatmap_wrapper2(fcs=output$fcs,
                                 color_clusters = clustercolors,
                                 cell_clustering = factor(output$cell_clustering1m, levels=clusterlevels), 
                                 subtype_markers=output$subtype_markers,
                                 fileName = 'clusteringheatmap_merged_new.pdf')
dev.off()

#==========DIFFERENTIAL ABUNDANCE PLOTS===========

counts_table <- table(output$cell_clustering1m, output$sample_ids)
props_table <- t(t(counts_table) / colSums(counts_table)) * 100
counts <- as.data.frame.matrix(counts_table)
props <- as.data.frame.matrix(props_table)

ggdf <- melt(data.frame(cluster = rownames(props), props, check.names = FALSE),
             id.vars = "cluster", value.name = "proportion", 
             variable.name = "sample_id")
ggdf$sample_id <- factor(ggdf$sample_id)
ggdf$cluster <- factor(ggdf$cluster, levels=clusterlevels)
ggdf$timepoint <- factor(output$meta_data$timepoint[match(ggdf$sample_id,output$meta_data$sample_id)], levels = timepointlevels)
ggdf$pubID <- factor(output$meta_data$pubID[match(ggdf$sample_id,output$meta_data$sample_id)], levels = ptlevels)
ggdf$study <- "J13108"

counts_tableb <- table(output_j14113$cell_clustering1m, output_j14113$sample_ids)
props_tableb <- t(t(counts_tableb) / colSums(counts_tableb)) * 100
countsb <- as.data.frame.matrix(counts_tableb)
propsb <- as.data.frame.matrix(props_tableb)

ggdfb <- melt(data.frame(cluster = rownames(propsb), propsb, check.names = FALSE),
             id.vars = "cluster", value.name = "proportion", 
             variable.name = "sample_id")
ggdfb$sample_id <- factor(ggdfb$sample_id)
ggdfb$cluster <- factor(ggdfb$cluster, levels=clusterlevels2)
ggdfb$timepoint <- factor(output_j14113$meta_data$timepoint[match(ggdfb$sample_id,output_j14113$meta_data$sample_id)], levels = c('C1D1','C3D1','C4D1'))
ggdfb$pubID <- factor(output_j14113$meta_data$pubID[match(ggdfb$sample_id,output_j14113$meta_data$sample_id)], levels = unique(output_j14113$meta_data$pubID))
ggdfb$condition <- factor(output_j14113$meta_data$condition[match(ggdfb$sample_id,output_j14113$meta_data$sample_id)], levels = c('A','B'))
ggdfb$study <- "J14113"
ggdfb<-ggdfb[ggdfb$timepoint!="C4D1",]
ggdfb<-ggdfb[ggdfb$pubID!="J14113_P86",]
ggdfb<-ggdfb[ggdfb$condition=="A",]

ggdf_wide<-ggdf[,c("cluster",'pubID','timepoint','proportion')]
ggdf_wide2<-ggdf_wide %>% pivot_wider(names_from=timepoint, values_from=proportion) #pivot the df by timepoint
ggdf_wide2<-ggdf_wide2 %>% mutate(FC = `7wk`/`0wk`) #take fold changes
ggdf_wide2<-ggdf_wide2[!is.na(ggdf_wide2$FC),] #remove all NAs (unpaired timepoints)
ggdf_wide2[is.infinite(ggdf_wide2$FC),]$FC<-100 #replace infinite fold changes (when there was 0 in C1D1) with '100'
ggdf_wide2$study <- "J13108"
colnames(ggdf_wide2)[3:4]<-c('C1D1','C3D1')

ggdf_wideb<-ggdfb[,c("cluster",'pubID','timepoint','proportion')]
ggdf_wide2b<-ggdf_wideb %>% pivot_wider(names_from=timepoint, values_from=proportion) #pivot the df by timepoint
ggdf_wide2b$C1D1<-as.numeric(unlist(ggdf_wide2b$C1D1))
ggdf_wide2b<-ggdf_wide2b %>% mutate(FC = C3D1/C1D1) #take fold changes
ggdf_wide2b<-ggdf_wide2b[!is.na(ggdf_wide2b$FC),] #remove all NAs (unpaired timepoints)
ggdf_wide2b[is.infinite(ggdf_wide2b$FC),]$FC<-100 #replace infinite fold changes (when there was 0 in C1D1) with '100'
ggdf_wide2b$study <- "J14113"

includeclusters<-intersect(ggdf_wide2$cluster,ggdf_wide2b$cluster)

ggdf_widecombined <- rbind(ggdf_wide2,ggdf_wide2b)
ggdf_widecombined <- ggdf_widecombined[ggdf_widecombined$cluster %in% includeclusters,]

lineplot_r <- ggplot(ggdf, aes(x=timepoint, y=proportion))+
  facet_wrap(~cluster, scales='free',ncol=5)+
  scale_shape_manual(values=c(0:46,0,1))+
  geom_point(aes(shape=pubID))+
  geom_line(show.legend = T, linewidth=0.25, aes(group=pubID))+
  ggtitle('Changes in abundance')+
  theme(axis.text.x = element_text(size=8, angle = 45, hjust=1, vjust=1, color="black"),
        axis.title.x = element_blank(),
        axis.ticks=element_line(color="black",size=0.25),
        axis.line.x = element_line(size=0.25),
        axis.line.y = element_line(size=0.25),
        axis.text.y = element_text(color="black"),
        axis.title.y = element_blank(),
        strip.background = element_rect(fill=NA),
        strip.text = element_text(size=10),
        panel.background = element_rect(fill="white"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.title = element_blank(),
        legend.key.size = unit(1,'lines'),
        legend.text = element_text(size=10),
        legend.key = element_rect(fill="white"))+
  scale_y_sqrt()

#for T cell oriented
pdf("abundanceLine_col_tp.pdf", width=7, height=8.5);lineplot_r;dev.off()
pdf("abundanceLine_col_tp_wstats.pdf", width=7, height=8.5)
lineplot_r+stat_compare_means(paired=T, size=1)
dev.off()

ggfc<-ggplot(ggdf_wide2, aes(x=reorder(paste(pubID,cluster), FC), y=FC))+
  facet_wrap(~cluster, scales='free',ncol=5)+
  geom_bar(stat='identity', aes(fill=FC))+
  geom_hline(yintercept=1, linetype='dotted', linewidth=0.25)+
  coord_cartesian(ylim=c(0,2))+
  theme(strip.background = element_rect(fill=NA),
        panel.grid = element_blank(),
        panel.border = element_rect(color='black', fill=NA),
        panel.background = element_rect(fill='white'),
        axis.line = element_blank(),
        axis.text.y = element_text(color='black'),
        axis.text.x = element_blank(),
        axis.ticks.y = element_line(color = 'black'),
        axis.ticks.x = element_blank())+
  xlab('pubID')+
  scale_fill_gradient2(low=muted('blue'),mid='lightgray',high=muted('red'),midpoint=1,limits=c(0,2),oob=squish)+
  ggtitle('Fold changes in abundance')

#for T cell oriented
pdf('abundanceFCbars.pdf',width=7,height=7);ggfc;dev.off()

#statistic for fold change differences between Arms A and B
ggdf_wide2_p<- ggdf_widecombined %>% 
  dplyr::select(cluster, pubID, study, FC) %>% 
  group_by(study, cluster) %>% 
  summarise(value=list(FC)) %>%
  spread(study, value) %>%
  group_by(cluster) %>%
  mutate(p_value = wilcox.test(unlist(J13108),unlist(J14113))$p.value,
         diffJ13108J14113 = mean(unlist(J13108))-mean(unlist(J14113))) #Arm A - Arm B = what does nivo do in addition to GVAX

write.csv(ggdf_wide2_p[,c('cluster','p_value','diffJ13108J14113')],'abundanceFC_J13108vJ14113.csv')


#==========DIFFERENTIAL STATE PLOTS=============

exprmntbl <- data.frame(fsApply(output$fcs,exprs)[, c(output$functional_markers)],
                        sample_id = output$sample_ids, 
                        cluster = output$cell_clustering1m) %>%
  group_by(sample_id, cluster) %>%
  summarize_all(funs(mean))
exprmntbl$cluster <- factor(exprmntbl$cluster, levels=clusterlevels)
exprmntbl$pubID <- factor(output$meta_data$pubID[match(exprmntbl$sample_id,output$meta_data$sample_id)], levels = ptlevels)
exprmntbl$timepoint <- factor(output$meta_data$timepoint[match(exprmntbl$sample_id,output$meta_data$sample_id)], levels = timepointlevels)
ggdf2<-melt(exprmntbl, id.var=c("cluster","sample_id","pubID","timepoint"))
#ggdf2<-ggdf2[ggdf2$timepoint!="C4D1",]
ggdf2$study<-"J13108"

ggdf2_wide<-ggdf2[,c("cluster",'pubID','timepoint','variable','value','study')]
ggdf2_wide2<-ggdf2_wide %>% pivot_wider(names_from=timepoint, values_from=value) #pivot the df by timepoint
ggdf2_wide2<-ggdf2_wide2 %>% mutate(FC = `7wk`/`0wk`) #take fold changes
ggdf2_wide2<-ggdf2_wide2[!is.na(ggdf2_wide2$FC),] #remove all NAs (unpaired timepoints)
ggdf2_wide2[is.infinite(ggdf2_wide2$FC),]$FC<-100 #replace infinite fold changes (when there was 0 in C1D1) with '100'
colnames(ggdf2_wide2)[5:6]<-c('C1D1','C3D1')

exprmntblb <- data.frame(fsApply(output_j14113$fcs,exprs)[, c(output$functional_markers)], #need to use the same markers as J13108 so that statistical tests will run properly
                        sample_id = output_j14113$sample_ids, 
                        cluster = output_j14113$cell_clustering1m) %>%
  group_by(sample_id, cluster) %>%
  summarize_all(funs(mean))
exprmntblb$cluster <- factor(exprmntblb$cluster, levels=clusterlevels2)
exprmntblb$pubID <- factor(output_j14113$meta_data$pubID[match(exprmntblb$sample_id,output_j14113$meta_data$sample_id)], levels = unique(output_j14113$meta_data$pubID))
exprmntblb$timepoint <- factor(output_j14113$meta_data$timepoint[match(exprmntblb$sample_id,output_j14113$meta_data$sample_id)], levels = c('C1D1','C3D1','C4D1'))
exprmntblb$condition <- factor(output_j14113$meta_data$condition[match(exprmntblb$sample_id,output_j14113$meta_data$sample_id)], levels = c('A','B'))

ggdf2b<-melt(exprmntblb, id.var=c("cluster","sample_id","pubID","timepoint","condition"))
ggdf2b<-ggdf2b[ggdf2b$condition=="A",]
ggdf2b<-ggdf2b[ggdf2b$timepoint!="C4D1",]
ggdf2b<-ggdf2b[ggdf2b$pubID!="J14113_P86",]
ggdf2b$study<-"J14113"

ggdf2_wideb<-ggdf2b[,c("cluster",'pubID','timepoint','variable','value','study')]
ggdf2_wide2b<-ggdf2_wideb %>% pivot_wider(names_from=timepoint, values_from=value) #pivot the df by timepoint
ggdf2_wide2b<-ggdf2_wide2b %>% mutate(FC = C3D1/C1D1) #take fold changes
ggdf2_wide2b<-ggdf2_wide2b[!is.na(ggdf2_wide2b$FC),] #remove all NAs (unpaired timepoints)
ggdf2_wide2b[is.infinite(ggdf2_wide2b$FC),]$FC<-100 #replace infinite fold changes (when there was 0 in C1D1) with '100'

includeclusters2<-intersect(ggdf2_wide2$cluster,ggdf2_wide2b$cluster)

ggdf_widecombined2 <- rbind(ggdf2_wide2,ggdf2_wide2b)
ggdf_widecombined2 <- ggdf_widecombined2[ggdf_widecombined2$cluster %in% includeclusters2,]




fmlistplot <- c(output$functional_markers)
#fmlistplot[which(fmlistplot=="41BB")]<-"X41BB"
dev.off()

#for lymphoid
pdf('functionalLineplots.pdf',width=8.5,height=8.5)

for(i in 1:length(fmlistplot)){
  ggp <- ggplot(ggdf2[ggdf2$variable==fmlistplot[i],], aes(x=timepoint, y=value))+
    ggtitle(fmlistplot[i])+
    facet_wrap(~cluster, scales='free',ncol=6)+
    scale_shape_manual(values=c(0:46,0,1))+
    geom_point(aes(shape=pubID))+
    geom_line(show.legend = F, aes(group=pubID), linewidth=0.25)+
    theme(axis.text.x = element_text(size=8, angle = 45, hjust=1, vjust=1, color="black"),
          axis.title.x = element_blank(),
          axis.ticks=element_line(color="black"),
          axis.line.x = element_line(size=0.25),
          axis.line.y = element_line(size=0.25),
          axis.text.y = element_text(color="black"),
          axis.title.y = element_blank(),
          strip.background = element_rect(fill=NA),
          strip.text = element_text(size=10),
          panel.background = element_rect(fill="white"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.title = element_blank(),
          legend.key.size = unit(1,'lines'),
          legend.text = element_text(size=10),
          legend.key = element_rect(fill="white"))
  print(ggp)
  print(ggp+stat_compare_means(paired=T, size=1))
}
dev.off()


#for lymphoid
pdf('functionalFCplots.pdf',width=7,height=7)

for(i in 1:length(fmlistplot)){
  ggfc2<-ggplot(ggdf2_wide2[ggdf2_wide2$variable==fmlistplot[i],], aes(x=reorder(paste(pubID,cluster), FC), y=FC))+
    facet_wrap(~cluster, scales='free',ncol=5)+
    geom_bar(stat='identity', aes(fill=FC))+
    geom_hline(yintercept=1, linetype='dotted', linewidth=0.25)+
    coord_cartesian(ylim=c(0,2))+
    theme(strip.background = element_rect(fill=NA),
          panel.grid = element_blank(),
          panel.border = element_rect(color='black', fill=NA),
          panel.background = element_rect(fill='white'),
          axis.line = element_blank(),
          axis.text.y = element_text(color='black'),
          axis.text.x = element_blank(),
          axis.ticks.y = element_line(color = 'black'),
          axis.ticks.x = element_blank())+
    xlab('patients')+
    scale_fill_gradient2(low=muted('blue'),mid='lightgray',high=muted('red'),midpoint=1,limits=c(0,2),oob=squish)+
    ggtitle(paste(fmlistplot[i],'Fold changes in marker expression',sep=": "))
  print(ggfc2)
}
dev.off()

ggdf2_wide2_p<- ggdf_widecombined2 %>% 
  dplyr::select(cluster, pubID, variable, study, FC) %>% 
  group_by(study, cluster, variable) %>% 
  summarise(value=list(FC), .groups = 'drop') %>%
  spread(study, value) %>%
  group_by(cluster, variable) %>%
  mutate(p_value = wilcox.test(unlist(J13108),unlist(J14113))$p.value,
         diffJ13108J14113 = mean(unlist(J13108))-mean(unlist(J14113))) #Arm A - Arm B = what does nivo do in addition to GVAX

write.csv(ggdf2_wide2_p[,c('cluster','variable','p_value','diffJ13108J14113')],'functionalFC_J13108vJ14113.csv')
