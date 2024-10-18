## Carter
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

####READ and CLUSTER FUNCTIONS####
returnfcs <- function(FDR_cutoff=.05,
                      metaDataFile='~/',
                      panelDataFile='~/',
                      dataDirectory='~/',
                      shape_conditions=NULL,
                      color_conditions=NULL){
  #This function generates an fcs file, subtype_markers, colors and shapes for clustering 
  require(scales);require(readxl);require(dplyr);require(flowCore)
  ##directory and metadatafile checking
  if(!dir.exists(dataDirectory)) {stop('ERR: cannot find data directory')}
  if(!file.exists(metaDataFile)) {stop('ERR: cannot find metadata.xlsx or .csv file')}
  ##readin metadata and clean
  ifelse(grepl(metaDataFile,pattern='.xls'),md <- read_excel(metaDataFile),md <- read.csv(metaDataFile,header = TRUE))#must be in xl format or csv
  md$timepoint <- factor(md$timepoint)
  md$batch <- factor(md$batch)
  rownames(md) = md$sample_id;md$sample_id <- md$sample_id
  #Make sure all files in metadata present in datadirectory
  if(!all(md$file_name %in% list.files(dataDirectory)[grep(list.files(dataDirectory),pattern = '.fcs')])){
    print(paste('ERR: not all filenames in metadata present in data folder - missing',md$file_name[!which(md$file_name %in% list.files(dataDirectory)[grep(list.files(dataDirectory),pattern = '.fcs')])],'Subsetting...'))
    md <- md[-c(!which(md$file_name %in% list.files(dataDirectory)[grep(list.files(dataDirectory),pattern = '.fcs')])),]
  }
  ##Define shapes for conditions
  if(is.null(shape_conditions)){shape_conditions <- c(0:25)[1:length(levels(md$timepoint))]}#can specify as long as number is same
  if(length(shape_conditions)!=length(levels(md$timepoint))){stop(paste0('ERR no. shapes specified is less than no. of conditions (',length(levels(md$condition)),')'))}
  names(shape_conditions) <- levels(md$timepoint)
  ## Define colors for the conditions
  if(is.null(color_conditions)){color_conditions <- hue_pal()(length(levels(md$timepoint)))}#can specify as long as number is same
  if(length(color_conditions)!=length(levels(md$timepoint))){stop(paste0('ERR no. shapes specified is less than no. of conditions (',length(levels(md$condition)),')'))}
  ## read fcs
  fcs_raw <- read.flowSet(md$file_name, path = dataDirectory, transformation = FALSE, truncate_max_range = FALSE)
  sample_ids <- rep(md$sample_id, fsApply(fcs_raw, nrow))
  panel <- read_excel(panelDataFile)
  head(data.frame(panel))
  ## Replace problematic characters
  panel$Metal <- gsub('-', '_', panel$Metal)
  panel_fcs <- pData(parameters(fcs_raw[[1]]))
  panel_fcs$desc <- gsub('-', '_', panel_fcs$desc)
  panel_fcs$desc[is.na(panel_fcs$desc)] <- paste0('NA_',which(is.na(panel_fcs$desc))) #was labelled 'todo'(mapping based on isotope for now just getting rid of NA keeping rownum) 
  # use panel$Antigen to fix description in panel_fcs
  # use metal+isotope as mapping between panel from xlsx and panel from the fcs files
  rownames(panel_fcs) = panel_fcs$name
  panel_fcs[paste0(panel$Metal,panel$Isotope,'Di'),2] <- panel$Antigen
  ## Replace paramater data in flowSet
  pData(parameters(fcs_raw[[1]])) <- panel_fcs
  ## Define variables indicating marker types
  subtype_markers <- panel$Antigen[panel$Subtype == 1]
  functional_markers <- panel$Antigen[panel$Functional == 1]
  if(!all(subtype_markers %in% panel_fcs$desc)){stop('ERR: Not all subtype_markers in panel_fcs$desc (isotopes)')}
  if(!all(functional_markers %in% panel_fcs$desc)){stop('ERR: Not all functional_markers in panel_fcs$desc (isotopes)')}
  ## arcsinh transformation and column subsetting
  fcs <- fsApply(fcs_raw, function(x, cofactor = 5){
    colnames(x) <- panel_fcs$desc
    expr <- exprs(x)
    expr <- asinh(expr[, union(subtype_markers,functional_markers)] / cofactor)
    exprs(x) <- expr
    x
  })
  return(list('fcs'=fcs,
              'subtype_markers'=subtype_markers,
              'functional_markers'=functional_markers,
              'shape_conditions'=shape_conditions,
              'color_conditions'=color_conditions,
              'sample_ids'=sample_ids,
              'meta_data'=md))
}

clusterfcs <- function(fcs=output$fcs,
                       subtype_markers = output$subtype_markers,
                       seed=1234,plottitle='consensus_plots',
                       numclusters=40){
  ## Cell population identification with FlowSOM and ConsensusClusterPlus
  require(dplyr);require(FlowSOM);require(ConsensusClusterPlus)
  set.seed(seed)
  som <- ReadInput(fcs, transform = FALSE, scale = FALSE) %>% BuildSOM(colsToUse = subtype_markers)
  ## Get the cell clustering into 100 SOM codes
  cell_clustering_som <- som$map$mapping[,1]
  ## Metaclustering into numclusters with ConsensusClusterPlus
  codes <- som$map$codes
  mc <- ConsensusClusterPlus(t(codes), maxK = numclusters, reps = 100,
                             pItem = 0.9, pFeature = 1, title = plottitle, 
                             plot = "png", clusterAlg = "hc", 
                             innerLinkage = "average", finalLinkage = "average",
                             distance = "euclidean", seed = 1234)
  
  ## Get cluster ids for each cell
  code_clustering <- mc[[numclusters]]$consensusClass#metaclusters consensus
  cell_clustering <- code_clustering[cell_clustering_som]#cell clustering from som
  return(list('code_clustering'=code_clustering,'cell_clustering'=cell_clustering,'metaclusters'=mc))
}

#Read output
setwd("E:/WJHLab/cytofAtlas/updated_data/J14113_PBMC/Myeloid/Myeloid")

workd<-getwd()

#remove patient info -> replace it with pub ID
pubID <- read_xlsx("E:/WJHLab/cytofAtlas/updated_data/J14113_PBMC/J14113_PubIDs_20240419.xlsx")
metadata <- read_xlsx(paste0(workd,"/Config","/metadata.xlsx"))

# Extract last three digits of pubID in df2
pubID$last_three_digits <- substr(pubID$`Subject ID`, start = 15, stop = nchar(pubID$`Subject ID`[1]))
pubID$last_three_digits <- ifelse(nchar(pubID$last_three_digits) == 2, paste0(pubID$last_three_digits, "0"), pubID$last_three_digits)
pubID <- pubID[-c(53:nrow(pubID)),]
pubid <- c()
for (i in 1:nrow(metadata)){
  pubid <- c(pubid,pubID$`Pub ID`[which(metadata$patient_id[i]==pubID$last_three_digits)])
}
metadata$pubID <- pubid
metadata <- metadata[,!names(metadata) %in% c("patient_id","last_three_digits")]
write.xlsx(metadata,paste0(workd,"/Config/J14113_metadataM_wPubID.xlsx"))

output <- readRDS(paste0(workd,"/backup_output.rds"))
output$meta_data <- metadata

#set levels
clusterMergeFile = paste0(workd,"/Config",'/merge.xlsx')

cluster_merging <- read_xlsx(clusterMergeFile)
cm<-read_xlsx(clusterMergeFile)
cluster_merging <- read_xlsx(clusterMergeFile)

clusterlevels = c("Mono_I",
                  "Mono_II",
                  "Mono_III",
                  "T",
                  "B",
                  "DC_I",
                  "DC_II",
                  "Gran")

samplevels <- c("batch1_BC1","batch1_BC2","batch1_BC3","batch1_BC4","batch1_BC5","batch1_BC6","batch1_BC7",
                "batch1_BC8","batch2_BC1","batch2_BC2","batch2_BC3","batch2_BC4","batch2_BC5","batch2_BC6",
                "batch3_BC1","batch3_BC2","batch3_BC3","batch3_BC4","batch3_BC5","batch3_BC6","batch3_BC7",
                "batch3_BC8","batch3_BC9","batch4_BC1","batch4_BC2","batch4_BC3","batch4_BC4","batch4_BC5",
                "batch4_BC6","batch4_BC7","batch4_BC8","batch5_BC1","batch5_BC2","batch5_BC3","batch5_BC4",
                "batch5_BC5","batch5_BC6","batch5_BC7","batch5_BC8","batch6_BC1","batch6_BC2","batch6_BC3",
                "batch6_BC4","batch6_BC5","batch6_BC6","batch6_BC7","batch6_BC8","batch7_BC1","batch7_BC2",
                "batch7_BC3","batch7_BC4","batch7_BC5","batch7_BC6","batch7_BC7","batch7_BC8","batch7_BC9",
                "batch8_BC1","batch8_BC2","batch8_BC3","batch8_BC4","batch8_BC5","batch8_BC6","batch8_BC7",
                "batch8_BC8","batch8_BC9","batch9_BC1","batch9_BC2","batch9_BC3","batch9_BC4","batch9_BC5",
                "batch9_BC6","batch9_BC7","batch9_BC8","batch9_BC9")

ptlevels <- unique(metadata$pubID)
  
conditionlevels=c("A","B")

timepointlevels=c("C1D1","C3D1","C4D1")

responselevels=c("PD", "SD", "PR")

survivallevels=c("2.6", "4.47", "4.73", "4.93", "5.03", "5.16", "5.88", "6.67",
                 "7.1", "8.64", "8.71", "9.76", "10.32", "12.71", "13.5", "18.17",
                 "20.3", "20.73", "21.29", "24.97", "NA")

batchlevels=c("1","2","3","4","5","6","7","8","9")



clustercolors <- as.character(c(cols25(n=25),alphabet(n=19)))

mm1 <- match(output$cell_clustering, cluster_merging$original_cluster)

cell_clustering1m <- cluster_merging$new_cluster[mm1]

output$cell_clustering1m <- cell_clustering1m
saveRDS(output,paste0(workd,"/FinalRun/output_wPubid.rds"))

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
                                 fileName = paste0(workd,'/FinalRun/clusteringheatmap_merged_new.pdf'))
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
ggdf$condition <- factor(output$meta_data$condition[match(ggdf$sample_id,output$meta_data$sample_id)], levels = conditionlevels)

ggdf_A<-ggdf[ggdf$condition=="A",]
#ggdf_A<-ggdf_A[ggdf_A$timepoint!="C4D1",]
#ggdf_A<-ggdf_A[ggdf_A$pubID!="038",] #remove unpaired sample

ggdf_B<-ggdf[ggdf$condition=="B",]
#ggdf_B<-ggdf_B[ggdf_B$timepoint!="C4D1",]

ggdf_wide<-ggdf[,c("cluster",'pubID','timepoint','proportion','condition')]
ggdf_wide2<-ggdf_wide %>% pivot_wider(names_from=timepoint, values_from=proportion) #pivot the df by timepoint
ggdf_wide2<-ggdf_wide2 %>% mutate(FC = C3D1/C1D1) #take fold changes "C2D1/C1D1 or C3D1/C1D1"
ggdf_wide2<-ggdf_wide2[!is.na(ggdf_wide2$FC),] #remove all NAs (unpaired timepoints)
ggdf_wide2[is.infinite(ggdf_wide2$FC),]$FC<-100 #replace infinite fold changes (when there was 0 in C1D1) with '100'
ggdf_wide2_A<-ggdf_wide2[ggdf_wide2$condition=="A",]
ggdf_wide2_B<-ggdf_wide2[ggdf_wide2$condition=="B",]

lineplot_r_A <- ggplot(ggdf_A, aes(x=timepoint, y=proportion))+
  facet_wrap(~cluster, scales='free',ncol=5)+
  scale_shape_manual(values=c(0:46,0,1))+
  geom_point(aes(shape=pubID))+
  geom_line(show.legend = T, linewidth=0.25, aes(group=pubID))+
  ggtitle('Changes in abundance - Arm A')+
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
lineplot_r_B <- ggplot(ggdf_B, aes(x=timepoint, y=proportion))+
  facet_wrap(~cluster, scales='free',ncol=5)+
  scale_shape_manual(values=c(0:46,0,1))+
  geom_point(aes(shape=pubID))+
  geom_line(show.legend = T, linewidth=0.25, aes(group=pubID))+
  ggtitle('Changes in abundance - Arm B')+
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
        legend.title = element_blank(),legend.key.size = unit(1,'lines'),
        legend.text = element_text(size=10),
        legend.key = element_rect(fill="white"))+
  scale_y_sqrt()

#for T cell oriented
pdf(paste0(workd,"/FinalRun/abundanceLine_col_tp_AvB.pdf"), width=7, height=8.5);lineplot_r_A; lineplot_r_B;dev.off()
pdf(paste0(workd,"/FinalRun/abundanceLine_col_tp_AvB_wstats.pdf"), width=7, height=8.5)
lineplot_r_A+stat_compare_means(paired=T, size=1)
lineplot_r_B+stat_compare_means(paired=T, size=1)
dev.off()

#for myeloid oriented
pdf(paste0(workd,"/FinalRun/abundanceLine_col_tp_AvB.pdf"), width=7, height=5);lineplot_r_A; lineplot_r_B;dev.off()
pdf(paste0(workd,"/FinalRun/abundanceLine_col_tp_AvB_wstats.pdf"), width=7, height=5)
lineplot_r_A+stat_compare_means(paired=T, size=1)
lineplot_r_B+stat_compare_means(paired=T, size=1)
dev.off()

###without separating out the arms
lineplot_r_all <- ggplot(ggdf, aes(x=timepoint, y=proportion))+
  facet_wrap(~cluster, scales='free',ncol=5)+
  scale_shape_manual(values=c(0:46,0,1))+
  geom_point(aes(shape=pubID))+
  geom_line(show.legend = T, linewidth=0.25, aes(group=pubID))+
  ggtitle('Changes in abundance - both Arms')+
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

#for lymphoid
pdf(paste0(workd,"/FinalRun/abundanceLine_col_tp_both.pdf"), width=7, height=8.5);lineplot_r_all;dev.off()
pdf(paste0(workd,"/FinalRun/abundanceLine_col_tp_both_wstats.pdf"), width=7, height=8.5)
lineplot_r_all+stat_compare_means(paired=T, size=1)
dev.off()

#for myeloid
pdf(paste0(workd,"/FinalRun/abundanceLine_col_tp_both.pdf"), width=7, height=5);lineplot_r_all;dev.off()
pdf(paste0(workd,"/FinalRun/abundanceLine_col_tp_both_wstats.pdf"), width=7, height=5)
lineplot_r_all+stat_compare_means(paired=T, size=1)
dev.off()


ggfcA<-ggplot(ggdf_wide2_A, aes(x=reorder(paste(pubID,cluster), FC), y=FC))+
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
  ggtitle('Fold changes in abundance - Arm A')

ggfcB<-ggplot(ggdf_wide2_B, aes(x=reorder(paste(pubID,cluster), FC), y=FC))+
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
  ggtitle('Fold changes in abundance - Arm B')

dev.off()

#for T cell oriented
pdf(paste0(workd,'/FinalRun/abundanceFCbars.pdf'),width=7,height=7);ggfcA;ggfcB;dev.off()

#for myeloid oriented
pdf(paste0(workd,'/FinalRun/abundanceFCbars.pdf'),width=7,height=5);ggfcA;ggfcB;dev.off()

#statistic for fold change differences between Arms A and B
ggdf_wide2_p<- ggdf_wide2 %>% 
  dplyr::select(cluster, pubID, condition, FC) %>% 
  group_by(condition, cluster) %>% 
  summarise(value=list(FC)) %>%
  spread(condition, value) %>%
  group_by(cluster) %>%
  mutate(p_value = wilcox.test(unlist(A),unlist(B))$p.value,
         diffAB = mean(unlist(A))-mean(unlist(B))) #Arm A - Arm B = what does nivo do in addition to GVAX

write.csv(ggdf_wide2_p[,c('cluster','p_value','diffAB')],paste0(workd,'/FinalRun/abundanceFC_AvB.csv'))


#==========DIFFERENTIAL STATE PLOTS=============

exprmntbl <- data.frame(fsApply(output$fcs,exprs)[, c(output$functional_markers)],
                        sample_id = output$sample_ids, 
                        cluster = output$cell_clustering1m) %>%
  group_by(sample_id, cluster) %>%
  summarize_all(funs(mean))
exprmntbl$cluster <- factor(exprmntbl$cluster, levels=clusterlevels)
exprmntbl$pubID <- factor(output$meta_data$pubID[match(exprmntbl$sample_id,output$meta_data$sample_id)], levels = ptlevels)
exprmntbl$timepoint <- factor(output$meta_data$timepoint[match(exprmntbl$sample_id,output$meta_data$sample_id)], levels = timepointlevels)
exprmntbl$condition <- factor(output$meta_data$condition[match(exprmntbl$sample_id,output$meta_data$sample_id)], levels = conditionlevels)

ggdf2<-melt(exprmntbl, id.var=c("cluster","condition","sample_id","pubID","timepoint"))
#ggdf2<-ggdf2[ggdf2$timepoint!="C4D1",]
ggdf2_A<-ggdf2[ggdf2$condition=="A",]
#ggdf2_A<-ggdf2_A[ggdf2_A$pubID!="038",]
ggdf2_B<-ggdf2[ggdf2$condition=="B",]

ggdf2_wide<-ggdf2[,c("cluster",'pubID','timepoint','variable','value','condition')]
ggdf2_wide2<-ggdf2_wide %>% pivot_wider(names_from=timepoint, values_from=value) #pivot the df by timepoint
ggdf2_wide2<-ggdf2_wide2 %>% mutate(FC = C3D1/C1D1) #take fold changes
ggdf2_wide2<-ggdf2_wide2[!is.na(ggdf2_wide2$FC),] #remove all NAs (unpaired timepoints)
ggdf2_wide2[is.infinite(ggdf2_wide2$FC),]$FC<-100 #replace infinite fold changes (when there was 0 in C1D1) with '100'
ggdf2_wide2_A<-ggdf2_wide2[ggdf2_wide2$condition=="A",]
ggdf2_wide2_B<-ggdf2_wide2[ggdf2_wide2$condition=="B",]

fmlistplot <- c(output$functional_markers)
fmlistplot[which(fmlistplot=="41BB")]<-"X41BB"
fmlistplot[which(fmlistplot=="Arginase I")]<-"Arginase.I"
fmlistplot[which(fmlistplot=="HLA-DR")]<-"HLA.DR"

dev.off()

#for lymphoid
pdf('functionalLineplots_A.pdf',width=8.5,height=8.5)
#for myeloid
pdf(paste0(workd,'/FinalRun/functionalLineplots_A.pdf'),width=8.5,height=5)

for(i in 1:length(fmlistplot)){
  ggp <- ggplot(ggdf2_A[ggdf2_A$variable==fmlistplot[i],], aes(x=timepoint, y=value))+
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
pdf('functionalLineplots_B.pdf',width=8.5,height=8.5)
#for myeloid
pdf(paste0(workd,'/FinalRun/functionalLineplots_B.pdf'),width=8.5,height=5)

for(i in 1:length(fmlistplot)){
  ggp <- ggplot(ggdf2_B[ggdf2_B$variable==fmlistplot[i],], aes(x=timepoint, y=value))+
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

###for both arms

#for lymphoid
pdf('functionalLineplots_all.pdf',width=8.5,height=8.5)
#for myeloid
pdf(paste0(workd,'/FinalRun/functionalLineplots_all.pdf'),width=8.5,height=5)

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
pdf('functionalFCplots_AvB.pdf',width=7,height=7)
#for myeloid
pdf(paste0(workd,'/FinalRun/functionalFCplots_AvB.pdf'),width=7,height=5)

for(i in 1:length(fmlistplot)){
  ggfc2A<-ggplot(ggdf2_wide2_A[ggdf2_wide2_A$variable==fmlistplot[i],], aes(x=reorder(paste(pubID,cluster), FC), y=FC))+
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
    ggtitle(paste(fmlistplot[i],'Fold changes in marker expression - Arm A',sep=": "))
  print(ggfc2A)
  
  ggfc2B<-ggplot(ggdf2_wide2_B[ggdf2_wide2_B$variable==fmlistplot[i],], aes(x=reorder(paste(pubID,cluster), FC), y=FC))+
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
    ggtitle(paste(fmlistplot[i],'Fold changes in marker expression - Arm B',sep=": "))
  print(ggfc2B) 
}
dev.off()

ggdf2_wide2_p<- ggdf2_wide2 %>% 
  dplyr::select(cluster, pubID, variable, condition, FC) %>% 
  group_by(condition, cluster, variable) %>% 
  summarise(value=list(FC), .groups = 'drop') %>%
  spread(condition, value) %>%
  group_by(cluster, variable) %>%
  mutate(p_value = wilcox.test(unlist(A),unlist(B))$p.value,
         diffAB = mean(unlist(A))-mean(unlist(B))) #Arm A - Arm B = what does nivo do in addition to GVAX

write.csv(ggdf2_wide2_p[,c('cluster','variable','p_value','diffAB')],paste0(workd,'/FinalRun/functionalFC_AvB.csv'))


saveRDS(output,'backup_output.rds')
