library(Seurat)
library(flowCore)
library(readxl)
library(SeuratObject)
library(SeuratDisk)
rm(list = ls())

### 14113M
setwd("E:/WJHLab/cytofAtlas/updated_data/J14113_PBMC/Myeloid/Myeloid")
output <- readRDS("backup_output_wjh.rds")
md <- read_xlsx("E:/WJHLab/cytofAtlas/updated_data/J14113_PBMC/Myeloid/Myeloid/Config/J14113_metadataM_wPubID.xlsx")

mm1 <- match(output$sample_ids, md$sample_id)
batch_id <- md$batch[mm1]
condition <- md$condition[mm1]
timepoint <- md$timepoint[mm1]
#patient_id <- md$patient_id[mm1]
response <- md$response[mm1]
survival <- md$survival[mm1]
pubID <- md$pubID[mm1]
output$timepoint <- timepoint
output$batch_id <- batch_id
output$condition <- condition
output$pubID <- pubID
output$response <- response
output$survival <- survival
output$meta_data <- md
saveRDS(output,paste0(getwd(),"/cellxgene_withPubID/output_wPubid.rds"))
output <- readRDS(paste0(getwd(),"/cellxgene_withPubID/output14113M_wPubid.rds"))
#aggr_flowframe <- fsApply(output$fcs, exprs)
#saveRDS(aggr_flowframe,"aggr_flowframe_14113M.rds")
#aggr_flowframe <- readRDS("aggr_flowframe_14113M.rds")
aggr_flowframe14113M <- fsApply(output$fcs, exprs)
rownames(aggr_flowframe14113M) <- c(1:nrow(aggr_flowframe14113M))
seurat_objb <- CreateSeuratObject(counts = t(aggr_flowframe14113M))
seurat_obj <- seurat_objb
seurat_obj@meta.data$celltype <- output$cell_clustering1m
seurat_obj@meta.data$sample_id <- output$sample_ids
seurat_obj@meta.data$timepoint <- output$timepoint
seurat_obj@meta.data$batch_id <- output$batch_id
seurat_obj@meta.data$condition <- output$condition
seurat_obj@meta.data$pubID <- output$pubID
seurat_obj@meta.data$response <- output$response
seurat_obj@meta.data$survival <- output$survival
seurat_obj <- seurat_obj[, sample(colnames(seurat_obj), size =140000, replace=F)]
saveRDS(seurat_obj,)

seurat_obj <- NormalizeData(seurat_obj)
seurat_obj <- ScaleData(seurat_obj)
seurat_obj <- FindVariableFeatures(seurat_obj)
saveRDS(seurat_obj,paste0(getwd(),"/cellxgene_withPubID/sobj14113Mpubid_14w_nUmap.rds"))
seurat_obj <- RunPCA(seurat_obj)
seurat_obj <- RunUMAP(object = seurat_obj, dims = 1:20)
saveRDS(seurat_obj,paste0(getwd(),"/cellxgene_withPubID/sobj14113Mpubid_14w_wUmap.rds"))


###14113 T STIM
setwd("E:/WJHLab/cytofAtlas/updated_data/J14113_PBMC/Tcell_stim/Tcell_stim")
output <- readRDS("backup_output_wjh.rds")
md <- read_xlsx("E:/WJHLab/cytofAtlas/updated_data/J14113_PBMC/Tcell_stim/Tcell_stim/Config/J14113_metadataTS_wPubID.xlsx")

mm1 <- match(output$sample_ids, md$sample_id)
batch_id <- md$batch[mm1]
condition <- md$condition[mm1]
timepoint <- md$timepoint[mm1]
pubID <- md$pubID[mm1]
response <- md$response[mm1]
survival <- md$survival[mm1]
output$timepoint <- timepoint
output$batch_id <- batch_id
output$condition <- condition
output$pubID <- pubID
output$response <- response
output$survival <- survival
output$meta_data <- md
saveRDS(output,paste0(getwd(),"/cellxgene_withPubID/output14113TS_wPubID.rds"))
output <- readRDS(paste0(getwd(),"/cellxgene_withPubID/output14113TS_wPubID.rds"))

aggr_flowframe14113TS <- fsApply(output$fcs, exprs)
rownames(aggr_flowframe14113TS) <- c(1:nrow(aggr_flowframe14113TS))
seurat_obj <- CreateSeuratObject(counts = t(aggr_flowframe14113TS))
seurat_obj@meta.data$celltype <- output$cell_clustering1m
seurat_obj@meta.data$sample_id <- output$sample_ids
seurat_obj@meta.data$timepoint <- output$timepoint
seurat_obj@meta.data$batch_id <- output$batch_id
seurat_obj@meta.data$condition <- output$condition
seurat_obj@meta.data$pubID <- output$pubID
seurat_obj@meta.data$response <- output$response
seurat_obj@meta.data$survival <- output$survival
seurat_obj <- seurat_obj[, sample(colnames(seurat_obj), size =140000, replace=F)]
seurat_obj <- NormalizeData(seurat_obj)
seurat_obj <- ScaleData(seurat_obj)
seurat_obj <- FindVariableFeatures(seurat_obj)
saveRDS(seurat_obj,paste0(getwd(),"/cellxgene_withPubID/sobj14113TS_full_subset140w.rds"))
sobj <- readRDS(paste0(getwd(),"/cellxgene_withPubID/sobj14113TS_full_subset140w.rds"))

seurat_obj <- RunPCA(seurat_obj)
seurat_obj <- RunUMAP(object = seurat_obj, dims = 1:20)
saveRDS(seurat_obj,paste0(getwd(),"/cellxgene_withPubID/sobj14113TSpubID_14w_wUmap.rds"))


###1790 M
rm(list=ls())
setwd("E:/WJHLab/cytofAtlas/updated_data/J1790_PBMC/Myeloid/Myeloid")
output <- readRDS("finalrun_revised/outputMyeloid_carter_fixed.rds")
md <- read_xlsx("E:/WJHLab/cytofAtlas/updated_data/J1790_PBMC/Myeloid/Myeloid/Config/J1790_metadataM_wPubID.xlsx")

mm1 <- match(output$sample_ids, md$sample_id)
batch_id <- md$batch[mm1]
condition <- md$condition[mm1]
timepoint <- md$timepoint[mm1]
pubID <- md$pubID[mm1]
response <- md$response[mm1]
survival <- md$survival[mm1]
output$timepoint <- timepoint
output$batch_id <- batch_id
output$condition <- condition
output$pubID <- pubID
output$response <- response
output$survival <- survival
saveRDS(output,paste0(getwd(),"/finalrun_revised/output1790M_wPubid.rds"))

aggr_flowframe1790M <- fsApply(output$fcs, exprs)
rownames(aggr_flowframe1790M) <- c(1:nrow(aggr_flowframe1790M))
seurat_obj <- CreateSeuratObject(counts = t(aggr_flowframe1790M))
seurat_obj@meta.data$celltype <- output$cell_clustering1m
seurat_obj@meta.data$sample_id <- output$sample_ids
seurat_obj@meta.data$timepoint <- output$timepoint
seurat_obj@meta.data$batch_id <- output$batch_id
seurat_obj@meta.data$condition <- output$condition
seurat_obj@meta.data$pubID <- output$pubID
seurat_obj@meta.data$response <- output$response
seurat_obj@meta.data$survival <- output$survival
seurat_obj <- seurat_obj[, sample(colnames(seurat_obj), size =40000, replace=F)]
seurat_obj <- NormalizeData(seurat_obj)
seurat_obj <- ScaleData(seurat_obj)
seurat_obj <- FindVariableFeatures(seurat_obj)
seurat_obj <- RunPCA(seurat_obj)
seurat_obj <- RunUMAP(object = seurat_obj, dims = 1:20)
saveRDS(seurat_obj,paste0(getwd(),"/finalrun_revised/sobj1790MpubID_4w_wUmap.rds"))
seurat_obj <- readRDS("E:/WJHLab/cytofAtlas/updated_data/J1790_PBMC/Myeloid/Myeloid/cellxgene_withPubID/sobj1790MpubID_4w_wUmap.rds")

###1790 T STIM
rm(list=ls())
setwd("E:/WJHLab/cytofAtlas/updated_data/J1790_PBMC/Tcell_stim/Tcell_stim/finalrun_revised")
output <- readRDS("outputTcellStim_carter_fixed.rds")
md <- read_xlsx("E:/WJHLab/cytofAtlas/updated_data/J1790_PBMC/Tcell_stim/Tcell_stim/Config/J1790_metadataM_wPubID.xlsx")

mm1 <- match(output$sample_ids, md$sample_id)
sample_id <- md$sample_id[mm1]
batch_id <- md$batch[mm1]
condition <- md$condition[mm1]
timepoint <- md$timepoint[mm1]
pubID <- md$pubID[mm1]
response <- md$response[mm1]
survival <- md$survival[mm1]
output$timepoint <- timepoint
output$batch_id <- batch_id
output$condition <- condition
output$pubID <- pubID
output$sample_id <- sample_id
output$response <- response
output$survival <- survival
saveRDS(output,"output1790TS_wPubid.rds")

aggr_flowframe1790TS <- fsApply(output$fcs, exprs)
rownames(aggr_flowframe1790TS) <- c(1:nrow(aggr_flowframe1790TS))
seurat_obj <- CreateSeuratObject(counts = t(aggr_flowframe1790TS))
seurat_obj@meta.data$celltype <- output$cell_clustering1m
seurat_obj@meta.data$sample_id <- output$sample_ids
seurat_obj@meta.data$timepoint <- output$timepoint
seurat_obj@meta.data$batch_id <- output$batch_id
seurat_obj@meta.data$condition <- output$condition
seurat_obj@meta.data$pubID <- output$pubID
seurat_obj@meta.data$response <- output$response
seurat_obj@meta.data$survival <- output$survival
seurat_obj <- seurat_obj[, sample(colnames(seurat_obj), size =40000, replace=F)]
seurat_obj <- NormalizeData(seurat_obj)
seurat_obj <- ScaleData(seurat_obj)
seurat_obj <- FindVariableFeatures(seurat_obj)
seurat_obj <- RunPCA(seurat_obj)
seurat_obj <- RunUMAP(object = seurat_obj, dims = 1:20)
saveRDS(seurat_obj,"sobj1790TSpubID_4w_wUmap.rds")

###J13108
rm(list=ls())
setwd("E:/WJHLab/cytofAtlas/updated_data/J13108/finalrun_revised")
output <- readRDS("output113108_wPubid.rds")
md <- read_xlsx("E:/WJHLab/cytofAtlas/updated_data/J13108/Config/metadata_wPubID.xlsx")

mm1 <- match(output$sample_ids, md$sample_id)
sample_id <- md$sample_id[mm1]
batch_id <- md$batch[mm1]
condition <- md$condition[mm1]
timepoint <- md$timepoint[mm1]
pubID <- md$Pub_ID[mm1]
output$timepoint <- timepoint
output$batch_id <- batch_id
output$condition <- condition
output$pubID <- pubID
output$sample_id <- sample_id
saveRDS(output,"output113108_wPubid.rds")

aggr_flowframe1790TS <- fsApply(output$fcs, exprs)
rownames(aggr_flowframe1790TS) <- c(1:nrow(aggr_flowframe1790TS))
seurat_obj <- CreateSeuratObject(counts = t(aggr_flowframe1790TS))
seurat_obj@meta.data$celltype <- output$cell_clustering1m
seurat_obj@meta.data$sample_id <- output$sample_ids
seurat_obj@meta.data$timepoint <- output$timepoint
seurat_obj@meta.data$batch_id <- output$batch_id
seurat_obj@meta.data$condition <- output$condition
seurat_obj@meta.data$pubID <- output$pubID
seurat_obj@meta.data$response <- output$response
seurat_obj@meta.data$survival <- output$survival
seurat_obj <- seurat_obj[, sample(colnames(seurat_obj), size =60000, replace=F)]
seurat_obj <- NormalizeData(seurat_obj)
seurat_obj <- ScaleData(seurat_obj)
seurat_obj <- FindVariableFeatures(seurat_obj)
seurat_obj <- RunPCA(seurat_obj)
seurat_obj <- RunUMAP(object = seurat_obj, dims = 1:20)
saveRDS(seurat_obj,"sobj13108pubID_6w_wUmap.rds")



wjh <- readRDS("E:/WJHLab/cytofAtlas/updated_data/J13108/backup_output.rds")


