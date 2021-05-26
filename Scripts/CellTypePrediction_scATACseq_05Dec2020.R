#Vikas Bansal 5 Dec 2020 #scATACseq
#Foundin project. Integrating with scRNA-seq data. Classify scATAC-seq cells based on an scRNA-seq experiment.

set.seed(786)
setwd("/data/vikas/FOUNDIN_scATAC/")


library(Signac)
library(Seurat)
library(GenomicRanges)
library(future)
library(harmony)
library(readxl)
library(stringi)
library(GenomicFeatures)
library(ballgown)
library(ggplot2)

plan("multiprocess", workers = 40)
options(future.globals.maxSize = 250000 * 1024^2) 

pathto.outPlots <- "/data/vikas/FOUNDIN_scATAC/OutputPlots/"
pathto.outData <- "/data/vikas/FOUNDIN_scATAC/OutputData/"
outName <- "iPSCsDopa_scATACseq"



hm.integrated <- readRDS("OutputData/iPSCsDopa_scATACseq_hm.integratedBeforeScaling.RDS")

# Load the pre-processed scRNA-seq data
iPSCsDopa.integrated <- readRDS("../FOUNDIN_scRNA/OutputData/iPSCsDopaALL/iPSCsDopaALL_integratedAfterBroadCellType.RDS")

iPSCsDopa.integratedSub <- subset(iPSCsDopa.integrated, downsample = 10000)

#saveRDS(iPSCsDopa.integratedSub,file = paste0(pathto.outData,"iPSCsDopa.integratedSub.RDS"))
iPSCsDopa.integratedSub <- readRDS("OutputData/iPSCsDopa.integratedSub.RDS")

comm_features <- rownames(iPSCsDopa.integratedSub)[(rownames(iPSCsDopa.integratedSub)%in%(rownames(hm.integrated)))]
transfer.anchors <- FindTransferAnchors(
  reference = iPSCsDopa.integratedSub,
  query = hm.integrated,
  reduction = 'cca', features = comm_features
)

#saveRDS(transfer.anchors,file = paste0(pathto.outData,outName,"_transfer.anchors.RDS"))
transfer.anchors <- readRDS("OutputData/iPSCsDopa_scATACseq_transfer.anchors.RDS")


##Harmony reduction for prediction
predicted.labels <- TransferData(
  anchorset = transfer.anchors,
  refdata = iPSCsDopa.integratedSub$CellType,
  weight.reduction = hm.integrated[['harmony']],
  dims = 2:30
)

hm.integrated <- AddMetaData(object = hm.integrated, metadata = predicted.labels)

Idents(hm.integrated) <- hm.integrated$predicted.id
hm.integrated@meta.data$predictedCellType <- hm.integrated@meta.data$predicted.id



png(paste0(pathto.outPlots,"UMAP_scATAC_harmony_predicted.id.png"), width=5000, height=5000, res = 300)
DimPlot(object = hm.integrated, label = TRUE) 
dev.off()

## Harmony reduction and broad cell type for prediction
iPSCsDopa.integratedSub@meta.data$BroadCellType <- Idents(iPSCsDopa.integratedSub)

predicted.labels <- TransferData(
  anchorset = transfer.anchors,
  refdata = iPSCsDopa.integratedSub$BroadCellType,
  weight.reduction = hm.integrated[['harmony']],
  dims = 2:30
)

hm.integrated <- AddMetaData(object = hm.integrated, metadata = predicted.labels)

Idents(hm.integrated) <- hm.integrated$predicted.id

hm.integrated@meta.data$predictedBroadCellType <- hm.integrated@meta.data$predicted.id


png(paste0(pathto.outPlots,"UMAP_scATAC_harmony_BroadCellType_predicted.id.png"), width=5000, height=5000, res = 300)
DimPlot(object = hm.integrated, label = TRUE) 
dev.off()



##Make broad cell type and cell type meta data
Idents(hm.integrated) <- hm.integrated$ATAC_snn_res.0.1
tail(sort(table(paste0(hm.integrated$ATAC_snn_res.0.1,";",hm.integrated$predicted.id))))

5;Dopaminergic Neurons                  2;Ependymal-like Cells
3346                                    3752
6;Late neuron Progenitor                2;Late neuron Progenitor
3820                                    4162
0;Neuroepithelial-like Cells                  4;Dopaminergic Neurons
6577                                    6806
1;Proliferating Floor Plate Progenitors 2;Proliferating Floor Plate Progenitors
7397                                    9459
3;Dopaminergic Neurons               1;Early neuron Progenitor
16808                                   17673
0;Immature Dopaminergic Neurons
41703


hm.integrated <- RenameIdents(object = hm.integrated, '0' = 'Immature Dopaminergic Neurons', '1' = 'Early neuron Progenitor', '2' = 'Proliferating Floor Plate Progenitors', '3' = 'Dopaminergic Neurons', '4' = 'Dopaminergic Neurons', '5' = 'Dopaminergic Neurons', '6' = 'Late neuron Progenitor', '7' = 'Immature Dopaminergic Neurons', '8' = 'Proliferating Floor Plate Progenitors')
hm.integrated@meta.data$BroadCellType <- Idents(hm.integrated)



png(paste0(pathto.outPlots,"UMAP_scATAC_harmony_BroadCellType.png"), width=5000, height=5000, res = 300)
DimPlot(object = hm.integrated, label = TRUE) 
dev.off()




##Make Neuroepithelial-like Cells in BroadCellType using resolution 0.2
hm.integrated@meta.data$BroadCellTypeV2 <- as.character(hm.integrated@meta.data$BroadCellType)

hm.integrated@meta.data[hm.integrated@meta.data$ATAC_snn_res.0.2 == 8,"BroadCellTypeV2"] <- "NE"

Idents(hm.integrated) <- hm.integrated$BroadCellTypeV2


saveRDS(hm.integratedBackup,file = paste0(pathto.outData,"iPSCsDopa_scATACseq_hm.integratedAfterBroadCellTypeV2.RDS"))





