#Vikas Bansal 20 Nov 2020 #scATACseq
#Foundin project.

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

dirs <- list.dirs(path="/data/vikas/FOUNDIN_scATAC/CountsCopy_scATACseq/", recursive = F)

# read in peak sets
Test_peaks <- read.table(
  file = paste0(dirs[1],"/outs/peaks.bed"),
  col.names = c("chr", "start", "end")
)

All_joined_countMatrixTemp <- (Test_peaks)
combined.peaks <- makeGRangesFromDataFrame(All_joined_countMatrixTemp)

for (i in 2:length(dirs)){
  cat("loop", i, "\n")
  
  
  # read in peak sets
  Test_peaks <- read.table(
    file = paste0(dirs[i],"/outs/peaks.bed"),
    col.names = c("chr", "start", "end")
  )
  
  gr.Test_peaks <- makeGRangesFromDataFrame(Test_peaks)
  
  gr.All_joined_countMatrixTemp <- reduce(x = c(gr.Test_peaks, combined.peaks))
  
  combined.peaks <- gr.All_joined_countMatrixTemp

}

# Filter out bad peaks based on length
peakwidths <- width(combined.peaks)
combined.peaks <- combined.peaks[peakwidths  < 10000 & peakwidths > 20]
combined.peaks




# load metadata
metaDataTest <- read.table(
  file = paste0(dirs[1],"/outs/singlecell.csv"),
  stringsAsFactors = FALSE,
  sep = ",",
  header = TRUE,
  row.names = 1
)[-1, ] # remove the first row

metaDataTest.bar <- read.delim(
  file = paste0(dirs[1],"/outs/filtered_peak_bc_matrix/barcodes.tsv"),
  stringsAsFactors = FALSE, header=F)

metaDataTest <- metaDataTest[metaDataTest.bar[,1], ]

Sample_name <- unlist(strsplit(dirs[1],"//"))[[2]]

metaDataTest$SampleID <- Sample_name


# create fragment objects
frags.test <- CreateFragmentObject(
  path = paste0(dirs[1],"/outs/fragments.tsv.gz"),
  cells = rownames(metaDataTest)
)

test.counts <- FeatureMatrix(
  fragments = frags.test,
  features = combined.peaks,
  cells = rownames(metaDataTest),process_n=20000)

colnames(test.counts) <- paste0(colnames(test.counts),"_",Sample_name)

rownames(metaDataTest) <- paste0(rownames(metaDataTest),"_",Sample_name)

ALL_Frag <- test.counts
ALL_metaData <- metaDataTest

for (i in 2:length(dirs)){
  cat("loop", i, "\n")
  
  
  # load metadata
  metaDataTest <- read.table(
    file = paste0(dirs[i],"/outs/singlecell.csv"),
    stringsAsFactors = FALSE,
    sep = ",",
    header = TRUE,
    row.names = 1
  )[-1, ] # remove the first row
  
  metaDataTest.bar <- read.delim(
    file = paste0(dirs[i],"/outs/filtered_peak_bc_matrix/barcodes.tsv"),
    stringsAsFactors = FALSE, header=F)
  
  metaDataTest <- metaDataTest[metaDataTest.bar[,1], ]
  
  Sample_name <- unlist(strsplit(dirs[i],"//"))[[2]]
  
  metaDataTest$SampleID <- Sample_name
  
  
  # create fragment objects
  frags.test <- CreateFragmentObject(
    path = paste0(dirs[i],"/outs/fragments.tsv.gz"),
    cells = rownames(metaDataTest)
  )
  
  test.counts <- FeatureMatrix(
    fragments = frags.test,
    features = combined.peaks,
    cells = rownames(metaDataTest),process_n=20000)
  
  colnames(test.counts) <- paste0(colnames(test.counts),"_",Sample_name)
  
  rownames(metaDataTest) <- paste0(rownames(metaDataTest),"_",Sample_name)
  
  #if(all.equal(rownames(test.counts),rownames(ALL_Frag))){
    Temp_FragAll <- cbind(test.counts,ALL_Frag)
    ALL_Frag <- Temp_FragAll
  #}else {
  #  stop()
  #}
  
  
  Temp_metaAll <- rbind(metaDataTest,ALL_metaData)
  ALL_metaData <- Temp_metaAll
  
}

saveRDS(ALL_Frag,file = paste0(pathto.outData,outName,"_ALL_Frag.RDS"))
saveRDS(ALL_metaData,file = paste0(pathto.outData,outName,"_ALL_metaData.RDS"))

#ALL_Frag <- readRDS("OutputData/iPSCsDopa_scATACseq_ALL_Frag.RDS")
#ALL_metaData <- readRDS("OutputData/iPSCsDopa_scATACseq_ALL_metaData.RDS")



#CreateChromatinAssay and filter min cell 100 and features 1000 as in scRNA-seq
ALL_Frag_assay <- CreateChromatinAssay(ALL_Frag, min.cells = 100,
                                       min.features = 1000)
ALL_Frag_assay



ALL_Frag_seurat <- CreateSeuratObject(ALL_Frag_assay, assay = "ATAC", meta.data = ALL_metaData)
ALL_Frag_seurat[["ATAC"]]



# compute TSS % per cell
ALL_Frag_seurat$pct.TSSfrag <- ALL_Frag_seurat$TSS_fragments/ALL_Frag_seurat$passed_filters * 100


ALL_Frag_seurat$pct_reads_in_peaks <- ALL_Frag_seurat$peak_region_fragments / ALL_Frag_seurat$passed_filters * 100

ALL_Frag_seurat$pct_mito <- ALL_Frag_seurat$mitochondrial / ALL_Frag_seurat$total * 100

png(paste0(pathto.outPlots,"nCountPercentMT_Vln.png"), width=2000, height=1000)
VlnPlot(
  object = ALL_Frag_seurat,
  features = c('pct_reads_in_peaks', 'peak_region_fragments',
               'pct.TSSfrag', 'pct_mito'),
  pt.size = 0.1,
  ncol = 2
)
dev.off()


#filter low quality cells
ALL_Frag_seurat <- subset(
  x = ALL_Frag_seurat,
  subset = peak_region_fragments > 3000 &
    peak_region_fragments < 50000 &
    pct_reads_in_peaks > 20 &
    pct.TSSfrag > 20 &
    pct_mito < 20
)
ALL_Frag_seurat

png(paste0(pathto.outPlots,"nCountPercentMT_VlnAfterFilt.png"), width=2000, height=1000)
VlnPlot(
  object = ALL_Frag_seurat,
  features = c('pct_reads_in_peaks', 'peak_region_fragments',
               'pct.TSSfrag', 'pct_mito'),
  pt.size = 0.1,
  ncol = 2
)
dev.off()



ALL_Frag_seurat <- RunTFIDF(ALL_Frag_seurat)
ALL_Frag_seurat <- FindTopFeatures(ALL_Frag_seurat, min.cutoff = 'q0')
ALL_Frag_seurat <- RunSVD(ALL_Frag_seurat)
ALL_Frag_seurat <- RunUMAP(ALL_Frag_seurat, dims = 2:50, reduction = 'lsi')
p1 <- DimPlot(ALL_Frag_seurat, group.by = 'SampleID', pt.size = 0.1)

png(paste0(pathto.outPlots,"UMAP_scATAC_unintegrated.png"), width=2000, height=1000)
DimPlot(ALL_Frag_seurat, group.by = 'SampleID', pt.size = 0.1) 
dev.off()

#integration using Harmony
hm.integrated <- RunHarmony(
  object = ALL_Frag_seurat,
  group.by.vars = 'SampleID',
  reduction = 'lsi',
  assay.use = 'ATAC',
  project.dim = FALSE
)

# re-compute the UMAP using corrected LSI embeddings
hm.integrated <- RunUMAP(hm.integrated, dims = 2:30, reduction = 'harmony')
p2 <- DimPlot(hm.integrated, group.by = 'SampleID', pt.size = 0.1) 

png(paste0(pathto.outPlots,"UMAP_scATAC.png"), width=2000, height=1000)
CombinePlots(plots = list(p1, p2))
dev.off()

png(paste0(pathto.outPlots,"UMAP_scATAC_Harmony.png"), width=2000, height=1000)
DimPlot(hm.integrated, group.by = 'SampleID', pt.size = 0.1) 
dev.off()

#saveRDS(hm.integrated,file = paste0(pathto.outData,outName,"_hm.integrated.RDS"))

hm.integrated <- FindNeighbors(object = hm.integrated, reduction = 'harmony', dims = 2:30)
hm.integrated <- FindClusters(object = hm.integrated, verbose = TRUE, algorithm = 1,resolution = 0.2)

png(paste0(pathto.outPlots,"UMAP_scATAC_Harmony_ClustRes0.2.png"), width=2000, height=1000)
DimPlot(object = hm.integrated, label = TRUE) 
dev.off()

hm.integrated <- FindClusters(object = hm.integrated, verbose = TRUE, algorithm = 1,resolution = 0.1)

png(paste0(pathto.outPlots,"UMAP_scATAC_Harmony_ClustRes0.1.png"), width=2000, height=1000)
DimPlot(object = hm.integrated, label = TRUE) 
dev.off()


saveRDS(hm.integrated,file = paste0(pathto.outData,outName,"_hm.integratedAlgo1.RDS"))

#hm.integrated <- readRDS("OutputData/iPSCsDopa_scATACseq_hm.integratedAlgo1.RDS")

##Include in-house meta info
Petertable <- read_excel("/data/vikas/FOUNDIN_scRNA/MetaInfo/Samples_names_codes_final_vikas.xlsx", sheet = 1)

Petertablev2 <- (as.data.frame(Petertable))
Petertablev3 <- Petertablev2[grep("CDI",Petertablev2[,"Barcode_DZNE"]),]


Petertablev3$Barcode_last4 <- (stri_sub(Petertablev3$Barcode_DZNE,-4,-1))

Petertablev3$SampleID <- (paste0("SCAT_PPMI",Petertablev3$PPMI_ID, "_", Petertablev3$Barcode_last4,"_da65"))

Petertablev4 <- Petertablev3
Petertablev4[c(nrow(Petertablev3)+1,nrow(Petertablev3)+2,nrow(Petertablev3)+3),] <- Petertablev3[grep("SCAT_PPMI3966_2813_da65",Petertablev3$SampleID),]

Petertablev4[c(nrow(Petertablev3)+1,nrow(Petertablev3)+2,nrow(Petertablev3)+3),"SampleID"] <- c('SCAT_PPMI3966B3_2813_da65', 'SCAT_PPMI3966E6_2813_da65', 'SCAT_PPMI3966E8_2813_da65')
Petertablev4[grep("SCAT_PPMI3966B3_2813_da65",Petertablev4$SampleID),"BATCH"] <- 3
Petertablev4[grep("SCAT_PPMI3966E6_2813_da65",Petertablev4$SampleID),"BATCH"] <- 5
Petertablev4[grep("SCAT_PPMI3966E8_2813_da65",Petertablev4$SampleID),"BATCH"] <- 5

PetertableToPut <- Petertablev4[(match(hm.integrated@meta.data$SampleID,Petertablev4$SampleID)),]



###******************************* Change the number of columns accordingly in future***************
hm.integrated@meta.data[,28:41] <- PetertableToPut[,1:14]

###*************************************************************************************************

##################################Add annotation gene wise counts####################################################################

gtf_file <- "/data/vikas/FOUNDIN_scRNA/MetaInfo/genes.gtf"

# read in FOUNDIN GTF as GRanges:
annotations = gffReadGR(gtf_file)
colnames(annotations@elementMetadata)[6] <- "gene_biotype"


Idents(hm.integrated) <- hm.integrated$ATAC_snn_res.0.1
Annotation(hm.integrated) <- annotations



extend <- function(x, upstream=0, downstream=0)     
{
  if (any(strand(x) == "*"))
    warning("'*' ranges were treated as '+'")
  on_plus <- strand(x) == "+" | strand(x) == "*"
  new_start <- start(x) - ifelse(on_plus, upstream, downstream)
  new_end <- end(x) + ifelse(on_plus, downstream, upstream)
  ranges(x) <- IRanges(new_start, new_end)
  trim(x)
}

annotationsExtend <- extend(annotations, 2000, 0)



combined.peaks <-  (annotationsExtend[annotationsExtend$type == "gene"])



# load metadata
metaDataTest <- read.table(
  file = paste0(dirs[1],"/outs/singlecell.csv"),
  stringsAsFactors = FALSE,
  sep = ",",
  header = TRUE,
  row.names = 1
)[-1, ] # remove the first row

metaDataTest.bar <- read.delim(
  file = paste0(dirs[1],"/outs/filtered_peak_bc_matrix/barcodes.tsv"),
  stringsAsFactors = FALSE, header=F)

metaDataTest <- metaDataTest[metaDataTest.bar[,1], ]

Sample_name <- unlist(strsplit(dirs[1],"//"))[[2]]

metaDataTest$SampleID <- Sample_name


# create fragment objects
frags.test <- CreateFragmentObject(
  path = paste0(dirs[1],"/outs/fragments.tsv.gz"),
  cells = rownames(metaDataTest)
)

test.counts <- FeatureMatrix(
  fragments = frags.test,
  features = combined.peaks,
  cells = rownames(metaDataTest),process_n=2500)

colnames(test.counts) <- paste0(colnames(test.counts),"_",Sample_name)

rownames(metaDataTest) <- paste0(rownames(metaDataTest),"_",Sample_name)

ALL_Frag <- test.counts
ALL_metaData <- metaDataTest

for (i in 2:length(dirs)){
  cat("loop", i, "\n")
  
  
  # load metadata
  metaDataTest <- read.table(
    file = paste0(dirs[i],"/outs/singlecell.csv"),
    stringsAsFactors = FALSE,
    sep = ",",
    header = TRUE,
    row.names = 1
  )[-1, ] # remove the first row
  
  metaDataTest.bar <- read.delim(
    file = paste0(dirs[i],"/outs/filtered_peak_bc_matrix/barcodes.tsv"),
    stringsAsFactors = FALSE, header=F)
  
  metaDataTest <- metaDataTest[metaDataTest.bar[,1], ]
  
  Sample_name <- unlist(strsplit(dirs[i],"//"))[[2]]
  
  metaDataTest$SampleID <- Sample_name
  
  
  # create fragment objects
  frags.test <- CreateFragmentObject(
    path = paste0(dirs[i],"/outs/fragments.tsv.gz"),
    cells = rownames(metaDataTest)
  )
  
  test.counts <- FeatureMatrix(
    fragments = frags.test,
    features = combined.peaks,
    cells = rownames(metaDataTest),process_n=2500)
  
  colnames(test.counts) <- paste0(colnames(test.counts),"_",Sample_name)
  
  rownames(metaDataTest) <- paste0(rownames(metaDataTest),"_",Sample_name)
  
  #if(all.equal(rownames(test.counts),rownames(ALL_Frag))){
  Temp_FragAll <- cbind(test.counts,ALL_Frag)
  ALL_Frag <- Temp_FragAll
  #}else {
  #  stop()
  #}
  
  Temp_metaAll <- rbind(metaDataTest,ALL_metaData)
  ALL_metaData <- Temp_metaAll
  
}


saveRDS(ALL_Frag,file = paste0(pathto.outData,outName,"_Gene_Frag.RDS"))

#ALL_Frag <- readRDS("OutputData/iPSCsDopa_scATACseq_Gene_Frag.RDS")

jikl <- (paste0(combined.peaks@seqnames,"-",combined.peaks@ranges,";",combined.peaks$gene_name))

jikl2 <- (jikl[-(grep("chrM",jikl))])
jikl3 <- (sapply(strsplit(jikl2,";"),"[[",2))
rownames(ALL_Frag) <- jikl3


Idents(hm.integrated) <- hm.integrated$ATAC_snn_res.0.1

##Create seurat object using gene/RNA counts
Gene_Frag_seurat <- CreateSeuratObject(ALL_Frag, assay = "RNA", min.cells = 100)
Gene_Frag_seuratV2 <- subset(Gene_Frag_seurat, cells=colnames(hm.integrated))

saveRDS(Gene_Frag_seuratV2,file = paste0(pathto.outData,outName,"_Gene_Frag_seuratV2.RDS"))


Gene_frag_filt <- GetAssayData(Gene_Frag_seuratV2,slot="counts", assay="RNA")


#Add RNA counts as RNA asay in harmony integrated object
hm.integrated[['RNA']]  <- CreateAssayObject(counts = Gene_frag_filt)

hm.integrated <- NormalizeData(
  object = hm.integrated,
  assay = 'RNA',
  normalization.method = 'LogNormalize',
  scale.factor = median(hm.integrated$nCount_RNA)
)

DefaultAssay(hm.integrated) <- 'RNA'

saveRDS(hm.integrated,file = paste0(pathto.outData,outName,"_hm.integratedBeforeScaling.RDS"))



hm.integrated <- ScaleData(hm.integrated)


saveRDS(hm.integrated,file = paste0(pathto.outData,outName,"_hm.integratedFINAL.RDS"))

sessionInfo()



R version 3.6.0 (2019-04-26)
Platform: x86_64-pc-linux-gnu (64-bit)
Running under: Debian GNU/Linux 9 (stretch)

Matrix products: default
BLAS/LAPACK: /usr/lib/libopenblasp-r0.2.19.so

locale:
 [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8     LC_MONETARY=en_US.UTF-8    LC_MESSAGES=C             
 [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                  LC_ADDRESS=C               LC_TELEPHONE=C             LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       

attached base packages:
[1] parallel  stats4    stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] ggplot2_3.3.2          ballgown_2.5.3         GenomicFeatures_1.36.4 AnnotationDbi_1.46.0   Biobase_2.44.0         stringi_1.4.3          readxl_1.3.1          
 [8] harmony_1.0            Rcpp_1.0.1             future_1.14.0          GenomicRanges_1.36.0   GenomeInfoDb_1.20.0    IRanges_2.18.1         S4Vectors_0.22.1      
[15] BiocGenerics_0.30.0    Seurat_3.2.0.9015      Signac_1.1.0          

loaded via a namespace (and not attached):
  [1] reticulate_1.12             tidyselect_1.1.0            RSQLite_2.1.1               htmlwidgets_1.3             grid_3.6.0                 
  [6] BiocParallel_1.18.0         Rtsne_0.15                  munsell_0.5.0               codetools_0.2-16            ica_1.0-2                  
 [11] miniUI_0.1.1.1              withr_2.1.2                 colorspace_1.4-1            OrganismDbi_1.26.0          knitr_1.23                 
 [16] rstudioapi_0.10             ROCR_1.0-7                  tensor_1.5                  listenv_0.7.0               GenomeInfoDbData_1.2.1     
 [21] polyclip_1.10-0             bit64_0.9-7                 farver_2.0.3                vctrs_0.3.5                 generics_0.0.2             
 [26] xfun_0.8                    biovizBase_1.32.0           lsa_0.73.2                  ggseqlogo_0.1               R6_2.4.0                   
 [31] rsvd_1.0.1                  AnnotationFilter_1.8.0      bitops_1.0-6                spatstat.utils_1.13-0       reshape_0.8.8              
 [36] DelayedArray_0.10.0         assertthat_0.2.1            promises_1.0.1              scales_1.0.0                nnet_7.3-12                
 [41] gtable_0.3.0                npsurv_0.4-0                sva_3.32.1                  globals_0.12.4              goftest_1.1-1              
 [46] ggbio_1.32.0                ensembldb_2.8.0             rlang_0.4.9                 genefilter_1.66.0           RcppRoll_0.3.0             
 [51] splines_3.6.0               rtracklayer_1.44.0          lazyeval_0.2.2              acepack_1.4.1               dichromat_2.0-0            
 [56] checkmate_1.9.4             BiocManager_1.30.4          reshape2_1.4.3              abind_1.4-5                 backports_1.1.4            
 [61] httpuv_1.5.1                Hmisc_4.2-0                 RBGL_1.60.0                 tools_3.6.0                 ellipsis_0.2.0.1           
 [66] gplots_3.0.4                RColorBrewer_1.1-2          ggridges_0.5.1              plyr_1.8.6                  base64enc_0.1-3            
 [71] progress_1.2.2              zlibbioc_1.30.0             purrr_0.3.4                 RCurl_1.95-4.12             prettyunits_1.0.2          
 [76] rpart_4.1-15                deldir_0.1-22               pbapply_1.4-0               cowplot_1.0.0               zoo_1.8-6                  
 [81] SummarizedExperiment_1.14.1 ggrepel_0.8.1               cluster_2.0.8               magrittr_1.5                data.table_1.12.8          
 [86] lmtest_0.9-37               RANN_2.6.1                  SnowballC_0.7.0             ProtGenerics_1.16.0         fitdistrplus_1.0-14        
 [91] matrixStats_0.56.0          hms_0.5.3                   patchwork_1.0.1.9000        lsei_1.2-0                  mime_0.7                   
 [96] xtable_1.8-4                XML_3.98-1.20               gridExtra_2.3               compiler_3.6.0              biomaRt_2.40.3             
[101] tibble_3.0.1                KernSmooth_2.23-15          crayon_1.3.4                htmltools_0.3.6             mgcv_1.8-28                
[106] later_0.8.0                 Formula_1.2-3               tidyr_1.1.0                 DBI_1.0.0                   tweenr_1.0.1               
[111] MASS_7.3-51.4               Matrix_1.2-17               gdata_2.18.0                igraph_1.2.4.1              pkgconfig_2.0.2            
[116] GenomicAlignments_1.20.1    foreign_0.8-71              plotly_4.9.0                annotate_1.62.0             XVector_0.24.0             
[121] stringr_1.4.0               VariantAnnotation_1.30.1    digest_0.6.20               sctransform_0.2.0           RcppAnnoy_0.0.12           
[126] graph_1.62.0                spatstat.data_1.4-0         Biostrings_2.52.0           cellranger_1.1.0            leiden_0.3.3               
[131] fastmatch_1.1-0             htmlTable_1.13.1            uwot_0.1.8                  curl_3.3                    shiny_1.3.2                
[136] Rsamtools_2.0.0             gtools_3.8.1                lifecycle_0.2.0             nlme_3.1-139                jsonlite_1.6               
[141] limma_3.40.6                viridisLite_0.3.0           BSgenome_1.52.0             pillar_1.4.4                lattice_0.20-38            
[146] GGally_1.4.0                httr_1.4.1                  survival_2.44-1.1           glue_1.4.2                  spatstat_1.60-1            
[151] png_0.1-7                   bit_1.1-14                  ggforce_0.3.1               blob_1.2.1                  latticeExtra_0.6-28        
[156] caTools_1.17.1.2            memoise_1.1.0               dplyr_1.0.2                 irlba_2.3.3                 future.apply_1.3.0         
[161] ape_5.3                    










