#!/usr/bin/Rscript

## damID normalization and detection of significant fusion regions
suppressPackageStartupMessages(library(HMMt))
suppressPackageStartupMessages(library(GenomicRanges))

#scripts='/Users/dinazielinski/Dropbox/curie/heard/analysis/LAD/scripts/'

## total bam counts
inPath='~/Dropbox/curie/heard/analysis/LAD/damid/results/counts_10kb/'

## merge reps and call LADs
outPath=paste0(inPath,"sLADs_merged/")
dir.create(outPath, recursive = TRUE)
## per replicate sLADs:
outPath.2=paste0(inPath,"sLADs/")
dir.create(outPath.2, recursive = TRUE)

## load supporting functions
source(paste0(scripts, "damid_hmmt_functions.R"))

## min number of combined fusion+dam reads
counts.in<-list.files(path=inPath, pattern=".txt", full.names=TRUE, recursive=FALSE)
pseudo=1

## call LADs on merged replicates
lapply(counts.in, function(file) {
  sampID=gsub("_merged_replicates.txt", "", basename(file))
  countTable <- read.table(file, header=TRUE, sep="\t")
  
  countTable<-countTable[countTable$FUSION_1 + countTable$DAM_1 > min_reads & countTable$FUSION_2 + countTable$DAM_2 > min_reads,]
  ## normalization
  ## replicate 1
  countTable$FUSION_NORM_1<-countTable$FUSION_1/sum(countTable$FUSION_1) * 1e6 + pseudo
  countTable$DAM_NORM_1<-countTable$DAM_1/sum(countTable$DAM_1) * 1e6 + pseudo
  countTable$logNORM_1<-log2(countTable$FUSION_NORM_1/countTable$DAM_NORM_1)
  
  ## replicate 2
  countTable$FUSION_NORM_2<-countTable$FUSION_2/sum(countTable$FUSION_2) * 1e6 + pseudo
  countTable$DAM_NORM_2<-countTable$DAM_2/sum(countTable$DAM_2) * 1e6 + pseudo
  countTable$logNORM_2<-log2(countTable$FUSION_NORM_2/countTable$DAM_NORM_2)
  
  ## avg log2(fusion/dam)
  countTable$avgRATIO=(countTable$logNORM_1 + countTable$logNORM_2) / 2
  
  normData <- countTable[,c("CHR","START","END","avgRATIO")]
  colnames(normData) = c("chr", "start", "end", "score")
  normData = normData[order(normData$chr, normData$start),]
  
  output_hmm <- file.path(outPath, paste0(sampID, "_HMM.bed"))
  output_ad <- file.path(outPath, paste0(sampID, "_AD.bed"))
  
  ## run HMM on combined replicates and select "associated domains"
  hmm_calls <- HMM(normData, na_solution = "NA")
  hmm_ranges <- getHMMRanges(as(hmm_calls, "GRanges"), score = "AD")
  
  ## Write data tables
  write.table(hmm_calls, file = output_hmm, 
              quote = F, sep = "\t", row.names = F, col.names = F)
  
  write.table(as(hmm_ranges, "data.frame")[, 1:3],
              file = output_ad, 
              quote = F, sep = "\t", row.names = F, col.names = F)
})


## call LADs on ind. replicates
lapply(counts.in, function(file) {
  #sampID=gsub("_merged_replicates.txt", "", basename(file))
  sampID=gsub("_merged_counts.txt", "", basename(file))
  countTable <- read.table(file, header=TRUE, sep="\t")
  
  #countTable<-countTable[countTable$FUSION_1 + countTable$DAM_1 > min_reads & countTable$FUSION_2 + countTable$DAM_2 > min_reads,]
  countTable<-countTable[countTable$FUSION_1 + countTable$DAM_1 > min_reads,]
  
  ## replicate 1
  ## normalization
  countTable$FUSION_NORM_1<-countTable$FUSION_1/sum(countTable$FUSION_1) * 1e6 + pseudo
  countTable$DAM_NORM_1<-countTable$DAM_1/sum(countTable$DAM_1) * 1e6 + pseudo
  countTable$logNORM_1<-log2(countTable$FUSION_NORM_1/countTable$DAM_NORM_1)
  
  ## call LADs
  normData.1 <- countTable[,c("CHR","START","END","logNORM_1")]
  colnames(normData.1) = c("chr", "start", "end", "score")
  normData.1 = normData.1[order(normData.1$chr, normData.1$start),]
  
  output_hmm.1 <- file.path(outPath.2, paste0(sampID, "_1_HMM.bed"))
  output_ad.1 <- file.path(outPath.2, paste0(sampID, "_1_AD.bed"))
  
  ## run HMM
  hmm_calls.1 <- HMM(normData.1, na_solution = "NA")
  hmm_ranges.1 <- getHMMRanges(as(hmm_calls.1, "GRanges"), score = "AD")
  
  ## Write data tables
  write.table(hmm_calls.1, file = output_hmm.1, 
              quote = F, sep = "\t", row.names = F, col.names = F)
  
  write.table(as(hmm_ranges.1, "data.frame")[, 1:3],
              file = output_ad.1, 
              quote = F, sep = "\t", row.names = F, col.names = F)
  ## replicate 2
  ## normalization
  countTable$FUSION_NORM_2<-countTable$FUSION_2/sum(countTable$FUSION_2) * 1e6 + pseudo
  countTable$DAM_NORM_2<-countTable$DAM_2/sum(countTable$DAM_2) * 1e6 + pseudo
  countTable$logNORM_2<-log2(countTable$FUSION_NORM_2/countTable$DAM_NORM_2)
  ## call LADs
  normData.2 <- countTable[,c("CHR","START","END","logNORM_2")]
  colnames(normData.2) = c("chr", "start", "end", "score")
  normData.2 = normData.2[order(normData.2$chr, normData.2$start),]
  
  output_hmm.2 <- file.path(outPath.2, paste0(sampID, "_2_HMM.bed"))
  output_ad.2 <- file.path(outPath.2, paste0(sampID, "_2_AD.bed"))
  
  ## run HMM
  hmm_calls.2 <- HMM(normData.2, na_solution = "NA")
  hmm_ranges.2 <- getHMMRanges(as(hmm_calls.2, "GRanges"), score = "AD")
  
  ## Write data tables
  write.table(hmm_calls.2, file = output_hmm.2, quote = F, sep = "\t", row.names = F, col.names = F)
  
  write.table(as(hmm_ranges.2, "data.frame")[, 1:3], file = output_ad.2, quote = F, sep = "\t", row.names = F, col.names = F)
})

