library(matrixStats)
library(plyranges)
library(pbmcapply)
library(cowplot)
library(Matrix)
library(tidyverse)
library(magrittr)
library(data.table)
library(fastseg)
library(tibble)
library(PCAtools)

examp.betas <- fread("../METAFORA_output/Population_methylation.tissue_Blood/Population_methylation.tissue_Blood.chrom_chr1.1.betas.mat.gz")
examp.depths <- fread("../METAFORA_output/Population_methylation.tissue_Blood/Population_methylation.tissue_Blood.chrom_chr1.1.coverage.mat.gz")

# calculate population mean beta value across all samples 
beta.mat <- as.matrix(examp.betas[,4:ncol(examp.betas)])
depth.mat <- as.matrix(examp.depths[,4:ncol(examp.depths)])

pop_mean <- data.table(chrom=examp.betas$chromosome, start=examp.betas$start, end=examp.betas$end)
pop_mean$total_depth <- rowSums(depth.mat, na.rm=T)
pop_mean$mean_beta <- rowSums(beta.mat*depth.mat, na.rm=T)/(pop_mean$total_depth) 
pop_mean %<>% mutate(M=log(mean_beta/(1-mean_beta)))
summary(pop_mean$mean_beta)
#segment mean profile
cat('segmenting population mean profile . . . \n')
#segment mean profile
breakup_large_segments <- function(segs, threshold=200, smaller_seg_size=100) {
    smaller_segments <- Reduce(rbind, lapply(1:nrow(segs), function(i) { 
            if (segs$num.mark[i] < threshold) {return(segs[i,])}
            N = round(segs$num.mark[i] / smaller_seg_size)
            stepsize = round(segs$num.mark[i] / N)
            starts <- seq(segs$start[i],segs$end[i], stepsize)[1:N]
            ends <- c(starts[2:length(starts)] - 1, segs$end[i])
            tmp <- map_dfr(seq_len(N), ~segs[i,])
            tmp$start <- starts
            tmp$end <- ends
            tmp$num.mark <- tmp$end - tmp$start + 1
            return(tmp)
    }))
    return(smaller_segments)
}
default_segment_blocks <- function(pop_mean, chrom_block, alpha = 0.01, minSeg = 10) {
  index <- c(1,which(diff(pop_mean$start) > 1000))
  last_start <- index[length(index)]
  last_end <- nrow(pop_mean)
  block_start <- index[-length(index)]
  block_end <- index[-1] - 1
  block_start <- c(block_start, last_start)
  block_end <- c(block_end, last_end)
  blocks = data.frame(start=block_start, end=block_end)
  
  segments <- Reduce(rbind,lapply(1:nrow(blocks), function(i) { 
    block_beta <- pop_mean[blocks$start[i]:blocks$end[i],]
    segs <- as.data.frame(fastseg(block_beta$mean_beta, alpha=alpha, minSeg=minSeg, segMedianT=c(.65,.35)))
    #segs <- breakup_large_segments(segs)
    segs$start <- block_beta$start[segs$start]
    segs$end <- block_beta$end[segs$end]
    segs$width <- segs$end - segs$start
    segs$seqnames <- gsub(chrom_block, pattern="(\\w+)\\.\\d+",replacement="\\1")
    segs$seg_id <- paste0(chrom_block,"_",i,"_",1:nrow(segs)) 
	segs
    }))
  makeGRangesFromDataFrame(segments,keep.extra.columns = T)
}
default_segments=default_segment_blocks(pop_mean,"chr1.1",minSeg=20)
default_seg_n<-length(default_segments)
betas.gr <- makeGRangesFromDataFrame(examp.betas)
ol <- findOverlaps(default_segments, betas.gr)
# create Segment x CpG identity matrix to indicate which cpgs belong to which segment
CpG_Identity <- sparseMatrix(i = queryHits(ol), j = subjectHits(ol), dims=c(length(default_segments),nrow(examp.betas)), x=1)
beta.mat[is.na(beta.mat)] <- 0
depth.mat[is.na(depth.mat)] <- 0
tmp.beta.mat <- ((beta.mat*depth.mat) + 1) / (depth.mat + 2)
tmp.depth.mat <- depth.mat + 2

# aggregate betas across each segment
seg_beta <- as.matrix(((CpG_Identity %*% (tmp.beta.mat*tmp.depth.mat))) / ((CpG_Identity %*% tmp.depth.mat)))
seg_depth <- as.matrix(CpG_Identity %*% depth.mat / rowSums(CpG_Identity))
M <- log(seg_beta/(1-seg_beta))
mpcs <- PCAtools::pca(M)
default_pcs <- data.frame(mpcs$rotated[,1:5])
default_pcs

## test segmentation parameters
alphas_test <- c(0.01,.05,0.1,0.25)
segMeds_test <- c(0.5,1.0,1.5,2.0)
minSeg_test <- c(5,10,20,40)
params_to_test <- rbind( 
                    expand.grid(alpha=alphas_test, segMedianT=segMeds_test[2], minSeg=minSeg_test[3])%>%mutate(vary="alpha"),
                    expand.grid(alpha=alphas_test[2], segMedianT=segMeds_test, minSeg=minSeg_test[3])%>%mutate(vary="segMedT"),
                    expand.grid(alpha=alphas_test[2], segMedianT=segMeds_test[2], minSeg=minSeg_test)%>%mutate(vary="minSeg")
                )
segment_blocks <- function(pop_mean, chrom_block, alpha = 0.01, minSeg = 10, segMedianT=1) {
  index <- c(1,which(diff(pop_mean$start) > 1000))
  last_start <- index[length(index)]
  last_end <- nrow(pop_mean)
  block_start <- index[-length(index)]
  block_end <- index[-1] - 1
  block_start <- c(block_start, last_start)
  block_end <- c(block_end, last_end)
  blocks = data.frame(start=block_start, end=block_end)
  
  segments <- Reduce(rbind,lapply(1:nrow(blocks), function(i) { 
    block_beta <- pop_mean[blocks$start[i]:blocks$end[i],]
    segs <- as.data.frame(fastseg(block_beta$M, alpha=alpha, minSeg=minSeg, segMedianT=c(segMedianT,-segMedianT)))
    segs$start <- block_beta$start[segs$start]
    segs$end <- block_beta$end[segs$end]
    segs$width <- segs$end - segs$start
    segs$seqnames <- gsub(chrom_block, pattern="(\\w+)\\.\\d+",replacement="\\1")
    segs$seg_id <- paste0(chrom_block,"_",i,"_",1:nrow(segs)) 
	segs
    }))
  makeGRangesFromDataFrame(segments,keep.extra.columns = T)
}
segment_results <- pbmclapply(1:nrow(params_to_test), function(i) segment_blocks(pop_mean,"chr1.1",alpha=params_to_test[i,"alpha"],params_to_test[i,"minSeg"],segMedianT=params_to_test[i,"segMedianT"]), mc.cores=6)
seg_summary <- lapply(segment_results, function(x) x %>% as.data.frame %>% summarize(n=n(), mean_cpgs=mean(num.mark), median_length=median(width), variable_seg_n=sum(abs(seg.mean)<1,na.rm=T), variable_seg_length=median(width[which(abs(seg.mean)<1)]))) %>% bind_rows
seg_summary <- cbind(params_to_test,seg_summary)
seg_summary

alph.plt<-ggplot(seg_summary%>%filter(vary=="alpha"), aes(alpha,n,fill=mean_cpgs)) + geom_line()  + 
    geom_point(shape=21,alpha=.8,size=5) +  
    scale_fill_continuous(name="Mean CpGs per Seg", limits=c(80,260)) +
    geom_hline(yintercept=default_seg_n,linetype="dashed",color="red4") +
    theme_minimal() + theme(axis.line=element_line()) + ylim(0,NA)  +
    labs(y="Number of segments")
segMed.plt <- ggplot(seg_summary%>%filter(vary=="segMedT"), aes(segMedianT,n,fill=mean_cpgs)) + geom_line()  + 
    geom_point(shape=21,alpha=.8,size=5) +  
    scale_fill_continuous(name="Mean CpGs per Seg", limits=c(80,260)) +
    geom_hline(yintercept=default_seg_n,linetype="dashed",color="red4") +
    theme_minimal() + theme(axis.line=element_line()) + ylim(0,NA)  +
    labs(y="Number of segments")
minSeg.plt <- ggplot(seg_summary%>%filter(vary=="minSeg"), aes(minSeg,n,fill=mean_cpgs)) + geom_line()  + 
    geom_point(shape=21,alpha=.8,size=5) +  
    scale_fill_continuous(name="Mean CpGs per Seg", limits=c(80,260)) +
    geom_hline(yintercept=default_seg_n,linetype="dashed",color="red4") +
    theme_minimal() + theme(axis.line=element_line()) + ylim(0,NA) +
    labs(y="Number of segments",x="minSeg Size")
number_segment.plt <- plot_grid(alph.plt,segMed.plt,minSeg.plt,ncol=1)

segment_beta_PC_corr <- function(idx, pcs_test=5) {
    meth_segs<-segment_results[[idx]]
    ol <- findOverlaps(meth_segs, betas.gr)
    # create Segment x CpG identity matrix to indicate which cpgs belong to which segment
    CpG_Identity <- sparseMatrix(i = queryHits(ol), j = subjectHits(ol), dims=c(length(meth_segs),nrow(examp.betas)), x=1)
    # aggregate betas across each segment
    seg_beta <- as.matrix(((CpG_Identity %*% (tmp.beta.mat*tmp.depth.mat))) / ((CpG_Identity %*% tmp.depth.mat)))
    seg_depth <- as.matrix(CpG_Identity %*% depth.mat / rowSums(CpG_Identity))
    M <- log(seg_beta/(1-seg_beta))
    mpcs <- PCAtools::pca(M)
    PCcorrelations <- sapply(1:pcs_test, function(i) cor(mpcs$rotated[,i], default_pcs[,i]))
    data.frame(index=idx,params_to_test[idx,], PC=paste0("PC",1:pcs_test), cor=PCcorrelations)
}
correlations.df <- pbmclapply(1:length(segment_results), segment_beta_PC_corr) %>% bind_rows()

alph.plt <-ggplot(correlations.df%>%filter(segMedianT==1,minSeg==20), aes(factor(alpha),cor**2,fill=PC)) + geom_col(color="black",position="dodge",width=.5) +
    theme_minimal() + theme(axis.line=element_line()) + ylim(0,1) +
    scale_fill_brewer(palette="Set2") +
    labs(y = expression("PC concordance (" * R^2 * ")"), x="alpha")
segMed.plt <-ggplot(correlations.df%>%filter(alpha==0.05,minSeg==20), aes(factor(segMedianT),cor**2,fill=PC)) + geom_col(color="black",position="dodge",width=.5) +
    theme_minimal() + theme(axis.line=element_line()) + ylim(0,1) +
    scale_fill_brewer(palette="Set2") +
    labs(y = expression("PC concordance (" * R^2 * ")"),x="segMedianT")
minSeg.plt <-ggplot(correlations.df%>%filter(alpha==.05,segMedianT==1), aes(factor(minSeg),cor**2,fill=PC)) + geom_col(color="black",position="dodge",width=.5) +
    theme_minimal() + theme(axis.line=element_line()) + ylim(0,1) +
    scale_fill_brewer(palette="Set2") +
    labs(y = expression("PC concordance (" * R^2 * ")"),x="minSeg size")
pc_cor.plt <- plot_grid(alph.plt,segMed.plt,minSeg.plt,ncol=1)
final.plt <- plot_grid(number_segment.plt,pc_cor.plt,ncol=2)
ggsave("segmentation_parameter.sensitivity_analysis.plots.pdf", width=12)


