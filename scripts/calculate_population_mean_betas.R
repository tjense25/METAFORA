library(argparser)
library(Matrix)
library(tidyverse)
library(magrittr)
library(data.table)
library(plyranges)
library(pbmcapply)
library(fastseg)
library(GenomicRanges)
library(tibble)
library(bedr)
 

parser <- arg_parser("Aggregate individual beta profiles to calculate population average")
parser <- add_argument(parser, "--filelist", help="new line separated list of beta profiles file paths")
parser <- add_argument(parser, "--chrom", help="chromosome to compute betas from")
parser <- add_argument(parser, "--cpgs", help="bed file of all reference genome cpgs (only on autosomes)")
parser <- add_argument(parser, "--min_segment_cpgs", help="minimum number of cpgs in segment for segmentation", type="integer", default=10)
parser <- add_argument(parser, "--beta_mat", help="where to write beta matrix output")
parser <- add_argument(parser, "--depth_mat", help="where to write depth matrix output")
parser <- add_argument(parser, "--pop_mean", help="where to write population mean beta")
parser <- add_argument(parser, "--segment_beta", help="where to write mean porfile segments")
parser <- add_argument(parser, "--segment_depth", help="where to write average depth over segments")
argv <- parse_args(parser)

this_chrom <- argv$chrom
methylation.dat <- fread(argv$cpgs)
colnames(methylation.dat) <- c("chromosome", "start","end")
methylation.dat <- methylation.dat[methylation.dat$chromosome == this_chrom,]
max_start<- max(methylation.dat$start) + 1
this_region=paste0(this_chrom,":",1,"-",max_start)
min_cpg_number <- argv$min_segment_cpgs
this_region

#read in files and merge depth and betas to cpg data.table 
file_list=read.table(argv$filelist,col.names="file")$file
for (f in file_list) {
    sample_id = str_split(basename(f), pattern='\\.',simplify =T)[1]
    cat(paste0('reading methylation profile for sample ', sample_id, ' from file ', f,'\n'))
    tmp.cpg <- tabix(region=this_region, f, verbose=F)
    colnames(tmp.cpg) <- c("chromosome","start","end","depth","beta")
    tmp.cpg$start %<>% as.integer()
    tmp.cpg$depth %<>% as.integer()
    tmp.cpg$beta %<>% as.numeric()
    cols=paste0(sample_id,c(".depth", ".beta"))
    methylation.dat[tmp.cpg, on=.(chromosome,start), (cols) := mget(c("i.depth", "i.beta"))]
}
cat("\nCalculating mean beta and total depth across all samples\n\n")
# create matrix of depths per cpg per sample
depth_cols <- c(1,2,3,seq(4,ncol(methylation.dat),2))
depth <- methylation.dat[,..depth_cols]
names(depth) <- gsub(names(depth), pattern="(\\w+)\\.?.*", replacement="\\1")

# create matrix of beta values per cpg per sample
beta_cols <- c(1,2,3,seq(5,ncol(methylation.dat),2))
beta <- methylation.dat[,..beta_cols]
names(beta) <- gsub(names(beta), pattern="(\\w+)\\.?.*", replacement="\\1")

depth.mat <- as.matrix(depth[,4:ncol(depth)])
beta.mat <- as.matrix(beta[,4:ncol(beta)])

#subset to cpgs that have coverage across samples
cat('subsetting to cpgs with reasonable coverage across most samples . . .\n')
median_depth <- apply(depth.mat,1,function(x) median(x,na.rm=T))
keep_cpgs <- (!is.na(median_depth)) & median_depth >= 5
depth.mat <- depth.mat[keep_cpgs,]
depth <- depth[keep_cpgs,]
beta.mat <- beta.mat[keep_cpgs,]
beta <- beta[keep_cpgs,]
methylation.dat <- methylation.dat[keep_cpgs, ]

# calculate population mean beta value across all samples 
cat('calculting population mean beta . . .\n')
pop_mean <- data.table(chrom=methylation.dat$chrom, start=methylation.dat$start, end=methylation.dat$end)
pop_mean$total_depth <- rowSums(depth.mat, na.rm=T)
pop_mean$mean_beta <- rowSums(beta.mat*depth.mat, na.rm=T)/(pop_mean$total_depth) 

#segment mean profile
cat('segmenting population mean profile . . . \n')
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
    } ))
    return(smaller_segments)
}

segment_blocks <- function(pop_mean, chrom, alpha = 0.01, minSeg = 10) {
  index <- c(1,which(diff(pop_mean$start) > 1000))
  last_start <- index[length(index)]
  last_end <- nrow(pop_mean)
  block_start <- index[-length(index)]
  block_end <- index[-1] - 1
  block_start <- c(block_start, last_start)
  block_end <- c(block_end, last_end)
  blocks = data.frame(start=block_start, end=block_end)
  blocks <- blocks[(blocks$end - blocks$start) > minSeg,]
  
  segments <- Reduce(rbind,pbmclapply(1:nrow(blocks), function(i) { 
    block_beta <- pop_mean[blocks$start[i]:blocks$end[i],]
    segs <- as.data.frame(fastseg(block_beta$mean_beta, alpha=alpha, minSeg=minSeg, segMedianT=c(.65,.35)))
    segs <- breakup_large_segments(segs)
    segs$start <- block_beta$start[segs$start]
    segs$end <- block_beta$end[segs$end]
    segs$width <- segs$end - segs$start
    segs$seqnames <- chrom
    segs$seg_id <- paste0(chrom,"_",i,"_",1:nrow(segs)) 
	segs
    }))
  makeGRangesFromDataFrame(segments,keep.extra.columns = T)
}
betas.gr <- makeGRangesFromDataFrame(beta)
# cap depth at 30 so all samples have similar influence on mean (higher coverage samples will not influence more)
#depth.mat[depth.mat > 20] <- 20
  
#pop_mean = rowSums(betas.mat*depth.mat, na.rm=T) / rowSums(depth.mat, na.rm=T)
meth_segs <- segment_blocks(pop_mean, this_chrom, minSeg = min_cpg_number)
  
ol <- findOverlaps(meth_segs, betas.gr)
# create Segment x CpG identity matrix to indicate which cpgs belong to which segment
CpG_Identity <- sparseMatrix(i = queryHits(ol), j = subjectHits(ol), dims=c(length(meth_segs),nrow(beta)), x=1)
  
beta.mat[is.na(beta.mat)] <- 0
depth.mat[is.na(depth.mat)] <- 0
beta.mat <- ((beta.mat*depth.mat) + 1) / (depth.mat + 2)
depth.mat <- depth.mat + 2

# aggregate betas across each segment
seg_beta <- as.matrix(((CpG_Identity %*% (beta.mat*depth.mat))) / ((CpG_Identity %*% depth.mat)))
seg_depth <- as.matrix(CpG_Identity %*% depth.mat / rowSums(CpG_Identity))
segments.df <- data.frame(chrom=this_chrom, start=start(meth_segs), end=end(meth_segs), seg_id=meth_segs$seg_id,
               num_cpg = meth_segs$num.mark, width=width(meth_segs), 
               cpg_density = meth_segs$num.mark/width(meth_segs), as.data.frame(seg_beta))

#save data
cat("Saving files...\n\n") 
fwrite(depth, file=argv$depth_mat, sep="\t", quote=F, scipen=999)
fwrite(beta, file=argv$beta_mat, sep="\t", quote=F, scipen=999)
fwrite(pop_mean, file=argv$pop_mean, sep="\t", quote=F, scipen=999)
fwrite(segments.df, file=argv$segment_beta, sep="\t", quote=F, scipen=999)
fwrite(as.data.frame(seg_depth), file=argv$segment_depth, sep="\t", quote=F, scipen=999)

cat("done!\n")
