library(argparser)
library(methylKit)
library(tidyverse)
library(magrittr)
library(data.table)
library(Biostrings)
library(plyranges)
library(bedr)
 
parser <- arg_parser("Aggregate individual beta profiles to calculate population average")
parser <- add_argument(parser, "--filelist", help="new line separated list of beta profiles file paths")
parser <- add_argument(parser, "--chrom", help="chromosome to compute betas from")
parser <- add_argument(parser, "--cpgs", help="bed file of all reference genome cpgs (only on autosomes)")
parser <- add_argument(parser, "--beta_mat", help="where to write beta matrix output")
parser <- add_argument(parser, "--depth_mat", help="where to write depth matrix output")
parser <- add_argument(parser, "--pop_mean", help="where to write population mean beta")

argv <- parse_args(parser)
#valid_chrom = paste0("chr", seq(1,22))
#ref = readDNAStringSet("/oak/stanford/groups/euan/annotations/homo_sapiens/GRCH38/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta")
#chrs <- names(ref)[1:22]
#cpgs <- lapply(chrs, function(x) start(matchPattern("CG", ref[[x]])))
#cpgr <- do.call(c, lapply(1:22, function(x) GRanges(valid_chrom[x], IRanges(cpgs[[x]],width=2))))
#cpgr %<>% as.data.frame()

this_chrom <- argv$chrom
methylation.dat <- fread(argv$cpgs)
colnames(methylation.dat) <- c("chromosome", "start","end")
methylation.dat <- methylation.dat[methylation.dat$chromosome == this_chrom,]
max_start<- max(methylation.dat$start) + 1
this_region=paste0(this_chrom,":",1,"-",max_start)

#read in files and merge depth and betas to cpg data.table 
file_list=read.table(argv$filelist,col.names="file")$file
for (f in file_list) {
    sample_id = str_split(basename(f), pattern='\\.',simplify =T)[1]
    cat(paste0('reading methylation profile for sample ', sample_id, ' from file ', f,'\n'))
    tmp.cpg <- tabix(region=this_region,f, verbose=F)
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

#save data
cat("Saving files...\n\n")
fwrite(depth, file=argv$depth_mat, sep="\t", quote=F, scipen=999)
fwrite(beta, file=argv$beta_mat, sep="\t", quote=F, scipen=999)
fwrite(pop_mean, file=argv$pop_mean, sep="\t", quote=F, scipen=999)

cat("done!\n")
