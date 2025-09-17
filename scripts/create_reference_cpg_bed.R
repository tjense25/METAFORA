library(tidyverse)
library(magrittr)
library(Biostrings)
library(GenomicRanges)
library(argparser)
library(data.table)

parser <- arg_parser("Read a reference genome fasta and output a bed of all CpG sites")

parser <- add_argument(parser, "--reference", help="Reference Genome Fasta File")
parser <- add_argument(parser, "--valid_chroms", help="comma-separated list of valid chromosomes to call methylation")
parser <- add_argument(parser, "--cpg_bed_out", help="where to write bed file of CpG sites")
parser <- add_argument(parser, "--block_size", help="how many cpgs to include in block")
parser <- add_argument(parser, "--block_bed_out", help="where to write bed file defining blocks")

argv <- parse_args(parser)

valid_chroms = unlist(strsplit(argv$valid_chroms, ","))
ref = readDNAStringSet(argv$reference)
names(ref) <- gsub(names(ref), pattern="(^\\w+)\\s.*$", replacement="\\1")
cpgs <- lapply(valid_chroms, function(x) start(matchPattern("CG", ref[[x]])))
cpgr <- do.call(c, lapply(1:length(valid_chroms), function(x) GRanges(valid_chroms[x], IRanges(cpgs[[x]],width=1))))
cpgr %<>% as.data.frame()
cpgr %<>% select(seqnames,start,end) %>% mutate(start = start -1, end=end-1)
fwrite(cpgr, argv$cpg_bed_out, col.names =F, row.names=F, sep="\t", scipen=999)

colnames(cpgr) <- c("seqnames", "start", "end")

block_size<-as.integer(argv$block_size)

cpgr$large_diff <- c(1,abs(diff(cpgr$start)) > 1000) #find cpgs that are more than 1kb away from next adjacent cpg as natural break points 
cpgr %<>% group_by(seqnames) %>% mutate(large_diff=ifelse(start==min(start)|end==max(end), yes=1, no=large_diff), cpg_num=1:n())
cpg_blocks <- cpgr %>% filter(large_diff==1) %>% group_by(seqnames) %>% 
    mutate(block=cpg_num%/%block_size+1) %>% group_by(seqnames,block) %>% 
    summarize(start=min(start), end=max(end), start_cpg_num=min(cpg_num), end_cpg_num=max(cpg_num)) %>% 
    mutate(num_cpgs=end_cpg_num-start_cpg_num) %>% 
    mutate(block=ifelse(num_cpgs<block_size/3,yes=max(1,block-1),no=block)) %>% 
    group_by(seqnames,block) %>% summarize(start=min(start), end=max(end), num_cpgs=sum(num_cpgs)) %>%
    mutate(block=paste0(seqnames,".",block)) %>%
    select(seqnames, start, end, block, num_cpgs)
cpg_blocks

fwrite(cpg_blocks, argv$cpg_bed_out, col.names=T, row.names=T, sep="\t")
