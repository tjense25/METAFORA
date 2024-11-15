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

argv <- parse_args(parser)

valid_chroms = unlist(strsplit(argv$valid_chroms, ","))
ref = readDNAStringSet(argv$reference)
names(ref) <- gsub(names(ref), pattern="(^\\w+)\\s.*$", replacement="\\1")
cpgs <- lapply(valid_chroms, function(x) start(matchPattern("CG", ref[[x]])))
cpgr <- do.call(c, lapply(1:length(valid_chroms), function(x) GRanges(valid_chroms[x], IRanges(cpgs[[x]],width=1))))
cpgr %<>% as.data.frame()
cpgr %<>% select(seqnames,start,end) %>% mutate(start = start -1, end=end-1)
fwrite(cpgr, argv$cpg_bed_out, col.names =F, row.names=F, sep="\t", scipen=999)
