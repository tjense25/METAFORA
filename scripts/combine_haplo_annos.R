library(argparser)
library(tidyverse)
library(magrittr)
library(data.table)

parser <- arg_parser("Combine haplotype annotated outliers into master file")
parser <- add_argument(parser, "--outlier_bed", help="outlier bed file containing outliers for all samples")
parser <- add_argument(parser, "--haplo_anno_list", help="list of sample-level outlier files annotated with phasing info, one file per line, SKIP if no samples were phased")
parser <- add_argument(parser, "--haplo_anno_out", help="where to write combined haplotype annotated outlier bed file")
argv <- parse_args(parser)

combined_outliers <- fread(argv$outlier_bed)

if (argv$haplo_anno_list == "SKIP") {
  fwrite(combined_outliers, argv$haplo_anno_out, sep="\t")
  quit(save="no")
}

file_list=read.table(argv$haplo_anno_list,col.names="file")$file
file_list <- file_list[file.info(file_list)$size > 1]
haplo_annos <- lapply(file_list, fread) %>% bind_rows

#filter samples that are phased from combined outliers so file contains only NOT-phased samples
combined_outliers <- combined_outliers[!(combined_outliers$seg_id %in% haplo_annos$seg_id),]

#set phase specific columns of remaining file to NA
if (nrow(combined_outliers)>0) {
    extra_cols = colnames(haplo_annos)[!(colnames(haplo_annos) %in% colnames(combined_outliers))]
    for (this_col in extra_cols) {
      combined_outliers[[this_col]] <- NA
    }
    haplo_annos <- rbind(haplo_annos, combined_outliers)
}

haplo_annos %<>% arrange(seqnames, start)

fwrite(haplo_annos, argv$haplo_anno_out, sep="\t")
