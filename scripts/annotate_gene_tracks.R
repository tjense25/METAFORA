library(argparser)
library(tidyverse)
library(plyranges)
library(data.table)
library(rtracklayer)

parser <- arg_parser("Annotate combined outlier file with gene info and supplied anno tracks")
parser <- add_argument(parser, "--outlier_bed", help="combined outlier bed file")
parser <- add_argument(parser, "--gene_model", help="gtf or gff of gene models")
parser <- add_argument(parser, "--annotation_tsv", help="tsv of additional user suppleid annotation tracks")
parser <- add_argument(parser, "--annotated_out", help="where to write annotated outlier file")
argv <- parse_args(parser)

outliers <- fread(argv$outlier_bed)
outliers.gr <- makeGRangesFromDataFrame(outliers, keep.extra.columns=T)

if (argv$annotation_tsv != "SKIP") {
    annos <- fread(argv$annotation_tsv)
    add_annotation_row <- function(outliers.gr, anno_idx) {
        cat(paste0("Annotating outliers with user-supplied track ",annos$name[anno_idx]," . . . \n"))
        tmp.outliers <- as.data.frame(outliers.gr)
        tmp.anno <- fread(annos$url[anno_idx])
        tmp.anno <- tmp.anno[,1:min(ncol(tmp.anno),4)]
        colnames(tmp.anno)[1:3] <- c("seqnames", "start","end")
        if (ncol(tmp.anno)==4) colnames(tmp.anno)[4] <- "Names"
        anno.gr <- makeGRangesFromDataFrame(tmp.anno,keep.extra.columns =T)
        tmp.outliers[[annos$name[anno_idx]]] <- overlapsAny(outliers.gr,anno.gr)
        tmp.overlap <- join_overlap_left(outliers.gr, anno.gr) %>% as.data.frame %>% group_by(seqnames,start,end,seg_id) %>% summarize(Names=paste(Names,collapse=",")) %>% as.data.frame
        rownames(tmp.overlap) <- tmp.overlap$seg_id
        tmp.outliers[[paste0(annos$name[anno_idx],"_names")]] <- tmp.overlap[tmp.outliers$seg_id,]$Names
        makeGRangesFromDataFrame(tmp.outliers, keep.extra.columns=T)
    }
    for (i in 1:nrow(annos)) {
        outliers.gr <- add_annotation_row(outliers.gr, i)
    }
}

if(argv$gene_model != "SKIP") {
    cat("Annotating gene model and promoters . . . \n")
    gff <- readGFF(argv$gene_model)
    genes <- gff[gff$type=="gene",]
    genes <- genes %>% as_tibble %>% select(seqid, start, end, strand, gene_type, gene_name, gene_id)

    genes.gr <- makeGRangesFromDataFrame(genes, keep.extra.columns=T)
    promoters <- flank_upstream(genes.gr, width=1000) %>% stretch(2000)
    genes.gr <- stretch(anchor_3p(genes.gr), 2000) #add promoter seq to gene

    outliers.gr$Promoter <- overlapsAny(outliers.gr, promoters)
    outliers.gr <- outliers.gr %>% join_overlap_left(genes.gr)
}
annotated_out <- outliers.gr %>% as.data.frame
fwrite(annotated_out, file=argv$annotated_out, sep="\t")
