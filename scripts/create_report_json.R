library(argparser)
library(jsonlite)
library(data.table)
library(dplyr)

parser <- arg_parser("Create json track config file for igv-reports from Metafora output")
parser <- add_argument(parser, "--sample", help="sample name")
parser <- add_argument(parser, "--outlier_bed", help="bed file of outliers for sample")
parser <- add_argument(parser, "--outlier_z_mat", help="matrix of zscores across all samples")
parser <- add_argument(parser, "--covariates", help="covariates matrix of z scores across outliers")
parser <- add_argument(parser, "--sample_table", help="sample table with bam paths to create paths")
parser <- add_argument(parser, "--annos", help="additional annotation track annotations to add to report")
parser <- add_argument(parser, "--output_bed", help="where to write tmp outlier bed")
parser <- add_argument(parser, "--other_out_bed", help="where to write tmp bed file for cohort outliers")
parser <- add_argument(parser, "--output_tsv", help="where to write tmp outlier tsv")
parser <- add_argument(parser, "--output_json", help="where to write output track json config file")
parser <- add_argument(parser, "--min_prio_score", help="minimum prioritization score to keep outlier in reprot, default=1", default=1)
parser <- add_argument(parser, "--flank_length", help="number of context base pairs flanking outlier to include in igv report", default=4000)
parser <- add_argument(parser, "--min_seg_size", help="minimum number of cpgs in outlier to include in report", default=30)
parser <- add_argument(parser, "--max_width", help="maxmimum width of outliers in bp to include in report", default=2e4)
parser <- add_argument(parser, "--top_outlier_n", help="maxmimum number of outliers to plot in report", default=20)
parser <- add_argument(parser, "--number_controls", help="number of control samples for comparison to include in igv report", default=20)
argv <- parse_args(parser)

this_sample <- argv$sample
outliers <- fread(argv$outlier_bed) %>% group_by(seqnames,start,end,seg_id) %>% summarize(gene_name=paste(gene_name,collapse="|"),
                                                                       gene_type=paste(gene_type,collapse="|"),
                                                                       gene_id=paste(gene_id,collapse="|"),
                                                                       across(everything(), first))
outlier_z_mat <- read.table(argv$outlier_z_mat, row.names=1)
covariates <- fread(argv$covariates)
sample_table <- fread(argv$sample_table)
annos <- fread(argv$annos)
MIN_PRIO_SCORE <- as.numeric(argv$min_prio_score)
FLANK_LENGTH <- as.numeric(argv$flank_length)
MIN_SEG_SIZE <- as.numeric(argv$min_seg_size)
MAX_WIDTH <- as.numeric(argv$max_width)
TOP_OUTLIER_N <- as.numeric(argv$top_outlier_n)
NUMBER_CONTROLS <- as.numeric(argv$number_controls)

INPUT_TYPE=sample_table$Technology[sample_table$Sample_name==this_sample]
if (!INPUT_TYPE %in% c("ONT","PacBio")) { #if sample is not ONT or PacBio . . . bams are not listed in sample table so do not plot IGV alignment as tracks
   NUMBER_CONTROLS<-0 
}
other_outs <- outliers[outliers$ID != this_sample,]
outliers <- outliers[outliers$ID == this_sample,]

#filter to prioritized outliers
outliers <- filter(outliers, 
                   PRIORITIZATION_SCORE >= MIN_PRIO_SCORE,
                   `num.mark` >= MIN_SEG_SIZE,
                   width < MAX_WIDTH)
outlier_z_mat <- outlier_z_mat[rownames(outlier_z_mat) %in% outliers$seg_id,]

outliers$coordinates <- paste0(outliers$seqnames,":",outliers$start,"-",outliers$end)
out_bed <- unique(outliers[,c("seqnames","start","end","seg_id","zscore")])
other_out_bed <- unique(other_outs[,c("seqnames","start","end","seg_id","zscore")])
outliers$start <- outliers$start - FLANK_LENGTH
outliers$end <- outliers$end + FLANK_LENGTH

#only make report for top N outliers
if(nrow(outlier_z_mat) >= TOP_OUTLIER_N) {
    top_outliers <- outliers$seg_id[order(10*outliers$PRIORITIZATION_SCORE+abs(outliers$zscore), decreasing=T)[1:TOP_OUTLIER_N]]
    outlier_z_mat <- outlier_z_mat[top_outliers,]
    outliers <- outliers[outliers$seg_id %in% top_outliers,]
    out_bed <- out_bed[out_bed$seg_id %in% top_outliers,]
}

outliers <- arrange(outliers, desc(10*PRIORITIZATION_SCORE+abs(zscore))) #sort outliers by prioritization score then outlier magnitude
fwrite(outliers,file=argv$output_tsv, sep="\t")
fwrite(out_bed, file=argv$output_bed, sep="\t", col.names=F)
fwrite(other_out_bed, file=argv$other_out_bed, sep="\t", col.names=F)

# Select control samples for comaprison to see  what "normal" methylation looks like
## choose samples from same batch / technology 
control_samples <- NULL
if (NUMBER_CONTROLS>0) {
    matched.outlier_z_mat <- outlier_z_mat
    if ("Batch" %in% colnames(covariates)) {
        this_batch <- covariates[covariates$Sample_name == this_sample,]$Batch
        matched.outlier_z_mat <- matched.outlier_z_mat[,colnames(matched.outlier_z_mat) %in% covariates$Sample_name[covariates$Batch==this_batch]]
    }
    if ("sex" %in% colnames(covariates)) {
        this_sex <- covariates[covariates$Sample_name == this_sample,]$sex
        matched.outlier_z_mat <- matched.outlier_z_mat[,colnames(matched.outlier_z_mat) %in% covariates$Sample_name[covariates$sex==this_sex]]
    }

    # if not enough Batch/Sex matched samples to pool from, pool from all samples instead 
    if (is.null(ncol(matched.outlier_z_mat)) || ncol(matched.outlier_z_mat) <= 3) {
       matched.outlier_z_mat <- outlier_z_mat 
    }

    #choose samples that are not themselves outliers over defined regions
    abs_maxes <- apply(abs(matched.outlier_z_mat), 2, max)
    if (sum(abs_maxes < 2) >= NUMBER_CONTROLS) {
        matched.outlier_z_mat <- matched.outlier_z_mat[,abs_maxes < 2]
    }

    # if this condition isnt satisfied optimize based on minimum abs sum of all zscores
    # trying to find the most "normal" samples over these regions, the ones closest to the median expected methylation
    min_abs_sums <- colSums(abs(matched.outlier_z_mat))
    control_samples <- names(min_abs_sums[order(min_abs_sums)][1:NUMBER_CONTROLS])
}

# make json track output
json_out <- list( 
    list( #outlier table
        name = paste0(this_sample, " Methylation Outliers"),
        url = normalizePath(argv$output_bed),
        color="red",
        height = 20
    ),
    list( #outlier in other samples
         name="Cohort Outliers",
         url = normalizePath(argv$other_out_bed),
         color="purple",
         height=40
    ))
if (INPUT_TYPE %in% c("ONT", "PacBio")) {
    json_out <- append(json_out, list( 
        list( #outlier sample bam
             name=paste0("Outlier Sample (", this_sample,")"),
             type="alignment",
             colorBy="basemod2:m",
             groupBy="tag:HP",
             alignmentRowHeight=12,
             hideSmallIndels="true",
             indelSizeThreshold="10",
             url=normalizePath(sample_table$Methylation_input[sample_table$Sample_name==this_sample]),
             height=400
    )))
}
for (this_control in control_samples) {
    json_out <- append(json_out, 
        list(list( #control sample
             name=paste0("Control (", this_control,")"),
             type="alignment",
             colorBy="basemod2:m",
             groupBy="tag:HP",
             alignmentRowHeight=12,
             hideSmallIndels="true",
             indelSizeThreshold="10",
             url=normalizePath(sample_table$Methylation_input[sample_table$Sample_name==this_control]),
             height=250
        )))
}
if (is.null(annos$color)) annos$color  <- rainbow(nrow(annos))
if (is.null(annos$height)) annos$height <- 50
for (i in 1:nrow(annos)) {
    json_out <- append(json_out, 
        list(list(
            name=annos$name[i],
            url=normalizePath(annos$url[i]),
            color=annos$color[i],
            height=annos$height[i]
        )))
}

json_out <- toJSON(json_out, pretty=TRUE, auto_unbox=T)
write(json_out, file=argv$output_json)
