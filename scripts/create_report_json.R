library(argparser)
library(jsonlite)
library(data.table)

parser <- arg_parser("Create json track config file for igv-reports from Metafora output")
parser <- add_argument(parser, "--sample", help="sample name")
parser <- add_argument(parser, "--outlier_bed", help="bed file of outliers for sample")
parser <- add_argument(parser, "--outlier_z_mat", help="matrix of zscores across all samples")
parser <- add_argument(parser, "--covariates", help="covariates matrix of z scores across outliers")
parser <- add_argument(parser, "--sample_table", help="sample table with bam paths to create paths")
parser <- add_argument(parser, "--annos", help="additional annotation track annotations to add to report")
parser <- add_argument(parser, "--output_bed", help="where to write tmp outlier bed")
parser <- add_argument(parser, "--output_tsv", help="where to write tmp outlier tsv")
parser <- add_argument(parser, "--output_json", help="where to write output track json config file")
argv <- parse_args(parser)

this_sample <- argv$sample

if(file.info(argv$outlier_bed)$size == 0) {
    fwrite(NULL, argv$output_bed)
    fwrite(NULL, argv$output_tsv)
    fwrite(NULL, argv$output_json)
    quit(save = "no")
}

outliers <- fread(argv$outlier_bed)
outlier_z_mat <- read.table(argv$outlier_z_mat, row.names=1)
covariates <- fread(argv$covariates)
sample_table <- fread(argv$sample_table)
annos <- fread(argv$annos)

if(is.null(outliers$hap_delta)) {
  outliers$haplotype_coverage_bias <- NA
  outliers$hap_delta <- NA
}

outliers$coordinates <- paste0(outliers$seqnames,":",outliers$start,"-",outliers$end)
out_bed <- outliers[,c("seqnames","start","end","seg_id","zscore")]
outliers$start <- outliers$start - 4000
outliers$end <- outliers$end + 4000

#only make report for top 100 outliers if more than a hundo
if(nrow(outlier_z_mat) > 100) {
    top_outliers <- order(abs(outlier_z_mat[,this_sample]), decreasing=T)[1:100]
    outlier_z_mat <- outlier_z_mat[1:nrow(outlier_z_mat) %in% top_outliers,]
    outliers <- outliers[1:nrow(outliers) %in% top_outliers,]
    out_bed <- out_bed[1:nrow(out_bed) %in% top_outliers,]
}

# Select 3 control samples to compare what "normal" methylation looks like
## choose samples from same batch / technology 
if ("Batch" %in% colnames(covariates)) {
    this_batch <- covariates[covariates$Sample_name == this_sample,]$Batch
    outlier_z_mat <- outlier_z_mat[,colnames(outlier_z_mat) %in% covariates$Sample_name[covariates$Batch==this_batch]]
}
if ("sex" %in% colnames(covariates)) {
    this_sex <- covariates[covariates$Sample_name == this_sample,]$sex
    outlier_z_mat <- outlier_z_mat[,colnames(outlier_z_mat) %in% covariates$Sample_name[covariates$sex==this_sex]]
}

# if not enough Batch/Sex matched samples to pool from, pool from all samples instead 
if (is.null(ncol(outlier_z_mat)) || ncol(outlier_z_mat) <= 4) {
   outlier_z_mat <- read.table(argv$outlier_z_mat, row.names=1) 
   outlier_z_mat <- outlier_z_mat[outliers$seg_id,]
}

#choose samples that are not themselves outliers over defined regions
abs_maxes <- apply(abs(outlier_z_mat), 2, max)
if (sum(abs_maxes < 2) >= 3) {
    outlier_z_mat <- outlier_z_mat[,abs_maxes < 2]
}

# if this condition isnt satisfied optimize based on minimum abs sum of all zscores
# trying to find the most "normal" samples over these regions, the ones closest to the median expected methylation
min_abs_sums <- colSums(abs(outlier_z_mat))
control_samples <- names(min_abs_sums[order(min_abs_sums)][1:3])

# make json track output
json_out <- list( 
    list( #outlier table
        name = "Methylation Outliers",
        url = argv$output_bed,
        height = 40
    ),
    list( #outlier sample bam
         name=paste0("Outlier Sample (", this_sample,")"),
         type="alignment",
         colorBy="basemod2",
         groupBy="tag:HP",
         url=sample_table$Methylation_input[sample_table$Sample_name==this_sample],
         height=250
    ))
for (this_control in control_samples) {
    json_out <- append(json_out, 
        list(list( #control sample
             name=paste0("Control (", this_control,")"),
             type="alignment",
             colorBy="basemod2",
             groupBy="tag:HP",
             url=sample_table$Methylation_input[sample_table$Sample_name==this_control],
             height=150
        )))
}
if (is.null(annos$color)) annos$color  <- rainbow(nrow(annos))
if (is.null(annos$height)) annos$height <- 50
for (i in 1:nrow(annos)) {
    json_out <- append(json_out, 
        list(list(
            name=annos$name[i],
            url=annos$url[i],
            color=annos$color[i],
            height=annos$height[i]
        )))
}

json_out <- toJSON(json_out, pretty=TRUE, auto_unbox=T)
write(json_out, file=argv$output_json)
fwrite(out_bed, file=argv$output_bed, sep="\t", col.names=F)
fwrite(outliers,file=argv$output_tsv, sep="\t")
