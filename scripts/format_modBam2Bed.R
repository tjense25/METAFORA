library(argparser)
library(data.table)

parser <- arg_parser("Format modBam2Bed bed files so they are consistent with METAFORA formatting")
parser <- add_argument(parser, "--input", help="input bed file from modBam2Bed")
parser <- add_argument(parser, "--output", help="output bed file of properly formatted bed")

argv <- parse_args(parser)


meth <- fread(argv$input, col.names=c("chromosome", "start", "end", "5mC", "score", "strand", "start_again", "end_again","truple","depth","beta"))
meth$start[meth$strand == "-"] <- meth$start[meth$strand=="-"] - 1
meth$end <- meth$start

meth$score <- meth$score / 1000
meth$beta <- meth$beta / 100

meth$depth <- round(meth$score * meth$depth) #correct depth for no_call and other filtered reads
meth$beta[is.na(meth$beta)] <- 0 #if beta is NA this is because effective depth is 0, (all reads filtered or no call) set beta to 0 to work with updated beta calcualtion below

meth <- meth[, .(depth=sum(depth), beta = sum(depth*beta)/sum(depth)), by=c("chromosome", "start", "end")]
fwrite(meth, file=argv$output, sep="\t", scipen=999)
