library(tidyverse)
library(magrittr)
library(data.table)
library(ggrepel)
library(cowplot)
library(argparser)


parser <- arg_parser("Plot Summary Figures")
parser <- add_argument(parser, "--outlier_files", help="comma-separated list of outlier summary files from different tissues")
parser <- add_argument(parser, "--plot_dir_out", help="where to save summary plots")
parser <- add_argument(parser, "--summary_out", help="where to write summary tsv file of outlier counts per sample")
parser <- add_argument(parser, "--covariates_files", help="optional comma-separated list of covariates to plot data stratified by Batch", default=NULL)
parser <- add_argument(parser, "--min_abs_delta", help="MIN_ABS_DELTA threshold", type="numeric", default=.25) 
parser <- add_argument(parser, "--min_abs_zscore", help="MIN_ABS_ZSCORE threhsold", type="numeric", default=3)
#argv <- NULL
#argv$outlier_files <- "../METAFORA_output/METAFORA.tissue_Blood.methylation_outliers.combined.tsv"
#argv$plot_dir_out <- "../METAFORA_output/summary_figures"
#argv$covariates_files <- "../METAFORA_output/Global_Methylation_PCA_tissue_Blood/PCA_covariates.txt"
#argv$min_abs_delta <- 0.25
#argv$min_abs_zscore <- 3
argv <- parse_args(parser)





outlier_files <- unlist(strsplit(argv$outlier_files, ","))
covariates_files <- NULL
if(!is.null(argv$covariates_files)) {
    covariates_files <- unlist(strsplit(argv$covariates_files,","))
}
MIN_ABS_DELTA <- argv$min_abs_delta
MIN_ABS_ZSCORE <- argv$min_abs_zscore

number_tissue = length(outlier_files)

plot_count_per_samp <- function(count_per_samp, label, color) {
    ggplot(count_per_samp, aes(rank,n, label=ifelse(outlier,Sample_name,""))) + geom_point(color=color) + theme_minimal() + 
    geom_hline(aes(yintercept=mean(n)),color="red") + 
    annotate("text", x=2,y=mean(count_per_samp$n)+1,label="Mean", color="red") +
    geom_text_repel(color="red") + 
    geom_hline(aes(yintercept=median(n)), color="red", linewidth=2) + 
    annotate("text", x=2, y=median(count_per_samp$n)+1,label="Median", color="red")+
    geom_hline(aes(yintercept=mean(n)+sd(n)), color="red", linetype="dashed") + 
    annotate("text", x=2,y=mean(count_per_samp$n)+sd(count_per_samp$n)+1,label="Mean + 1 s.d.", color="red") +
    geom_hline(aes(yintercept=mean(n)+2*sd(n)), color="red", linetype="dashed") + 
    annotate("text", x=2,y=mean(count_per_samp$n)+2*sd(count_per_samp$n)+1,label="Mean + 2 s.d.", color="red") +
    geom_hline(aes(yintercept=mean(n)+3*sd(n)), color="red", linetype="dashed") + 
    annotate("text", x=2,y=mean(count_per_samp$n)+3*sd(count_per_samp$n)+1,label="Mean + 3 s.d.", color="red") +
    theme(panel.border=element_rect(color="black", fill=NA,linewidth=2)) +
    xlab("rank-order samples by number of outliers") +  ylab("count of outlier per genome") + 
    ggtitle(label)
}

combined_outliers <- NULL
for (i in 1:number_tissue) {
    f <- outlier_files[i]
    this_tissue <- gsub(x=f, pattern=".+\\.tissue_(\\w+)\\..+", replacement="\\1")
    outliers  <- fread(f)
    outliers %<>% dplyr::rename(Sample_name=ID) 
    covariates <- NULL
    outliers$Batch <- "Full Cohort"
    outliers$sex <- "NA"
    if(!is.null(covariates_files)) {
        covariates <- fread(covariates_files[i])
        if ("Batch" %in% colnames(covariates)) {
            outliers$Batch <- covariates$Batch[match(outliers$Sample_name, covariates$Sample_name)]
        }
    }

    outliers %>% filter(CHROM_TYPE=="AUTOSOME") %>% group_by(Sample_name, Batch,CHROM_TYPE) %>% summarize(hyper=sum(delta > 0), hypo=sum(delta < 0)) %>%
        pivot_longer(c(hyper,hypo), names_to="direction", values_to="count") %>% 
        ggplot(aes(Batch, count, fill=direction)) + geom_boxplot(alpha=.8) + theme_minimal() + facet_wrap(~CHROM_TYPE) + 
        theme(text=element_text(size=16), panel.border=element_rect(color="black",fill=NA,size=2), axis.text.x = element_text(angle=25,hjust=1))  + 
        ylab("number of Methylation Outliers per genome") + 
        scale_fill_manual(values=c("brown4", "cadetblue4"))
    ggsave(paste0(argv$plot_dir_out, "/Autosomes.number_outliers_per_genome.box_plot.tissue_", this_tissue,".pdf"))
    outliers %>% filter(CHROM_TYPE=="AUTOSOME") %>% group_by(Sample_name, Batch,CHROM_TYPE) %>% summarize(hyper=sum(delta > 0), hypo=sum(delta < 0)) %>%
        pivot_longer(c(hyper,hypo), names_to="direction", values_to="count") %>% 
        ggplot(aes(Batch, count, fill=direction)) + geom_boxplot(alpha=.8) + theme_minimal() + facet_wrap(~CHROM_TYPE) + 
        theme(text=element_text(size=16), panel.border=element_rect(color="black",fill=NA,size=2), axis.text.x = element_text(angle=25,hjust=1))  + 
        ylab("number of Methylation Outliers per genome") + scale_y_log10()
        scale_fill_manual(values=c("brown4", "cadetblue4"))
    ggsave(paste0(argv$plot_dir_out, "/Autosomes.number_outliers_per_genome.box_plot.tissue_", this_tissue,"logscale_y.pdf"))

    #if("sex" %in% colnames(covariates)) {
    #    outliers$sex <- covariates$sex[match(outliers$Sample_name, covariates$Sample_name)]
    #    outliers %>% mutate(sex=ifelse(sex,"XY", "XX")) %>% 
    #        filter(CHROM_TYPE=="SEX_CHROM") %>% 
    #        group_by(seqnames, Sample_name, sex) %>% summarize(hyper=sum(delta>0), hypo=sum(delta < 0)) %>% 
    #        pivot_longer(c(hyper,hypo), names_to="direction", values_to="count") %>% 
    #        ggplot(aes(sex, count, fill=direction)) + geom_boxplot() + theme_minimal() +
    #        facet_grid(~seqnames) + 
    #        theme(panel.border=element_rect(fill=NA)) + ggtitle("Sex Chromosomes") +
    #        scale_fill_manual(values=c("brown4", "cadetblue4"))
    #    ggsave(paste0(argv$plot_dir_out,"/Sex_chromosomes.number_outliers_stratified_by_sex.box_plot.tissue_",this_tissue,".pdf"))
    #}

    ## Outlier counts per sample
    count_per_samp <- outliers %>% group_by(Sample_name, Batch) %>% summarize(n=n()) %>% ungroup() %>% arrange(n) %>% mutate(rank=1:length(n)) %>% 
        mutate(outlier = n > mean(n)+3*sd(n))
    combined <- plot_count_per_samp(count_per_samp, "Total Number of Methylation Outlier", "black")
    count_per_samp <- outliers %>% group_by(Sample_name, Batch) %>% summarize(n=sum(delta<0)) %>% ungroup() %>% arrange(n) %>% mutate(rank=1:length(n)) %>% 
        mutate(outlier = n > mean(n)+3*sd(n))
    hypo <- plot_count_per_samp(count_per_samp, "Number of HYPOmethylation Outliers", "cadetblue4")
    count_per_samp <- outliers %>% group_by(Sample_name, Batch) %>% summarize(n=sum(delta>0)) %>% ungroup() %>% arrange(n) %>% mutate(rank=1:length(n)) %>% 
        mutate(outlier = n > mean(n)+3*sd(n))
    hyper <- plot_count_per_samp(count_per_samp, "Number of HYPERmethylation Outliers", "brown4")
    plot_grid(combined, hypo, hyper, ncol=1) 
    ggplot2::ggsave(paste0(argv$plot_dir_out, "/Rank_order_number_outliers_per_genome.dot_plot.tissue_",this_tissue,".pdf"), height=21)

    ggplot(outliers, aes(abs(zscore),fill=Batch)) + geom_density(alpha=.1) + theme_minimal()
    ggsave(paste0(argv$plot_dir_out,"/zscore_distribution_across_batch.tissue_",this_tissue,".pdf"))

    outlierlength <- outliers %>% mutate(length=end-start) %>% 
        ggplot(aes(length)) + geom_density(fill="chartreuse", alpha=.2) + scale_x_log10() + theme_minimal() +
        xlab("Length of outlier segment (bp)") +
        theme(panel.border=element_rect(color="black",fill=NA,size=2),
                text = element_text(size=16))
    outlier_nummark <- outliers %>% 
        ggplot(aes(num.mark)) + geom_density(fill="goldenrod2", alpha=.2) + scale_x_log10() + theme_minimal() +
        xlab("Number of CpG in outlier segment") +
        theme(panel.border=element_rect(color="black",fill=NA,size=2),
                text = element_text(size=16))
    cpgdensity <- outliers %>% mutate(length=end-start) %>% mutate(cpg_density=num.mark/length) %>% 
        ggplot(aes(cpg_density)) + geom_density(fill="dodgerblue2", alpha=.2) + theme_minimal() +
        xlab("CpG Density (number_CpG/length)") +
        theme(panel.border=element_rect(color="black",fill=NA,size=2),
                text = element_text(size=16))
    outlier_popmedian <- outliers %>% 
        ggplot(aes(pop_median)) + geom_density(fill="deeppink", alpha=.2) + theme_minimal() +
        xlab("Population Median Methylation of Outlier Segments") +
        theme(panel.border=element_rect(color="black",fill=NA,size=2),
                text = element_text(size=16))
    plot_grid(outlierlength,outlier_nummark,cpgdensity,outlier_popmedian, ncol=2)
    ggplot2::ggsave(paste0(argv$plot_dir_out, "/Outlier_summary_stats.density_plots.tissue_",this_tissue,".pdf"), width=12,height=12)

    length_by_nummark <- outliers %>% mutate(length=end-start) %>% 
        ggplot(aes(length, num.mark)) + geom_point() + scale_x_log10() + theme_minimal() + 
        xlab("Length of outlier segment (bp)") + ylab("Number of CpG in outlier segment") +
        theme(panel.border=element_rect(color="black",fill=NA,size=2),
                text = element_text(size=16))
    cpgdensity_by_popmedian <- outliers %>% mutate(length=end-start) %>% mutate(cpg_density=num.mark/length) %>% 
        ggplot(aes(cpg_density, pop_median)) + geom_point() + theme_minimal() +
        xlab("CpG Density in outlier segment") + ylab("Population Median Methylation") + 
        theme(panel.border=element_rect(color="black",fill=NA,size=2),
                text = element_text(size=16))
    plot_grid(length_by_nummark,cpgdensity_by_popmedian)
    ggplot2::ggsave(paste0(argv$plot_dir_out, "/Outlier_summary_stats.scatter_plots.tissue_",this_tissue,".pdf"), width=12)

    combined_outliers <- rbind(combined_outliers, outliers)
}

## Volcano Plot
ggplot(combined_outliers, aes(delta, abs(zscore), shape=CHROM_TYPE, color=Tissue)) + geom_point() + 
    theme_minimal() + xlim(-1,1) + scale_color_brewer(palette ="Set3") +
    facet_wrap(~Tissue) + theme(panel.border=element_rect(color="black", fill=NA)) + 
    geom_vline(xintercept=-MIN_ABS_DELTA, linetype="dashed", color="grey30") +
    geom_vline(xintercept=MIN_ABS_DELTA, linetype="dashed", color="grey30") +
    geom_hline(yintercept=MIN_ABS_ZSCORE, color="grey50") +
    theme(text=element_text(size=16))
ggsave(paste0(argv$plot_dir_out, "/Outlier_delta_zscore.volcano_plot.pdf"), width=number_tissue*6)

outlier_counts <- combined_outliers %>% group_by(Sample_name, CHROM_TYPE, Tissue) %>% summarize(total_outliers=n(), 
                                                                                       number_hypo_outliers=sum(delta < 0 ),
                                                                                       number_hyper_outliers=sum(delta > 0))

ggplot(outlier_counts, aes(Tissue, total_outliers, fill=Tissue)) +
        geom_boxplot(alpha=.8) + theme_minimal() + facet_wrap(~CHROM_TYPE) + 
        theme(text=element_text(size=16), panel.border=element_rect(color="black",fill=NA,size=2))  + 
        ylab("number of Methylation Outliers per genome") + scale_y_log10() + 
        scale_fill_brewer(palette="Set3")
ggsave(paste0(argv$plot_dir_out, "/Number_outliers_per_genome.box_plot.across_tissues.pdf"))
fwrite(outlier_counts, argv$summary_out, col.names=T, row.names=F, sep="\t")
