library(tidyverse)
library(magrittr)
library(data.table)

metafora.outs <- fread("./GREGoR_U08.METAFORA_outlier_regions.combined_tissues.gene_annotated.haplotype_annotated.ImprintDisord_Column.txt")
head(metafora.outs)
colnames(cand.patients)
cand.patients <- metafora.outs %>% filter(ImprintDisordDMR) %>% filter(pop_median < .65, pop_median > .35)
fwrite(cand.patients, "GREGoR_U08.imprinting_disorder_outliers.csv")
#added manual curation based on visualization of outliers and phenotypes
cand.patients <- fread("./GREGoR_U08.imprinting_disorder_outliers.manual_curation_annot.csv")
table(cand.patients$ImprintDisordDMR_names)

cand.patients$Sample_name <- gsub(cand.patients$seg_id, pattern="chr\\d+.\\d+_(\\w+)_\\d+", replacement="\\1")
cand.patients %<>% left_join(samp.tab  %>% select(Sample_name,Technology))

fisher.test(table(cand.patients$Technology, cand.patients$HAP_BIAS))

cand.patients.barplt.df <- rbind( 
     cand.patients %>% select(ImprintDisordDMR_names, Technology, affected_status, Tissue) %>% mutate(filter="OVERLAP_DMR"),
     cand.patients %>% filter(HIGH_COV) %>% select(ImprintDisordDMR_names, Technology, affected_status, Tissue) %>% mutate(filter="OVERLAP_DMR+HIGH_COV"),
     cand.patients %>% filter(HIGH_COV,!HAP_BIAS) %>% select(ImprintDisordDMR_names, Technology, affected_status, Tissue) %>% mutate(filter="OVERLAP_DMR+HIGH_COV+NO_HAP_BIAS"),
     cand.patients %>% filter(HIGH_COV,!HAP_BIAS,PHENO_MATCH) %>% select(ImprintDisordDMR_names, Technology, affected_status, Tissue) %>% mutate(filter="OVERLAP_DMR+HIGH_COV+NO_HAP_BIAS+PHENO_MATCH")
)

cand.patients.barplt.df$ImprintDisordDMR_names %<>% factor()
cand.patients.barplt.df$filter %<>% factor(levels=c("OVERLAP_DMR+HIGH_COV+NO_HAP_BIAS+PHENO_MATCH","OVERLAP_DMR+HIGH_COV+NO_HAP_BIAS","OVERLAP_DMR+HIGH_COV","OVERLAP_DMR"))
cand.patients.barplt.df %>% group_by(ImprintDisordDMR_names, filter, .drop=F) %>% summarize(n=n()) %>%
    mutate(color_group=factor(case_when(n==0~"0",n==1 ~ "1", n<=3 ~ "2-3", n<=6 ~ "4-6", n <= 10 ~ "7-9"), levels=c("0","1","2-3","4-6","7-9"))) %>% 
    ggplot(aes(ImprintDisordDMR_names,filter,fill=color_group)) + 
    geom_tile(color="black") + 
    geom_text(aes(label = ifelse(n != 0, n, "")), color = "white") + 
    theme_classic() + scale_fill_manual(values=list("0"="white","1"="lightpink","2-3"="firebrick1","4-6"="firebrick3","7-9"="firebrick4"))  +
    theme(legend.position="none", axis.text.x=element_text(angle=45,vjust=1,hjust=1)) 
cand.patients.barplt.df %>% group_by(ImprintDisordDMR_names, filter, .drop=F) %>% summarize(n=n()) %>%
    ggplot(aes(ImprintDisordDMR_names,filter,fill=log(n+1))) + 
    geom_tile(color="black") + 
    geom_text(aes(label = ifelse(n != 0, n, "")), color = "white", fontface="bold") + 
    scale_fill_gradient(low = "grey95", high = "firebrick3") + 
    theme(legend.position="none", axis.text.x=element_text(angle=45,vjust=1,hjust=1)) 
ggsave("ImprintDMR_heatmap.by_filters.pdf",height=4, width=8)

cand.patients.barplt.df$affected_status[cand.patients.barplt.df$affected_status=="Unknown"]="Unaffected"
ggplot(cand.patients.barplt.df, aes(filter,fill=affected_status)) + geom_bar(width=.5,color="black") + theme_minimal() + coord_flip() +  
    scale_fill_manual(values=c("firebrick3","grey60")) +
    scale_alpha_manual(values=c(.2,1))
ggsave("ImprintDMR_barplot.by_filters.pdf",height=3)

cand.patients$combined_depth
cand.patients$haplotype_coverage_bias
table(cand.patients$affected_status)
cand.patients %<>% mutate(hap_biased = abs(haplotype_coverage_bias)>1)
biased_outliers_remove <- cand.patients  %>% filter(hap_biased & affected_status!="Affected") %>% pull(seg_id)
cand.patients %<>% filter(!seg_id %in% biased_outliers_remove)

colnames(cand.patients)
cand.patients.pretty <- cand.patients %>% mutate(region=paste0(seqnames,":",start,'-',end)) %>% select(region,participant_id,seg.mean,seg_id,pop_median,delta,zscore,Tissue,haplotype_coverage_bias,hap_delta,ImprintDisordDMR_names,affected_status,hap_biased)
cand.patients.pretty
fwrite(cand.patients.pretty, "GREGoR_U08.candidate_imprinting_disorder_patients.METAFORA_outliers.csv")
