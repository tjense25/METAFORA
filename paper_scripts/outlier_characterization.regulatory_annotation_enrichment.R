library(data.table)
library(tidyverse)
library(plyranges)
library(matrixStats)
library(magrittr)
library(cowplot)
library(rtracklayer)

outliers <- fread("../METAFORA_output/METAFORA_methylation_outlier_regions.tissue_Blood.ALL_CHROM_COMBINED.haplotype_annotated.gene_track_annotated.bed") %>% filter(CHROM_TYPE=="AUTOSOME")
out.segs <- outliers %>% select(seqnames,start,end,seg_id,pop_median) %>% mutate(class="outlier") %>% dplyr::rename(chrom=seqnames)

pop_segment <- fread("../METAFORA_output/Population_methylation.tissue_Blood/Meth_segments.tissue_Blood.segment_betas.bed")  %>% filter(!chrom %in% c("chrX","chrY"))
pop_segment.mat <- as.matrix(pop_segment[,8:ncol(pop_segment)])
pop_segment$pop_median <- rowMedians(pop_segment.mat)
bg.segs <- pop_segment %>% select(chrom,start,end,seg_id,pop_median) %>% mutate(class="background")

#sample 
set.seed(123)
seg.gr <- makeGRangesFromDataFrame(rbind(out.segs, bg.segs[sample(1:nrow(bg.segs),size=1e5)]), keep.extra.columns=T)

gnocchi.gr <- fread("../references/gnocchi_1kb_window.scores.bed.gz")
colnames(gnocchi.gr) <- c("chr","start","end","gnocchi_z")
gnocchi.gr %<>% as_granges(seqnames=chr)

seg.gnocchi <- seg.gr %>% join_overlap_left(gnocchi.gr) %>% group_by(seg_id,class) %>% summarize(gnocchi_z=mean(gnocchi_z,na.rm=T))
fit <- glm(data=seg.gnocchi %>% as.data.frame%>%filter(!is.na(gnocchi_z)), formula=(class=="outlier") ~ gnocchi_z)
summary(fit)$coefficients

length.plt <- seg.gr %>% as.data.frame %>%
    ggplot(aes(width,fill=class)) + geom_density(alpha=.4) + scale_x_log10(limits=c(50,1e5)) + scale_fill_manual(values=c("black", "red")) + theme_minimal() +
    theme(legend.position="none")
pop_med.plt <- seg.gr %>% as.data.frame %>%
    ggplot(aes(pop_median,fill=class)) + geom_density(alpha=.4) + scale_fill_manual(values=c("black","red")) + theme_minimal() +
    theme(legend.position="none")
noncoding_constraint.plt <- seg.gnocchi %>% as.data.frame %>%
    ggplot(aes(gnocchi_z, fill=class)) + geom_density(alpha=.4) + xlim(-10,10) +  scale_fill_manual(values=c("black","red")) + theme_minimal() +
    theme(legend.position="none")
plot_grid(length.plt, pop_med.plt, noncoding_constraint.plt, ncol=1)
ggsave('../paper_plots/outlier_summary.length_popmed_noncoding_constraint.pdf', height=4,width=4)

seg.gr %<>% mutate(variable_methylated=abs(pop_median-.5)<.25, hypomethylated=pop_median<.25, hypermethylated=pop_median>.75) #segment variably methylated in popmedian beta between 0.25-0.75
fisher.test(table(seg.gr$variable_methylated, seg.gr$class=="outlier"))
fisher.test(table(seg.gr$hypomethylated, seg.gr$class=="outlier"))
fisher.test(table(seg.gr$hypermethylated, seg.gr$class=="outlier"))

#regulatory annotation enrichments

#collate granges coordinates for variaous regulatory annotations
#CpG islands downloaded from UCSC genome browser
autosomes <- paste0("chr",1:22)
cpg.islands <- fread("../references/cpgIslandExt.GRCh38.bed.gz") %>% filter(V1 %in% autosomes)
islands.gr <- makeGRangesFromDataFrame(cpg.islands, seqnames.field="V1", start.field="V2", end.field="V3")
shores.gr <- c( #cpg shores are 2kb up/downstream of islands
                   islands.gr %>% flank_upstream(2e3),
                   islands.gr %>% flank_downstream(2e3)
            )
shelfs.gr <- c( #cpg shores are 2-4kb up/downstream of islands,)
                   islands.gr %>% flank_upstream(2e3) %>% flank_upstream(2e3),
                   islands.gr %>% flank_downstream(2e3) %>% flank_downstream(2e3)
                   )

# ABC elements downloaded from Engreitz website, lifted over to GRCh38, and merged across all celltypes
ABC <- fread("../references/ABC_enhancers.merged.GRCh38.bed.gz") %>% filter(V1 %in% autosomes)
abc.gr <- makeGRangesFromDataFrame(ABC,seqnames.field="V1", start.field="V2", end.field="V3")
abc.gr$element_class <- ABC$V5
table(abc.gr$element_class)
abc.prom.gr <- abc.gr %>% filter(element_class=="promoter")
abc.genic.gr <- abc.gr %>% filter(element_class=="genic")
abc.intergenic.gr <- abc.gr %>% filter(element_class=="intergenic")

#ATAC peak joint peak set generated from GSS ATAC peaks following encode protocol for ATAC
ATAC <- fread("../references/GSS_PBMC_ATACpeaks.open_chromatin.GRCh38.bed.gz") %>% filter(V1 %in% autosomes)
pbmc.atac.gr <- makeGRangesFromDataFrame(ATAC,seqnames.field="V1",start.field="V2",end.field="V3")

#ENCODE cCRES downloaded from UCSC genome browser
ccre.gr <- import.bb("../references/encodeCcreCombined.bb")
ccre.gr <- ccre.gr[seqnames(ccre.gr) %in% autosomes]
prom.gr <- ccre.gr %>% filter(ucscLabel=="prom")
enhP.gr <- ccre.gr %>% filter(ucscLabel=="enhP")
enhD.gr <- ccre.gr %>% filter(ucscLabel=="enhD")
ctcf.gr <- ccre.gr %>% filter(ucscLabel=="CTCF")
k4m3.gr <- ccre.gr %>% filter(ucscLabel=="K4m3")

enrichment <- function(seg.gr, anno.gr, label=deparse(substitute(anno.gr))) {
    seg.gr$anno_overlap <- seg.gr %>% overlapsAny(anno.gr)
    fisher_result <- fisher.test(table(seg.gr$class=="outlier",seg.gr$anno_overlap))
    c(label,fisher_result$estimate, fisher_result$p.value, fisher_result$conf.int)
}

enrichment_results <- data.frame(rbind( 
    enrichment(seg.gr, islands.gr),
    enrichment(seg.gr, shores.gr),
    enrichment(seg.gr, shelfs.gr),
    enrichment(seg.gr,abc.gr),
    enrichment(seg.gr,abc.prom.gr),
    enrichment(seg.gr,abc.genic.gr),
    enrichment(seg.gr,abc.intergenic.gr),
    enrichment(seg.gr,pbmc.atac.gr),
    enrichment(seg.gr,ccre.gr),
    enrichment(seg.gr,prom.gr),
    enrichment(seg.gr,enhP.gr),
    enrichment(seg.gr,enhD.gr),
    enrichment(seg.gr,ctcf.gr),
    enrichment(seg.gr,k4m3.gr)))
colnames(enrichment_results) <- c("label", "odds_ratio", "p_value", "lower", "upper")
enrichment_results$odds_ratio %<>% as.numeric()
enrichment_results$lower %<>% as.numeric()
enrichment_results$upper %<>% as.numeric()
enrichment_results$label %<>% factor(levels=rev(enrichment_results$label))
fwrite(enrichment_results, "outlier_functional_annotations.enrichment.csv")
#enrichment_results <- fread("./outlier_functional_annotations.enrichment.csv")

ggplot(enrichment_results, aes(odds_ratio, label, xmin=lower,xmax=upper)) + geom_pointrange(size=.1) + theme_minimal() +
    geom_vline(xintercept=1,linetype="dashed",color="grey") + 
    theme(legend.position="none",text=element_text(size=18), panel.border=element_rect(fill=NA)) + 
    xlab("odds ratio (enrichment)") + ylab("Regulatory annotation") 
ggsave("../paper_plots/outlier_enrichment.regulatory_annotations.pdf", width=5,height=6)


multiomics.enrichment <- fread("./methylation_U08.odds_ratio.omics_overlap.txt") %>% select(-ome1,-assay1,-assay2)
colnames(multiomics.enrichment) <- c("label", "odds_ratio", "p_value", "lower", "upper")
combined_enrich <- rbind( 
    enrichment_results %>% mutate(enrich_type="segment_level"),
    multiomics.enrichment %>% mutate(enrich_type="gene_level_multiomics")
)
combined_enrich$label %<>% factor(levels=rev(combined_enrich$label))
ggplot(combined_enrich, aes(odds_ratio,label,xmin=lower,xmax=upper,color=enrich_type)) + geom_pointrange(size=0) + theme_minimal() +
    geom_vline(xintercept=1,linetype="dashed",color="grey") +
    theme(legend.position="none",text=element_text(size=11),panel.border=element_rect(fill=NA)) +
    xlab("odds ratio (enrichment)") + ylab("Regulatory annotation") +
    scale_color_manual(values=c("firebrick2","black")) + scale_x_log10()
ggsave("../paper_plots/combined_enrichment.functional_annos_and_multiomics.pointrange.pdf", width=4,height=4) 

    
combined_enrich



