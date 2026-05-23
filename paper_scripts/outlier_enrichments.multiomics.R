library(tidyverse)
library(data.table)
library(magrittr)
library(plyranges)
library(rtracklayer)
library(epitools)
library(viridis)
library(RColorBrewer)
library(cowplot)

## Extract and Flank Gene Annotation ##
gff <- readGFF("../references/gencode.v32.annotation.gff3.gz")
genes <- gff[gff$type=="gene",]

#look in 10kb window around gene for outliers
genes.gr <- as_granges(genes) %>% select(gene_id, gene_type, gene_name)
start(genes.gr) <- start(genes.gr) - 1e4
end(genes.gr) <- end(genes.gr) + 1e4

#### Genome table -- genes intersecting or nearby rare SVs ####
rare_svs <- fread("../../LongRead_Methylation/output/MERGED_SVs/U08.SV_rare_genotypes.long_table.tsv")
rare_svs$sample <- gsub(rare_svs$sample,pattern="\\d+_(\\w+)_sniffles", replacement="\\1")
rare_svs$SVTYPE <- gsub(rare_svs$SVTYPE,pattern="\\d+_(\\w+)",replacement="\\1")
rare_svs$end[rare_svs$end < rare_svs$start] <- rare_svs$start[rare_svs$end < rare_svs$start]
SV_outliers <- rare_svs %>% group_by(SVTYPE, sample) %>% summarize(n=n()) %>% ungroup %>% group_by(SVTYPE) %>% mutate(burden_z=scale(n)) %>% filter(burden_z>5) %>% pull(sample)
rare_svs %<>% filter(!sample %in% SV_outliers)
rare_svs %<>% filter(rare_svs$svlen > 50)

rare_svs.gr <- as_granges(rare_svs, seqnames=rare_svs$chrom)
rare_svs.gene_intersect <- rare_svs.gr %>% join_overlap_left(genes.gr)
rare_svs.gene_intersect <- rare_svs.gene_intersect[!is.na(rare_svs.gene_intersect$gene_id),]
rare_svs.gene_intersect$gene_id %<>% gsub(pattern="(ENSG\\d+)\\.\\d+", replacement="\\1")

rare_svs.gene_intersect <- rare_svs.gene_intersect[rare_svs.gene_intersect$gene_type=="protein_coding",]
rare_sv_genes <- unique(rare_svs.gene_intersect$gene_id)
rare_sv_samples <- unique(rare_svs.gene_intersect$sample)

rare_svs.gene_collapsed <- rare_svs.gene_intersect %>% group_by(sample, gene_id) %>% summarize(has_rare_sv=n()>0) %>% as.data.frame
genome_table <- expand.grid(sample=rare_sv_samples, gene_id=rare_sv_genes) %>% as.data.frame %>% left_join(rare_svs.gene_collapsed)
genome_table$has_rare_sv[is.na(genome_table$has_rare_sv)] <- 0

genome_table %>% group_by(sample) %>% summarize(n = sum(has_rare_sv)) %>% pull(n) %>% median

#### methylation table -- genes nearby or overlapplying methylation outliers -- ####
meth_outlier_regions <- fread("../METAFORA_output/METAFORA_methylation_outlier_regions.tissue_Blood.ALL_CHROM_COMBINED.haplotype_annotated.gene_track_annotated.bed")
global_outliers <- meth_outlier_regions %>% group_by(ID) %>% summarize(n=n()) %>% mutate(outlier_burden_z=scale(n)) %>% filter(outlier_burden_z>3) %>% pull(ID)
length(global_outliers)
meth_outlier_regions <- meth_outlier_regions[meth_outlier_regions$gene_type=="protein_coding",] #only protein coding genes

##remove global outliers
meth_outlier_regions <- meth_outlier_regions[!(meth_outlier_regions$ID %in% global_outliers),]
meth_out_genes <- unique(meth_outlier_regions$gene_id)
meth_out_samples <- unique(meth_outlier_regions$ID)

length(meth_out_genes)
length(meth_out_samples)

meth_outlier_regions %<>% dplyr::rename(sample=ID)
meth_out.gene_collapsed <- meth_outlier_regions %>% group_by(sample, gene_id) %>% summarize(has_methylation_outlier=n()>0) %>% as.data.frame
methylome_table <- expand.grid(sample=meth_out_samples, gene_id=meth_out_genes) %>% as.data.frame %>% left_join(meth_out.gene_collapsed)
methylome_table$has_methylation_outlier[is.na(methylome_table$has_methylation_outlier)] <- 0
methylome_table$gene_id <- gsub(methylome_table$gene_id, pattern="(ENSG\\d+)\\.\\d+", replacement="\\1")

methylome_table %>% group_by(sample) %>% summarize(n = sum(has_methylation_outlier)) %>% pull(n) %>% median

#### NEW ATAC (EPIOUT Data) With all 3 batches #### 
epiout <- fread("../../EpiOut/all_batches_shared_peak_set/all_batches_joint_called_reference_outliers.tsv")
epiout$padj %<>% pmax(1e-15)
epiout$zscore <- sign(epiout$l2fc) * qnorm(1-epiout$padj)

epiout_global_outliers <- epiout %>% group_by(sample) %>% summarize(n=dplyr::n()) %>% ungroup %>% mutate(globalZ=scale(n)) %>% filter(globalZ > 4) %>% pull(sample)
epiout %<>% filter(!sample %in% epiout_global_outliers)
epiout %<>% separate(col = peak, into=c("seqnames", "start", "end"), sep=":|-")
epiout$start %<>% as.numeric
epiout$end %<>% as.numeric
epiout$sample <- gsub(epiout$sample, pattern="(GSS\\d+)_.+", replacement="\\1")
length(unique(epiout$sample))

epiout_gr <- as_granges(epiout)
epiout.gene_intersect <- epiout_gr %>% join_overlap_left(genes.gr)

epiout.gene_intersect <- epiout.gene_intersect[!is.na(epiout.gene_intersect$gene_id),]
epiout.gene_intersect$gene_id %<>% gsub(pattern="(ENSG\\d+)\\.\\d+", replacement="\\1")

epiout.gene_collapsed <- epiout.gene_intersect %>% group_by(sample,gene_id) %>% summarize(has_outlier_peak=n() > 0) %>% as.data.frame
atac_out_genes <- unique(epiout.gene_intersect$gene_id)
atac_out_samples <- unique(epiout.gene_intersect$sample)
ATAC_table <- expand.grid(sample=atac_out_samples, gene_id=atac_out_genes) %>% as.data.frame %>% left_join(epiout.gene_collapsed)
ATAC_table$has_outlier_peak[is.na(ATAC_table$has_outlier_peak)] <- 0
ATAC_table %>% group_by(sample) %>% summarize(n = sum(has_outlier_peak)) %>% pull(n) %>% median

ATAC_table$sample <- gsub(ATAC_table$sample, pattern="(GSS\\d+)_.+", replacement="\\1")

#### RNA-SEQ ####
rna.zscores <- fread("../../gss_prospective_rnaseq/aggregate/expression/eoutliers_withUDN_blood.txt.gz")
length(unique(rna.zscores$sample))
rna_global_outliers <- rna.zscores %>% group_by(sample) %>% summarize(n= sum(abs(zscore) > 2), na.rm=T) %>% mutate(globalZ=scale(n)) %>% ungroup %>% filter(globalZ > 2) %>% pull(sample)
rna.zscores %<>% filter(!sample %in% rna_global_outliers)
rna.zscores %<>% dplyr::rename(gene_id=gene)

protein_coding_genes <- gsub(pattern="(ENSG\\d+)\\.\\d+", replacement="\\1", genes$gene_id[genes$gene_type=="protein_coding"])
rna.zscores %<>% filter(gene_id %in% protein_coding_genes)

expression_out_genes <- unique(rna.zscores %>% filter(abs(zscore) > 2) %>% pull(gene_id))
expression_out_samples <- unique(rna.zscores %>% filter(abs(zscore) > 2) %>% pull(sample))

length(expression_out_genes)
length(expression_out_samples)

rna.zscores.summarized <- rna.zscores %>% group_by(sample, gene_id) %>% summarize(has_expression_outlier=any(abs(zscore) > 2))

RNA_table <- expand.grid(sample=expression_out_samples, gene_id=expression_out_genes) %>% as.data.frame %>% left_join(rna.zscores.summarized)
RNA_table$has_expression_outlier[is.na(RNA_table$has_expression_outlier)] <- 0

RNA_table %>% group_by(sample) %>% summarize(n = sum(has_expression_outlier)) %>% pull(n) %>% median

### Protein Z-scores
prot.zscores <- fread("../../GSS_Multiomics/Olink_proteomics/GSS_Olink_proteomics.zscores.tsv")
length(unique(prot.zscores$SampleID))
length(unique(prot.zscores$OlinkID))
prot.zscores %<>% dplyr::rename(sample=SampleID)
prot.zscores %<>% filter(!sample %in% prot_global_outliers)

prot_out_genes <- unique(prot.zscores %>% filter(abs(zscore) > 2) %>% pull(gene_id))
prot_out_samples <- unique(prot.zscores %>% filter(abs(zscore) > 2) %>% pull(sample))
prot.zscores

head(prot.zscores)
prot.zscores.summarized <- prot.zscores %>% group_by(sample, gene_id) %>% summarize(has_protein_outlier=any(abs(zscore) > 2))

PROT_table <- expand.grid(sample=prot_out_samples, gene_id=prot_out_genes) %>% as.data.frame %>% left_join(prot.zscores.summarized)
PROT_table$has_protein_outlier[is.na(PROT_table$has_protein_outlier)] <- 0

PROT_table %>% group_by(sample) %>% summarize(n = sum(has_protein_outlier)) %>% pull(n) %>% median

length(prot_out_genes)
length(prot_out_samples)

# convert UDN IDs to GSS IDs if they are retrospective cases
retros <- fread("../../metadata/GSS_UDN_MAP.csv")

correct_ids <- function(ids) { 
    ids_chopped <- gsub(x=ids, pattern="_Blood|_Fibro", replacement="")
    ids_chopped[ids_chopped %in% retros$UDN] <- retros$GSS[match(ids_chopped[ids_chopped %in% retros$UDN], retros$UDN)]
    ids_chopped
}
genome_table$sample <- correct_ids(genome_table$sample)
methylome_table$sample <- correct_ids(methylome_table$sample)
ATAC_table$sample <- correct_ids(ATAC_table$sample)
RNA_table$sample <- correct_ids(RNA_table$sample)
PROT_table$sample <- correct_ids(PROT_table$sample)

# compute pair-wise cross ome enrichments!!
enrichment <- function(table1, table2) {
    sample_intersect <- intersect(unique(table1$sample), unique(table2$sample))
    gene_intersect <- intersect(unique(table1$gene_id), unique(table2$gene_id))
    cat(paste0("number samples overlap: ", length(sample_intersect), ", number geens overlap: ", length(gene_intersect), "\n"))
    table1 %<>% filter(sample %in% sample_intersect, gene_id %in% gene_intersect)
    table2 %<>% filter(sample %in% sample_intersect, gene_id %in% gene_intersect)

    overlap <- left_join(table1, table2)
    fisher_result <- fisher.test(table(overlap[,3], overlap[,4]))
    c(colnames(overlap)[3:4], fisher_result$estimate, fisher_result$p.value, fisher_result$conf.int)
}

enrichment_results <- data.frame(rbind( 
    enrichment(methylome_table, genome_table),
    enrichment(methylome_table, ATAC_table),
    enrichment(methylome_table, RNA_table),
    enrichment(methylome_table, PROT_table)))
colnames(enrichment_results) <- c("ome1", "ome2", "odds_ratio", "p_value", "lower", "upper")
assay_map <- list("has_expression_outlier" = "RNA", "has_methylation_outlier"="Methylation", "has_outlier_peak"="ATAC", "has_rare_sv"="Genome", "has_protein_outlier"="Protein")
enrichment_results$assay1 <- sapply(enrichment_results$ome1, function(x) assay_map[[x]])
enrichment_results$assay2 <- sapply(enrichment_results$ome2, function(x) assay_map[[x]])
enrichment_results$odds_ratio %<>% as.numeric()
enrichment_results$lower %<>% as.numeric()
enrichment_results$upper %<>% as.numeric()
enrichment_results$assay1 %<>% factor(levels=c("Genome", "Methylation", "ATAC", "RNA"))
enrichment_results$assay2 %<>% factor(levels=c("Protein", "RNA", "ATAC", "Genome"))
fwrite(enrichment_results,sep="\t",file="methylation_U08.odds_ratio.omics_overlap.txt")

ggplot(enrichment_results, aes(odds_ratio, assay2, color=assay2, xmin=lower, xmax=upper)) + geom_pointrange(size=.5, linewidth=2) + theme_minimal() + geom_vline(xintercept=1, linetype="dashed", color="grey") + theme(legend.position="none", text= element_text(size=18), panel.border=element_rect(fill=NA))  + xlab("odds ratio (enrichment)") + ylab("Molecular Outlier")  + scale_x_log10()
ggsave('../paper_plots/methylation_outlier_enrichment_results.point_ranges.pdf', width=5, height=5)

ATAC_meth_examples <- high_hap.gr %>% join_overlap_left(epiout_gr) %>% filter(ID==sample)
ATAC_meth_examples %>% filter(gene_name=="FAM193B")
epiout_gr %>% filter(sample=="GSS175297")
epiout_full <- fread("../../EpiOut/all_batches_shared_peak_set/all_batches_joint_called_reference_all_results.tsv") 
epiout_full$sample <- gsub(epiout_full$sample, pattern="(\\w+)_S.+", replacement="\\1")
set.seed(0)
atac.plt <- epiout_full %>% filter(peak=="chr5:177553743-177554946") %>% 
    ggplot(aes(l2fc, -log10(pval), shape=outlier, fill=sample=="GSS175297")) + geom_point(size=2,color="black",alpha=.7) + theme_minimal() + 
    scale_fill_manual(values=c("grey", "red")) + 
    scale_shape_manual(values=c(21, 23)) + 
    xlab("Peak log2 Fold Change") + ylab("-log10(p-value)") + theme(legend.position="none")
fam193b_rna <- fread("./FAM193b.expression_zscore.gregor.txt")  %>% melt %>% set_colnames(c("gene_id","sample","zscore")) %>% mutate(outlier=sample=="GSS175297")
rna.plt <- ggplot(fam193b_rna,aes(zscore,fill=outlier)) + geom_histogram(alpha=.7) + theme_minimal() + scale_fill_manual(values=c("grey", "red")) + theme(legend.position="none") +
     geom_vline(xintercept=fam193b_rna$zscore[fam193b_rna$sample=="GSS175297"],linetype="dashed",color="red") + xlab("expression z-score")

plot_grid(atac.plt, rna.plt, ncol=2)
ggsave("FAM193.mutiomics_outliers.pdf", width=5, height=3)
