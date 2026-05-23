library(tidyverse)
library(data.table)
library(limma)
library(edgeR)
library(matrixStats)

autosomes <- paste0("chr",1:22)
segment_blocks <- fread("./METAFORA_output/Chromosome_block.paralleliztion.bed")
segment_blocks$block

meth_segs <- fread("./METAFORA_output/Population_methylation.tissue_Blood/Meth_segments.tissue_Blood.segment_betas.bed") 
meth_depths <- fread("./METAFORA_output/Population_methylation.tissue_Blood/Meth_segments.tissue_Blood.segment_coverage.mat")
meth_depths <- meth_depths[meth_segs$chrom %in% autosomes,]
meth_segs <- meth_segs[meth_segs$chrom %in% autosomes,]
depth.mat <- as.matrix(meth_depths)
depths <- colMedians(depth.mat)
meth_segs.mat <- as.matrix(meth_segs[,8:ncol(meth_segs)])
rownames(meth_segs.mat) <- meth_segs$seg_id

covar <- as.data.frame(fread("./U08_GREGoR.covariates.txt"))
rownames(covar) <- covar$Sample_name

covar <- covar[colnames(meth_segs.mat),]
covar$depth <- depths

meth_segs.pacbio_vs_r10 <- meth_segs.mat[,covar$Batch!="ONT_R9"]
covar.pacbio_vs_r10 <- covar[covar$Batch!="ONT_R9",]

model.matrix <- as.matrix(cbind(as.numeric(covar.pacbio_vs_r10$Batch=="PacBio"), covar.pacbio_vs_r10$depth))
colnames(model.matrix) <- c("Technology", "Depth")

eps <- 1e-6
beta_clipped <- pmin(pmax(meth_segs.pacbio_vs_r10, eps), 1 - eps)
M <- log2(beta_clipped / (1 - beta_clipped))

lm.results <- lmFit(M, model.matrix)
fit <- eBayes(lm.results)
tech_effect  <- topTable(fit, coef="Technology", number=Inf)
tech_effect$segment_type = ifelse(beta.ont>.5, yes="hyper", no="hypo")
    

beta.ont<-rowMedians(meth_segs.pacbio_vs_r10[,covar.pacbio_vs_r10$Batch=="ONT_R10"],na.rm=T)
beta.pb<-rowMedians(meth_segs.pacbio_vs_r10[,covar.pacbio_vs_r10$Batch=="PacBio"],na.rm=T)
tech_effect$delta_beta_hat <- beta.pb-beta.ont
tech_effect$segment_type <- factor(ifelse(beta.ont>.8, yes="hyper", ifelse(beta.ont<.2, yes= "hypo", no="variable")), levels=c("hypo", "variable","hyper"))

tech_effect$sig <- abs(tech_effect$delta_beta_hat) > .1
ggplot(tech_effect, aes(delta_beta_hat, -log10(P.Value), color=sig)) + geom_point(alpha=.7) + facet_wrap(~segment_type, ncol=1)  + ylim(6,NA) + xlim(-.5,.5) + theme_minimal() + scale_color_manual(values=c("black", "red")) + 
    geom_vline(xintercept=-.1,color="grey20",linetype='dashed') + geom_vline(xintercept=.1, color="grey20",linetype="dashed") + theme(legend.position="none")
ggsave("tech_effect.volcano.pdf", height=6, width=5)

tmp.meth.seg <- data.frame(Batch=covar.pacbio_vs_r10$Batch, depth=covar.pacbio_vs_r10$depth, seg_beta=as.numeric(meth_segs.pacbio_vs_r10[500,])) +xlim(-.3,3)
ggplot(tmp.meth.seg, aes(Batch, depth,fill=Batch)) + geom_violin()
ggplot(tmp.meth.seg, aes(Batch, seg_beta,fill=Batch)) + geom_violin()
ggsave("tmp.violin.tech_effect.pdf")

volcanoplot(fit)
summary(fit$p.value)
lm.results$coefficients
