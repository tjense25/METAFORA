library(data.table)
library(limma)
library(plyranges)
library(ggrepel)
library(tibble)
library(fastseg)
library(GenomicRanges)
library(tidyverse)
library(PCAtools)
library(pbmcapply)
library(Matrix)
library(magrittr)
library(cowplot)
library(matrixStats)
library(parallelly)
ncores <- availableCores()
options(mc.cores=ncores)
options(dplyr.summarise.inform = FALSE, dplyr.join.inform = FALSE)

segment_chrom <- function(meth.sample, this_chrom, segment_alpha=.01, min_seg_size = 10, median_seg_z = 2) {
  # segment deviance score profile
  meth.sample <- meth.sample[meth.sample$chrom == this_chrom,]
  
  index <- c(1,which(diff(meth.sample$start) > 1000)) #break into blocks of contiguous cpgs in at least 1000bp window
  last_start <- index[length(index)]
  last_end <- nrow(meth.sample) #number of cpgs on chrom
  block_start <- index[-length(index)]
  block_end <- index[-1] - 1
  block_start <- c(block_start, last_start)
  block_end <- c(block_end, last_end)

  blocks <- data.frame(start = block_start, end = block_end)
  blocks <- blocks[blocks$end - blocks$start >= min_seg_size,] #do not segment blocks with less than minimum seg_size contiguous cpgs

  cand.segs <-Reduce(rbind,lapply(1:nrow(blocks), function(i) {  
    tmp.meth <- meth.sample[blocks$start[i]:blocks$end[i],]
    as.data.frame(fastseg(tmp.meth$deviance_score, alpha = segment_alpha, 
            minSeg = min_seg_size, segMedianT = c(median_seg_z,-median_seg_z))) %>%
      mutate(start = tmp.meth$start[start], end=tmp.meth$start[end],
             seqnames=this_chrom)}))
  return(cand.segs)
}


MIN_Z_THRESH <- function(D) {0.76983 + 0.02709*D} #intercept optimized by simulation benchmark experiments
MAX_Z_THRESH <- function(D) {4.112 + 0.1395*D} #intercept and slope optimized by simulation benchmark experiments

segment_candidate_outliers <- function(pop_mean, betas, depths, this_sample, this_chrom, segment_alpha=.01, min_seg_size = 10, MAX_DEPTH=100) {
  betas.sample <- betas %>% select("chromosome","start",all_of(this_sample))
  colnames(betas.sample) <- c("chromosome","start","sample_beta")
  
  depth.sample <- depths %>% select("chromosome","start",all_of(this_sample))
  colnames(depth.sample) <- c("chromosome", "start", "sample_depth")
  D=median(depth.sample$sample_depth)
  
  meth.sample <- pop_mean %>% left_join(betas.sample,join_by(chromosome,start)) %>% left_join(depth.sample,join_by(chromosome,start)) %>% mutate(cpg_num = 1:n())
  meth.sample$sample_depth %<>% pmin(MAX_DEPTH) #cap depth at specified value to not inflate p-values for high-coverage samples
  MIN_Z=max(1,MIN_Z_THRESH(D))
  MAX_Z=MAX_Z_THRESH(D)

  #beta correction on sample
  meth.sample$sample_beta <- (meth.sample$sample_beta*meth.sample$sample_depth+1)/(meth.sample$sample_depth+2)
  
  # calculate deviance score
  meth.sample %<>% mutate(deviance_score = qnorm(pbeta(q=sample_beta, mean_beta*sample_depth, (1-mean_beta)*sample_depth)))
  meth.sample$deviance_score %<>% pmax(-MAX_Z)
  meth.sample$deviance_score %<>% pmin(MAX_Z)
  meth.sample <- meth.sample[!is.na(meth.sample$sample_beta),]

  # segment zscore profile
  cand.segs <- Reduce(rbind,lapply(c(this_chrom), function(x) {  
        segs <- segment_chrom(meth.sample, x, segment_alpha, min_seg_size, MIN_Z)
        if (nrow(segs) == 0) {return(NULL)}
        return(segs %>% mutate(ID=this_sample, seg_id=paste0(x,"_",this_sample,"_",1:nrow(segs))))
      }))
  cand.segs <- cand.segs[!is.na(cand.segs$`seg.mean`),] #remove cand segments that are NA over all cpgs
  cand.segs <- cand.segs[abs(cand.segs$seg.mean) > MIN_Z,]
  return(list("meth.sample"=meth.sample,"cand.segs"=cand.segs))
}

call_outliers <-function(cand.segs, betas, depth, sample_id, MIN_ABS_ZSCORE = 3, covariates=NULL) {
  betas.gr <- makeGRangesFromDataFrame(betas)
  cands.gr <- makeGRangesFromDataFrame(cand.segs, keep.extra.columns = T)
  
  betas.mat <-  as.matrix(betas[,4:ncol(betas)])
  depth.mat <- as.matrix(depth[,4:ncol(depth)])
  
  ol <- findOverlaps(cands.gr, betas.gr)
  # create Segment x CpG identity matrix to indicate which cpgs belong to which segment
  CpG_Identity <- sparseMatrix(i = queryHits(ol), j = subjectHits(ol), dims=c(length(cands.gr),nrow(betas)), x=1)
  
  betas.mat[is.na(betas.mat)] <- 0
  depth.mat[is.na(depth.mat)] <- 0
  # aggregate betas across each segment
  region_beta <- as.matrix(((CpG_Identity %*% (betas.mat*depth.mat))+1) / ((CpG_Identity %*% depth.mat)+2))
  cand.segs$pop_median <- rowMedians(region_beta,na.rm=T)

  colnames(region_beta) <- colnames(betas.mat)
  if (!is.null(covariates) && "Batch" %in% colnames(covariates)) {
      this_batch <- covariates$Batch[rownames(covariates) == sample_id]

      #create batch specific profile only if batch group has more than 10 samples
      if (sum(covariates$Batch == this_batch) > 10){
          cand.segs$pop_median <- rowMedians(matrix(region_beta[,covariates$Batch==this_batch],nrow=dim(region_beta)[1]),na.rm=T)
      }
  }
  cand.segs$delta <- region_beta[,sample_id] - cand.segs$pop_median

  B <- region_beta
  M <- log(B/(1-B))
  #M.scaled <- t(scale(t(M),center=T,scale=T))
  Batch=NULL
  if(!is.null(covariates) && "Batch" %in% colnames(covariates)) {
      Batch=covariates$Batch
      covariates <- covariates %>% select(-Batch)

      #if Batch only has 1 group or any Batch group has less than 4 individuals 
      #(not big enough sample size to reliable correct batch effects) set batch to NULL
      if((length(table(Batch)) < 2) || any(table(Batch) <= 5)) {
          Batch=NULL
      }
  }
  M.corrected <- removeBatchEffect(M, batch=Batch, covariates=covariates)
  zscores <- t(scale(t(M.corrected)))
  cand.segs$zscore <- zscores[,sample_id]

  zscores <- zscores[abs(cand.segs$zscore) > MIN_ABS_ZSCORE,]
  cand.segs <- cand.segs[abs(cand.segs$zscore) > MIN_ABS_ZSCORE,]
  if (nrow(cand.segs)==0) { return(list("outlier.segs"=NULL, "z.mat"=NULL)) }

  zscores <- matrix(zscores, nrow=nrow(cand.segs), ncol=ncol(betas.mat))
  colnames(zscores) <- colnames(betas.mat)
  rownames(zscores) <- cand.segs$seg_id
  return(list("outlier.segs"=cand.segs,"z.mat"=zscores))
}

outlier_pipeline <- function(betas,depths,this_chrom="chrT",MIN_ABS_ZSCORE=2,MIN_SEG_SIZE=5,MIN_ABS_DELTA=0.1,MAX_DEPTH=100) {
    # Read in data
    covariates <- NULL

    ## Correct pop mean for specific batch
    beta.mat <- as.matrix(betas[,4:ncol(betas)])
    depth.mat <- as.matrix(depths[,4:ncol(depths)])
    tmp.bmat <- beta.mat
    tmp.dmat <- depth.mat
    tmp.bmat[is.na(tmp.bmat)] <- 0
    tmp.dmat[is.na(tmp.dmat)] <- 0
    tmp.dmat[tmp.dmat > MAX_DEPTH] <- MAX_DEPTH

    #make population mean specific to Batch to make less susceptible to batch effects:
    pop_mean <- (rowSums(tmp.bmat*tmp.dmat) + rowSums(tmp.dmat > 0)) / (rowSums(tmp.dmat) + 2*rowSums(tmp.dmat > 0))
    pop_sd <- rowSds(tmp.bmat, na.rm=T)
    total_depth <- rowSums(tmp.dmat, na.rm=T)

    pop_mean <- data.table(betas[,1:3], total_depth, mean_beta=pop_mean, sd_beta=pop_sd)

    cand.outliers <- segment_candidate_outliers(pop_mean, betas, depths,this_sample="outlier_sample", this_chrom=this_chrom, min_seg_size=MIN_SEG_SIZE, MAX_DEPTH=MAX_DEPTH)
    meth.sample = cand.outliers[["meth.sample"]]
    cand.segs <- cand.outliers[["cand.segs"]]
    
    cand.segs
    if(nrow(cand.segs)==0||is.null(nrow(cand.segs))) {
        return(NULL)
    }
    # calculate region aggregated M values and call zscores across samples
    outliers <- call_outliers(cand.segs, betas, depths, sample_id="outlier_sample", MIN_ABS_ZSCORE=MIN_ABS_ZSCORE, covariates=covariates)
    outlier.segs <- outliers[["outlier.segs"]]
    return(outlier.segs)
}
 
simulate_outlier <- function(pop_mean_beta, N=30, D=30, M=20, pop_median=0.5, delta=0.3, noise=0.01, context_length=1500, plot=F) {
    random_starts <- sample(1:(nrow(pop_mean_beta)-context_length), 2)
    fake_profile <- c(pop_mean_beta$mean_beta[random_starts[1]:(random_starts[1]+context_length-1)], #left context
                      pop_median + noise*rnorm(M), #outlier region
                      pop_mean_beta$mean_beta[random_starts[2]:(random_starts[2]+context_length-1)]) # right context
    fake_profile[fake_profile < (1/D)] <- 1/D
    fake_profile[fake_profile > (1-(1/D))] <- 1-(1/D)
    cpgs=M+2*context_length
    depths <- matrix(data=rpois(N*cpgs, D), nrow=cpgs, ncol=N)
    betas <- apply(depths, 2, function(depth_col) rbeta(cpgs,shape1=fake_profile*depth_col,shape2=(1-fake_profile)*depth_col)+noise*rnorm(cpgs))
    
    outlier_profile <- fake_profile
    outlier_profile[(context_length+1):(context_length+M)] <- outlier_profile[(context_length+1):(context_length+M)] + delta
    outlier_profile[outlier_profile < (1/D)] <- 1/D
    outlier_profile[outlier_profile > (1-(1/D))] <- 1-(1/D)
    
    betas[,1] <- rbeta(cpgs,shape1=outlier_profile*depths[,1],shape2=(1-outlier_profile)*depths[,1])+noise*rnorm(cpgs)
    betas[betas < 0] <- 0
    betas[betas > 1] <- 1
    #  scale_color_manual(values=c("grey20", "red")) + annotate("rect",xmin=context_length+2,xmax=context_length+M+1,ymin=0,ymax=1, fill="red", alpha=.25)
    
    prefix <- data.frame(chromosome="chrT", start=1:cpgs, end=1:cpgs)
    colnames(betas) <- c("outlier_sample", paste0("control", 1:(N-1)))
    colnames(depths) <- colnames(betas)
    
    betas <- as.data.frame(cbind(prefix,betas))
    depths <- as.data.frame(cbind(prefix,depths))
    if(plot) {
        betas.df <- betas %>% pivot_longer(cols=-c("chromosome","start","end")) %>% mutate(outlier=name=="outlier_sample")
        ggplot(betas.df, aes(start,value,group=name,color=outlier,alpha=outlier,lineype=outlier)) + geom_line() + theme_minimal() +
          scale_alpha_manual(values=c(.2,1)) + scale_color_manual(values=c("black","red")) + scale_linetype_manual(values=c("dashed","solid")) +
          annotate("rect",xmin=context_length+2,xmax=context_length+M+1,ymin=0,ymax=1,fill="red",alpha=.25) +
          xlim(context_length-50,context_length+M+50)
        ggsave("outlier_plot.pdf")
    }
    return(list("betas"=betas,"depths"=depths))
}

evaluate_outlier <- function(outlier, context_length=1500, M) {
    if(is.null(outlier)||nrow(outlier)==0) {
        return(data.frame(tp=F,overlap=NA,segmean=NA,zscore=NA,fps=0,delta_hat=NA,pop_median_hat=NA,max_FP_size=NA,max_FP_delta=NA))
    }
    true_outlier <- makeGRangesFromDataFrame(data.frame(seqnames="chrT",start=context_length+1,end=context_length+M))
    o.gr <- makeGRangesFromDataFrame(outlier, keep.extra.columns = T)
    if(!any(overlapsAny(o.gr,true_outlier))) {
       return(data.frame(tp=F,overlap=NA,segmean=NA,zscore=NA,fps=nrow(outlier),delta_hat=NA,pop_median_hat=NA,max_FP_size=max(o.gr$`num.mark`),max_FP_delta=max(abs(o.gr$delta))))
    }
    overlap_props <- o.gr %>% join_overlap_intersect(true_outlier) %>% as.data.frame %>% 
        mutate(propA=width/`num.mark`, propB=width/M, 
               prop=pmin(propA,propB))  %>%
        select(seqnames,start,end,prop) %>% 
        makeGRangesFromDataFrame(keep.extra.columns =T)
    o.gr %<>% join_overlap_left(overlap_props)
    fps <- o.gr[is.na(o.gr$prop)]
    eval.out <- data.frame(tp=any(o.gr$prop >= .5,na.rm=T), 
               overlap=max(o.gr$prop,na.rm=T), 
               fps=sum(is.na(o.gr$prop)), 
               segmean=o.gr$`seg.mean`[which.max(o.gr$prop)],
               zscore=o.gr$zscore[which.max(o.gr$prop)],
               delta_hat=o.gr$delta[which.max(o.gr$prop)],
               pop_median_hat=o.gr$pop_median[which.max(o.gr$prop)],
               max_FP_size=NA,
               max_FP_delta=NA)
    if(eval.out$fps>0) {
        eval.out$max_FP_size=max(fps$`num.mark`)
        eval.out$max_FP_delta=max(abs(fps$delta))
    }
    return(eval.out)
}

pop_mean_beta <- fread("../../1000G_ONT_Consortium/METAFORA/METAFORA_output/Population_methylation.tissue_LCL/Population_methylation.tissue_LCL.chrom_chr1.population_mean_betas.tsv.gz")

# set up parameter grid to test
input_deltas <- c(.1,.2,.3,.4,.5,.6,.7,.8,.9)
N_inputs=c(15,30,50,100)
M_inputs=c(5,10,20,50,100)
D_inputs=c(5,10,15,20,25,30)
noise_inputs=c(0,0.01,0.025,0.05,0.1)

parameter_grid <- expand.grid(pop_median=input_deltas, delta=input_deltas, N=N_inputs, M=M_inputs, D=D_inputs, noise=noise_inputs)
parameter_grid$delta <- round((1-parameter_grid$delta) - parameter_grid$pop_median,1)
parameter_grid <- parameter_grid[parameter_grid$delta != 0,] 
parameter_grid$config_id <- 1:nrow(parameter_grid)

R=10 #how many replicates per parameter
C=1000
i=1
parameter_grid[i,]
eval_results <- pbmclapply(parameter_grid$config_id, function(i) {
        pop_median=parameter_grid[i,"pop_median"] 
        delta=parameter_grid[i,"delta"]
        N=parameter_grid[i,"N"]
        M=parameter_grid[i,"M"]
        D=parameter_grid[i,"D"]
        noise=parameter_grid[i,"noise"]
        set.seed(2)
        tmp.eval.df <- replicate(R, {
                simulated_out <- simulate_outlier(pop_mean_beta, pop_median=pop_median, delta=delta, N=N, M=M, D=D, noise=noise, context_length=C, plot=F)
                MIN_SEG_SIZE=min(10,M)
                outlier <- outlier_pipeline(simulated_out$betas,simulated_out$depths, MIN_SEG_SIZE=MIN_SEG_SIZE)
                evaluate_outlier(outlier, M=M, context_length=C)
        }, simplify=F)  %>% bind_rows
        tmp.eval.df$config_id <- i
        tmp.eval.df
   }, mc.cores=ncores, mc.preschedule=T)
saveRDS(eval_results, "eval_results.rds")
eval.df <- eval_results %>% bind_rows %>% left_join(parameter_grid)
fwrite(eval.df, file="Outlier_simulation.evaluation_metrics.txt", sep="\t")
