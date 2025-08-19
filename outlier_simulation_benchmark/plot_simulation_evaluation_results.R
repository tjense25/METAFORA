library(tidyverse)
library(data.table)
library(magrittr)
library(ggExtra)

eval.df <- fread("./Outlier_simulation.evaluation_metrics.txt")

power.df <- eval.df %>% group_by(abs(delta),D,noise) %>% summarize(power=mean(tp))
colnames(power.df) <- c("delta", "depth","noise", "power")
power.df$depth %<>% factor()

# power by delta and depth
ggplot(power.df%>%filter(noise==0.01), aes(delta,power,color=depth)) + geom_line() + geom_point(color="black",alpha=.3) + theme_minimal() + scale_color_manual(values=c("violet","violetred1","violetred2","violetred3","violetred4","red4")) + xlab("|delta| outlier deviance from population median") + ylab("power")
ggsave("./power_plot.pdf")

# power by noise and outlier length?
noise.df <- eval.df %>% filter(D==30) %>% group_by(M,noise) %>% summarize(power=mean(tp))
ggplot(noise.df, aes(factor(M),power,fill=factor(noise))) + geom_col(position="dodge", width=.8, color="black") + theme_minimal() + 
    scale_fill_brewer() +
    xlab("noise") + ylab("power")
ggsave("./noise_M_power_plot.pdf")

#look at power, false positives, z_score / segmeans
fp.df <- eval.df %>% filter(M>5,N==15) %>% group_by(D,noise) %>% summarize(fps=mean(fps)/2)
ggplot(fp.df, aes(factor(noise),fps,fill=factor(D))) + geom_col(position="dodge",width=.8, color="black") + theme_minimal() +
    scale_fill_manual(values=c("violet","violetred1","violetred2","violetred3","violetred4","red4")) +
    xlab("noise") + ylab("false positive rate per 1000CpGs")
ggsave("./false_positive_rate.barplot.pdf")

#eval.df %>% filter(tp) %>% 
#    ggplot(aes(factor(abs(delta)), abs(segmean), fill=factor(D))) + geom_violin(scale="width") + theme_minimal() 
gg <- eval.df %>% filter(fps > 0) %>% filter(noise==0)  %>% filter(config_id==1) %>% 
    ggplot(aes(max_FP_size, max_FP_delta)) + geom_jitter(alpha=.7) + theme_minimal() +
    geom_hline(yintercept=0.2,color="red", linetype="dashed") + 
    geom_vline(xintercept=10.5,color="red", linetype="dashed") + 
    xlab("false positive size (number of CpGs)") + 
    ylab("false positive effect size (outlier delta)")
ggsave(plot=gg,"./false_positive.statistics.pdf")

fps <- eval.df %>% filter(fps > 0) %>% group_by(max_FP_size>10,max_FP_delta>.2) %>% summarize(n=dplyr::n())
fps$prop <- fps$n/sum(fps$n)
fps

#look at z-score as a function of N
eval.df %>% filter(tp) %>% filter(abs(delta) %in% c(.2,.4,.6,.8))  %>%
    ggplot(aes(factor(abs(delta)), abs(zscore), fill=factor(N))) + geom_violin(scale="width") + theme_minimal()  + scale_fill_brewer(palette=2)
ggsave("./cohort_size.zscores.pdf")

eval.df %>% filter(abs(delta)==.1) %>% group_by(N,D,abs(delta)) %>% summarize(power=mean(tp), n=sum(tp)) %>% filter(D==30)
eval.df %>% filter(abs(delta)==.2) %>% group_by(N,D,abs(delta)) %>% summarize(power=mean(tp), n=sum(tp)) %>% filter(D==30)
eval.df %>% filter(abs(delta)==.3) %>% group_by(N,D,abs(delta)) %>% summarize(power=mean(tp), n=sum(tp)) %>% filter(D==30)
