# JLW - 2026 - Permafrost Warming Experiments

################################################################################
# Packages ---------------------------------------------------------------------
################################################################################

library(dplyr)
library(data.table)
library(SPRING)
library(stringr) 
library(mclust)
library(gRodon)
library(ggplot2)
library(ggpubr)
library(ggridges)
library(ggcorrplot)
library(corrplot)

merge.easy <- function(df1,df2,key){
  df1 <- data.table(df1,key=key)
  df2 <- data.table(df2,key=key)
  return(unique(merge(df1,df2,all.x=TRUE,by=.EACHI,allow.cartesian=TRUE)))
}

################################################################################
# Load Data (see warming_analysis_prep_data.R) ---------------------------------
################################################################################

setwd("~/../../Documents/permafrost_warming/data")
load("warming_data_for_analysis.rda")

bs <- 14 #font base size


################################################################################
################################################################################
# Main Text Figures  -----------------------------------------------------------
################################################################################
################################################################################

################################################################################
# MAG Figure  ------------------------------------------------------------------
################################################################################

## IOC Clusters ----

p <- (colSums(mC$z)/sum(mC$z))
x_seq <- seq(.1,10,0.01)
cl_df <- data.frame(x=x_seq,
                    cl1=dnorm(x_seq,
                              mean=mC$parameters$mean[1],
                              sd=sqrt(mC$p$variance$sigmasq[1]))*p[1],
                    cl2=dnorm(x_seq,
                              mean=mC$parameters$mean[2],
                              sd=sqrt(mC$p$variance$sigmasq[1]))*p[2])
pC <- ggplot() + 
  geom_density(data=mag_df,aes(x=IoC),lwd=2,color=alpha("black", .4)) + 
  geom_density(data=permafrost,aes(x=IoC),lwd=2) + 
  geom_polygon(data=cl_df,aes(x=x, y=cl1, fill="Oligotroph Cluster"),
               alpha=0.5, color="black") +
  geom_polygon(data=cl_df,aes(x=x, y=cl2, fill="Copiotroph Cluster"), 
               alpha=0.5, color="black") +
  theme_pubclean(base_size = bs) + 
  xlab("Index of Copiotrophy")  + 
  ylab("Density") +
  theme(legend.position = "right") + 
  labs(fill="") +
  geom_vline(xintercept=c_h,lty=2,color="#1B9E77") +
  geom_vline(xintercept=c_l,lty=2,color="#7570B3") +
  annotate("text",x=c_h+0.25,
           y=0.3,
           angle=90,
           fontface = "bold",
           label="P(Copiotroph)>0.95",
           color="#1B9E77") +
  annotate("text",x=c_l-0.25,
           y=0.15,
           angle=90,
           fontface = "bold",
           label="P(Oligotroph)>0.95",
           color="#7570B3") +
  scale_fill_manual(values=c("#1B9E77","#7570B3")) +
  theme(legend.position = c(0.8,0.8))
pC

## IoC By Phylum ----

phyplot <- ggplot(mag_df,
                  aes(x=IoC, y=reorder(Phylum,IoC,FUN=median))) + 
  stat_density_ridges(
    geom = "density_ridges_gradient", calc_ecdf = TRUE,
    quantiles = 2, quantile_lines = F,
    jittered_points = F,
    position = position_points_jitter(width = 0.05, height = 0),
    scale=0.5,
    fill="darkgray",
    point_shape = '|', point_size = 3, point_alpha = 1, alpha = 0.5
  ) +
  geom_jitter(data=mag_df,
              aes(x=IoC, y=reorder(Phylum,IoC,FUN=median),fill=IoCc),
              alpha=0.5,size=3,pch=21,
              height=0.1,width=0) +
  scale_fill_brewer(palette="Dark2",direction=-1) +
  labs(fill="") + 
  theme_pubclean(base_size = bs) + 
  geom_vline(xintercept=c_h,lty=2,color="#1B9E77") +
  geom_vline(xintercept=c_l,lty=2,color="#7570B3") +
  xlab("Index of Copiotrophy") + 
  ylab("") +
  theme(legend.position = "none",
        legend.justification="left",
        legend.direction='vertical',
        legend.box.background = element_rect(colour = "black")) 


## Relative Abundance ----

ram_g$IoCc <- "Oligotroph"
ram_g$IoCc[ram_g$IoC>c_l] <- "Undefined"
ram_g$IoCc[ram_g$IoC>c_h] <- "Copiotroph"
ram_g$IoCc <- factor(ram_g$IoCc,
                     levels=c("Oligotroph",
                              "Undefined",
                              "Copiotroph"))
ram_g_permafrost$IoCc <- "Oligotroph"
ram_g_permafrost$IoCc[ram_g_permafrost$IoC>c_l] <- "Undefined"
ram_g_permafrost$IoCc[ram_g_permafrost$IoC>c_h] <- "Copiotroph"
ram_g_permafrost$IoCc <- factor(ram_g_permafrost$IoCc,
                     levels=c("Oligotroph",
                              "Undefined",
                              "Copiotroph"))
relab <- ggplot(ram_g %>% subset(value.pre+value.post > 0),
       aes(x=IoC,y=(value.post-value.pre),fill=IoCc)) +
  geom_point(alpha=.15,size=3,pch=21) +
  geom_point(data=ram_g_permafrost %>% subset(value.pre+value.post > 0),
             alpha=1,size=3,pch=21) +
  geom_hline(yintercept = 0,lty=2) +
  geom_smooth(method="gam",color=alpha("black", .4),fill="black",alpha=0.2) +
  geom_smooth(data=ram_g_permafrost %>% subset(value.pre+value.post > 0),
              method="gam",color="black",fill="black") +
  scale_fill_brewer(palette="Dark2",direction=-1) +
  theme_pubclean(base_size = bs)+
  labs(fill="") +
  theme(legend.position = c(0.8,0.1),
        legend.background = element_blank(),
        legend.key = element_blank(),
        legend.box.background = element_rect(colour = "black"),
        legend.title = element_blank()) +
  ylab("Avg. Change in Transformed Rel. Abundance") +
  xlab("Index of Copiotrophy") +
  geom_vline(xintercept=c_h,lty=2,color="#1B9E77") +
  geom_vline(xintercept=c_l,lty=2,color="#7570B3") +
  stat_cor(data=ram_g_permafrost %>% subset(value.pre+value.post > 0),
           label.y=1.95,
           inherit.aes=F,aes(x=IoC,y=(value.post-value.pre))) +
  stat_cor(alpha=0.5,label.y=1.85,
           inherit.aes=F,aes(x=IoC,y=(value.post-value.pre))) +
  guides(fill = guide_legend(override.aes = list(size = 5,alpha=1) ) )
relab

cor.test(ram_g_permafrost$delta[ram_g_permafrost$value.pre+ram_g_permafrost$value.post > 0],
         ram_g_permafrost$IoC[ram_g_permafrost$value.pre+ram_g_permafrost$value.post > 0])

cor.test(ram_g$delta[ram_g$value.pre+ram_g$value.post > 0],
         ram_g$IoC[ram_g$value.pre+ram_g$value.post > 0])

## Absolute Abundance Barplot ----

abam_g$IoCc <- "Oligotroph"
abam_g$IoCc[abam_g$IoC>c_l] <- "Undefined"
abam_g$IoCc[abam_g$IoC>c_h] <- "Copiotroph"
abam_g$IoCc <- factor(abam_g$IoCc,
                      levels=c("Oligotroph",
                               "Undefined",
                               "Copiotroph"))
binom_df <- as.data.frame(table(abam_g$IoCc,Increased=abam_g$delta>0)) %>%
  group_by(Var1) %>%
  summarise(n=sum(Freq),
            xI=sum(Freq[Increased=="TRUE"]),
            xD=sum(Freq[Increased=="FALSE"])) %>%
  mutate(pI=xI/n,pD=xD/n)
binom_df$LowerI <- mapply(FUN=binom.test,binom_df$xI,binom_df$n)["conf.int",] %>% lapply("[",1) %>% unlist()
binom_df$UpperI <- mapply(FUN=binom.test,binom_df$xI,binom_df$n)["conf.int",] %>% lapply("[",2) %>% unlist()
binom_df$LowerD <- mapply(FUN=binom.test,binom_df$xD,binom_df$n)["conf.int",] %>% lapply("[",1) %>% unlist()
binom_df$UpperD <- mapply(FUN=binom.test,binom_df$xD,binom_df$n)["conf.int",] %>% lapply("[",2) %>% unlist()

pAbabBar <- ggarrange(ggplot(binom_df, 
                             aes(x=Var1, y = pI, fill=Var1)) + 
                        geom_bar(stat="identity") +
                        geom_errorbar(aes(ymin=LowerI,ymax=UpperI,width=0.2)) +
                        theme_pubclean(base_size = bs) +
                        theme(legend.position = "none",
                              axis.text.x = element_text(angle = 60,hjust=1)) +
                        ylab("Proportion Increased Absolute Abundance") +
                        xlab("") +
                        ylim(0,0.9) +
                        scale_fill_brewer(palette="Dark2",direction=-1),
                      ggplot(binom_df, 
                             aes(x=Var1, y = pD, fill=Var1)) + 
                        geom_bar(stat="identity") +
                        geom_errorbar(aes(ymin=LowerD,ymax=UpperD,width=0.2)) +
                        theme_pubclean(base_size = bs) +
                        theme(legend.position = "none",
                              axis.text.x = element_text(angle = 60,hjust=1)) +
                        ylab("Proportion Decreased Absolute Abundance") +
                        xlab("") +
                        ylim(0,0.9) +
                        scale_fill_brewer(palette="Dark2",direction=-1))

## Put it Together ----

png("../figures/Figure_MAGs.png",width=400,height=250,units="mm",res=1000)
ggarrange(ggarrange(pC,
                    ggarrange(relab,
                              pAbabBar,
                              ncol=2,
                              widths=c(5,3),
                              hjust=0,
                              labels=c("(b)","(c)")),
                    nrow=2,
                    heights=c(2,3),
                    labels=c("(a)","")),
          phyplot,
          ncol=2,
          labels=c("","(d)"))
dev.off()

chisq.test(table(abam_g$IoCc,Increased=abam_g$delta>0))
chisq.test(table(ram_g$IoCc,Increased=ram_g$delta>0))

table(abam_g$IoCc,Increased=abam_g$delta>0)
table(ram_g$IoCc,Increased=ram_g$delta>0)

################################################################################
# Metagenome Figure  -----------------------------------------------------------
################################################################################


b <- runif(nrow(metagenome_df), -0.1, 0.1) #same jitter for both violin plots

## IoC Metagenomes ----


pIoC <- ggplot() +
  geom_violin(data=metagenome_df,
              aes(y=IoC,
                  x=as.numeric(pre_post_thaw=="post"),
                  group=pre_post_thaw),
              fill="gray") +
  geom_boxplot(data=metagenome_df,
               aes(y=IoC,
                   x=as.numeric(pre_post_thaw=="post"),
                   group=pre_post_thaw),
               width=0.2) +
  geom_point(data=metagenome_df,
             aes(y=IoC,
                 x=as.numeric(pre_post_thaw=="post")+b,
                 group=id,
                 size=gbp,
                 shape=gbp<1,
                 color=IoCc),
             alpha=0.75) +
  geom_line(data=metagenome_df,aes(y=IoC,x=as.numeric(pre_post_thaw=="post")+b,
                       group=id,color=IoCc),lwd=1.5) +
  theme_pubclean(base_size = 14) +
  theme(legend.position = "right") + 
  scale_x_continuous(breaks = c(0,1),
                     labels = c("Pre-Thaw", "Post-Thaw"),
                     expand = c(0,0))+
  scale_size_continuous(breaks=c(1,10),limits = c(0.1,35),range=c(3,10)) +
  scale_shape_discrete(solid=T) +
  xlab("") +
  labs(color="",fill="",size="Gbp",shape="<1 Gbp") +
  ylab("Avg. Index of Copiotrophy") +
  theme(legend.position = "right",
        axis.text.x = element_text(angle = 60,hjust=1)) +
  geom_hline(yintercept=c_h,lty=2,color="#1B9E77") +
  geom_hline(yintercept=c_l,lty=2,color="#7570B3") +
  scale_fill_brewer(palette="Dark2",direction=-1) +
  scale_color_brewer(palette="Dark2",direction=-1)
pIoC

## IoC Barplot ----

ioc_metagenomes <- as.data.frame(table(metagenome_df$IoCc,ppt=metagenome_df$pre_post_thaw))
ioc_metagenomes$ppt <- factor(ioc_metagenomes$ppt,
                              levels=c("pre",
                                       "post"))
pIoCbar <- ggplot(ioc_metagenomes, 
       aes(x=ppt, y = Freq, fill=forcats::fct_rev(Var1))) + 
  geom_bar(stat="identity") +
  theme_pubclean(base_size = 14) +
  theme(legend.position = "none") +
  ylab("Number of Samples") +
  xlab("") +
  scale_fill_brewer(palette="Dark2",direction=) +
  scale_x_discrete(labels=c("Pre-Thaw","Post-Thaw")) +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 60,hjust=1))

## MuMax Metagenome ----

pMuMax <- ggplot() +
  geom_violin(data=metagenome_df,aes(y=log(2)/d,x=as.numeric(pre_post_thaw=="post"),
                         group=pre_post_thaw),fill="gray") +
  geom_boxplot(data=metagenome_df,aes(y=log(2)/d,x=as.numeric(pre_post_thaw=="post"),
                          group=pre_post_thaw),
               width=0.2) +
  geom_point(data=metagenome_df,aes(y=log(2)/d,
                                    x=as.numeric(pre_post_thaw=="post")+b,
                                    group=id,
                                    size=gbp,
                                    shape=gbp<1,
                                    color=core),
             alpha=0.75) +
  geom_line(data=metagenome_df,aes(y=log(2)/d,x=as.numeric(pre_post_thaw=="post")+b,
                       group=id,color=core),lwd=1.5) +
  theme_pubclean(base_size = 14) +
  scale_color_brewer(palette = "Set1") +
  scale_fill_brewer(palette = "Set1") +
  scale_size_continuous(breaks=c(1,10),limits = c(0.1,35),range=c(3,10)) +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 60,hjust=1)) + 
  scale_x_continuous(breaks = c(0,1),
                     labels = c("Pre-Thaw", "Post-Thaw"),
                     expand = c(0,0))+
  xlab("") +
  labs(color="Permafrost\nCore",fill="Permafrost\nCore") +
  ylab("Avg. Maximum Growth Rate (1/Hours)") 

## Carbon Plot 

carbon <- data.frame(core=c("AT3",
                            "AT4",
                            "BEO1",
                            "BEO2B",
                            "FL1C4",
                            "FL2C3"),
                     carbon_percent=c(3.6,
                                      1.8,
                                      30.9,
                                      21.8,
                                      37.9,
                                      39.6))
x <- merge.easy(metagenome_df,carbon,key="core")                     

dupl.id <- x$id[duplicated(x$id)]
x.paired <- x[x$id %in% dupl.id,]
x.post <- x.paired %>% subset(pre_post_thaw=="post")
x.post$diff <- 100*(x.paired[x.paired$pre_post_thaw=="pre",]$d - 
  x.paired[x.paired$pre_post_thaw=="post",]$d)/x.paired[x.paired$pre_post_thaw=="pre",]$d

t.test(log(2)/x.paired[x.paired$pre_post_thaw=="pre",]$d,
       log(2)/x.paired[x.paired$pre_post_thaw=="post",]$d,
       paired=T)

t.test(x.paired[x.paired$pre_post_thaw=="pre",]$IoC,
       x.paired[x.paired$pre_post_thaw=="post",]$IoC,
       paired=T)

pCarbon <- ggplot(x.post,
                  aes(x=carbon_percent,y=-diff,fill=core,group=core)) +
  geom_boxplot(fill="white") +
  geom_point(size=5,pch=21,alpha=0.75) +
  scale_fill_brewer(palette = "Set1") +
  theme_pubclean(base_size = 12) +
  #scale_y_log10() +
  #scale_x_log10() +
  geom_smooth(fill="black",method="lm",color="black") +
  labs(fill="Permafrost\nCore") +
  ylab("% Change in Min. Doubling Time Over Experiment") +
  xlab("Organic C at Start of Experiment (%)") +
  theme(legend.position = "right") +
  geom_hline(yintercept=0,lty=2)

## Density IoC ----

pMdense <- ggplot(mag_df,
       aes(x=IoC)) +
  geom_density(lwd=2,alpha=0.25) +
  geom_jitter(aes(y=-0.035),width=0,height=0.03,alpha=0.25,size=0.75) +
  geom_density(data=metagenome_df,
               aes(color=pre_post_thaw,fill=pre_post_thaw),
               lwd=1.1,alpha=0.25) +
  geom_jitter(data=metagenome_df,
              aes(fill=pre_post_thaw,
                  y=-0.035),
              width=0,height=0.03,size=3,alpha=0.75,pch=21) +
  theme_pubclean(base_size = 14) +
  geom_vline(xintercept=c_h,lty=2,color="#1B9E77") +
  geom_vline(xintercept=c_l,lty=2,color="#7570B3") +
  theme(legend.position = "right") +
  scale_color_brewer(palette = "Accent") +
  scale_fill_brewer(palette = "Accent") +
  labs(fill="Thaw",color="Thaw") +
  xlab("Index of Copiotrophy") +
  ylab("Density")
pMdense 

## Put it Together ----

png("../figures/Figure_Metagenomes.png",width=250,height=250,units="mm",res=1000)
ggarrange(ggarrange(pMuMax,
                    pIoCbar,
                    pIoC,
                    widths=c(3,1,4),
                    ncol=3,
                    labels=c("(a)","(c)","(d)"),
                    hjust=0),
          ggarrange(pCarbon,
                    pMdense,
                    ncol=2,
                    widths=c(3,5),
                    labels=c("(b)","(e)")),
          nrow=2)
dev.off()


################################################################################
################################################################################
# Supplemental Figures  --------------------------------------------------------
################################################################################
################################################################################

################################################################################
# Coverage Histogram -----------------------------------------------------------
################################################################################

gbp1_cov <- ra_cov[substr(ra_cov$fastq_id,1,4) %in% substr(metagenome_df$fastq_id[metagenome_df$gbp>1],1,4),]

png("../figures/SFigure_Coverage.png",width=180,height=60,units="mm",res=1000)
ggarrange(ggplot(ra_cov,aes(x=RA)) +
  geom_histogram() +
  scale_x_log10(limits=c(0.01,0.7)) +
  theme_pubclean(base_size = 10) +
  xlab("Total Relative Abundance in MAGs") +
  geom_vline(xintercept=median(ra_cov$RA),lty=2,color="red"),
  ggplot(gbp1_cov,aes(x=RA)) +
    geom_histogram() +
    scale_x_log10(limits=c(0.01,0.7)) +
    theme_pubclean(base_size = 10) +
    xlab("Total Relative Abundance in MAGs") +
    geom_vline(xintercept=median(gbp1_cov$RA),lty=2,color="red"),
  labels=c("(a)","(b)"),
  hjust=0)
dev.off()

mean(ra_cov$RA)
mean(gbp1_cov$RA)





################################################################################
# Rel vs Absolute Abundance ----------------------------------------------------
################################################################################

# Relative Abundance ----

ram_g$IoCc <- "Oligotroph"
ram_g$IoCc[ram_g$IoC>c_l] <- "Undefined"
ram_g$IoCc[ram_g$IoC>c_h] <- "Copiotroph"
ram_g$IoCc <- factor(ram_g$IoCc,
                     levels=c("Oligotroph",
                              "Undefined",
                              "Copiotroph"))
relab <- ggplot(ram_g %>% subset(value.pre+value.post > 0),
                aes(x=IoC,y=(value.post-value.pre),fill=IoCc)) +
  geom_point(alpha=.25,size=3,pch=21) +
  geom_hline(yintercept = 0,lty=2) +
  geom_smooth(method="gam",color="black",fill="black") +
  scale_fill_brewer(palette="Dark2",direction=-1) +
  theme_pubclean(base_size = 12)+
  labs(fill="") +
  theme(legend.position = "left",
        legend.background = element_blank(),
        legend.key = element_blank(),
        legend.box.background = element_rect(colour = "black"),
        legend.title = element_blank()) +
  ylab("Avg. Change in Transformed Rel. Abundance") +
  xlab("Index of Copiotrophy") +
  geom_vline(xintercept=c_h,lty=2) +
  geom_vline(xintercept=c_l,lty=2,lwd=1.2) +
  guides(fill = guide_legend(override.aes = list(size = 5,alpha=1) ) )

## Relative Abundance Boxplots ----

relab_h <- ggplot(ram_g %>% subset(value.pre+value.post > 0),
                  aes(x=IoCc,y=(value.post-value.pre),fill=IoCc)) +
  geom_violin(fill="gray",width=1) +
  geom_boxplot(width=0.1,fill="white",aes(color=IoCc)) +
  # geom_jitter(pch=21,width=0.1,alpha=0.5,size=2) +
  scale_fill_brewer(palette="Dark2",direction=-1) +
  scale_color_brewer(palette="Dark2",direction=-1) +
  theme_pubclean(base_size = 12)+
  labs(fill="Predicted\nCopiotrophic?") +
  theme(legend.position = "none") +
  ylab(expression("Avg. "*Delta*" in Transformed Rel. Ab.")) +
  xlab("")
pb <- ggplot(ram_g,aes(x=IoCc,y=(value.post-value.pre))) +
  geom_boxplot(width=0.6)
dat_1 <- ggplot_build(pb)$data[[1]]
dat_2 <- ggplot_build(relab_h)$data[[2]]
dat <- data.frame(xmin=c(unique(dat_1$xmin)),
                  xmax=c(unique(dat_1$xmax)),
                  middle=c(dat_2$middle),
                  IoCc=factor(c("Oligotroph",
                                "Undefined",
                                "Copiotroph"),
                              levels=c("Oligotroph",
                                       "Undefined",
                                       "Copiotroph")))
relab_h2  <- relab_h + geom_segment(data=dat, aes(x=xmin, xend=xmax, 
                                                  y=middle, yend=middle,
                                                  color=IoCc), 
                                    size=1,lty=1) +
  theme(axis.text.x = element_text(angle = 45, hjust=1)) +
  geom_hline(yintercept=0,lty=2) + 
  stat_pwc() + 
  stat_anova_test(label.y=2.75,
                  label.x=1.1)

## Absolute Abundance ----

abam_g$IoCc <- "Oligotroph"
abam_g$IoCc[abam_g$IoC>c_l] <- "Undefined"
abam_g$IoCc[abam_g$IoC>c_h] <- "Copiotroph"
abam_g$IoCc <- factor(abam_g$IoCc,
                      levels=c("Oligotroph",
                               "Undefined",
                               "Copiotroph"))
abab <- ggplot(abam_g,
               aes(x=IoC,y=(deltarel+1)/2,fill=IoCc)) +
  geom_point(alpha=.5,size=5,pch=21) +
  geom_hline(yintercept = 0.5,lty=2) +
  geom_smooth(method="glm",method.args=list(family="binomial"),color="black",fill="black") +
  scale_fill_brewer(palette="Dark2",direction=-1) +
  theme_pubclean(base_size = 10)+
  labs(fill="") +
  theme(legend.position = "none") +
  ylab("Rel. Change in Absolute Abundance: ((Start-End)/(Start+End)+1)/2") +
  xlab("Index of Copiotrophy")

## Absolute Abundance Barplots ----

binom_df <- as.data.frame(table(abam_g$IoCc,Increased=abam_g$delta>0)) %>%
  group_by(Var1) %>%
  summarise(n=sum(Freq),
            xI=sum(Freq[Increased=="TRUE"]),
            xD=sum(Freq[Increased=="FALSE"])) %>%
  mutate(pI=xI/n,pD=xD/n)
binom_df$LowerI <- mapply(FUN=binom.test,binom_df$xI,binom_df$n)["conf.int",] %>% lapply("[",1) %>% unlist()
binom_df$UpperI <- mapply(FUN=binom.test,binom_df$xI,binom_df$n)["conf.int",] %>% lapply("[",2) %>% unlist()
binom_df$LowerD <- mapply(FUN=binom.test,binom_df$xD,binom_df$n)["conf.int",] %>% lapply("[",1) %>% unlist()
binom_df$UpperD <- mapply(FUN=binom.test,binom_df$xD,binom_df$n)["conf.int",] %>% lapply("[",2) %>% unlist()
pAbabBar <- ggarrange(ggplot(binom_df, 
                             aes(x=Var1, y = pI, fill=Var1)) + 
                        geom_bar(stat="identity") +
                        geom_errorbar(aes(ymin=LowerI,ymax=UpperI,width=0.2)) +
                        theme_pubclean(base_size = 12) +
                        theme(legend.position = "none",
                              axis.text.x = element_text(angle = 60,hjust=1)) +
                        ylab("Proportion Increased Absolute Abundance") +
                        xlab("") +
                        ylim(0,0.9) +
                        scale_fill_brewer(palette="Dark2",direction=-1),
                      ggplot(binom_df, 
                             aes(x=Var1, y = pD, fill=Var1)) + 
                        geom_bar(stat="identity") +
                        geom_errorbar(aes(ymin=LowerD,ymax=UpperD,width=0.2)) +
                        theme_pubclean(base_size = 12) +
                        theme(legend.position = "none",
                              axis.text.x = element_text(angle = 60,hjust=1)) +
                        ylab("Proportion Decreased Absolute Abundance") +
                        xlab("") +
                        ylim(0,0.9) +
                        scale_fill_brewer(palette="Dark2",direction=-1))

## Put it Together ----

png("../figures/SFigure_AbsoluteAb.png",width=200,height=250,units="mm",res=1000)
ggarrange(ggarrange(relab,
                    relab_h2,
                    ncol=2,
                    labels=c("(a)","(b)"),
                    hjust=0,
                    widths=c(3,1)),
          ggarrange(abab,
                    pAbabBar,
                    ncol=2,
                    labels=c("(c)","(d)"),
                    hjust=0,
                    widths=c(3,2)),
          nrow=2)
dev.off()

################################################################################
# MuMax  -----------------------------------------------------------------------
################################################################################


pMuRel <- ggplot(ram_g %>% subset(value.pre+value.post > 0),
       aes(x=log(2)/d,y=(value.post-value.pre))) +
  geom_point(alpha=.25,size=2,pch=19) +
  geom_hline(yintercept = 0,lty=2) +
  geom_smooth(method="lm",color="black",fill="black") +
  theme_pubclean(base_size = 8)+
  labs(fill="") +
  theme(legend.position = "left",
        legend.background = element_blank(),
        legend.key = element_blank(),
        legend.box.background = element_rect(colour = "black"),
        legend.title = element_blank()) +
  ylab("Avg. Change in Transformed Rel. Abundance") +
  xlab("Predicted Max. Growth Rate") +
  scale_x_log10()

pMuAb <- ggplot(abam_g,
       aes(x=log(2)/d,y=(deltarel+1)/2)) +
  geom_point(alpha=.5,size=2,pch=19) +
  geom_hline(yintercept = 0.5,lty=2) +
  geom_smooth(method="glm",method.args=list(family="binomial"),color="black",fill="black") +
  scale_fill_brewer(palette="Dark2",direction=-1) +
  theme_pubclean(base_size = 8)+
  labs(fill="") +
  theme(legend.position = "none") +
  ylab("Rel. Change in Absolute Abund.: ((Start-End)/(Start+End)+1)/2") +
  xlab("Predicted Max. Growth Rate") +
  scale_x_log10()

png("../figures/SFigure_mumax.png",width=200,height=100,units="mm",res=1000)
ggarrange(pMuRel,pMuAb,ncol=2,labels=c("(a)","(b)"))
dev.off()

################################################################################
# IoC Component ----------------------------------------------------------------
################################################################################

## IoC Dimensions ----

pIoC <- ggplot(mag_df,aes(x=-dCUB,y=nCAZy,fill=IoCc,size=nGenes)) +
  geom_point(alpha=.5,pch=21) +
  scale_fill_brewer(palette="Dark2",direction=-1) +
  theme_pubclean(base_size = 14)+
  labs(fill="",size="Number of\nGenes") +
  scale_y_log10() +
  theme(legend.position = "right") +
  geom_vline(xintercept = 0.08,lty=2) +
  ylab("Number of CAZymes") +
  xlab("Codon Usage Bias (-dCUB)") + 
  guides(fill = guide_legend(override.aes = list(size=5)))

## Source Breakdown ----

pSource <- ggplot(mag_df,aes(x=-dCUB,y=nCAZy,fill=Source,size=nGenes)) +
  geom_point(alpha=.5,pch=21) +
  scale_fill_brewer(palette="Set1",direction=1) +
  theme_pubclean(base_size = 14)+
  labs(fill="",size="Number of\nGenes") +
  scale_y_log10() +
  theme(legend.position = "right") +
  geom_vline(xintercept = 0.08,lty=2) +
  ylab("Number of CAZymes") +
  xlab("Codon Usage Bias (-dCUB)") + 
  guides(fill = guide_legend(override.aes = list(size=5))) +
  facet_wrap(~Source)

## Correlation between components ----

x <- mag_df %>%
  mutate(RelNumCAZymes=rCAZy,
         CUB=-dCUB,
         NumGenes=nGenes,
         IoC=IoC) %>%
  subset(select=c(RelNumCAZymes,CUB,NumGenes,IoC))
M <- cor(x)
Mp <- cor.mtest(x,conf.level=0.99)

pCor <- ggarrange(ggcorrplot(M, hc.order = TRUE,
                             type = "lower", 
                             p.mat = Mp$p,lab=T),
                  ggplot(x,
                         aes(x=NumGenes,y=RelNumCAZymes)) +
                    geom_point() +
                    scale_x_log10() +
                    scale_y_log10() +
                    geom_smooth(method="lm") +
                    ylab("Rel. # CAZy") +
                    xlab("# Genes") +
                    theme_pubclean(),
                  ncol=2,
                  labels=c("(c)","(d)"),
                  hjust = 0)

## Put it Together ----

png("../figures/SFigure_IoC_components.png",width=250,height=250,units="mm",res=1000)
ggarrange(pIoC,
          pSource,
          pCor,
          nrow=3,
          heights=c(3,3,2),labels=c("(a)","(b)",""))
dev.off()


################################################################################
# IoC Components Metagenome ----------------------------------------------------
################################################################################

png("../figures/SFigure_IoC_components_metagenome.png",width=200,height=150,units="mm",res=1000)
ggplot(mag_df %>% subset(Source=="Present Study"),
       aes(x=-dCUB,y=nCAZy,size=nGenes)) +
  geom_point(alpha=.5,pch=19) +
  geom_point(data=metagenome_df,aes(fill=IoCc),alpha=0.5,pch=22) +
  theme_pubclean(base_size = 14)+
  labs(fill="",size="Number of\nGenes\n(Per Genome)") +
  scale_y_log10() +
  theme(legend.position = "right") +
  geom_vline(xintercept = 0.08,lty=2) +
  scale_fill_brewer(palette="Dark2",direction=-1) +
  ylab("Number of CAZymes (Per Genome)") +
  xlab("Codon Usage Bias (-dCUB)") + 
  guides(fill = guide_legend(override.aes = list(size=5)))
dev.off()

################################################################################
# Completeness -----------------------------------------------------------------
################################################################################

complete_50 <- mag_df$Name[mag_df$Completeness>50]
complete_90 <- mag_df$Name[mag_df$Completeness>90]

## 50% Completeness ----

rel50 <- ggplot(ram_g %>% 
         subset(value.pre+value.post > 0) %>%
         subset(Genome %in% complete_50),
       aes(x=IoC,y=(value.post-value.pre),fill=IoCc)) +
  geom_point(alpha=.25,size=3,pch=21) +
  geom_hline(yintercept = 0,lty=2) +
  geom_smooth(method="gam",color="black",fill="black") +
  scale_fill_brewer(palette="Dark2",direction=-1) +
  theme_pubclean(base_size = 12)+
  labs(fill="") +
  theme(legend.position = c(0.8,0.1),
        legend.background = element_blank(),
        legend.key = element_blank(),
        legend.box.background = element_rect(colour = "black"),
        legend.title = element_blank()) +
  ylab("Avg. Change in Transformed Rel. Abundance") +
  xlab("Index of Copiotrophy") +
  geom_vline(xintercept=c_h,lty=2) +
  geom_vline(xintercept=c_l,lty=2) +
  guides(fill = guide_legend(override.aes = list(size = 5,alpha=1) ) ) +
  stat_cor(inherit.aes=F,aes(x=IoC,y=(value.post-value.pre))) +
  ggtitle(">50% Completeness") 

binom50_df <- as.data.frame(table(abam_g$IoCc[abam_g$Genome %in% complete_50],
                                Increased=abam_g$delta[abam_g$Genome %in% complete_50]>0)) %>%
  group_by(Var1) %>%
  summarise(n=sum(Freq),
            xI=sum(Freq[Increased=="TRUE"]),
            xD=sum(Freq[Increased=="FALSE"])) %>%
  mutate(pI=xI/n,pD=xD/n)
binom50_df$LowerI <- mapply(FUN=binom.test,binom50_df$xI,binom50_df$n)["conf.int",] %>% lapply("[",1) %>% unlist()
binom50_df$UpperI <- mapply(FUN=binom.test,binom50_df$xI,binom50_df$n)["conf.int",] %>% lapply("[",2) %>% unlist()
binom50_df$LowerD <- mapply(FUN=binom.test,binom50_df$xD,binom50_df$n)["conf.int",] %>% lapply("[",1) %>% unlist()
binom50_df$UpperD <- mapply(FUN=binom.test,binom50_df$xD,binom50_df$n)["conf.int",] %>% lapply("[",2) %>% unlist()
ab50 <- ggarrange(ggplot(binom50_df, 
                             aes(x=Var1, y = pI, fill=Var1)) + 
                        geom_bar(stat="identity") +
                        geom_errorbar(aes(ymin=LowerI,ymax=UpperI,width=0.2)) +
                        theme_pubclean(base_size = 12) +
                        theme(legend.position = "none",
                              axis.text.x = element_text(angle = 60,hjust=1)) +
                        ylab("Proportion Increased Absolute Abundance") +
                        xlab("") +
                        ylim(0,1) +
                        scale_fill_brewer(palette="Dark2",direction=-1),
                      ggplot(binom_df, 
                             aes(x=Var1, y = pD, fill=Var1)) + 
                        geom_bar(stat="identity") +
                        geom_errorbar(aes(ymin=LowerD,ymax=UpperD,width=0.2)) +
                        theme_pubclean(base_size = 12) +
                        theme(legend.position = "none",
                              axis.text.x = element_text(angle = 60,hjust=1)) +
                        ylab("Proportion Decreased Absolute Abundance") +
                        xlab("") +
                        ylim(0,1) +
                        scale_fill_brewer(palette="Dark2",direction=-1))


## 90% Completeness ----

rel90 <- ggplot(ram_g %>% 
                  subset(value.pre+value.post > 0) %>%
                  subset(Genome %in% complete_90),
                aes(x=IoC,y=(value.post-value.pre),fill=IoCc)) +
  geom_point(alpha=.25,size=3,pch=21) +
  geom_hline(yintercept = 0,lty=2) +
  geom_smooth(method="gam",color="black",fill="black") +
  scale_fill_brewer(palette="Dark2",direction=-1) +
  theme_pubclean(base_size = 12)+
  labs(fill="") +
  theme(legend.position = "none",
        legend.background = element_blank(),
        legend.key = element_blank(),
        legend.box.background = element_rect(colour = "black"),
        legend.title = element_blank()) +
  ylab("Avg. Change in Transformed Rel. Abundance") +
  xlab("Index of Copiotrophy") +
  geom_vline(xintercept=c_h,lty=2) +
  geom_vline(xintercept=c_l,lty=2) +
  guides(fill = guide_legend(override.aes = list(size = 5,alpha=1) ) )+
  stat_cor(inherit.aes=F,aes(x=IoC,y=(value.post-value.pre))) +
  ggtitle(">90% Completeness")

binom90_df <- as.data.frame(table(abam_g$IoCc[abam_g$Genome %in% complete_90],
                                Increased=abam_g$delta[abam_g$Genome %in% complete_90]>0)) %>%
  group_by(Var1) %>%
  summarise(n=sum(Freq),
            xI=sum(Freq[Increased=="TRUE"]),
            xD=sum(Freq[Increased=="FALSE"])) %>%
  mutate(pI=xI/n,pD=xD/n)
binom90_df$LowerI <- mapply(FUN=binom.test,binom90_df$xI,binom90_df$n)["conf.int",] %>% lapply("[",1) %>% unlist()
binom90_df$UpperI <- mapply(FUN=binom.test,binom90_df$xI,binom90_df$n)["conf.int",] %>% lapply("[",2) %>% unlist()
binom90_df$LowerD <- mapply(FUN=binom.test,binom90_df$xD,binom90_df$n)["conf.int",] %>% lapply("[",1) %>% unlist()
binom90_df$UpperD <- mapply(FUN=binom.test,binom90_df$xD,binom90_df$n)["conf.int",] %>% lapply("[",2) %>% unlist()
ab90 <- ggarrange(ggplot(binom90_df, 
                         aes(x=Var1, y = pI, fill=Var1)) + 
                    geom_bar(stat="identity") +
                    geom_errorbar(aes(ymin=LowerI,ymax=UpperI,width=0.2)) +
                    theme_pubclean(base_size = 12) +
                    theme(legend.position = "none",
                          axis.text.x = element_text(angle = 60,hjust=1)) +
                    ylab("Proportion Increased Absolute Abundance") +
                    xlab("") +
                    ylim(0,1) +
                    scale_fill_brewer(palette="Dark2",direction=-1),
                  ggplot(binom_df, 
                         aes(x=Var1, y = pD, fill=Var1)) + 
                    geom_bar(stat="identity") +
                    geom_errorbar(aes(ymin=LowerD,ymax=UpperD,width=0.2)) +
                    theme_pubclean(base_size = 12) +
                    theme(legend.position = "none",
                          axis.text.x = element_text(angle = 60,hjust=1)) +
                    ylab("Proportion Decreased Absolute Abundance") +
                    xlab("") +
                    ylim(0,1) +
                    scale_fill_brewer(palette="Dark2",direction=-1))

## Put it Together ----

png("../figures/SFigure_completeness.png",width=200,height=250,units="mm",res=1000)
ggarrange(ggarrange(rel50,
                    ab50,
                    ncol=2,
                    widths=c(2,2),
                    hjust=0,
                    labels=c("(a)","(b)")),
          ggarrange(rel90,
                    ab90,
                    ncol=2,
                    widths=c(2,2),
                    hjust=0,
                    labels=c("(c)","(d)")),
          nrow=2)
dev.off()

################################################################################
# Carbon Vs MuMax/IoC
################################################################################

carbon <- data.frame(core=c("AT3",
                            "AT4",
                            "BEO1",
                            "BEO2B",
                            "FL1C4",
                            "FL2C3"),
                     carbon_percent=c(3.6,
                                      1.8,
                                      30.9,
                                      21.8,
                                      37.9,
                                      39.6))
x <- merge.easy(metagenome_df,carbon,key="core") 

pCioc <- ggplot(x,
       aes(x=carbon_percent,
           y=IoC,
           color=core,
           fill=core,
           group=core,
           shape=pre_post_thaw)) +
  geom_boxplot(color="black") +
  geom_point(size=5,color="black",alpha=0.75) +
  geom_point(size=3,alpha=0.75) +
  scale_color_brewer(palette = "Set1") +
  scale_fill_brewer(palette = "Set1") +
  scale_shape(solid=T) +
  theme_pubclean(base_size=12) +
  theme(legend.position = "right") +
  labs(shape="Thaw",color="Permafrost\nCore",fill="Permafrost\nCore") +
  ylab("Avg. Index of Copiotrophy") +
  xlab("Organic C at Start of Experiment (%)")

pCmu <- ggplot(x,
       aes(x=carbon_percent,
           y=log(2)/d,
           color=core,
           fill=core,
           group=core,
           shape=pre_post_thaw)) +
  geom_boxplot(color="black") +
  geom_point(size=5,color="black",alpha=0.75) +
  geom_point(size=3,alpha=0.75) +
  scale_color_brewer(palette = "Set1") +
  scale_fill_brewer(palette = "Set1") +
  scale_shape(solid=T) +
  theme_pubclean(base_size=12) +
  theme(legend.position = "none") +
  labs(shape="Thaw",color="Permafrost\nCore",fill="Permafrost\nCore") +
  ylab("Index of Copiotrophy") +
  xlab("Organic C at Start of Experiment (%)") +
  ylab("Predicted Avg. Maximum Growth Rate")

png("../figures/SFigure_carbon.png",width=220,height=100,units="mm",res=1000)
ggarrange(pCmu,pCioc,ncol=2,widths=c(3,4),labels=c("(a)","(b)"),hjust=0)
dev.off()

################################################################################
# MAG Quality
################################################################################


png("../figures/SFigure_mag_qual.png",width=200,height=150,units="mm",res=1000)
ggplot(mag_df,
       aes(x=Completeness,y=nHE,color=Contamination)) +
  geom_point(alpha=0.75,size=5) +
  scale_color_viridis_c() +
  theme_pubclean() +
  theme(legend.position = "right") +
  ylab("Number of Ribosomal Proteins")
dev.off()



################################################################################
# RelAb/AbAb Per Experiment
################################################################################

ram$IoCc <- "Oligotroph"
ram$IoCc[ram$IoC>c_l] <- "Undefined"
ram$IoCc[ram$IoC>c_h] <- "Copiotroph"
ram$IoCc <- factor(ram$IoCc,
                     levels=c("Oligotroph",
                              "Undefined",
                              "Copiotroph"))
png("../figures/SFigure_by_id.png",width=200,height=250,units="mm",res=1000)
ggplot(ram %>% subset(delta!=0),
       aes(x=IoCc,
           y=delta,
           fill=IoCc)) + 
  geom_boxplot(aes(color=IoCc),fill="white") +
  geom_jitter(alpha=.25,size=3,pch=21,width=0.1) +
  geom_hline(yintercept = 0,lty=2) +
  scale_fill_brewer(palette="Dark2",direction=-1) +
  scale_color_brewer(palette="Dark2",direction=-1) +
  theme_bw()+
  labs(fill="",color="") +
  theme(legend.position = "none",
        legend.background = element_blank(),
        legend.key = element_blank(),
        legend.box.background = element_rect(colour = "black"),
        legend.title = element_blank()) +
  ylab("Avg. Change in Transformed Rel. Abundance") +
  xlab("") +
  guides(fill = guide_legend(override.aes = list(size = 5,alpha=1) ) ) +
  facet_wrap(~id) +
  stat_anova_test() +
  theme(legend.position = "right",
        axis.text.x = element_text(angle = 60,hjust=1))
dev.off()


png("../figures/SFigure_by_id_relab.png",width=200,height=250,units="mm",res=1000)
ggplot(ram %>% subset(delta!=0),
       aes(x=value.pre,
           y=value.post,
           fill=IoCc,
           size=IoC)) + 
  geom_point(alpha=.5,pch=21,width=0.1) +
  scale_fill_brewer(palette="Dark2",direction=-1) +
  scale_color_brewer(palette="Dark2",direction=-1) +
  theme_pubclean(base_size = 12)+
  labs(fill="",color="") +
  theme(legend.position = "none",
        legend.background = element_blank(),
        legend.key = element_blank(),
        legend.box.background = element_rect(colour = "black"),
        legend.title = element_blank()) +
  guides(fill = guide_legend(override.aes = list(size = 5,alpha=1) ) ) +
  facet_wrap(~id) +
  geom_abline(slope=1,intercept=0,lty=2) +
  xlab("Transformed Relative Abundance Pre Thaw") +
  ylab("Transformed Relative Abundance Post Thaw") +
  theme_bw() +
  stat_anova_test() +
  theme(legend.position = "right",
        axis.text.x = element_text(angle = 60,hjust=1))
dev.off()

abam$IoCc <- "Oligotroph"
abam$IoCc[abam$IoC>c_l] <- "Undefined"
abam$IoCc[abam$IoC>c_h] <- "Copiotroph"
abam$IoCc <- factor(abam$IoCc,
                   levels=c("Oligotroph",
                            "Undefined",
                            "Copiotroph"))
png("../figures/SFigure_by_id_abab.png",width=200,height=250,units="mm",res=1000)
ggplot(abam %>% subset(delta!=0),
       aes(x=Absolute_16S_copy_gsoil.pre,
           y=Absolute_16S_copy_gsoil.post,
           fill=IoCc,
           size=IoC)) + 
  geom_point(alpha=.5,pch=21,width=0.1) +
  scale_fill_brewer(palette="Dark2",direction=-1) +
  scale_color_brewer(palette="Dark2",direction=-1) +
  theme_pubclean(base_size = 12)+
  labs(fill="",color="") +
  theme(legend.position = "none",
        legend.background = element_blank(),
        legend.key = element_blank(),
        legend.box.background = element_rect(colour = "black"),
        legend.title = element_blank()) +
  guides(fill = guide_legend(override.aes = list(size = 5,alpha=1) ) ) +
  facet_wrap(~id) +
  geom_abline(slope=1,intercept=0,lty=2) +
  scale_x_log10() +
  scale_y_log10() +
  xlab("Absolute Abundance Pre Thaw") +
  ylab("Absolute Abundance Post Thaw") +
  theme_bw() +
  stat_anova_test() +
  theme(legend.position = "right",
        axis.text.x = element_text(angle = 60,hjust=1))
dev.off()


################################################################################
# RelAb/AbAb Per Phylum
################################################################################

phy <- mag_df %>% 
  mutate(Genome=Name) %>%
  subset(select=c(Genome,Phylum))
ram_g <- merge.easy(ram_g,phy,key="Genome")
phy_keep <- names(table(ram_g$Phylum[ram_g$IoCc=="Copiotroph"]))[table(ram_g$Phylum[ram_g$IoCc=="Copiotroph"])>10]

ram_gp <- ram_g %>% subset(Phylum %in% phy_keep)

png("../figures/SFigure_by_phy.png",width=200,height=250,units="mm",res=1000)
ggplot(ram_gp %>% subset(delta!=0),
       aes(x=IoCc,
           y=delta,
           fill=IoCc)) + 
  geom_boxplot(aes(color=IoCc),fill="white") +
  geom_jitter(alpha=.25,size=3,pch=21,width=0.1) +
  geom_hline(yintercept = 0,lty=2) +
  scale_fill_brewer(palette="Dark2",direction=-1) +
  scale_color_brewer(palette="Dark2",direction=-1) +
  theme_pubclean(base_size = 12)+
  labs(fill="",color="") +
  theme(legend.position = "none",
        legend.background = element_blank(),
        legend.key = element_blank(),
        legend.box.background = element_rect(colour = "black"),
        legend.title = element_blank()) +
  ylab("Avg. Change in Transformed Rel. Abundance") +
  xlab("") +
  guides(fill = guide_legend(override.aes = list(size = 5,alpha=1) ) ) +
  facet_wrap(~Phylum) +
  theme_bw() +
  stat_anova_test() +
  theme(legend.position = "right",
        axis.text.x = element_text(angle = 60,hjust=1))
dev.off()


################################################################################
# Change Vs. Initial Abundance
################################################################################

png("../figures/SFigure_initial_vs_changes.png",
    width=180,height=120,units="mm",res=1000)
ggplot(ram_g %>% subset(value.pre+value.post > 0),
       aes(x=value.pre,y=(value.post-value.pre),fill=IoCc,size=IoC)) +
  geom_point(alpha=.75,pch=21) +
  geom_hline(yintercept = 0,lty=2) +
  scale_fill_brewer(palette="Dark2",direction=-1) +
  theme_pubclean(base_size = 12)+
  labs(fill="",size="IoC") +
  theme(legend.position = "right") +
  ylab("Avg. Change in Transformed Rel. Abundance") +
  xlab("Transformed Pre-Thaw Rel. Abundance") 
dev.off()

################################################################################
# Distribution Across Samples Pre/Post Thaw
################################################################################

mag_dist <- ram %>%
  mutate(is.pre=value.pre>0,
         is.post=value.post>0) %>%
  group_by(Genome) %>%
  summarise(pre_thaw=sum(is.pre),
            post_thaw=sum(is.post),
            IoCc=unique(IoCc))

mag_dist$bloomer <- (mag_dist$pre_thaw==0)
mag_dist$bloomer[mag_dist$post_thaw==0] <- -1

bloom_df <- mag_dist %>%
  group_by(IoCc,bloomer) %>%
  count() %>%
  group_by(IoCc) %>%
  summarise(nTot=sum(n),
            xB=sum(n[bloomer==1]),
            xO=sum(n[bloomer==0]),
            xC=sum(n[bloomer==-1])) %>%
  mutate(pB=xB/nTot,
         pO=xO/nTot,
         pC=xC/nTot,
         n=nTot)
bloom_df$LowerB <- mapply(FUN=binom.test,bloom_df$xB,bloom_df$n)["conf.int",] %>% lapply("[",1) %>% unlist()
bloom_df$UpperB <- mapply(FUN=binom.test,bloom_df$xB,bloom_df$n)["conf.int",] %>% lapply("[",2) %>% unlist()
bloom_df$LowerO <- mapply(FUN=binom.test,bloom_df$xO,bloom_df$n)["conf.int",] %>% lapply("[",1) %>% unlist()
bloom_df$UpperO <- mapply(FUN=binom.test,bloom_df$xO,bloom_df$n)["conf.int",] %>% lapply("[",2) %>% unlist()
bloom_df$LowerC <- mapply(FUN=binom.test,bloom_df$xC,bloom_df$n)["conf.int",] %>% lapply("[",1) %>% unlist()
bloom_df$UpperC <- mapply(FUN=binom.test,bloom_df$xC,bloom_df$n)["conf.int",] %>% lapply("[",2) %>% unlist()



pBloom <- ggarrange(ggplot(bloom_df, 
                 aes(x=IoCc, y = pB, fill=IoCc)) + 
            geom_bar(stat="identity",color="#1AFF1A",lwd=2) +
            geom_errorbar(aes(ymin=LowerB,ymax=UpperB,width=0.2)) +
            theme_pubclean(base_size = 10) +
            theme(legend.position = "none",
                  axis.text.x = element_text(angle = 60,hjust=1)) +
            ylab("Proportion \"Bloomers\" (Undetected Pre-Thaw)") +
            xlab("") +
            ylim(0,0.75) +
            scale_fill_brewer(palette="Dark2",direction=-1),
          ggplot(bloom_df, 
                 aes(x=IoCc, y = pC, fill=IoCc)) + 
            geom_bar(stat="identity",color="#4B0092",lwd=2) +
            geom_errorbar(aes(ymin=LowerC,ymax=UpperC,width=0.2)) +
            theme_pubclean(base_size = 10) +
            theme(legend.position = "none",
                  axis.text.x = element_text(angle = 60,hjust=1)) +
            ylab("Proportion \"Crashers\" (Undetected Post-Thaw)") +
            xlab("") +
            ylim(0,0.75) +
            scale_fill_brewer(palette="Dark2",direction=-1))


pNumDist <- ggplot(mag_dist,
       aes(x=pre_thaw+post_thaw,
           color=IoCc)) + 
  geom_density(lwd=1.2) +
  labs(color="") +
  scale_color_brewer(palette="Dark2",direction=-1) +
  theme_pubclean(base_size = 12) +
  theme(legend.position = "right") +
  xlab("No. Samples With Non-Zero Abundance") 


pNumBar <- ggplot(bloom_df, 
       aes(x=IoCc, y = n, fill=factor(bloomer))) + 
  geom_bar(stat="identity",position="fill") +
  scale_fill_manual(values=c("#4B0092","lightgray","#1AFF1A"),labels=c("Crash","Other","Bloom")) +
  labs(fill="") +
  theme_pubclean(base_size = 12) +
  theme(legend.position = "right",
        axis.text.x = element_text(angle = 60,hjust=1)) +
  xlab("") +
  ylab("Proportion")

pNumScatter <- ggplot(mag_dist,
       aes(x=pre_thaw,
           y=post_thaw,
           fill=IoCc)) +
  geom_jitter(pch=21,alpha=0.75,size=3,width=0.25,height=0.25) +
  scale_fill_brewer(palette="Dark2",direction=-1) +
  theme_pubclean(base_size = 12) +
  annotate("rect", ymin = 0.25, ymax = 9.75, xmin = -0.75, xmax = 0.75,
           alpha = 0, color= "#1AFF1A") +
  annotate("rect", xmin = 0.25, xmax = 9.75, ymin = -0.75, ymax = 0.75,
           alpha = 0, color= "#4B0092") +
  theme(legend.position = "none") +
  geom_abline(slope=1,intercept=0,lty=2) +
  xlab("No. Pre-Thaw Samples With Non-Zero Abundance") +
  ylab("No. Post-Thaw Samples With Non-Zero Abundance") +
  labs(fill="") +
  facet_wrap(~IoCc)

png("../figures/SFigure_distribution_across_samples.png",
    width=250,height=200,units="mm",res=1000)
ggarrange(pNumDist,
          ggarrange(pNumScatter,
                    pBloom,
                    ncol=2,
                    widths=c(2,1),
                    labels=c("(b)","(c)"),
                    hjust=0),
          nrow=2,
          heights=c(2,3)
          ,labels=c("(a)",""),
          hjust=0)
dev.off()

################################################################################
# High Coverage Metagenomes ----------------------------------------------------
################################################################################

metagenome_df <- metagenome_df %>% subset(gbp>=1)
b <- runif(nrow(metagenome_df), -0.1, 0.1) #same jitter for both violin plots

## IoC Metagenomes ----


pIoC <- ggplot() +
  geom_violin(data=metagenome_df,
              aes(y=IoC,
                  x=as.numeric(pre_post_thaw=="post"),
                  group=pre_post_thaw),
              fill="gray") +
  geom_boxplot(data=metagenome_df,
               aes(y=IoC,
                   x=as.numeric(pre_post_thaw=="post"),
                   group=pre_post_thaw),
               width=0.2) +
  geom_point(data=metagenome_df,
             aes(y=IoC,
                 x=as.numeric(pre_post_thaw=="post")+b,
                 group=id,
                 size=gbp,
                 shape=gbp<5,
                 color=IoCc),
             alpha=0.75) +
  geom_line(data=metagenome_df,aes(y=IoC,x=as.numeric(pre_post_thaw=="post")+b,
                                   group=id,color=IoCc),lwd=1.5) +
  theme_pubclean(base_size = 14) +
  theme(legend.position = "right") + 
  scale_x_continuous(breaks = c(0,1),
                     labels = c("Pre-Thaw", "Post-Thaw"),
                     expand = c(0,0))+
  scale_size_continuous(breaks=c(1,10),limits = c(0.1,35),range=c(3,10)) +
  scale_shape_discrete(solid=T) +
  xlab("") +
  labs(color="",fill="",size="Gbp",shape="<5 Gbp") +
  ylab("Avg. Index of Copiotrophy") +
  theme(legend.position = "right",
        axis.text.x = element_text(angle = 60,hjust=1)) +
  geom_hline(yintercept=c_h,lty=2,color="#1B9E77") +
  geom_hline(yintercept=c_l,lty=2,color="#7570B3") +
  scale_fill_brewer(palette="Dark2",direction=-1) +
  scale_color_brewer(palette="Dark2",direction=-1)
pIoC

## IoC Barplot ----

ioc_metagenomes <- as.data.frame(table(metagenome_df$IoCc,ppt=metagenome_df$pre_post_thaw))
ioc_metagenomes$ppt <- factor(ioc_metagenomes$ppt,
                              levels=c("pre",
                                       "post"))
pIoCbar <- ggplot(ioc_metagenomes, 
                  aes(x=ppt, y = Freq, fill=forcats::fct_rev(Var1))) + 
  geom_bar(stat="identity") +
  theme_pubclean(base_size = 14) +
  theme(legend.position = "none") +
  ylab("Number of Samples") +
  xlab("") +
  scale_fill_brewer(palette="Dark2",direction=) +
  scale_x_discrete(labels=c("Pre-Thaw","Post-Thaw")) +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 60,hjust=1))
pIoCbar

## MuMax Metagenome ----

pMuMax <- ggplot() +
  geom_violin(data=metagenome_df,aes(y=log(2)/d,x=as.numeric(pre_post_thaw=="post"),
                                     group=pre_post_thaw),fill="gray") +
  geom_boxplot(data=metagenome_df,aes(y=log(2)/d,x=as.numeric(pre_post_thaw=="post"),
                                      group=pre_post_thaw),
               width=0.2) +
  geom_point(data=metagenome_df,aes(y=log(2)/d,
                                    x=as.numeric(pre_post_thaw=="post")+b,
                                    group=id,
                                    size=gbp,
                                    shape=gbp<5,
                                    color=core),
             alpha=0.75) +
  geom_line(data=metagenome_df,aes(y=log(2)/d,x=as.numeric(pre_post_thaw=="post")+b,
                                   group=id,color=core),lwd=1.5) +
  theme_pubclean(base_size = 14) +
  scale_color_brewer(palette = "Set1") +
  scale_fill_brewer(palette = "Set1") +
  scale_size_continuous(breaks=c(1,10),limits = c(0.1,35),range=c(3,10)) +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 60,hjust=1)) + 
  scale_x_continuous(breaks = c(0,1),
                     labels = c("Pre-Thaw", "Post-Thaw"),
                     expand = c(0,0))+
  xlab("") +
  labs(color="Permafrost\nCore",fill="Permafrost\nCore") +
  ylab("Avg. Maximum Growth Rate (1/Hours)") 
pMuMax

## Carbon Plot 

carbon <- data.frame(core=c("AT3",
                            "AT4",
                            "BEO1",
                            "BEO2B",
                            "FL1C4",
                            "FL2C3"),
                     carbon_percent=c(3.6,
                                      1.8,
                                      30.9,
                                      21.8,
                                      37.9,
                                      39.6))
x <- merge.easy(metagenome_df,carbon,key="core")                     

dupl.id <- x$id[duplicated(x$id)]
x.paired <- x[x$id %in% dupl.id,]
x.post <- x.paired %>% subset(pre_post_thaw=="post")
x.post$diff <- 100*(x.paired[x.paired$pre_post_thaw=="pre",]$d - 
                      x.paired[x.paired$pre_post_thaw=="post",]$d)/x.paired[x.paired$pre_post_thaw=="pre",]$d

t.test(log(2)/x.paired[x.paired$pre_post_thaw=="pre",]$d,
       log(2)/x.paired[x.paired$pre_post_thaw=="post",]$d,
       paired=T)

t.test(x.paired[x.paired$pre_post_thaw=="pre",]$IoC,
       x.paired[x.paired$pre_post_thaw=="post",]$IoC,
       paired=T)

pCarbon <- ggplot(x.post,
                  aes(x=carbon_percent,y=-diff,fill=core,group=core)) +
  geom_boxplot(fill="white") +
  geom_point(size=5,pch=21,alpha=0.75) +
  scale_fill_brewer(palette = "Set1") +
  theme_pubclean(base_size = 12) +
  #scale_y_log10() +
  #scale_x_log10() +
  geom_smooth(fill="black",method="lm",color="black") +
  labs(fill="Permafrost\nCore") +
  ylab("% Change in Min. Doubling Time Over Experiment") +
  xlab("Organic C at Start of Experiment (%)") +
  theme(legend.position = "right") +
  geom_hline(yintercept=0,lty=2)

## Density IoC ----

pMdense <- ggplot(mag_df,
                  aes(x=IoC)) +
  geom_density(lwd=2,alpha=0.25) +
  geom_jitter(aes(y=-0.035),width=0,height=0.03,alpha=0.25,size=0.75) +
  geom_density(data=metagenome_df,
               aes(color=pre_post_thaw,fill=pre_post_thaw),
               lwd=1.1,alpha=0.25) +
  geom_jitter(data=metagenome_df,
              aes(fill=pre_post_thaw,
                  y=-0.035),
              width=0,height=0.03,size=3,alpha=0.75,pch=21) +
  theme_pubclean(base_size = 14) +
  geom_vline(xintercept=c_h,lty=2,color="#1B9E77") +
  geom_vline(xintercept=c_l,lty=2,color="#7570B3") +
  theme(legend.position = "right") +
  scale_color_brewer(palette = "Accent") +
  scale_fill_brewer(palette = "Accent") +
  labs(fill="Thaw",color="Thaw") +
  xlab("Index of Copiotrophy") +
  ylab("Density")
pMdense 

## Put it Together ----

png("../figures/SFigure_Metagenomes_highcov.png",width=250,height=250,units="mm",res=1000)
ggarrange(ggarrange(pMuMax,
                    pIoCbar,
                    pIoC,
                    widths=c(3,1,4),
                    ncol=3,
                    labels=c("(a)","(c)","(d)"),
                    hjust=0),
          ggarrange(pCarbon,
                    pMdense,
                    ncol=2,
                    widths=c(3,5),
                    labels=c("(b)","(e)")),
          nrow=2)
dev.off()


















