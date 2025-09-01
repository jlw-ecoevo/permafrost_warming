# JLW - 2025 - Permafrost Warming Experiments

# Packages ---------------------------------------------------------------------

library(dplyr)
library(ggplot2)
library(ggpubr)
library(data.table)
library(SPRING)
library(stringr) 
library(corrplot)
library(ggfortify)
library(mclust)
library(ape)
library(phytools)
library(rcartocolor)
library(ggridges)


merge.easy <- function(df1,df2,key){
  df1 <- data.table(df1,key=key)
  df2 <- data.table(df2,key=key)
  return(unique(merge(df1,df2,all.x=TRUE,by=.EACHI,allow.cartesian=TRUE)))
}

# Load Data --------------------------------------------------------------------

setwd("C:/Users/jlwei/Documents/permafrost_warming/data")

load("warming_gRodon_predictions_ms_scaffolds.rda")
pred_scaff <- pred_df %>%
  mutate(fastq_id=basename(Genome) %>%
           gsub(pattern="*.fasta.ffn",replace=""))

arch_names <- read.delim("gtdbtk.ar53.summary.tsv")$user_genome
bac_names <- read.delim("gtdbtk.bac120.summary.tsv")$user_genome
load("warming_gRodon_predictions_ms.rda")
pred_df <- pred_df %>%
  subset(grepl("prokka",Genome))
pred_df$arch <- grepl("arch",pred_df$Genome)
pred_df$Name <- basename(pred_df$Genome) %>%
  gsub(pattern="_ms.*",replace="")
pred_arch <- pred_df %>%
  subset(arch==1 & Name %in% c(arch_names))
pred_bac <- pred_df %>%
  subset(arch==0 & Name %in% c(bac_names))
pred_df <- rbind(pred_bac,pred_arch)

load("cog_counts.rda")
cog_df$arch <- grepl("arch",cog_df$Genome)
cog_df$Name <- basename(cog_df$Genome) %>%
  gsub(pattern="_ms.*",replace="")
cog_arch <- cog_df %>%
  subset(arch==1 & Name %in% c(arch_names))
cog_bac <- cog_df %>%
  subset(arch==0 & Name %in% c(bac_names))
cog_df <- rbind(cog_bac,cog_arch)


cazy_bac <- read.delim("dbcan_counts_warming_ms.tsv",header=F) %>%
  mutate(Name=basename(V1),
         nCAZy=V2) %>%
  subset(select=c(Name,nCAZy))
cazy_arch <- read.delim("dbcan_counts_warming_ms_arch.tsv",header=F) %>%
  mutate(Name=basename(V1),
         nCAZy=V2) %>%
  subset(select=c(Name,nCAZy))
cazy_df <- rbind(cazy_bac %>% subset(Name %in% bac_names),
                 cazy_arch %>% subset(Name %in% arch_names))
pred_df <- merge.easy(pred_df,cazy_df,key="Name") %>%
  subset(nHE>=10 & !is.na(dCUB))


# Run PCA, take PC1
pred_df$rCAZy <- pred_df$nCAZy/pred_df$nGenes
x <- prcomp(pred_df[,c("rCAZy","dCUB","nGenes")],
            scale=T)
autoplot(x,loadings=T,loadings.label=T,loadings.label.vjust = 1.3)
pred_df$PC1 <- x$x[,1]
#Make it positive for plotting
pred_df$IoC <- pred_df$PC1-min(pred_df$PC1)+1
save(x,file="PCA_model_warming.rda")

meta <- read.delim("warming_metadata.tsv") %>%
  mutate(fastq_id=str_sub(fastq_id, end = -2))
qual <- read.delim("quality_report.tsv")
qual$Genome <- paste0(qual$Name,".ffn")
coverage <- read.delim("bin_coverage_coverm_metaspades.tsv")

tax <- rbind(
  read.delim("gtdbtk.bac120.summary.tsv",sep="\t") %>%
    mutate(Name=user_genome,
           Taxonomy=classification) %>%
    subset(select=c(Name,Taxonomy)),
  read.delim("gtdbtk.ar53.summary.tsv",sep="\t") %>%
    mutate(Name=user_genome,
           Taxonomy=classification) %>%
    subset(select=c(Name,Taxonomy))
)
pred_df <- merge.easy(pred_df,tax,key="Name")
pred_df$Phylum <- pred_df$Taxonomy %>%
  gsub(pattern=";c__.*",replace="") %>%
  gsub(pattern=".*;p__",replace="")

qual <- merge.easy(qual,pred_df,key="Name")
table(Contam=qual$Contamination<5,nHE=qual$nHE>=10)

pred_df$fastq_id <- basename(pred_df$Genome) %>% 
  gsub(pattern="_bin.*",replacement="")

x <- merge.easy(meta,pred_scaff,key="fastq_id") %>% 
  subset(nHE>=10) %>%
  subset(!is.na(core))
x$id <- x$fastq_id %>% 
  gsub(pattern="_.*",replace="") %>%
  gsub(pattern="p",replace="") 

qual$nHE[is.na(qual$nHE)] <- 0

mags <- pred_df[grep("bin",pred_df$Genome),] %>% 
  subset(nHE>=10) %>%
  subset(Name %in% qual$Name[qual$Contamination<5])
mags$pre_post <- !grepl("p",mags$fastq_id)
mags$id <- mags$fastq_id %>% 
  gsub(pattern="_.*",replace="") %>%
  gsub(pattern="p",replace="")


ra <- coverage[,c(grep("Relative.Abundance",names(coverage)))] 
ra <- t(t(ra)/rowSums(t(ra)))
ra[is.na(ra)] <- 0
ra <- ra %>% 
  t() %>%
  mclr() %>%
  t() %>%
  as.data.frame() %>%
  mutate(Genome=coverage$Genome) %>%
  as.data.table() %>%
  melt(id.vars = "Genome")
ra$fastq_id <- ra$variable %>%
  gsub(pattern="_R1_001.clean.fastq.Relative.Abundance....",replace="") %>%
  gsub(pattern="X",replace="") %>%
  str_sub(., end = -2)
ra <- merge.easy(ra,x,key="fastq_id") %>%
  mutate(Genome=Genome.x) %>%
  subset(select=c(Genome,value,core,pre_post_thaw,OGT,dCUB,d,nHE,nGenes,id)) %>%
  subset(!is.na(pre_post_thaw)) %>%
  reshape(direction = "wide",timevar = "pre_post_thaw",
          idvar = c("Genome","core","id"))
m <- mags %>%
  mutate(Genome=Name) %>%
  subset(select=c(Genome,dCUB,IoC,nCAZy,nGenes,d,nHE,pre_post,Taxonomy))
ram <- merge.easy(ra,m,key="Genome") %>%
  subset(!is.na(d) & !is.na(value.pre) & !is.na(value.post))
ram$delta <- (ram$value.post - ram$value.pre)/((ram$value.post+ram$value.pre))
ram$source <- "Pre-Thaw"
ram$source[ram$pre_post==T] <- "Post-Thaw"
ram_g <- ram %>% 
  group_by(Genome) %>% 
  summarise_all(mean)


ra_cov <- coverage[,c(grep("Relative.Abundance",names(coverage)))] 
ra_cov <- t(t(ra_cov)/rowSums(t(ra_cov)))
ra_cov[is.na(ra_cov)] <- 0
ra_cov <- ra_cov %>% 
  t() %>%
  t() %>%
  as.data.frame() %>%
  mutate(Genome=coverage$Genome) %>%
  as.data.table() %>% 
  melt(id.vars = "Genome")
ra_cov$fastq_id <- ra_cov$variable %>%
  gsub(pattern="_R1_001.clean.fastq.Relative.Abundance....",replace="") %>%
  gsub(pattern="X",replace="") %>%
  str_sub(., end = -2)
ra_cov <- ra_cov %>% 
  subset(Genome %in% pred_df$Name) %>%
  subset(select=c(fastq_id,value)) %>%
  group_by(fastq_id) %>%
  summarise(RA=sum(value))


# Mag Qual Plot ----------------------------------------------------------------


pQual <- ggplot(qual,aes(x=Completeness,y=Contamination,fill=(nHE>=10))) + 
  geom_point(pch=21,size=2,alpha=0.75) + 
  scale_y_log10() +
  theme_pubclean(base_size = 14) +
  scale_fill_manual(values=c("white","black")) +
  geom_hline(yintercept = 5,lty=2,color="red") +
  theme(legend.position = "right") +
  labs(fill="At least 10\nRibosomal Proteins?") +
  xlab("Completeness (%)") +
  ylab("Contamination (%)")

b2 <- runif(nrow(mags), -0.1, 0.1)
pSource <- ggplot() +
  geom_violin(data=mags,aes(y=d,x=as.numeric(pre_post),
                            group=pre_post),fill="gray") +
  geom_boxplot(data=mags,aes(y=d,x=as.numeric(pre_post),
                             group=pre_post),
               width=0.2) +
  geom_point(data=mags,aes(y=d,x=as.numeric(pre_post)+b2),
             size=1,alpha=0.75) +
  theme_pubclean(base_size = 14) +
  scale_color_brewer(palette = "Set1") +
  theme(legend.position = "right",
        axis.text.x = element_text(angle = 60, hjust=1)) + 
  scale_x_continuous(breaks = c(0,1),
                     labels = c("Pre-Thaw", "Post-Thaw"),
                     expand = c(0,0))+
  xlab("Assembly Source") +
  ylab("Predicted Min. Doubling Time (Hours)")

setwd("../figures/")

png("warming_experiments_mags_qual.png",width=8,height=5,units="in",res=500)
ggarrange(pQual,pSource,widths=c(3,1),labels=c("(a)","(b)"))
dev.off()

median(ra_cov$RA)

t.test(d~pre_post,data=mags)

range(mags$d)

# Metagenome Plot --------------------------------------------------------------

b <- runif(nrow(x), -0.1, 0.1)
pMeta <- ggplot() +
  geom_violin(data=x,aes(y=d,x=as.numeric(pre_post_thaw=="post"),
                         group=pre_post_thaw),fill="gray") +
  geom_boxplot(data=x,aes(y=d,x=as.numeric(pre_post_thaw=="post"),
                     group=pre_post_thaw),
               width=0.2) +
  geom_point(data=x,aes(y=d,x=as.numeric(pre_post_thaw=="post")+b,
                    group=id,fill=core),
              size=5,alpha=0.75,pch=21) +
  geom_line(data=x,aes(y=d,x=as.numeric(pre_post_thaw=="post")+b,
                         group=id,color=core),lwd=1.5) +
  theme_pubclean(base_size = 14) +
  scale_color_brewer(palette = "Set1") +
  scale_fill_brewer(palette = "Set1") +
  theme(legend.position = "right") + 
  scale_x_continuous(breaks = c(0,1),
                     labels = c("Pre-Thaw", "Post-Thaw"),
                     expand = c(0,0))+
  xlab("") +
  labs(color="Permafrost\nCore") +
  ylab("Predicted Avg. Min. Doubling Time (Hours)") +
  theme(legend.position = "none")



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
x <- merge.easy(x,carbon,key="core")                     




dupl.id <- x$id[duplicated(x$id)]
x.paired <- x[x$id %in% dupl.id,]
x.post <- x.paired %>% subset(pre_post_thaw=="post")
x.post$diff <- x.paired[x.paired$pre_post_thaw=="post",]$d - 
  x.paired[x.paired$pre_post_thaw=="pre",]$d

pCarbon <- ggplot(x.post,
       aes(x=carbon_percent,y=-diff,fill=core,group=core)) +
  geom_boxplot(fill="white") +
  geom_point(size=5,pch=21,alpha=0.75) +
  scale_fill_brewer(palette = "Set1") +
  theme_pubclean(base_size = 14) +
  #scale_y_log10() +
  #scale_x_log10() +
  geom_smooth(fill="black",method="lm",color="black") +
  labs(fill="Permafrost\nCore") +
  ylab("Decrease in Predicted Min. Doubling Time Over Experiment (Hours)") +
  xlab("Organic C at Start of Experiment (%)") +
  theme(legend.position = "right") +
  geom_hline(yintercept=0,lty=2)
pCarbon


png("warming_experiments_metagenomes.png",width=7,height=7,units="in",res=500)
ggarrange(pMeta,pCarbon,labels = c("(a)","(b)"),hjust = 0)
dev.off()


t.test(x.paired[x.paired$pre_post_thaw=="pre",]$d,
       x.paired[x.paired$pre_post_thaw=="post",]$d,
       paired=T)


mean(x.paired[x.paired$pre_post_thaw=="pre",]$d - 
     x.paired[x.paired$pre_post_thaw=="post",]$d)


cor.test(x.post$diff,x.post$carbon_percent)


png("warming_experiments_S1.png",width=5,height=5,units="in",res=500)
ggplot(x,
       aes(x=carbon_percent,y=d,color=core,group=core,shape=pre_post_thaw)) +
  geom_point(size=7,alpha=0.75,color="black") +
  geom_point(size=5,alpha=0.75) +
  scale_color_brewer(palette = "Set1",guide = "legend",) +
  theme_pubclean(base_size = 14) +
  # scale_shape_manual(guide = "legend") +
  labs(color="Permafrost\nCore",shape="") +
  ylab("Predicted Min. Doubling Time (Hours)") +
  xlab("Organic C at Start of Experiment (%)") +
  theme(legend.position = "right") +
  geom_boxplot(fill=rgb(0,0,0,0.1),shape=NA,color="gray10") 
dev.off()


# IoC --------------------------------------------------------------------------

#Clustering
mC <- Mclust(log10(pred_df$IoC),
             verbose=T)
mC$parameters$mean
mC$parameters$variance$sigmasq
p <- (colSums(mC$z)/sum(mC$z))
x_seq <- seq(-0.2,1,0.01)
cl_df <- data.frame(x=10^x_seq,
                    cl1=dnorm(x_seq,
                              mean=mC$parameters$mean[1],
                              sd=sqrt(mC$p$variance$sigmasq[1]))*p[1],
                    cl2=dnorm(x_seq,
                              mean=mC$parameters$mean[2],
                              sd=sqrt(mC$p$variance$sigmasq[1]))*p[2])


#Confidence limits
plot(10^mC$data,mC$uncertainty)
points(10^mC$data,mC$classification/10,col="red")
abline(v=2.22)
abline(v=2.77)
abline(v=3.39)
abline(h=0.1)
c_l <- max(10^mC$data[mC$uncertainty<0.1 & 10^mC$data<2.5])
c_h <- min(10^mC$data[mC$uncertainty<0.1 & 10^mC$data>2.5])
c_m <- 10^mC$data[which.max(mC$uncertainty)]


#Plot: Cluster Plot for IoC
pC <- ggplot() + 
  scale_x_log10() +
  geom_density(data=pred_df,aes(x=IoC),lwd=2) + 
  geom_polygon(data=cl_df,aes(x=x, y=cl1, fill="Cluster 1"),
               alpha=0.75, color="black") +
  geom_polygon(data=cl_df,aes(x=x, y=cl2, fill="Cluster 2"), 
               alpha=0.75, color="black") +
  theme_pubclean(base_size = 14) + 
  xlab("Index of Copiotrophy")  + 
  ylab("Density") +
  theme(legend.position = "right") + 
  labs(fill="") +
  geom_vline(xintercept=c_l,lty=2) +
  geom_vline(xintercept=c_h,lty=2) +
  #geom_vline(xintercept=c_m,lty=2) +
  scale_fill_manual(values=c("#7570B3","#1B9E77")) +
  theme(legend.position = "bottom")
pC

pred_df$IoCc <- "Oligotroph"
pred_df$IoCc[pred_df$IoC>c_l] <- "Undefined"
pred_df$IoCc[pred_df$IoC>c_h] <- "Copiotroph"
pred_df$IoCc <- factor(pred_df$IoCc,
                       levels=c("Oligotroph",
                                "Undefined",
                                "Copiotroph"))
pIoC <- ggplot(pred_df,aes(x=-dCUB,y=nCAZy,fill=IoCc,size=nGenes)) +
  geom_point(alpha=.5,pch=21) +
  #geom_smooth(method="gam",color="black",fill="black") +
  scale_fill_brewer(palette="Dark2",direction=-1) +
  theme_pubclean(base_size = 14)+
  labs(fill="",size="Number of\nGenes") +
  scale_y_log10() +
  theme(legend.position = "right") +
  geom_hline(yintercept = 30,lty=2) +
  geom_vline(xintercept = 0.08,lty=2) +
  ylab("Number of CAZymes") +
  xlab("Codon Usage Bias (-dCUB)") + 
  guides(fill = guide_legend(override.aes = list(size=5)))


phyplot <- ggplot(pred_df,
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
  geom_jitter(data=pred_df,
             aes(x=IoC, y=reorder(Phylum,IoC,FUN=median),fill=IoCc),
             alpha=0.5,size=3,pch=21,
             height=0.1,width=0) +
  scale_fill_brewer(palette="Dark2",direction=-1) +
  labs(fill="Quartile") + 
  theme_pubclean(base_size = 14) + 
  scale_x_log10() +
  geom_vline(xintercept=c_l,lty=2) +
  geom_vline(xintercept=c_h,lty=2) +
  xlab("Index of Copiotrophy") + 
  ylab("") +
  theme(legend.position = "none")
phyplot


png("warming_experiments_mags_ioc.png",width=10,height=10,units="in",res=500)
ggarrange(ggarrange(pC,
                    pIoC,
                    labels=c("(a)","(b)"),
                    heights=c(2,3),
                    nrow=2,
                    hjust=0),
          phyplot,
          labels=c("","(c)"),
          ncol=2)
dev.off()

# Rel. Abundance ---------------------------------------------------------------

ram_g$IoCc <- "Oligotroph"
ram_g$IoCc[ram_g$IoC>c_l] <- "Undefined"
ram_g$IoCc[ram_g$IoC>c_h] <- "Copiotroph"
ram_g$IoCc <- factor(ram_g$IoCc,
                   levels=c("Oligotroph",
                            "Undefined",
                            "Copiotroph"))
p3c <- ggplot(ram_g %>% subset(value.pre+value.post > 0),
       aes(x=IoC,y=(value.post-value.pre),fill=IoCc)) +
  geom_point(alpha=.5,size=5,pch=21) +
  geom_hline(yintercept = 0,lty=2) +
  geom_smooth(method="gam",color="black",fill="black") +
  scale_fill_brewer(palette="Dark2",direction=-1) +
  theme_pubclean(base_size = 14)+
  labs(fill="") +
  theme(legend.position = c(0.8,0.1)) +
  ylab("Average Change in Relative Abundance") +
  xlab("Index of Copiotrophy")
p3c

p3d <- ggplot(ram_g %>% subset(value.pre+value.post > 0),
              aes(x=d,y=(value.post-value.pre))) +
  geom_point(alpha=.5,size=5,pch=21,fill="black") +
  geom_hline(yintercept = 0,lty=2) +
  geom_smooth(method="gam",color="black",fill="black") +
  scale_fill_brewer(palette="Dark2",direction=-1) +
  theme_pubclean(base_size = 14)+
  labs(fill="") +
  scale_x_log10() +
  theme(legend.position = c(0.8,0.1)) +
  ylab("Average Change in Relative Abundance") +
  xlab("Predicted Min. Doubling Time (Hours)")
p3d


p4b <- ggplot(ram_g %>% subset(value.pre+value.post > 0),
              aes(x=IoCc,y=(value.post-value.pre),fill=IoCc)) +
  geom_violin(fill="gray",width=1.5) +
  geom_boxplot(width=0.1,fill="white",aes(color=IoCc)) +
  geom_jitter(pch=21,width=0.1,alpha=0.5,size=2) +
  scale_fill_brewer(palette="Dark2",direction=-1) +
  scale_color_brewer(palette="Dark2",direction=-1) +
  theme_pubclean(base_size = 14)+
  labs(fill="Predicted\nCopiotrophic?") +
  theme(legend.position = "none") +
  ylab("Average Change in Relative Abundance") +
  xlab("")
pb <- ggplot(ram_g,aes(x=IoCc,y=(value.post-value.pre))) +
  geom_boxplot(width=0.6)
dat_1 <- ggplot_build(pb)$data[[1]]
dat_2 <- ggplot_build(p4b)$data[[2]]
dat <- data.frame(xmin=c(unique(dat_1$xmin)),
                  xmax=c(unique(dat_1$xmax)),
                  middle=c(dat_2$middle),
                  IoCc=factor(c("Oligotroph",
                         "Undefined",
                         "Copiotroph"),
                         levels=c("Oligotroph",
                                  "Undefined",
                                  "Copiotroph")))
p4b <- p4b + geom_segment(data=dat, aes(x=xmin, xend=xmax, 
                                      y=middle, yend=middle,
                                      color=IoCc), 
                        size=1,lty=1) +
  theme(axis.text.x = element_text(angle = 30, hjust=1)) +
  geom_hline(yintercept=0,lty=2)
p4b


y <- ram_g %>% subset(value.pre+value.post > 0)
y$del <- y$value.post-y$value.pre
cor.test(y$del,y$IoC)
cor.test(y$del,y$d)

ram_g$del0 <- (ram_g$value.post -ram_g$value.pre)>0
table(del0=ram_g$del0,IoC=ram_g$IoCc)


png("warming_experiments_mags_relab.png",width=7,height=10,units="in",res=500)
ggarrange(p3d,
          ggarrange(p3c,
                    p4b,
                    ncol=2,
                    widths=c(2,1),
                    labels=c("(b)","(c)")),
          nrow=2,
          labels=c("(a)",""))
dev.off()

# Panels -----------------------------------------------------------------------



png("warming_experiments_panels2.png",width=17,height=10,units="in",res=500)
ggarrange(ggarrange(ggarrange(pC,
                              pIoC,
                              labels=c("(a)","(b)"),
                              heights=c(2,3),
                              nrow=2,
                              hjust=0),
                    phyplot,
                    labels=c("","(c)"),
                    ncol=2),
          ggarrange(p3d,
                    ggarrange(p3c,
                              p4b,
                              ncol=2,
                              widths=c(2,1),
                              labels=c("(e)","(f)")),
                    nrow=2,
                    labels=c("(d)","")),
          ncol=2,
          widths=c(10,7))
dev.off()





