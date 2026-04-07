# JLW - 2026 - Permafrost Warming Experiments

################################################################################
# Packages ---------------------------------------------------------------------
################################################################################

library(dplyr)
library(data.table)
library(gRodon)
library(ggplot2)
library(ggpubr)
library(GGally)
library(ggrepel)
library(gridExtra)
library(cowplot)

merge.easy <- function(df1,df2,key){
  df1 <- data.table(df1,key=key)
  df2 <- data.table(df2,key=key)
  return(unique(merge(df1,df2,all.x=TRUE,by=.EACHI,allow.cartesian=TRUE)))
}


################################################################################
# Load Data (see warming_analysis_prep_data.R) ---------------------------------
################################################################################

setwd("~/../../Documents/permafrost_warming/data")

## Permafrost Data ----

load("warming_data_for_analysis.rda")
mag_df <- mag_df %>% subset(select=c(nGenes,nCAZy,dCUB)) %>%
  mutate(Environment="Permafrost")

## Pacific Data ----

load("P16_pacific.rda")
P16_df <- P16_df %>% subset(select=c(nGenes,nCAZy,dCUB)) %>%
  mutate(Environment="Pacific")

## Prefit IoC (permafrost from this paper and pacific from Zakem et al. 2025) --

load("warming_contained_mags_ioc.rda")
IoC_permafrost <- warming_contained_mags_ioc
load("PCA_model_P16.rda")
IoC_pacific <- x

## SMAG DB SOIL ----

load("smag_gRodon_out.rda")
pred_df$Name <- pred_df$Genome %>% 
  basename() %>%
  gsub(pattern=".ffn.gz",replace="")

cazy_df <- read.delim("dbcan_counts_warming_smag.tsv",header=F) %>%
  mutate(Name=gsub(".gz","",V1),
         nCAZy=V2) %>%
  subset(select=c(Name,nCAZy))

smag <- merge.easy(pred_df,cazy_df,key="Name") %>%
  subset(select=c(Name,dCUB,nHE,nGenes,nCAZy,d))

qual <- read.delim("quality_report.tsv")%>%
  mutate(Name=gsub(".fna","",Name))
smag <- merge.easy(smag,qual,key="Name") %>%
  subset(nHE>=10 & !is.na(nCAZy) & !is.na(Completeness))
smag$nGenes_corrected <- 100*smag$nGenes/smag$Completeness


## PASOLLI GUT ----

load("pasolli_gRodon_out.rda")
pred_df$Name <- pred_df$Genome %>% 
  basename() %>%
  gsub(pattern=".ffn",replace="")

cazy_df <- read.delim("dbcan_counts_passoli.tsv",header=F,sep=" ") %>%
  mutate(Name=gsub(".gz","",V2),
         nCAZy=V1) %>%
  subset(select=c(Name,nCAZy))

human <- merge.easy(pred_df,cazy_df,key="Name") %>%
  subset(select=c(Name,dCUB,nHE,nGenes,nCAZy,d)) 

qual <- read.csv("pasolli_quality.csv") %>%
  mutate(Name=gsub("[.]","_",Genome.Name))
human <- merge.easy(human,qual,key="Name") %>%
  subset(nHE>=10 & !is.na(nCAZy) & !is.na(Completeness))
human$nGenes_corrected <- 100*human$nGenes/human$Completeness

## Fit SMAG IoC ----

smag$rCAZy <- smag$nCAZy/smag$nGenes
smag$nGenes <- smag$nGenes_corrected
IoC_smag <- fitIoC(x=smag %>% as.data.frame())

## Fit PASOLLI IoC ----

human$rCAZy <- human$nCAZy/human$nGenes
human$nGenes <- human$nGenes_corrected
IoC_human <- fitIoC(x=human %>% as.data.frame())


################################################################################
# Compare IoC ------------------------------------------------------------------
################################################################################

smag$Soil <- IoC_smag$IoC
smag$Permafrost <- predictIoC(x=smag %>% as.data.frame(),model=IoC_permafrost$IoCmod)
smag$Pacific <- predictIoC(x=smag %>% as.data.frame(),model=IoC_pacific)
smag$Human <- predictIoC(x=smag %>% as.data.frame(),model=IoC_human$model)

loading_df <- rbind(abs(IoC_smag$model$rotation[,"PC1"])/sum(abs(IoC_smag$model$rotation[,"PC1"])),
                     abs(IoC_permafrost$IoCmod$model$rotation[,"PC1"])/sum(abs(IoC_permafrost$IoCmod$model$rotation[,"PC1"])),
                     abs(-IoC_pacific$rotation[,"PC1"])/sum(abs(-IoC_pacific$rotation[,"PC1"])),
                     abs(IoC_human$model$rotation[,"PC1"])/sum(abs(IoC_human$model$rotation[,"PC1"]))) %>%
  as.data.frame() %>%
  mutate(environment=c("Soil","Permafrost","Pacific","Human")) %>%
  melt(id.vars="environment")

pPairs <- ggpairs(smag,columns=20:23)

pLoad <- ggplot(loading_df, 
       aes(x=environment, y=variable, fill= (value))) + 
  geom_tile() +
  scale_fill_viridis_c() +
  theme_bw() +
  labs(fill="Relative\nLoading") +
  xlab("") +
  ylab("")



loading_mat <- rbind(IoC_smag$model$rotation[,"PC1"],
                     IoC_permafrost$IoCmod$model$rotation[,"PC1"],
                     -IoC_pacific$rotation[,"PC1"],
                     IoC_human$model$rotation[,"PC1"]) %>%
  as.data.frame() %>%
  mutate(environment=c("Soil","Permafrost","Pacific","Human"))


pVar <- ggplot(loading_mat,
       aes(x=dCUB,y=nGenes,fill=rCAZy,label=environment)) +
  geom_point(size=8,pch=21) +
  scale_fill_viridis_c(option="G") +
  theme_pubclean() +
  geom_text_repel(color="black",box.padding = 1, max.overlaps = Inf) +
  xlim(0,0.8) +
  ylim(-0.8,0) +
  theme(legend.position = c(0.25,0.7)) +
  xlab("Loading (dCUB)") +
  ylab("Loading (nGenes)") +
  labs(fill="Loading\n(rCAZy)")

pPC <- ggarrange(pVar,pLoad,ncol=2,labels=c("(a)","(b)"),widths=c(2,3))

png("../figures/SFigure_ioc_environment_comparison.png",width=200,height=200,units="mm",res=1000)
ggarrange(pPC,
          ggmatrix_gtable(pPairs),
          labels=c("","(c)"),
          hjust=0,
          nrow=2)
dev.off()


################################################################################
# Compare Data -----------------------------------------------------------------
################################################################################

smag <- smag %>% subset(select=c(nGenes,nCAZy,dCUB)) %>%
  mutate(Environment="Soil")

human <- human %>% subset(select=c(nGenes,nCAZy,dCUB)) %>%
  mutate(Environment="Human")

plot_df <- rbind(mag_df,
                 smag,
                 human,
                 P16_df) %>%
  mutate(rCAZy=nCAZy/nGenes)

pVars <- ggplot(plot_df,
  aes(x=dCUB,y=nGenes,fill=rCAZy)) +
  geom_point(data = subset(plot_df, select = -Environment), 
             fill = "black", pch = 19,size=1,alpha=0.1) +
  geom_point(pch=21,size=3,alpha=0.75,color="gray") +
  facet_wrap(~Environment) +
  scale_y_log10() +
  scale_fill_viridis_c(option="F") +
  theme_bw() +
  geom_smooth(method="lm") +
  stat_cor()

pG <- ggplot(plot_df,
       aes(x=Environment,group=Environment,y=nGenes)) +
  geom_violin(fill="gray") +
  geom_boxplot(width=0.1) +
  stat_anova_test(label.y=5.75) +
  stat_pwc() +
  theme_pubclean() +
  scale_y_log10() +
  xlab("")+
  theme(axis.text.x = element_text(angle = 60,hjust=1))

pD <- ggplot(plot_df,
       aes(x=Environment,group=Environment,y=dCUB)) +
  geom_violin(fill="gray") +
  geom_boxplot(width=0.1) +
  stat_anova_test(label.y=1.25) +
  stat_pwc() +
  theme_pubclean() +
  xlab("")+
  theme(axis.text.x = element_text(angle = 60,hjust=1))

pC <- ggplot(plot_df,
       aes(x=Environment,group=Environment,y=rCAZy)) +
  geom_violin(fill="gray") +
  geom_boxplot(width=0.1) +
  stat_anova_test(label.y=0.225) +
  stat_pwc() +
  theme_pubclean() +
  xlab("") +
  theme(axis.text.x = element_text(angle = 60,hjust=1))

png("../figures/SFigure_vars_environment_comparison.png",width=400,height=200,units="mm",res=1000)
ggarrange(pVars,pG,pC,pD,ncol=4,widths = c(3,1,1,1))
dev.off()

################################################################################
# Super IoC  -------------------------------------------------------------------
################################################################################

plot_df <- plot_df %>% 
  subset(Environment!="Soil") %>%
  group_by(Environment) %>%
  sample_n(size=500)

IoC_all <- fitIoC(x=plot_df %>% as.data.frame())
plot_df$Soil <- predictIoC(x=plot_df %>% as.data.frame(),model=IoC_smag$model)
plot_df$Permafrost <- predictIoC(x=plot_df %>% as.data.frame(),model=IoC_permafrost$IoCmod)
plot_df$Pacific <- predictIoC(x=plot_df %>% as.data.frame(),model=IoC_pacific)
plot_df$Human <- predictIoC(x=plot_df %>% as.data.frame(),model=IoC_human$model)
plot_df$All <- predictIoC(x=plot_df %>% as.data.frame(),model=IoC_all$model)

png("../figures/SFigure_ioc_all.png",width=200,height=200,units="mm",res=1000)
ggpairs(plot_df,columns=6:10)
dev.off()

IoC.soil <- IoC_smag$model
IoC.permafrost <- IoC_permafrost$IoCmod$model
IoC.pacific <- IoC_pacific
IoC.human <- IoC_human$model
IoC.all <- IoC_all$model

IoC_models <- list(Soil=IoC.soil,
                   Permafrost=IoC.permafrost,
                   Pacific=IoC.pacific,
                   Human=IoC.human,
                   All=IoC.all)

save(IoC_models,
     file="IoC_Models.rda")

