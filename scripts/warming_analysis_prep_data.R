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

merge.easy <- function(df1,df2,key){
  df1 <- data.table(df1,key=key)
  df2 <- data.table(df2,key=key)
  return(unique(merge(df1,df2,all.x=TRUE,by=.EACHI,allow.cartesian=TRUE)))
}

################################################################################
# MAG Data ---------------------------------------------------------------------
################################################################################

setwd("~/../../Documents/permafrost_warming/data")

## MAGMAX PERMAFROST ----

load("warming_gRodon_predictions_magmax.rda")
pred_df$Name <- pred_df$Genome %>% gsub(pattern=".ffn",replace="")

cazy_df <- read.delim("dbcan_counts_warming_magmax.tsv",header=F) %>%
  mutate(Name=V1,
         nCAZy=V2) %>%
  subset(select=c(Name,nCAZy))

permafrost <- merge.easy(pred_df,cazy_df,key="Name") %>%
  subset(select=c(Name,dCUB,nHE,nGenes,nCAZy,d))


## WOODCROFT PERMAFROST ----

load("woodcroft_gRodon_predictions.rda")
pred_df$Name <- pred_df$Genome %>% gsub(pattern=".ffn",replace="")

cazy_df <- read.delim("dbcan_counts_warming_woodcroft.tsv",header=F) %>%
  mutate(Name=V1,
         nCAZy=V2) %>%
  subset(select=c(Name,nCAZy))

woodcroft <- merge.easy(pred_df,cazy_df,key="Name") %>%
  subset(select=c(Name,dCUB,nHE,nGenes,nCAZy,d))


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


pred_df <- rbind(woodcroft,smag,permafrost)


contained_mags <- readLines("genomes_contained_cov.txt")
table(contained_mags %in% pred_df$Name)

pred_df <- pred_df %>%
  subset(Name %in% contained_mags) %>%
  subset(nHE>=10) %>%
  subset(!is.na(dCUB) & !is.na(nCAZy))

qual <- read.delim("quality_report_contained.tsv")
pred_df <- merge.easy(pred_df,qual,key="Name")
pred_df$nGenes_corrected <- 100*pred_df$nGenes/pred_df$Completeness

# TAXONOMY ---------------------------------------------------------------------

tax <- rbind(read.delim("gtdbtk.ar53.summary.contained.tsv"),
             read.delim("gtdbtk.bac120.summary.contained.tsv")) %>%
  mutate(Name=user_genome,
         Phylum=classification %>%
           gsub(pattern=";c__.*",replace="") %>%
           gsub(pattern=".*;p__",replace="")) %>%
  subset(select=c(Name,Phylum,classification))

pred_df <- merge.easy(pred_df,tax,key="Name")

# COVERAGE ---------------------------------------------------------------------

coverage <- read.delim("bin_coverage_coverm_magmax_contained.tsv")
meta <- read.delim("warming_metadata.tsv") %>%
  mutate(fastq_id=str_sub(fastq_id, end = -2)) %>%
  mutate(fastq_id=substr(fastq_id,1,4))
meta$id <- meta$fastq_id %>% 
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
  gsub(pattern="1.fastq.Relative.Abundance....",replace="") %>%
  gsub(pattern="X",replace="") %>%
  str_sub(., end = -2) %>%
  substr(.,1,4)
ra <- merge.easy(ra,meta,key="fastq_id") %>%
  subset(select=c(Genome,core,pre_post_thaw,id,value)) %>%
  subset(!is.na(pre_post_thaw)) %>%
  reshape(direction = "wide",timevar = "pre_post_thaw",
          idvar = c("Genome","core","id"))
ra[is.na(ra)] <- 0


# TOTAL COVERAGE ---------------------------------------------------------------

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
  gsub(pattern="1.fastq.Relative.Abundance....",replace="") %>%
  gsub(pattern="X",replace="") %>%
  str_sub(., end = -2)
ra_cov <- ra_cov %>% 
  subset(Genome %in% pred_df$Name) %>%
  subset(select=c(fastq_id,value)) %>%
  group_by(fastq_id) %>%
  summarise(RA=sum(value))

ra_cov$id <- ra_cov$fastq_id %>% 
  gsub(pattern="_.*",replace="") %>%
  gsub(pattern="p",replace="") 

ra_cutoff <- ra %>%
  mutate(tot=value.pre+value.post) %>%
  group_by(Genome) %>%
  summarise(tot=sum(tot))

non_zero_cov <- ra_cutoff$Genome[ra_cutoff$tot>0] %>% unique()
pred_df <- pred_df %>%
  subset(Name %in% non_zero_cov)

################################################################################
# Fit IoC ----------------------------------------------------------------------
################################################################################

pred_df$rCAZy <- pred_df$nCAZy/pred_df$nGenes
pred_df$nGenes <- pred_df$nGenes_corrected
IoC <- fitIoC(x=pred_df %>% as.data.frame())
pred_df$IoC <- IoC$IoC


# COVERAGE x IoC ---------------------------------------------------------------


m <- pred_df %>%
  mutate(Genome=Name) %>%
  subset(select=c(Genome,dCUB,d,IoC,nCAZy,nGenes))
ram <- merge.easy(ra,m,key="Genome") %>%
  subset(!is.na(dCUB) & !is.na(value.pre) & !is.na(value.post))
ram$delta <- (ram$value.post - ram$value.pre)
ram$source <- "Pre-Thaw"
ram$source[ram$pre_post==T] <- "Post-Thaw"
ram_g <- ram %>% 
  group_by(Genome) %>% 
  summarise_all(mean) %>%
  subset(delta!=0)

# ABSOLUTE ABUNDANCES ----------------------------------------------------------

abund <- read.csv("16Scopy_number.csv") %>%
  subset(select=c(Sample_ID,copy_numb_gsoil))
names(abund)[1] <- "sample_ID"
meta <- merge.easy(meta,abund,key="sample_ID")

aba <- coverage[,c(grep("Relative.Abundance",names(coverage)))] 
aba <- t(t(aba)/rowSums(t(aba)))
aba[is.na(aba)] <- 0
aba <- aba %>% 
  as.data.frame() %>%
  mutate(Genome=coverage$Genome) %>%
  as.data.table() %>%
  melt(id.vars = "Genome")
aba$fastq_id <- aba$variable %>%
  gsub(pattern="X",replace="") %>%
  str_sub(., end = -2) %>%
  substr(.,1,4)
aba <- merge.easy(aba,meta,key="fastq_id") 
aba$Absolute_16S_copy_gsoil <- aba$value*aba$copy_numb_gsoil

aba <- aba %>%
  subset(select=c(Genome,Absolute_16S_copy_gsoil,core,pre_post_thaw,id)) %>%
  subset(!is.na(pre_post_thaw)) %>%
  reshape(direction = "wide",timevar = "pre_post_thaw",
          idvar = c("Genome","core","id"))

m <- pred_df %>%
  mutate(Genome=Name) %>%
  subset(select=c(Genome,dCUB,d,IoC,nCAZy,nGenes))
abam <- merge.easy(aba,m,key="Genome") %>%
  subset(!is.na(dCUB) & !is.na(Absolute_16S_copy_gsoil.pre) & !is.na(Absolute_16S_copy_gsoil.post))
abam$delta <- (abam$Absolute_16S_copy_gsoil.post - abam$Absolute_16S_copy_gsoil.pre)
abam$source <- "Pre-Thaw"
abam$source[abam$pre_post==T] <- "Post-Thaw"
abam$deltarel <- (abam$Absolute_16S_copy_gsoil.post-abam$Absolute_16S_copy_gsoil.pre)/(abam$Absolute_16S_copy_gsoil.post+abam$Absolute_16S_copy_gsoil.pre)
abam_g <- abam %>% 
  group_by(Genome)%>%
  summarise_all(mean,na.rm=T)

#avoid overwriting
pred_df$Source <- "SMAG DB"
pred_df$Source[pred_df$Name %in% woodcroft$Name] <- "Woodcroft et al."
pred_df$Source[pred_df$Name %in% permafrost$Name] <- "Present Study"
mag_df <- pred_df

################################################################################
# Mags Binned From This Study Only ---------------------------------------------
################################################################################

qual <- read.delim("quality_report_contained.tsv")
permafrost <- merge.easy(permafrost,qual,key="Name") 
permafrost$nGenes_corrected <- 100*permafrost$nGenes/permafrost$Completeness
permafrost$nGenes <- permafrost$nGenes_corrected
permafrost$IoC <- predictIoC(x=permafrost %>% as.data.frame(),
                              model = IoC)
permafrost <- permafrost %>%
  subset(nHE>=10)


# Relative Abundance -----------------------------------------------------------

coverage <- read.delim("bin_coverage_coverm_magmax.tsv")

permafrost <- permafrost %>% subset(nHE>=10)

ra_cov_permafrost <- coverage[,c(grep("Relative.Abundance",names(coverage)))] 
ra_cov_permafrost <- t(t(ra_cov_permafrost)/rowSums(t(ra_cov_permafrost)))
ra_cov_permafrost[is.na(ra_cov_permafrost)] <- 0
ra_cov_permafrost <- ra_cov_permafrost %>% 
  t() %>%
  t() %>%
  as.data.frame() %>%
  mutate(Genome=coverage$Genome) %>%
  as.data.table() %>% 
  melt(id.vars = "Genome") %>%
  subset(Genome %in% permafrost$Name)
ra_cov_permafrost$fastq_id <- ra_cov_permafrost$variable %>%
  gsub(pattern="1.fastq.Relative.Abundance....",replace="") %>%
  gsub(pattern="X",replace="") %>%
  str_sub(., end = -2)
ra_cov_permafrost <- ra_cov_permafrost %>% 
  subset(Genome %in% permafrost$Name) %>%
  subset(select=c(fastq_id,value)) %>%
  group_by(fastq_id) %>%
  summarise(RA=sum(value))

ra_cov_permafrost$id <- ra_cov_permafrost$fastq_id %>% 
  gsub(pattern="_.*",replace="") %>%
  gsub(pattern="p",replace="") 

median(ra_cov$RA[gsub("_.*","",ra_cov$fastq_id) %in% metagenome_df$ID])
median(ra_cov_permafrost$RA[gsub("_.*","",ra_cov_permafrost$fastq_id) %in% metagenome_df$ID])

ra_permafrost <- coverage[,c(grep("Relative.Abundance",names(coverage)))] 
ra_permafrost <- t(t(ra_permafrost)/rowSums(t(ra_permafrost)))
ra_permafrost[is.na(ra_permafrost)] <- 0
ra_permafrost <- ra_permafrost %>% 
  t() %>%
  mclr() %>%
  t() %>%
  as.data.frame() %>%
  mutate(Genome=coverage$Genome) %>%
  as.data.table() %>%
  melt(id.vars = "Genome")
ra_permafrost$fastq_id <- ra_permafrost$variable %>%
  gsub(pattern="1.fastq.Relative.Abundance....",replace="") %>%
  gsub(pattern="X",replace="") %>%
  str_sub(., end = -2) %>%
  substr(.,1,4)
ra_permafrost <- merge.easy(ra_permafrost,meta,key="fastq_id") %>%
  subset(select=c(Genome,core,pre_post_thaw,id,value)) %>%
  subset(!is.na(pre_post_thaw)) %>%
  reshape(direction = "wide",timevar = "pre_post_thaw",
          idvar = c("Genome","core","id"))
ra_permafrost[is.na(ra_permafrost)] <- 0

m <- permafrost %>%
  mutate(Genome=Name) %>%
  subset(select=c(Genome,dCUB,d,IoC,nCAZy,nGenes))
ram_permafrost <- merge.easy(ra_permafrost,m,key="Genome") %>%
  subset(!is.na(dCUB) & !is.na(value.pre) & !is.na(value.post))
ram_permafrost$delta <- (ram_permafrost$value.post - ram_permafrost$value.pre)
ram_permafrost$source <- "Pre-Thaw"
ram_permafrost$source[ram_permafrost$pre_post==T] <- "Post-Thaw"
ram_g_permafrost <- ram_permafrost %>% 
  group_by(Genome) %>% 
  summarise_all(mean) %>%
  subset(delta!=0)

# Absolute Abundance -----------------------------------------------------------

aba_permafrost <- coverage[,c(grep("Relative.Abundance",names(coverage)))] 
aba_permafrost <- t(t(aba_permafrost)/rowSums(t(aba_permafrost)))
aba_permafrost[is.na(aba_permafrost)] <- 0
aba_permafrost <- aba_permafrost %>% 
  as.data.frame() %>%
  mutate(Genome=coverage$Genome) %>%
  as.data.table() %>%
  melt(id.vars = "Genome")
aba_permafrost$fastq_id <- aba_permafrost$variable %>%
  gsub(pattern="X",replace="") %>%
  str_sub(., end = -2) %>%
  substr(.,1,4)
aba_permafrost <- merge.easy(aba_permafrost,meta,key="fastq_id") 
aba_permafrost$Absolute_16S_copy_gsoil <- aba_permafrost$value*aba_permafrost$copy_numb_gsoil

aba_permafrost <- aba_permafrost %>%
  subset(select=c(Genome,Absolute_16S_copy_gsoil,core,pre_post_thaw,id)) %>%
  subset(!is.na(pre_post_thaw)) %>%
  reshape(direction = "wide",timevar = "pre_post_thaw",
          idvar = c("Genome","core","id"))

m <- permafrost %>%
  mutate(Genome=Name) %>%
  subset(select=c(Genome,dCUB,d,IoC,nCAZy,nGenes))
abam_permafrost <- merge.easy(aba_permafrost,m,key="Genome") %>%
  subset(!is.na(dCUB) & !is.na(Absolute_16S_copy_gsoil.pre) & !is.na(Absolute_16S_copy_gsoil.post))
abam_permafrost$delta <- (abam_permafrost$Absolute_16S_copy_gsoil.post - abam_permafrost$Absolute_16S_copy_gsoil.pre)
abam_permafrost$source <- "Pre-Thaw"
abam_permafrost$source[abam_permafrost$pre_post==T] <- "Post-Thaw"
abam_permafrost$deltarel <- (abam_permafrost$Absolute_16S_copy_gsoil.post-abam_permafrost$Absolute_16S_copy_gsoil.pre)/(abam_permafrost$Absolute_16S_copy_gsoil.post+abam_permafrost$Absolute_16S_copy_gsoil.pre)
abam_g_permafrost <- abam_permafrost %>% 
  group_by(Genome)%>%
  summarise_all(mean,na.rm=T)

# Cluster IoC ------------------------------------------------------------------

#Clustering
mC <- Mclust(permafrost$IoC,
             verbose=T)
p <- (colSums(mC$z)/sum(mC$z))

#Confidence limits
c_l <- max(mC$data[mC$uncertainty<0.05 & mC$data<4])
c_h <- min(mC$data[mC$uncertainty<0.05 & mC$data>4])
c_m <- mC$data[which.max(mC$uncertainty)]

permafrost$IoCc <- "Oligotroph"
permafrost$IoCc[permafrost$IoC>c_l] <- "Undefined"
permafrost$IoCc[permafrost$IoC>c_h] <- "Copiotroph"
permafrost$IoCc <- factor(permafrost$IoCc,
                          levels=c("Oligotroph",
                                   "Undefined",
                                   "Copiotroph"))

mag_df$IoCc <- "Oligotroph"
mag_df$IoCc[mag_df$IoC>c_l] <- "Undefined"
mag_df$IoCc[mag_df$IoC>c_h] <- "Copiotroph"
mag_df$IoCc <- factor(mag_df$IoCc,
                       levels=c("Oligotroph",
                                "Undefined",
                                "Copiotroph"))

warming_contained_mags_ioc <- list(IoCmod=IoC)
save(warming_contained_mags_ioc,file="warming_contained_mags_ioc.rda")

################################################################################
# Metagenome-Level Data --------------------------------------------------------
################################################################################

## Read Counts ----

reads <- read.delim("clean_read_counts.tsv",
                    header=F)
reads$ID <- reads$V1 %>%
  basename() %>%
  gsub(pattern="_.*",
       replace="") 
reads <- reads %>% 
  group_by(ID) %>%
  summarise(reads=sum(V2),
            bp=sum(V2)*150) %>%
  mutate(gbp=bp/1e9)


## Single Copy Marker Gene Normalized nGenes and nCAZy -----

setwd("scmg_analysis/")
id_vec <- list.files()[grep("cov",list.files())] %>% 
  gsub(pattern=".cov",replace="")
setwd("../")
scmg_df <- data_frame(fastq_id=character(length(id_vec)),
                      nGenes=numeric(length(id_vec)),
                      nCAZy=numeric(length(id_vec)),
                      rCAZy=numeric(length(id_vec)))
for(i in 1:length(id_vec)){
  print(id_vec[i])
  cazy <- readLines(paste0("scmg_analysis/",id_vec[i],".cazy"))
  
  if(length(cazy)>0){
    scmg <- read.delim(paste0("scmg_analysis/",id_vec[i],".scmg"),header = F)
    names(scmg) <- c("SCMG","Gene")
    coverage <- read.delim(paste0("scmg_analysis/",
                                  id_vec[i],
                                  ".cov"),
                           header = F) %>%
      mutate(Gene=V1,AvgCov=V7) %>%
      subset(select=c(Gene,AvgCov))
    coverage$CAZy <- 0
    coverage$CAZy[coverage$Gene %in% cazy] <- 1
    cov_scmg <- merge.easy(coverage,scmg,key="Gene")
    
    tot_cov <- sum(cov_scmg$AvgCov)
    tot_cazy <- sum(cov_scmg$AvgCov[cov_scmg$CAZy==1])
    med_scmg <- cov_scmg %>%
      subset(!is.na(SCMG)) %>%
      group_by(SCMG) %>%
      summarize(cov=sum(AvgCov)) %>%
      pull(cov) %>%
      median()
    num_genes <- tot_cov/med_scmg
    num_cazy <- tot_cazy/med_scmg
    
    scmg_df$fastq_id[i] <- id_vec[i]
    scmg_df$nGenes[i] <- tot_cov/med_scmg
    scmg_df$nCAZy[i] <- tot_cazy/med_scmg
    scmg_df$rCAZy[i] <- num_cazy/num_genes 
  } else {
    scmg_df$fastq_id[i] <- id_vec[i]
    scmg_df$nGenes[i] <- NA
    scmg_df$nCAZy[i] <- NA
    scmg_df$rCAZy[i] <- NA 
  }
}

## gRodon Output -----

load("warming_contained_mags_ioc.rda")
meta <- read.delim("warming_metadata.tsv") %>%
  mutate(fastq_id=str_sub(fastq_id, end = -2))
meta$ID <- meta$fastq_id %>% 
  gsub(pattern="_.*",replace="")
meta$id <- meta$fastq_id %>% 
  gsub(pattern="_.*",replace="") %>%
  gsub(pattern="p",replace="")
meta <- merge.easy(meta,reads,key="ID")
load("warming_gRodon_predictions_ms_scaffolds.rda")
pred_scaff <- pred_df %>%
  mutate(fastq_id=basename(Genome) %>%
           gsub(pattern="*.fasta.ffn",replace=""))

## Metagenome Stats ----

metagenome_stats <- merge.easy(pred_scaff,meta,key="fastq_id") %>%
  subset(select=c(fastq_id,sample_ID,site,core,pre_post_thaw,id,reads,gbp,Concentration..ng.ul.,nHE)) %>%
  mutate(Included.In.Analysis=nHE>=10) %>%
  mutate(fastq_id=gsub("_.*","",fastq_id))
write.csv(metagenome_stats,file="metagenome_stats.csv")

## Combine & Predict IoC ----

metagenome_df <- merge.easy(scmg_df,pred_scaff,key="fastq_id") %>%
  merge.easy(.,meta,key="fastq_id") %>%
  mutate(nGenes=nGenes.x) %>%
  subset(nHE>=10)
metagenome_df$IoC <- predictIoC(as.data.frame(metagenome_df),
                                model=warming_contained_mags_ioc$IoCmod)
metagenome_df$IoCc <- "Oligotroph"
metagenome_df$IoCc[metagenome_df$IoC>c_l] <- "Undefined"
metagenome_df$IoCc[metagenome_df$IoC>c_h] <- "Copiotroph"
metagenome_df$IoCc <- factor(metagenome_df$IoCc,
                 levels=c("Oligotroph",
                          "Undefined",
                          "Copiotroph"))

################################################################################
# Save Dataset for Analysis ----------------------------------------------------
################################################################################

save(mC,
     c_m,
     c_l,
     c_h,
     IoC,
     mag_df,
     ra_cov,
     abam_g,
     abam,
     ram_g,
     ram,
     permafrost,
     abam_g_permafrost,
     abam_g_permafrost,
     ram_g_permafrost,
     ram_permafrost,
     metagenome_df,
     file = "warming_data_for_analysis.rda")


