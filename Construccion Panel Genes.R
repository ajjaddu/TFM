#################################################################################
#################################################################################
#######################             TFM           ###############################
#######################              V3           ###############################
#################################################################################
#################################################################################
#################################################################################

################ Opciones ################
rm(list = ls())
Metodo = "TCGA"
cancer <- "breast"
setwd("D:/TFM/Qiime2/Compartida/Bowtie2")
###########################################

################ Librerias ################
library(dplyr)
library(data.table)
library(stringr)
library(vcfR)
library(tidyr)
library(bedr)
###########################################

################ CONSTRUCTOR ################
input <- list()
input[[1]] <- list.files(pattern = ".vcf.gz") 
names(input) <- c("vcf")
#############################################

############################## Panel BED TCGA #################################
if(Metodo=="TCGA"){
  my.vcf <- read.vcfR(input$vcf[2], verbose = F)
  my.vcf.df <- cbind(as.data.frame(getFIX(my.vcf)), INFO2df(my.vcf))
  
  panel <- my.vcf.df[,c(1:5,12)]
  panel <- cbind(panel,strsplit2matrix(panel$GENEINFO,split = "|",fixed = T))
  panel$GENEINFO <- NULL
  panel$CHROM <- paste("chr",panel$CHROM,sep = "")
  
  ## De este modo separamos los distintos alelos en varias filas una para cada uno.
  bedfile <- panel %>% separate_rows(ALT)
  
  names(bedfile) <- c("chr","start","rs","tipo","cambio","proteina","p2","p3")
  bedfile$start <- as.numeric(bedfile$start)
  save(bedfile,file = "panel_breast_htseq.RData")
  #write.csv(panel2,file = "panel_colorectal.csv")
}
###############################################################################

############################## Panel BED clinvar #################################
if(Metodo=="ClinVar"){
  my.vcf <- read.vcfR(input$vcf[[7]], verbose = F)
  my.vcf.df <- cbind(as.data.frame(getFIX(my.vcf)), INFO2df(my.vcf))
  
  panel <- my.vcf.df[,c(1:5,12,18,25)]
  cancer_specific <- grepl(pattern = "Pathogenic",panel$CLNSIG,ignore.case = TRUE)
  panel <- panel[which(cancer_specific==TRUE),]
  cancer_specific2 <- grepl(pattern = cancer,panel$CLNDN,ignore.case = TRUE)
  panel <- panel[which(cancer_specific2==TRUE),]
  panel$CLNDN <- NULL
  panel$CLNSIG <- NULL
  panel <- cbind(panel,strsplit2matrix(panel$GENEINFO,split = "|",fixed = T))
  panel$GENEINFO <- NULL
  panel$CHROM <- paste("chr",panel$CHROM,sep = "")
  
  ## De este modo separamos los distintos alelos en varias filas una para cada uno.
  bedfile <- panel %>% separate_rows(ALT)
  
  names(bedfile) <- c("chr","start","rs","tipo","cambio","proteina","p2","p3")
  bedfile$start <- as.numeric(bedfile$start)
  save(bedfile,file = "panel_breast_clinvar.RData")
  #write.csv(panel2,file = "panel_colorectal.csv")
}
##################################################################################