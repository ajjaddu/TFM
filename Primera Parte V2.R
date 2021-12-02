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
alineamiento <- "file.out"
setwd("D:/TFM/Qiime2/Compartida/GSE33328 (Datos RNA-Seq tejido tumoral)")

################ Librerias ################
library(readr)
library(data.table)
#library(CrispRVariants)
library(Biostrings)
library(dplyr)

#library(GenomicAlignments)
#library(DuffyTools)

#################### Funciones ############
seqsToAln <- function(cigar, dnaseq, target, del_char = "-", 
                      aln_start = NULL, reverse_complement = FALSE,
                      allow_partial = FALSE){
  # Additional trimming is necessary because deletion operations may overhang
  # boundaries of target region
  
  if (as.character(GenomicRanges::strand(target)) == "*") strand <- "+"
  
  result <- GenomicAlignments::sequenceLayer(dnaseq, cigar,
                                             D.letter = del_char,
                                             N.letter = "N")
  if (! is.null(aln_start)){
    shifts <- aln_start - GenomicAlignments::start(target)
    result <- Biostrings::stackStrings(result, shift = shifts,
                                       from = 1, to = width(target),
                                       Lpadding.letter = "+",
                                       Rpadding.letter = ".")
    if (! isTRUE(allow_partial) & isTRUE(any(grep("\\+|\\.", result)))){
      stop(paste("When allow_partial is FALSE, dnaseq to be",
                 "trimmed must span the target location"))
    }
  }
  
  if (isTRUE(reverse_complement)){
    result <- Biostrings::reverseComplement(result)
  }
  result <- as.character(result)
  result <- gsub("\\.", "<", gsub("\\+", ">", result))
  result
} # ------

###########################################

############## Cargar archivos ############
BAM_Normal_procesado <- read_delim("Secuencias_SRR1751231_STAR.txt", 
                                   "\t", escape_double = FALSE, col_names = FALSE, 
                                   trim_ws = TRUE)
BAM_Normal_procesado <- fread("Secuencias_SRR14562364_STAR.txt",sep = "\t")

colnames(BAM_Normal_procesado) <- c("name","rname","pos","cigar","algo3","seq","ascii")
BAM_Normal_procesado <- as.data.table(BAM_Normal_procesado)
cromosomas <- c(paste("chr",1:22,sep = ""),"chrX","chrY")
BAM_Normal_procesado <- BAM_Normal_procesado[rname %in% cromosomas]
cromosomas <- unique(BAM_Normal_procesado$rname)
gc()
for(n in 1:length(cromosomas)){
  print(n)
  BAM_Normal_procesado[rname==cromosomas[n]]$seq <- with(BAM_Normal_procesado[rname==cromosomas[n]],seqsToAln(cigar,dnaseq = DNAStringSet(seq),target = "+"))
  gc()
}
BAM_Normal_procesado$longitud_fragmento <- nchar(BAM_Normal_procesado$seq)
BAM_Normal_procesado$end <- BAM_Normal_procesado$pos + BAM_Normal_procesado$longitud_fragmento -1
BAM_Normal_procesado <- BAM_Normal_procesado[longitud_fragmento<=50000]
gc()
###########################################

################### TCGA ###################

if(Metodo=="TCGA"){
  #bedfile <- load(bedfile)
  genes <- c("MALAT1:378938",
             "AHNAK:79026",
             "GATA3:2625",
             "RNF213:57674",
             "MYH9:4627",
             "CDH1:999",
             "FOXA1:3169",
             "USP34:9736",
             "MACF1:23499",
             "PRKDC:5591",
             "ARID1A:8289",
             "HUWE1:10075",
             "PTEN:5728",
             "DYNC1H1:1778",
             "AKAP9:10142",
             "KMT2C:58508",
             "BIRC6:57448",
             "VPS13B:157680",
             "MAP3K1:4214",
             "KMT2D:8085")
  
  ### Hacemos un subset de los BAM teniendo en cuenta solo los cromosomas del BED y las regiones de 
  ### los BAM que se encuentren dentro de las regiones de los genes del BED. 
  reduccir_tamaño_bam <- function(bam,bed,genes){
    bam <- bam %>% filter(rname %in% bedfile$chr)
    desired_length <- length(genes) # or whatever length you want
    output <- vector(mode = "list", length = desired_length)
    for (n in 1:length(genes)){
      print(n)
      posicion_gen <- bed %>% filter(proteina == genes[n] | p2 == genes[n] | p3 == genes[n])
      rango_gen <-  (min(posicion_gen$start)-ceiling(mean(bam$longitud_fragmento))):max(posicion_gen$start)
      a <- bam %>% filter(rname == posicion_gen$chr[1],pos %in% rango_gen )
      output[[n]] <- bind_cols(a,rep(genes[n],nrow(a)))
      gc()
    }
    names(output) <- genes
    return(output)
  }
  
  ### Reducimos el tamaño de los BAM separandolos en DF distintos para cada gen.
  BAM_Normal_procesado_reducido_por_genes <- reduccir_tamaño_bam(BAM_Normal_procesado,bedfile,genes)
  
  #### Extraemos los nombres de los alineamiento reducidos para volver a reducir el bam desde linux. 
  BAM_Normal_procesado_reducido_por_genes2 <- rbindlist(BAM_Normal_procesado_reducido_por_genes)
  nombres <- as.data.frame(BAM_Normal_procesado_reducido_por_genes2$name)
  write.table(nombres,file = "nombres_alineamiento_positivo_TCGA.txt",row.names = FALSE,quote = FALSE,col.names = FALSE)
  save(BAM_Normal_procesado_reducido_por_genes,file = "BAM_Normal_procesado_reducido_por_genes_TCGA.RData")
}

############################################

################### ClinVar ###################

if (Metodo=="ClinVar"){
  #bedfile <- load(bedfile_clinvar)
  genes <- unique(bedfile$proteina)
  
  reduccir_tamaño_bam <- function(bam,bed,genes){
    bam <- bam %>% filter(rname %in% bedfile$chr)
    desired_length <- length(genes) # or whatever length you want
    output <- vector(mode = "list", length = desired_length)
    for (n in 1:length(genes)){
      posicion_gen <- bed %>% filter(proteina == genes[n] | p2 == genes[n] | p3 == genes[n])
      rango_gen <-  (min(posicion_gen$start)-ceiling(mean(bam$longitud_fragmento))):max(posicion_gen$start)
      a <- bam %>% filter(rname == posicion_gen$chr[1],pos %in% rango_gen )
      output[[n]] <- bind_cols(a,rep(genes[n],nrow(a)))
    }
    names(output) <- genes
    return(output)
  }
  
  ### Reducimos el tamaño de los BAM separandolos en DF distintos para cada gen.
  BAM_Normal_procesado_reducido_por_genes <- reduccir_tamaño_bam(BAM_Normal_procesado,bedfile,genes)
  
  #### Extraemos los nombres de los alineamiento reducidos para volver a reducir el bam desde linux. 
  BAM_Normal_procesado_reducido_por_genes2 <- rbindlist(BAM_Normal_procesado_reducido_por_genes)
  nombres <- as.data.frame(BAM_Normal_procesado_reducido_por_genes2$name)
  write.table(nombres,file = "nombres_alineamiento_negativp_clinvar.txt",row.names = FALSE,quote = FALSE,col.names = FALSE)
  save(BAM_Normal_procesado_reducido_por_genes,file = "BAM_Normal_procesado_reducido_por_genes.RData")
}

###############################################
