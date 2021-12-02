
################## Lectura #############
library(data.table)
BAM_Normal_procesado <- fread("Secuencias_SRR14562364_STAR.txt",sep = "\t",select = c(1:6))
colnames(BAM_Normal_procesado) <- c("name","rname","pos","cigar","algo3","seq")

library(readr)
rangos_genes_hg19 <- read_delim("C:/Users/anton/OneDrive/Documentos/rangos_genes_hg19.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
rangos_genes_hg19 <- rangos_genes_hg19[!duplicated(rangos_genes_hg19), ]
colnames(rangos_genes_hg19)[1] <- "chrom"

#########################################
################# Salida para sam2pairwise ##############

library(GenomicAlignments)
library(data.table)
insercciones <- grepl(pattern = "I",BAM_Normal_procesado$cigar)
insercciones <- which(insercciones==TRUE)
BAM_Normal_procesado_inserciones <- BAM_Normal_procesado[insercciones,]

nombres <- BAM_Normal_procesado_inserciones$name
write.table(nombres,file = "nombres_inserciones.txt",row.names = FALSE,quote = FALSE,col.names = FALSE)

#########################################################


############## Modificacion de sam2pairwise ############
file <- as.matrix(file)
matriz_modificada <- matrix(nrow = nrow(file)/4, ncol = 11)
matriz_modificada[1,1:9] <- file[1,]
matriz_modificada[1,10] <- file[2,1]
matriz_modificada[1,11] <- file[4,1]
numero1 <- 1
numero2 <- 2
numero3 <- 4
for (n in 1:((nrow(file)/4)-1)){
  numero1 <- numero1 + 4
  numero2 <- numero2 + 4
  numero3 <- numero3 + 4
  matriz_modificada[n+1,1:9] <- file[numero1,]
  matriz_modificada[n+1,10] <- file[numero2,1]
  matriz_modificada[n+1,11] <- file[numero3,1]
}
matriz_modificada <- as.data.frame(matriz_modificada)
colnames(matriz_modificada) <- c("name","algo","rname","pos","algo2","cigar","igual","pos2","algo3","seq","WT")
matriz_modificada$pos <- as.numeric(matriz_modificada$pos)
matriz_modificada$pos2 <- as.numeric(matriz_modificada$pos2)
matriz_modificada$algo <- as.numeric(matriz_modificada$algo)
matriz_modificada$algo2 <- as.numeric(matriz_modificada$algo2)
matriz_modificada$algo3 <- as.numeric(matriz_modificada$algo3)

matriz_modificada$longitud_fragmento <- nchar(matriz_modificada$seq)
matriz_modificada$end <- matriz_modificada$pos + nchar(matriz_modificada$longitud_fragmento)-1
matriz_modificada <- as.data.table(matriz_modificada)

insercciones <- grepl(pattern = "I",matriz_modificada$cigar)
insercciones <- which(insercciones==TRUE)
BAM_Normal_procesado_inserciones <- matriz_modificada[insercciones,]
########################################################


################## Inserciones totales ################## 



#nombres <- BAM_Normal_procesado_inserciones$name
#write.table(nombres,file = "nombres_inserciones.txt",row.names = FALSE,quote = FALSE,col.names = FALSE)
cigar_inserciones <- BAM_Normal_procesado_inserciones$cigar
cigar_inserciones <- cigarToRleList(cigar_inserciones)


#matriz_inserciones <- matrix(nrow = 500000, ncol=4)
#matriz_inserciones <- as.data.frame(matriz_inserciones)
#colnames(matriz_inserciones) <- c("longitud_inserccion","chr","pos_inicio","pos_final")


valor_insercciones2 <- runValue(cigar_inserciones)
valor_insercciones2 <- as.list(valor_insercciones2)
valor_insercciones2 <- as.data.frame(do.call(qpcR:::rbind.na, valor_insercciones2))

longitud2 <- runLength(cigar_inserciones)
longitud2 <- as.list(longitud2)
longitud2 <- as.data.frame(do.call(qpcR:::rbind.na, longitud2))

BAM_Normal_procesado_inserciones$longitud_fragmento <- apply(longitud2,1,sum,na.rm=T)
BAM_Normal_procesado_inserciones$end <- BAM_Normal_procesado_inserciones$pos + BAM_Normal_procesado_inserciones$longitud_fragmento -1

eliminar <- which(valor_insercciones2$V1=="S")
valor_insercciones2[eliminar,] <- cbind(valor_insercciones2[eliminar,2:ncol(valor_insercciones2)],NA)
longitud2[eliminar,] <- cbind(longitud2[eliminar,2:ncol(longitud2)],NA)

BAM_Normal_procesado_inserciones$longitud_fragmento <- apply(longitud2,1,sum,na.rm=T)
BAM_Normal_procesado_inserciones$end <- BAM_Normal_procesado_inserciones$pos + BAM_Normal_procesado_inserciones$longitud_fragmento -1

buscar_INDELS <- function(valor,longitud,tipo,bam,bam2,media_long,genes_hg19){
  valor_insercciones2 <- valor
  BAM_Normal_procesado_inserciones <- bam
  longitud2 <- longitud
  longitud2[is.na(longitud2)] <- 0
  longitud2_cumsum <- t(apply(longitud2,1,cumsum))
  
  test1 <- which(valor_insercciones2==tipo,arr.ind = T)
  test1 <- test1[order(test1[,1]),]
  
  longitud_inserccion <- longitud2[test1]
  chr <- BAM_Normal_procesado_inserciones$rname[test1[,1]]
  posicion <- BAM_Normal_procesado_inserciones$pos[test1[,1]]
  test2 <- test1 
  test2[,2] <- test2[,2]-1
  pos_inicio <- posicion + longitud2_cumsum[test2]
  pos_final <- posicion + longitud2_cumsum[test1]
  
  matriz_inserciones5 <- data.frame(cbind(longitud_inserccion,chr,pos_inicio,pos_final))
  matriz_inserciones5 <- na.omit(matriz_inserciones5)
  library(plyr)
  matriz_inserciones6 <- ddply(matriz_inserciones5,.(longitud_inserccion,chr,pos_inicio,pos_final),nrow)
  #matriz_inserciones6 <- matriz_inserciones6[-which(matriz_inserciones6$longitud_inserccion==0),]
  matriz_inserciones6$pos_final <- as.numeric(matriz_inserciones6$pos_final)
  matriz_inserciones6$pos_inicio <- as.numeric(matriz_inserciones6$pos_inicio)
  colnames(matriz_inserciones6)[5] <- "Numero_Mutadas"
  matriz_inserciones6 <- matriz_inserciones6[which(matriz_inserciones6$Numero_Mutadas>=5),]
  matriz_inserciones6 <- as.data.table(matriz_inserciones6)
  
  #setkey(bam2, rname, pos)
  #cromosoma <- matriz_inserciones6$chr
  #inicio <- matriz_inserciones6$pos_inicio
  output <- list()
  contador <- 0
  cromosomas <- unique(matriz_inserciones6$chr)
  for (z in 1:length(cromosomas)){
    bam <- BAM_Normal_procesado[rname==cromosomas[z]]
    bed <- matriz_inserciones6[chr==cromosomas[z]]
    for (n in 1:nrow(bed)){
      contador <- contador + 1
      print(contador)
      mini_df <- bam[end>=(bed$pos_final[n]-1) & pos<= bed$pos_inicio[n]]
      inserciones <- grepl(pattern = "I",mini_df$cigar)
      eliminar <- 0
      eliminar2 <- length(inserciones)
      for (w in 1:length(inserciones)){
        if (inserciones[w]==FALSE){
          eliminar <- eliminar+1
        }
        else (break)
      }
      for (y in length(inserciones):1){
        if (inserciones[y]==FALSE){
          eliminar2 <- eliminar2-1
        }
        else (break)
      }
      inserciones <- inserciones[(eliminar+1):eliminar2]
      output[[contador]] <- length(inserciones)
    }
  }
  matriz_inserciones6$Totales <- unlist(output)
  
  final <- matriz_inserciones6$pos_final
  salida <- list()
  for(n in 1:nrow(matriz_inserciones6)){
    print(n)
    a <- genes_hg19 %>% filter(chrom==cromosoma[n],txStart<=inicio[n],txEnd>=final[n])
    salida[[n]] <- a$name2
  }
  salida <- lapply(salida,paste,collapse=",")
  salida <- as.data.frame(do.call(qpcR:::rbind.na, salida))
  matriz_inserciones6 <- cbind(matriz_inserciones6,salida)
  return(matriz_inserciones6)
}

matriz_inserciones_SRR14562364 <- buscar_INDELS(valor_insercciones2,longitud2,"I",BAM_Normal_procesado_inserciones,BAM_Normal_procesado,200,rangos_genes_hg19)

######################################################################## 

################## Delecciones totales ################## 

library(GenomicAlignments)
library(data.table)
delecciones <- grepl(pattern = "D",BAM_Normal_procesado$cigar)
delecciones <- which(delecciones==TRUE)
BAM_Normal_procesado_delecciones <- BAM_Normal_procesado[delecciones,]
cigar_delecciones <- BAM_Normal_procesado_delecciones$cigar
cigar_delecciones <- cigarToRleList(cigar_delecciones)

valor_delecciones2 <- runValue(cigar_delecciones)
valor_delecciones2 <- as.list(valor_delecciones2)
valor_delecciones2 <- as.data.frame(do.call(qpcR:::rbind.na, valor_delecciones2))

longitud_delecciones2 <- runLength(cigar_delecciones)
longitud_delecciones2 <- as.list(longitud_delecciones2)
longitud_delecciones2 <- as.data.frame(do.call(qpcR:::rbind.na, longitud_delecciones2))

BAM_Normal_procesado_delecciones$longitud_fragmento <- apply(longitud_delecciones2,1,sum,na.rm=T)
BAM_Normal_procesado_delecciones$end <- BAM_Normal_procesado_delecciones$pos + BAM_Normal_procesado_delecciones$longitud_fragmento -1

eliminar <- which(valor_delecciones2$V1=="S")
valor_delecciones2[eliminar,] <- cbind(valor_delecciones2[eliminar,2:ncol(valor_delecciones2)],NA)
longitud_delecciones2[eliminar,] <- cbind(longitud_delecciones2[eliminar,2:ncol(longitud_delecciones2)],NA)

BAM_Normal_procesado_delecciones$longitud_fragmento <- apply(longitud_delecciones2,1,sum,na.rm=T)
BAM_Normal_procesado_delecciones$end <- BAM_Normal_procesado_delecciones$pos + BAM_Normal_procesado_delecciones$longitud_fragmento -1

matriz_deleciones_SRR1751229 <- buscar_INDELS(valor_delecciones2,longitud_delecciones2,"D",BAM_Normal_procesado_delecciones,BAM_Normal_procesado,200,rangos_genes_hg19)

######################################################################## 
prueba <- merge(matriz_inserciones2,matriz_inserciones_totales_SRR1751230,by=c("longitud_inserccion","chr","pos_inicio","pos_final"))


reduccir_tamaño_bam2 <- function(bam,bed){
  setkey(bam, rname, pos)
  cromosoma <- bed$chr
  inicio <- bed$pos_inicio
  output <- list()
    for (n in 1:nrow(bed)){
      print(n)
      output[[n]] <- nrow(na.omit(bam[J(cromosoma[n],(inicio[n]-200):inicio[n])]))
    }
  return(output)
}
bam_inserciones <- reduccir_tamaño_bam2(BAM_Normal_procesado,matriz_inserciones_SRR1751229)
bam_deleciones <- reduccir_tamaño_bam2(BAM_Normal_procesado,matriz_delecciones_SRR1751230)

matriz_inserciones_SRR1751229$Totales <- unlist(bam_inserciones)
matriz_delecciones_SRR1751230$Totales <- unlist(bam_deleciones)

################## Filtrar ####################
rangos_genes_hg19 <- rangos_genes_hg19[!duplicated(rangos_genes_hg19), ]
colnames(rangos_genes_hg19)[1] <- "chrom"




rapido <- function(x){
  x$porcentaje <- (x$Numero_Mutadas/x$Totales)*100
  return(x)
}
  
matriz_deleciones_SRR1751229 <- rapido(matriz_deleciones_SRR1751229)
matriz_deleciones_SRR1751230 <- rapido(matriz_deleciones_SRR1751230)
matriz_inserciones_SRR1751229 <- rapido(matriz_inserciones_SRR1751229)
matriz_inserciones_SRR1751230 <- rapido(matriz_inserciones_SRR1751230)


prueba <- matriz_deleciones_SRR1751230[which(matriz_deleciones_SRR1751230$porcentaje<=99 & matriz_deleciones_SRR1751230$porcentaje>=1),]

genes <- c(panel2$proteina,panel2$p2,panel2$p3)
genes <- unique(genes)
genes2 <- strsplit(genes, split = ":")
genes2 <- as.data.frame(do.call(qpcR:::rbind.na, genes2))
genes2 <- genes2$V1

pruebesita <- c()
for(n in 1:length(genes2)){
 pruebesita <- c(pruebesita,grep(genes2[n],matriz_inserciones_SRR1751229$V1))
 pruebesita <- unique(pruebesita)
}


matriz_deleciones_SRR1751230_malas <- matriz_deleciones_SRR1751230[pruebesita,]
matriz_deleciones_SRR1751229_malas <- matriz_deleciones_SRR1751229[pruebesita,]
matriz_inserciones_SRR1751230_malas <- matriz_inserciones_SRR1751230[pruebesita,]
matriz_inserciones_SRR1751229_malas <- matriz_inserciones_SRR1751229[pruebesita,]


prueba2 <- merge(matriz_inserciones_SRR1751229_malas,matriz_inserciones_SRR1751230_malas,by=c("chr","pos_inicio","pos_final"),all.x = T,all.y = T)
