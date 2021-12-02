valor,longitud,tipo,bam,bam2,media_long,genes_hg19

valor_insercciones2,longitud2,"I",BAM_Normal_procesado_inserciones,
BAM_Normal_procesado,200,rangos_genes_hg19

valor <- valor_insercciones2
longitud <- longitud2
tipo <- "I"
bam <- BAM_Normal_procesado_inserciones
bam2 <- BAM_Normal_procesado
media_long <- 200

########################################

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
pos_final <- posicion + longitud2_cumsum[test1] - 1

test4 <- as.data.frame(test1)
test3 <- ddply(test4,.(row),nrow)
BAM_Normal_procesado_inserciones$rep <- test3$V1
BAM_Normal_procesado_inserciones$row <- test3$row
probosita <- with(BAM_Normal_procesado_inserciones,BAM_Normal_procesado_inserciones[rep(row,rep),])
probosita$inicio <- longitud2_cumsum[test2] + 1 
probosita$final <- longitud2_cumsum[test1] 
bases <- with(probosita,substr(seq,(inicio),(final)))

matriz_inserciones5 <- data.frame(cbind(longitud_inserccion,chr,pos_inicio,pos_final,bases))
matriz_inserciones5 <- na.omit(matriz_inserciones5)
library(plyr)
matriz_inserciones6 <- ddply(matriz_inserciones5,.(longitud_inserccion,chr,pos_inicio,pos_final,bases),nrow)
#matriz_inserciones6 <- matriz_inserciones6[-which(matriz_inserciones6$longitud_inserccion==0),]
matriz_inserciones6$pos_final <- as.numeric(matriz_inserciones6$pos_final)
matriz_inserciones6$pos_inicio <- as.numeric(matriz_inserciones6$pos_inicio)
colnames(matriz_inserciones6)[5] <- "Numero_Mutadas"
matriz_inserciones6 <- matriz_inserciones6[which(matriz_inserciones6$Numero_Mutadas>=5),]
matriz_inserciones6 <- as.data.table(matriz_inserciones6)
matriz_inserciones6$pos_final <- matriz_inserciones6$pos_final-1


bam <- BAM_Normal_procesado[rname=="chr2"]
mini_df <- bam[end>=(55199565) & pos<= 55199563]

esto <- with(mini_df,substr(seq,(55199563-pos+1),(55199565-pos+1)))
table(esto)


test3 <- ddply(test1,.(row),nrow)
BAM_Normal_procesado_inserciones$rep <- test3$V1
BAM_Normal_procesado_inserciones$row <- test3$row
probosita <- with(BAM_Normal_procesado_inserciones,BAM_Normal_procesado_inserciones[rep(row,rep),])
probosita$inicio <- longitud2_cumsum[test2] + 1 
probosita$final <- longitud2_cumsum[test1] 
esto <- with(probosita,substr(seq,(inicio),(final)))

probosita <- list()
for (n in 1:nrow(test1)){
  df[rep(1, 5), ]
}



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



totales <- with(matriz_inserciones_SRR14562364,nrow(BAM_Normal_procesado[rname==chr & end>=pos_final & pos<= pos_inicio]))

totales <- with(BAM_Normal_procesado,nrow(BAM_Normal_procesado[rname==chr & end>=pos_final & pos<= pos_inicio]))
