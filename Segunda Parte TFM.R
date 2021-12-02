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
alineamiento <- "file_TCGA.out"
setwd("D:/TFM/Qiime2/Compartida/GSE33328 (Datos RNA-Seq tejido tumoral)")
###########################################

################ Librerias ################
library(dplyr)
library(data.table)
library(stringr)
library(readr)
library(DuffyTools)
###########################################
  genes <- unique(bedfile$proteina)
  genes <- c("KRAS:3845",
             "BRAF:673",
             "TP53:7157",
             "APC:324",
             "PIK3CA:5290",
             "FBXW7:55294",
             "SMAD4:4089",
             "LRP1B:53353",
             "CTNNB1:1499",
             "FAT4:79633",
             "NRAS:4893",
             "TCF7L2:6934",
             "KMT2C:58508",
             "ATM:472",
             "RNF43:54894",
             "KMT2D:8085",
             "PTEN:5728",
             "ARID1A:8289",
             "AMER1:139285")
  
  genes <- c("TP53:7157",
             "PIK3CA:5290",
             "TTN:7273",
             "MUC4:4585",
             "MUC16:94025",
             "CDH1:999",
             "GATA3:2625",
             "KMT2C:58508",
             "RYR2:6262",
             "MAP3K1:4214",
             "HMCN1:83872",
             "SYNE1:23345",
             "USH2A:7399",
             "FLG:2312",
             "SPTA1:6708",
             "PTEN:5728",
             "DST:667",
             "DMD:1756",
             "NEB:4703")
  
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
  
  BAM_Normal_procesado_reducido_por_genes2 <- rbindlist(BAM_Normal_procesado_reducido_por_genes)
  
  ##### Volvemos a incorporar el archivo con los alineamientos de sam3pairwise y lo modificamos
  file <- read_delim(alineamiento, 
                     "\t", escape_double = FALSE, col_names = FALSE, 
                     trim_ws = TRUE)
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
  
  prueba2 <- left_join(matriz_modificada,BAM_Normal_procesado_reducido_por_genes2,by=c("name","rname","pos","cigar","algo3"))
  prueba2 <- na.omit(prueba2)
  prueba2$seq.y <- prueba2$ascii
  prueba2$ascii <- NULL
  
  #### Modificar las ASCII para que encajen con los alineamientos
  prueba2$longitud_fragmento2 <- nchar(prueba2$seq.x)
  
  library(stringr)
  for (n in 1:nrow(prueba2)){
    print(n)
    if(prueba2$longitud_fragmento[n]>150){
      secuencia <- str_split(prueba2$seq.x[n],pattern = "")[[1]]
      longitud <- length(secuencia)
      llenas <- which(secuencia!=".")
      nuevo_ascii <- 1:longitud
      nuevo_ascii[llenas] <- str_split(prueba2$seq.y[n],pattern = "")[[1]]
      nuevo_ascii[-llenas] <- "!"
      prueba2$seq.y[n] <- paste(nuevo_ascii,collapse = "")
    }
  }
  
  ##### Recortar las secuencias para que encajen perfectamente con las coordenadas 
  ##### genomicas. 
  library(stringr)
  for(n in 1:nrow(prueba2)){
    print(n)
    a <- str_split(prueba2$WT[n],pattern = "")
    eliminar <- 0
    eliminar2 <- length(a[[1]])
    for (z in 1:length(a[[1]])){
      if (a[[1]][z]=="N"){
        eliminar <- eliminar+1
      }
      else (break)
    }
    for (y in length(a[[1]]):1){
      if (a[[1]][y]=="N"){
        eliminar2 <- eliminar2-1
      }
      else (break)
    }
    if(eliminar>0 | eliminar2<length(a[[1]])){
      prueba2$seq.x[n] <- substr(prueba2$seq.x[n],eliminar+1,eliminar2)
      prueba2$seq.y[n] <- substr(prueba2$seq.y[n],eliminar+1,eliminar2)
      prueba2$WT[n] <- substr(prueba2$WT[n],eliminar+1,eliminar2)}
  }
  
  prueba2$end <- prueba2$pos + nchar(prueba2$seq.x)-1
  prueba2$longitud_fragmento <- nchar(prueba2$seq.x)
  colnames(prueba2)[10] <- "seq"
  
  pruebesita <- str_split(prueba2$WT,pattern = "")
  for (n in 1:length(pruebesita)){
    recortables <- which(pruebesita[[n]]=="-")
    numero_recortables <- length(recortables)
    if (numero_recortables>=1){
      para_recortar <- str_split(prueba2$seq[n],pattern = "")[[1]]
      para_recortar2 <- str_split(prueba2$WT[n],pattern = "")[[1]]
      para_recortar3 <- str_split(prueba2$seq.y[n],pattern = "")[[1]]
      para_recortar <- para_recortar[-recortables]
      para_recortar2 <- para_recortar2[-recortables]
      para_recortar3 <- para_recortar3[-recortables]
      prueba2$seq[n] <- paste(para_recortar,collapse = "")
      prueba2$WT[n] <- paste(para_recortar2,collapse = "")
      prueba2$seq.y[n] <- paste(para_recortar3,collapse = "")
    }
  }
  
  prueba2$end <- prueba2$pos + nchar(prueba2$seq)-1
  prueba2$longitud_fragmento <- nchar(prueba2$seq)
  
  BAM_Normal_procesado_reducido_por_genes <- list()
  for (n in 1:length(genes)){
    BAM_Normal_procesado_reducido_por_genes[[n]] <- prueba2[which(prueba2$...10==genes[n]),]
  }

BAM_Normal_procesado_reducido_por_genes <- lapply(BAM_Normal_procesado_reducido_por_genes,as.data.table)

### Reducimos el tamaño del BED separandolo en DF distintos para cada gen.
if(Metodo=="TCGA"){
    
    #### V1 de la funcion. Todas las mutaciones de cada uno de los genes
    reduccir_tamaño_bed <- function(bed,genes){
      desired_length <- length(genes) # or whatever length you want
      output <- vector(mode = "list", length = desired_length)
      for (n in 1:length(genes)){
        posicion_gen <- bed %>% filter(proteina == genes[n] | p2 == genes[n] | p3 == genes[n])
        posicion_gen$proteina <- genes[n]
        output[[n]] <- posicion_gen
      }
      names(output) <- genes
      return(output)
    }
    
    Bed_reducido <- reduccir_tamaño_bed(bedfile,genes)
    Bed_reducido <- lapply(Bed_reducido,as.data.table)
}else{
  
  #### V1 de la funcion. Todas las mutaciones de cada uno de los genes 
  reduccir_tamaño_bed <- function(bed,genes){
    desired_length <- length(genes) # or whatever length you want
    output <- vector(mode = "list", length = desired_length)
    for (n in 1:length(genes)){
      posicion_gen <- bed %>% filter(proteina == genes[n] | p2 == genes[n])
      posicion_gen$proteina <- genes[n]
      output[[n]] <- posicion_gen
    }
    names(output) <- genes
    return(output)
  }
    
    Bed_reducido <- reduccir_tamaño_bed(bedfile,genes)
  }


############## V2 mucho mas rapida ################

funcion_ascii <- function(x,valor_ascii){
  return(any(x<valor_ascii))
}
funcion_comparacion_y_tabla <- function(bed,bam,valor_ascii){
  matriz_final <- matrix(ncol = 13,nrow = sum(unlist(lapply(bed,nrow))))
  colnames(matriz_final) <- c("Rs","WT","Mutada","Mutacion Desconocida",
                              "No_seq","Pos Mutacion Desconocida","Proteina",
                              "REF","ALT","Numero Secuencias","WT%","TUMOR%","Desc%")
  contador <- 0
  for (n in 1:length(bam)){
    print(n)
    bedfile <- bed[[n]]
    bamfile <- bam[[n]]
    for (z in 1:nrow(bedfile)){
      #print(z)
      contador <- contador +1 
      longitud_WT <- nchar(bedfile$tipo[z])
      longitud_Mutacion <- nchar(bedfile$cambio[z])
      marco <- max(longitud_WT,longitud_Mutacion)
      mas_grande <- which.max(c(longitud_WT,longitud_Mutacion))
      posicion_start_bed <- bedfile$start[z]
      mini_df <- bamfile[end>=(posicion_start_bed+marco-1) & pos<= posicion_start_bed]
      WT <- 0
      Mutada <- 0
      Nidea <- 0
      No_seq <- 0
      pos_nidea <- c()
      if (dim(mini_df)[1]>=1){
        wild_type <- with(mini_df,substr(seq,(posicion_start_bed-pos+1),(posicion_start_bed-pos+1+longitud_WT-1)))
        mutacion <- with(mini_df,substr(seq,(posicion_start_bed-pos+1),(posicion_start_bed-pos+1+longitud_Mutacion-1)))
        ascii_wt <- with(mini_df,phredScoreStringToInt(substr(seq.y,(posicion_start_bed-pos+1),(posicion_start_bed-pos+1+longitud_WT-1)),scoreType = "Phred33"))
        ascii_mutacion <- with(mini_df,phredScoreStringToInt(substr(seq.y,(posicion_start_bed-pos+1),(posicion_start_bed-pos+1+longitud_Mutacion-1)),scoreType = "Phred33"))
        if(longitud_WT == 1){ascii_wt <- t(ascii_wt)}; if(longitud_Mutacion == 1){ascii_mutacion <- t(ascii_mutacion)}
        ascii_combinado <- cbind(ascii_wt,ascii_mutacion)
        no_seq <- which(apply(ascii_combinado,1,funcion_ascii,valor_ascii)==TRUE)
        No_seq <- length(no_seq)
        if (No_seq>0) {wild_type <- wild_type[-no_seq]; mutacion <- mutacion[-no_seq]}
        Ambiguos <- which(wild_type==bedfile$tipo[z] & mutacion==bedfile$cambio[z])
        if(mas_grande==1){WT <- WT + length(Ambiguos)} else {Mutada <- Mutada + length(Ambiguos)}
        if (length(Ambiguos)>0){wild_type <- wild_type[-Ambiguos]; mutacion <- mutacion[-Ambiguos]}
        wt <- which(wild_type==bedfile$tipo[z])
        WT <- length(wt) + WT
        if(length(wt) > 0) {wild_type <- wild_type[-wt]; mutacion <- mutacion[-wt]}
        mutada <- which(mutacion==bedfile$cambio[z])
        Mutada <- Mutada + length(mutada)
        if (length(mutada) >0) {wild_type <- wild_type[-mutada]; mutacion <- mutacion[-mutada]}
        #enes <- which(grepl("N",wild_type)==TRUE)
        #if (length(enes)>0) {wild_type <- wild_type[-enes]}
        Nidea <- length(wild_type)
      }
      matriz_final[contador,2] <- WT
      matriz_final[contador,3] <- Mutada
      matriz_final[contador,4] <- Nidea
      matriz_final[contador,5] <- No_seq
      matriz_final[contador,6] <- toString(pos_nidea)
      matriz_final[contador,10] <- nrow(mini_df)
      matriz_final[contador,11] <- (WT/(nrow(mini_df)-No_seq))*100
      matriz_final[contador,12] <- (Mutada/(nrow(mini_df)-No_seq))*100
      matriz_final[contador,13] <- (Nidea/(nrow(mini_df)-No_seq))*100
    }
  }
  matriz_final <- as.data.frame(matriz_final)
  Bed_reducido2 <- bind_rows(bed)
  matriz_final$Rs <- Bed_reducido2$rs
  matriz_final$Proteina <- Bed_reducido2$proteina
  matriz_final$REF <- Bed_reducido2$tipo
  matriz_final$ALT <- Bed_reducido2$cambio
  matriz_final[,c(2,3,4,5,10:13)] <- apply(matriz_final[,c(2,3,4,5,10:13)],2,as.numeric)
  return(matriz_final)
}

mutaciones_BAM_normal <- funcion_comparacion_y_tabla(Bed_reducido,BAM_Normal_procesado_reducido_por_genes,30)
#mutaciones_BAM_tumor <- funcion_comparacion_y_tabla(bedfile_reducido,BAM_Tumor_procesado_reducido_por_genes,30)

#mutaciones_comunes <- merge(mutaciones_BAM_normal,mutaciones_BAM_negativo,mutaciones_BAM_positivo,by=c("Rs","Proteina","REF","ALT"))
#mutaciones_comunes[,c(5,6,7,10,11,12,13,14,15,18,19,20)] <- apply(mutaciones_comunes[,c(5,6,7,10,11,12,13,14,15,18,19,20)],2,as.numeric)

#df_list <- list(mutaciones_BAM_normal, mutaciones_BAM_negativo, mutaciones_BAM_positivo)
#mutaciones_comunes <- Reduce(function(x, y) merge(x, y, by=c("Rs","Proteina","REF","ALT")), df_list, accumulate=FALSE)

############################################################

###### Filtrar mutaciones en busca de micrometastasis ######
library(bedr)
resesitos <- read_delim("D:/TFM/Qiime2/Compartida/Canceres usados en el estudio/Cancer de mama con nodulo centinela/resesitos.txt", 
                        "\t", escape_double = FALSE, trim_ws = TRUE)
resesitos$func <- strsplit2matrix(resesitos$func,split = ",")[,1]

mutaciones_BAM_SRR1751228 <- as.data.table(mutaciones_BAM_SRR1751228)
mutaciones_BAM_SRR1751229 <- as.data.table(mutaciones_BAM_SRR1751229)
mutaciones_BAM_SRR1751230 <- as.data.table(mutaciones_BAM_SRR1751230)
mutaciones_BAM_SRR1751231 <- as.data.table(mutaciones_BAM_SRR1751231)
mutaciones_BAM_SRR1751232 <- as.data.table(mutaciones_BAM_SRR1751232)

mutaciones_BAM_SRR1751228$`Numero Secuencias` <- apply(mutaciones_BAM_SRR1751228[,c(2,3,4)],1,sum)
mutaciones_BAM_SRR1751229$`Numero Secuencias` <- apply(mutaciones_BAM_SRR1751229[,c(2,3,4)],1,sum)
mutaciones_BAM_SRR1751230$`Numero Secuencias` <- apply(mutaciones_BAM_SRR1751230[,c(2,3,4)],1,sum)
mutaciones_BAM_SRR1751231$`Numero Secuencias` <- apply(mutaciones_BAM_SRR1751231[,c(2,3,4)],1,sum)
mutaciones_BAM_SRR1751232$`Numero Secuencias` <- apply(mutaciones_BAM_SRR1751232[,c(2,3,4)],1,sum)

mutaciones_BAM_normal$`Numero Secuencias` <- apply(mutaciones_BAM_normal[,c(2,3,4)],1,sum)
mutaciones_BAM_SRR14562366$`Numero Secuencias` <- apply(mutaciones_BAM_SRR14562366[,c(2,3,4)],1,sum)
mutaciones_BAM_SRR14562367$`Numero Secuencias` <- apply(mutaciones_BAM_SRR14562367[,c(2,3,4)],1,sum)
mutaciones_BAM_SRR14562368$`Numero Secuencias` <- apply(mutaciones_BAM_SRR14562368[,c(2,3,4)],1,sum)
mutaciones_BAM_SRR14562369$`Numero Secuencias` <- apply(mutaciones_BAM_SRR14562369[,c(2,3,4)],1,sum)
mutaciones_BAM_SRR14562370$`Numero Secuencias` <- apply(mutaciones_BAM_SRR14562370[,c(2,3,4)],1,sum)
mutaciones_BAM_SRR14562371$`Numero Secuencias` <- apply(mutaciones_BAM_SRR14562371[,c(2,3,4)],1,sum)
mutaciones_BAM_SRR14562372$`Numero Secuencias` <- apply(mutaciones_BAM_SRR14562372[,c(2,3,4)],1,sum)
mutaciones_BAM_SRR14562373$`Numero Secuencias` <- apply(mutaciones_BAM_SRR14562373[,c(2,3,4)],1,sum)


min_max_norm <- function(x) {
  (x - min(x)) / (max(x) - min(x))
}

mutaciones_BAM_SRR1751228$`Numero Secuencias 2` <- min_max_norm(mutaciones_BAM_SRR1751228$`Numero Secuencias`)
mutaciones_BAM_SRR1751229$`Numero Secuencias 2` <- min_max_norm(mutaciones_BAM_SRR1751229$`Numero Secuencias`)
mutaciones_BAM_SRR1751230$`Numero Secuencias 2` <- min_max_norm(mutaciones_BAM_SRR1751230$`Numero Secuencias`)
mutaciones_BAM_SRR1751231$`Numero Secuencias 2` <- min_max_norm(mutaciones_BAM_SRR1751231$`Numero Secuencias`)
mutaciones_BAM_SRR1751232$`Numero Secuencias 2` <- min_max_norm(mutaciones_BAM_SRR1751232$`Numero Secuencias`)


calcular_hetero_homo <- function(x,columna,condicion,alternativa){
  total <- sum(as.numeric(x[2]),as.numeric(x[3]),as.numeric(x[4]))
  a <- binom.test(as.numeric(x[columna]),total,condicion,alternative = alternativa)
  return(a$p.value)
}

funcion_clasificacion_final <- function(x,min_seq,resesitos){
  tabla_porcentajes_finales <- data.table::data.table()
  tabla_porcentajes_finales$Totales <- nrow(x)
  tabla_porcentajes_finales$No_encontrados <- nrow(x[`TUMOR%`=="NaN"])
  tabla_porcentajes_finales$Poco_representados <- nrow(x[`Numero Secuencias`< min_seq & `TUMOR%`!="NaN"])
  
  x <- x[`Numero Secuencias`>= min_seq & `TUMOR%`!="NaN"]
  
  heterocigotos_1 <- apply(x,1,calcular_hetero_homo,2,0.5,"two.sided")
  heterocigotos_2 <- apply(x,1,calcular_hetero_homo,3,0.5,"two.sided")
  heterocigotos_3 <- apply(x,1,calcular_hetero_homo,4,0.5,"two.sided")
  a <- which(heterocigotos_1>=0.05)
  b <- which(heterocigotos_2>=0.05)
  c <- which(heterocigotos_3>=0.05)
  d <- unique(c(a,b,c))
  
  if(length(d)>0) {x <- x[-d,]}
  
  tabla_porcentajes_finales$Heterocigotos <- length(d)
  
  homocigotos_1 <- which(x$`WT%`==100)
  homocigotos_2 <- which(x$`TUMOR%`==100)
  homocigotos_3 <- which(x$`Desc%`==100)
  homocigotos_4 <- unique(c(homocigotos_1,homocigotos_2,homocigotos_3))
  
  tabla_porcentajes_finales$Homocigotos_WT <- length(homocigotos_1)
  tabla_porcentajes_finales$Homocigotos_Tumor <- length(homocigotos_2)
  tabla_porcentajes_finales$Homocigotos_Desc <- length(homocigotos_3)
  
  x <- x[-homocigotos_4,]
  
  x <- x[`TUMOR%` < 50 & `TUMOR%` > 0]
  
  descartables <- apply(x,1,calcular_hetero_homo,3,0.001,"greater")
  descartables2 <- which(descartables>=0.05)
  
  x <- x[-descartables2,]
  x <- merge(x,resesitos,by.x = "Rs",by.y = "name",all.x=TRUE)
  
  cambios_funcionales <- list()
  
  tabla_porcentajes_finales$`0-1` <- nrow(x[`TUMOR%` > 0 & `TUMOR%` <= 1]); cambios_funcionales[[1]] <- t(as.matrix(table(x[`TUMOR%` > 0 & `TUMOR%` <= 1]$func),row.names = TRUE))
  tabla_porcentajes_finales$`1-2` <- nrow(x[`TUMOR%` > 1 & `TUMOR%` <= 2]); cambios_funcionales[[2]] <- t(as.matrix(table(x[`TUMOR%` > 1 & `TUMOR%` <= 2]$func),row.names = TRUE))
  tabla_porcentajes_finales$`2-3` <- nrow(x[`TUMOR%` > 2 & `TUMOR%` <= 3]); cambios_funcionales[[3]] <- t(as.matrix(table(x[`TUMOR%` > 2 & `TUMOR%` <= 3]$func),row.names = TRUE))
  tabla_porcentajes_finales$`3-4` <- nrow(x[`TUMOR%` > 3 & `TUMOR%` <= 4]); cambios_funcionales[[4]] <- t(as.matrix(table(x[`TUMOR%` > 3 & `TUMOR%` <= 4]$func),row.names = TRUE))
  tabla_porcentajes_finales$`4-5` <- nrow(x[`TUMOR%` > 4 & `TUMOR%` <= 5]); cambios_funcionales[[5]] <- t(as.matrix(table(x[`TUMOR%` > 4 & `TUMOR%` <= 5]$func),row.names = TRUE))
  tabla_porcentajes_finales$`5-6` <- nrow(x[`TUMOR%` > 5 & `TUMOR%` <= 6]) ; cambios_funcionales[[6]] <- t(as.matrix(table(x[`TUMOR%` > 5 & `TUMOR%` <= 6]$func),row.names = TRUE))
  tabla_porcentajes_finales$`6-7` <- nrow(x[`TUMOR%` > 6 & `TUMOR%` <= 7]); cambios_funcionales[[7]] <- t(as.matrix(table(x[`TUMOR%` > 6 & `TUMOR%` <= 7]$func),row.names = TRUE))
  tabla_porcentajes_finales$`7-8` <- nrow(x[`TUMOR%` > 7 & `TUMOR%` <= 8]); cambios_funcionales[[8]] <- t(as.matrix(table(x[`TUMOR%` > 7 & `TUMOR%` <= 8]$func),row.names = TRUE))
  tabla_porcentajes_finales$`8-9` <- nrow(x[`TUMOR%` > 8 & `TUMOR%` <= 9]); cambios_funcionales[[9]] <- t(as.matrix(table(x[`TUMOR%` > 8 & `TUMOR%` <= 9]$func),row.names = TRUE))
  tabla_porcentajes_finales$`9-10` <- nrow(x[`TUMOR%` > 9 & `TUMOR%` <= 10]); cambios_funcionales[[10]] <- t(as.matrix(table(x[`TUMOR%` > 9 & `TUMOR%` <= 10]$func),row.names = TRUE))
  tabla_porcentajes_finales$`10-11` <- nrow(x[`TUMOR%` > 10 & `TUMOR%` <= 11]); cambios_funcionales[[11]] <- t(as.matrix(table(x[`TUMOR%` > 10 & `TUMOR%` <= 11]$func),row.names = TRUE))
  tabla_porcentajes_finales$`11-12` <- nrow(x[`TUMOR%` > 11 & `TUMOR%` <= 12]); cambios_funcionales[[12]] <- t(as.matrix(table(x[`TUMOR%` > 11 & `TUMOR%` <= 12]$func),row.names = TRUE))
  tabla_porcentajes_finales$`12-13` <- nrow(x[`TUMOR%` > 12 & `TUMOR%` <= 13]); cambios_funcionales[[13]] <- t(as.matrix(table(x[`TUMOR%` > 12 & `TUMOR%` <= 13]$func),row.names = TRUE))
  tabla_porcentajes_finales$`13-14` <- nrow(x[`TUMOR%` > 13 & `TUMOR%` <= 14]); cambios_funcionales[[14]] <- t(as.matrix(table(x[`TUMOR%` > 13 & `TUMOR%` <= 14]$func),row.names = TRUE))
  tabla_porcentajes_finales$`14-15` <- nrow(x[`TUMOR%` > 14 & `TUMOR%` <= 15]); cambios_funcionales[[15]] <- t(as.matrix(table(x[`TUMOR%` > 14 & `TUMOR%` <= 15]$func),row.names = TRUE))
  tabla_porcentajes_finales$`15-50` <- nrow(x[`TUMOR%` > 15 & `TUMOR%` <= 50]); cambios_funcionales[[16]] <- t(as.matrix(table(x[`TUMOR%` > 15 & `TUMOR%` <= 50]$func),row.names = TRUE))
  
  tabla_porcentajes_finales <- as.data.frame(t(tabla_porcentajes_finales))
  #rownames(cambios_funcionales) <- rownames(tabla_porcentajes_finales)[8:23]
  #cambios_funcionales <- as.data.frame(cambios_funcionales)
  cambios_funcionales <- lapply(cambios_funcionales,as.data.frame)
  cambitos <- lapply(cambios_funcionales,length); cambitos <- which(cambitos == 0)
  if (length(cambitos)>0){
    for (n in 1:length(cambitos)){
    cambios_funcionales[[cambitos[n]]] <- data.frame(ncRNA = 0)}
    }
  cambios_funcionales <- as.data.table(do.call(plyr::rbind.fill, cambios_funcionales))
  rownames(cambios_funcionales) <- rownames(tabla_porcentajes_finales)[8:23]
  cambios_funcionales[is.na(cambios_funcionales)] <- 0
  
  irrelevantes <- c("unknown","coding-synon","intron","near-gene-3","near-gene-5")
  relevantes <- c("nonsense","missense","stop-loss","frameshift","cds-indel","untranslated-3","untranslated-5","splice-3","splice-5")
  ncRNA <- c("ncRNA")
  
  irrelevantes <- irrelevantes[which(irrelevantes %in% colnames(cambios_funcionales) == TRUE)]
  relevantes <- relevantes[which(relevantes %in% colnames(cambios_funcionales) == TRUE)] 
  
  irrelevantes <- cambios_funcionales[,..irrelevantes]
  relevantes <- cambios_funcionales[,..relevantes]
  
  irrelevantes <- apply(irrelevantes,1,sum)
  relevantes <- apply(relevantes,1,sum)
  
  cambios_funcionales2 <- as.data.table(cbind(relevantes,irrelevantes,cambios_funcionales$ncRNA))
  return(list(tabla_porcentajes_finales,x,table(x$Proteina),cambios_funcionales,cambios_funcionales2))
}

porcentajes_finales_SRR1751228 <- funcion_clasificacion_final(mutaciones_BAM_SRR1751228,100,resesitos)
porcentajes_finales_SRR1751229 <- funcion_clasificacion_final(mutaciones_BAM_SRR1751229,100,resesitos)
porcentajes_finales_SRR1751230 <- funcion_clasificacion_final(mutaciones_BAM_SRR1751230,100,resesitos)
porcentajes_finales_SRR1751231 <- funcion_clasificacion_final(mutaciones_BAM_SRR1751231,100,resesitos)
porcentajes_finales_SRR1751232 <- funcion_clasificacion_final(mutaciones_BAM_SRR1751232,100,resesitos)

porcentajes_finales_normal <- funcion_clasificacion_final(mutaciones_BAM_normal,15,resesitos)
porcentajes_finales_SRR14562366 <- funcion_clasificacion_final(mutaciones_BAM_SRR14562366,15,resesitos)
porcentajes_finales_SRR14562367 <- funcion_clasificacion_final(mutaciones_BAM_SRR14562367,15,resesitos)
porcentajes_finales_SRR14562368 <- funcion_clasificacion_final(mutaciones_BAM_SRR14562368,15,resesitos)
porcentajes_finales_SRR14562369 <- funcion_clasificacion_final(mutaciones_BAM_SRR14562369,15,resesitos)
porcentajes_finales_SRR14562370 <- funcion_clasificacion_final(mutaciones_BAM_SRR14562370,15,resesitos)
porcentajes_finales_SRR14562371 <- funcion_clasificacion_final(mutaciones_BAM_SRR14562371,15,resesitos)
porcentajes_finales_SRR14562372 <- funcion_clasificacion_final(mutaciones_BAM_SRR14562372,15,resesitos)
porcentajes_finales_SRR14562373 <- funcion_clasificacion_final(mutaciones_BAM_SRR14562373,15,resesitos)

esto <- 5

cosa1 <- cbind(porcentajes_finales_SRR1751228[[esto]],porcentajes_finales_SRR1751229[[esto]],porcentajes_finales_SRR1751230[[esto]],porcentajes_finales_SRR1751231[[esto]],porcentajes_finales_SRR1751232[[esto]])
cosa1 <- cbind(porcentajes_finales_normal[[esto]],porcentajes_finales_SRR14562366[[esto]],porcentajes_finales_SRR14562371[[esto]])
cosa3 <- cbind(porcentajes_finales_normal[[esto]],porcentajes_finales_SRR14562366[[esto]],
               porcentajes_finales_SRR14562367[[esto]],porcentajes_finales_SRR14562368[[esto]],
               porcentajes_finales_SRR14562369[[esto]],porcentajes_finales_SRR14562370[[esto]],
               porcentajes_finales_SRR14562371[[esto]],porcentajes_finales_SRR14562372[[esto]],
               porcentajes_finales_SRR14562373[[esto]])


cosa2 <- c(porcentajes_finales_SRR1751228[[2]]$Rs,porcentajes_finales_SRR1751229[[2]]$Rs,porcentajes_finales_SRR1751230[[2]]$Rs,porcentajes_finales_SRR1751231[[2]]$Rs,porcentajes_finales_SRR1751232[[2]]$Rs)
cosa2 <- c(porcentajes_finales_normal[[2]]$Rs,porcentajes_finales_SRR14562366[[2]]$Rs,
           porcentajes_finales_SRR14562367[[2]]$Rs,porcentajes_finales_SRR14562368[[2]]$Rs,
           porcentajes_finales_SRR14562369[[2]]$Rs,porcentajes_finales_SRR14562370[[2]]$Rs,
           porcentajes_finales_SRR14562371[[2]]$Rs,porcentajes_finales_SRR14562372[[2]]$Rs,
           porcentajes_finales_SRR14562373[[2]]$Rs)


#### Sacar los rs
cosa2 <- unique(cosa2)
write.table(cosa2,file = "reses.txt",row.names = FALSE,quote = FALSE,col.names = FALSE)
############################################################

library(rsnps)

ereses <- unique(porcentajes_finales_SRR1751230[[2]]$Rs)
write.table(cosa2,file = "ereses.txt",row.names = FALSE,quote = FALSE,col.names = FALSE)


cosita <- merge(porcentajes_finales_SRR1751228[[2]],resesitos,by.x = "Rs",by.y = "name")

############################################################

creacion_plots <- function(x){
  genes <- table(x[[2]]$Proteina)
  genes <- genes[order(genes,decreasing = T)]
  barplot(genes,cex.names=0.5)
  
}

creacion_plots(porcentajes_finales_SRR1751230)
