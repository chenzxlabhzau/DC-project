---
title: "Untitled"
author: "qians"
date: "1/11/2021"
output: html_document
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

```{bash}
cd ~/Pseudo/Data/Ref/Human
grep -w gene Homo_sapiens.GRCh38.98.gtf |awk -F ";" '{print $1"\t"$3"\t"$5}' > gene.ENSG.name.txt
```

```{r}
rm(list = ls())
S="Human"
wd = "/home/qians/Pseudo/Data/Seqdata/CardosoKaessmann2019NatureRNAseq/"
a <- fread(file.path(wd,S,"counts.tsv")) %>% as.data.frame()
a[1:3,1:3]
a <- dplyr::filter(a,!grepl("_PAR_Y",V1))
a$V1 <- lapply(a$V1,function(x)strsplit(x,".",fixed=T)[[1]][1]) %>% unlist()
b <- read.csv(file.path("/home/qians/Pseudo/Data/Ref",S,paste0("gene",S,".bed")),
              header = FALSE,sep = "\t",stringsAsFactors = F)
a$type <- b[match(a$V1,b$V4),5]
a <- dplyr::filter(a,grepl("protein_coding",type,ignore.case = T))
row.names(a) <- a$V1
a <- a[,-c(1,ncol(a))]
a <- apply(a, 2, function(x){10^6 * x/sum(x)})
  
coldata <- read.csv(file = file.path("~/Pseudo/Result",S,"Savedata","coldata.csv"),
                    row.names = 1, sep = ",")
coldata$time <- lapply(coldata$stage2,function(x){substr(x,1,str_length(x)-3)}) %>% unlist() %>% as.numeric()
for (i in 1:nrow(coldata)) {
  if (coldata[i,4]=="week") {
    coldata[i,9] = coldata[i,8] * 7
  }
  if (coldata[i,4]=="day") {
    coldata[i,9] = coldata[i,8] * 1 + 280
  }
  if (coldata[i,4]=="month") {
    coldata[i,9] = coldata[i,8] * 30 + 280
  }
  if (coldata[i,4]=="year") {
    coldata[i,9] = coldata[i,8] * 365 + 280
  }
}
coldata$V9 <- log2(coldata$V9)
coldata <- coldata[order(coldata$V9),]

#for (T in unique(coldata$condition)) {
T = "Brain"
focal.tissue = dplyr::filter(coldata,condition== T)
focal.express = a[,focal.tissue$name] %>% as.data.frame()
focal.express = focal.express[apply(focal.express, 1, max)>0,]

cor.r = t(do.call(rbind, lapply(1:nrow(focal.express), function(i){
  cor.test(as.numeric(focal.express[i,]), focal.tissue$V9)$estimate})))
cor.p = t(do.call(rbind, lapply(1:nrow(focal.express), function(i){
  cor.test(as.numeric(focal.express[i,]), focal.tissue$V9)$p.value})))
focal.express$R = as.numeric(cor.r)
focal.express$P = as.numeric(cor.p)
focal.express$chr = b[match(row.names(focal.express),b$V4),1]
focal.express %<>% dplyr::filter(.,chr %in% c(1:100,"X"))
focal.express[focal.express$chr!="X","chr"]="A"
table(focal.express[focal.express$R>=0.8&focal.express$P<0.05,"chr"])
id = fread("~/Pseudo/Data/Ref/Human/gene.ENSG.name.txt") %>% as.data.frame()
id$ENSG = lapply(id$V9,function(x)strsplit(x,split = "\"",fixed = TRUE)[[1]][2]) %>% unlist()
id$Name = lapply(id$V10,function(x)strsplit(x,split = "\"",fixed = TRUE)[[1]][2]) %>% unlist()
ad = read.csv("~/MamDC/Data/Seqdata/HuWang2017AlzheimerRTgenelist/ADassociated.gene.csv",header = TRUE,sep = ",",stringsAsFactors = FALSE)
ad$ENSG = id[match(ad$Genesymbol,id$Name),"ENSG"]
ad[ad$Genesymbol=="DOPEY2","ENSG"] = "ENSG00000142197"
ad[ad$Genesymbol=="GSTT1","ENSG"] = "ENSG00000277656"
ad[ad$Genesymbol=="KIAA1033","ENSG"] = "ENSG00000136051"
#ad[ad$Genesymbol=="LOC642487","ENSG"] = 
ad[ad$Genesymbol=="PVRL2","ENSG"] = "ENSG00000130202"
ad[ad$Genesymbol=="SEPT3","ENSG"] = "ENSG00000100167"
ad = rbind(ad,c("KDM6A","KDM6A","ENSG00000147050"))
ad[,c("chr","type")] = b[match(ad$ENSG,b$V4),c(1,5)]
ad %<>% na.omit() %>% dplyr::filter(.,chr %in% c(1:100,"X") & type == " protein_coding")
ad[,c("R","P")] = focal.express[match(ad$ENSG,row.names(focal.express)),c("R","P")]


fisher.test(matrix(c(87,5172,404,18917),nrow = 2))

cutoff = 0.5
stagecor =  focal.express[abs(focal.express$R) >cutoff & focal.express$P <0.05, ]
c(sum(ad$ENSG %in% row.names(stagecor[stagecor$R>0,])), # overlap positive gene with ad gene
sum(stagecor$R>0), # number of stage positive gene
nrow(ad), # number of ad associated gene
nrow(focal.express)) %>% matrix(.,nrow = 2) %>% fisher.test(.,alternative = "greater") # number of gene tested

c(sum(ad$ENSG %in% row.names(stagecor[stagecor$R<0,])), # overlap negative gene with ad gene
sum(stagecor$R<0), # number of stage negative gene
nrow(ad), # number of ad associated gene
nrow(focal.express))  %>% matrix(.,nrow = 2) %>% fisher.test(.,alternative = "less") # number of gene tested

```

