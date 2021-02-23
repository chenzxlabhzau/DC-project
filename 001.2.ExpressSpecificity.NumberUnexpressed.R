{ #pdf version
  rm(list = ls())
  wd = "~/MamDC/Data/Seqdata/CardosoKaessmann2019NatureRNAseq/"
  library(dplyr)
  library(tidyr)
  library(ggpubr)
  for( Species in c("Human","Mouse")) {
    a <- fread(paste0(wd,Species,"/",Species,".allgene.fpkm.txt")) %>% as.data.frame()
    colnames(a) %<>% gsub("KidneyTestis_4w_Male","Kidney_4w_Male",.) %>%
      gsub("Forebrain","Brain",.) %>% gsub("Hindbrain","Cerebellum",.) %>%
      gsub("Ovary","Gonad",.) %>% gsub("Testis","Gonad",.)
    b <- read.csv(file.path("~/MamDC/Data/Ref",Species,paste0("gene",Species,".bed")),header = F,sep = "\t") 
    b$V5 %<>% gsub(" ","",.)
    a[,c("chr","type")] <- b[match(a$V1,b$V4),c(1,5)]
    a %<>% filter(chr %in% c(1:100,"X") & type %in% "protein_coding")
    a[a$chr!="X","chr"]="A"
    
    #cutoff = 0, 1
    for (cutoff in c(0,1)) {
      unexpress = t(apply(a[,-c(1,ncol(a)-1,ncol(a))], 2, function(x)tapply(x,a$chr,function(x)sum(x <=cutoff)))) %>% as.data.frame()
      unexpress$A = unexpress$A/sum(a$chr=="A") *100 # %
      unexpress$X = unexpress$X/sum(a$chr=="X") *100
      unexpress$tissue = lapply(row.names(unexpress),function(x)(strsplit(x,"_")[[1]][1])) %>% unlist()
      unexpress$sex = lapply(row.names(unexpress),function(x)(strsplit(x,"_")[[1]][3])) %>% unlist()
      mean(unexpress$X/unexpress$A)
      tapply(unexpress$X, unexpress$tissue,mean)
      pdf(paste0("/home/qians/MamDC/Result/CardosoKaessmann2019NatureRNAseq/",Species,"/Picture/",Species,
                 ".Frequency.unexpressed.",cutoff,"FPKM.tissue.pdf"),
          width = 12, height = 7)
      print(ggpaired(unexpress,cond1="A",cond2="X",
                     color = "condition",palette =c("#035782","#9c0b1f"),line.size = 0.05,
                     facet.by ="tissue")+
              stat_compare_means(paired = TRUE)+
              xlab("Chromosome")+
              ylab("Frequency of unexpressed genes"))
      dev.off()
      
      pdf(paste0("/home/qians/MamDC/Result/CardosoKaessmann2019NatureRNAseq/",Species,"/Picture/",Species,
                 ".Frequency.unexpressed.",cutoff,"FPKM.pdf"),
          width = 5, height = 4)
      print(ggpaired(unexpress,cond1="A",cond2="X",
                     color = "condition",palette =c("#035782","#9c0b1f"),line.size = 0.03)+
              stat_compare_means(paired = TRUE)+
              xlab("Chromosome")+
              ylab("Frequency of unexpressed genes"))
      dev.off()
    }
  }
}

{ #pptx version
  rm(list = ls())
  wd = "~/MamDC/Data/Seqdata/CardosoKaessmann2019NatureRNAseq/"
  library(dplyr)
  library(tidyr)
  library(ggpubr)
  for( Species in c("Human","Mouse")) {
    a <- fread(paste0(wd,Species,"/",Species,".allgene.fpkm.txt")) %>% as.data.frame()
    colnames(a) %<>% gsub("KidneyTestis_4w_Male","Kidney_4w_Male",.) %>%
      gsub("Forebrain","Brain",.) %>% gsub("Hindbrain","Cerebellum",.) #%>%
     # gsub("Ovary","Gonad",.) %>% gsub("Testis","Gonad",.)
    b <- read.csv(file.path("~/MamDC/Data/Ref",Species,paste0("gene",Species,".bed")),header = F,sep = "\t") 
    b$V5 %<>% gsub(" ","",.)
    a[,c("chr","type")] <- b[match(a$V1,b$V4),c(1,5)]
    a %<>% filter(chr %in% c(1:100,"X") & type %in% "protein_coding")
    a[a$chr!="X","chr"]="A"
    
    #cutoff = 0, 1
    for (cutoff in c(0,1)) {
      unexpress = t(apply(a[,-c(1,ncol(a)-1,ncol(a))], 2, function(x)tapply(x,a$chr,function(x)sum(x <=cutoff)))) %>% as.data.frame()
      unexpress$A = unexpress$A/sum(a$chr=="A") *100 # %
      unexpress$X = unexpress$X/sum(a$chr=="X") *100
      unexpress$tissue = lapply(row.names(unexpress),function(x)(strsplit(x,"_")[[1]][1])) %>% unlist()
      unexpress$sex = lapply(row.names(unexpress),function(x)(strsplit(x,"_")[[1]][3])) %>% unlist()
      p = ggpaired(unexpress,cond1="A",cond2="X",
                     color = "condition",palette =c("#035782","#9c0b1f"),line.size = 0.05,
                     facet.by ="tissue")+
              stat_compare_means(paired = TRUE)+
              xlab("Chromosome")+
              ylab("Frequency of unexpressed genes")
      topptx(p,paste0("/home/qians/MamDC/Result/CardosoKaessmann2019NatureRNAseq/",Species,"/Picture/",Species,
             ".Frequency.unexpressed.",cutoff,"FPKM.tissue.pptx"),width = 12, height = 7)
      
      
      pdf(paste0("/home/qians/MamDC/Result/CardosoKaessmann2019NatureRNAseq/",Species,"/Picture/",Species,
                 ".Frequency.unexpressed.",cutoff,"FPKM.pdf"),
          width = 5, height = 4)
      p = ggpaired(unexpress,cond1="A",cond2="X",
                     color = "condition",palette =c("#035782","#9c0b1f"),line.size = 0.01)+
              stat_compare_means(paired = TRUE)+
              xlab("Chromosome")+
              ylab("Frequency of unexpressed genes")
      topptx(p,paste0("/home/qians/MamDC/Result/CardosoKaessmann2019NatureRNAseq/",Species,"/Picture/",Species,
                      ".Frequency.unexpressed.",cutoff,"FPKM.pptx"),
             width = 5, height = 4)
    }
  }
}
