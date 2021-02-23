if (TRUE) { # Expression distribution, an example
  rm(list = ls());gc();rm(list = ls())
  library(dplyr)
  library(tidyr)
  library(ggplotify)
  library(eoffice)
  wd = "~/MamDC/Data/Seqdata/CardosoKaessmann2019NatureRNAseq/"
  for( Species in c("Human","Mouse")) {
    a <- fread(paste0(wd,Species,"/",Species,".allgene.fpkm.txt")) %>% as.data.frame()
    colnames(a) %<>% gsub("KidneyTestis_4w_Male","Kidney_4w_Male",.) %>%
      gsub("Forebrain","Brain",.) %>% gsub("Hindbrain","Cerebellum",.)
    b <- read.csv(file.path("~/MamDC/Data/Ref",Species,paste0("gene",Species,".bed")),header = F,sep = "\t") 
    b$V5 %<>% gsub(" ","",.)
    a[,c("chr","type")] <- b[match(a$V1,b$V4),c(1,5)]
    a %<>% filter(chr %in% c(1:100,"X") & type %in% "protein_coding")
    #a[a$chr!="X","chr"]="A"
    if (Species=="Human") {
      age = read.csv("~/MamDC/Data/Seqdata/Genorigin/Human/Homo_sapiens.csv",header = TRUE,sep = ",",stringsAsFactors = FALSE)
      age$gene_age %<>% gsub(">","",.) %>% as.numeric() 
      a$age = age[match(a$V1,age$ensembl_id),2]
      tissue = "Brain_4w_Male"
    }
    if (Species=="Mouse") {
      age = read.csv("~/MamDC/Data/Seqdata/Genorigin/Mouse/Mus_musculus.csv",header = TRUE,sep = ",",stringsAsFactors = FALSE)
      age$gene_age %<>% gsub(">","",.) %>% as.numeric() 
      a$age = age[match(a$V1,age$ensembl_id),2]
      tissue = "Brain_s11_Male"
    }
    
    focal.tissue = a %<>% dplyr::select(.,c(chr,age,starts_with(tissue)))
    #focal.tissue$age = 0-focal.tissue$age
    focal.tissue$express = apply(focal.tissue[,-(1:2)], 1, mean)
    focal.tissue[focal.tissue$chr!="X","chr"]="A"
    focal.tissue %<>% dplyr::filter(.,express >1)
    

    #focal.tissue$age %<>% factor()
    ymax = 0.23
    if (Species=="Mouse") {
      ymax = 0.20
    }
    p = ggplot(focal.tissue,aes(log2(express),group=chr,color=chr))+geom_density(adjust=0.5)+
      theme_classic()+xlab("Expression (log2)")+ylab(strsplit(colnames(focal.tissue)[ncol(focal.tissue)-1],split = "_rep")[[1]][1])+
      scale_color_manual(values = c("#035782","#9c0b1f"))+
      theme(axis.text.x = element_text(size=14),axis.title.x = element_text(size=16),
            axis.text.y = element_text(size=14),axis.title.y = element_text(size=16))+
      guides(color=FALSE)+
      annotate("text",x=7.5,y=ymax,
               label=paste0("Mean A =",round(log2(mean(focal.tissue[focal.tissue$chr=="A","express"])),2),",","\n",
                            paste0("Mean X =",round(log2(mean(focal.tissue[focal.tissue$chr=="X","express"])),2),",","\n",
                                   paste0("Median A =",round(log2(median(focal.tissue[focal.tissue$chr=="A","express"])),2),",","\n",
                                          paste0("Median X =",round(log2(median(focal.tissue[focal.tissue$chr=="X","express"])),2),
                                          size=5)))))
    ggsave(p,filename = file.path("/home/qians/MamDC/Result/Genorigin",Species,"Picture/Expression.distribution.AX.pdf"),
           width = 5,height = 4)
    topptx(p,filename = file.path("/home/qians/MamDC/Result/Genorigin",Species,"Picture/Expression.distribution.AX.pptx"),
           width = 5,height = 4)
    ggplot(focal.tissue,aes(log2(express),fill=chr))+geom_histogram(aes(y = ..density..),position='dodge',bins = 50)
  }
}