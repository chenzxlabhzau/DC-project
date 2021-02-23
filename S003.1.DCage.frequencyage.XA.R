if (TRUE) { #Frequency ~ age (X A)
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
    wil = wilcox.test(focal.tissue[focal.tissue$chr!="X","age"],focal.tissue[focal.tissue$chr=="X","age"])
    wil$p.value
    
    p = ggplot(focal.tissue,aes(age,group=chr,color=chr))+geom_density(adjust = 3.5)+theme_bw()+
      xlab("Million year (log10)")+ylab("Frequency")+
      theme(panel.grid.major =element_blank(), 
            panel.grid.minor = element_blank(),
            axis.text.x = element_text(size=14),axis.title.x = element_text(size=16),
            axis.text.y = element_blank(),axis.title.y = element_text(size=16),
            axis.ticks.y = element_blank()#,axis.ticks.x = element_blank()
      )+
      scale_color_manual(values = c(rep("#8db8d5",length(unique(focal.tissue$chr))-1),"#9c0b1f"))+
      geom_vline(xintercept = median(focal.tissue[focal.tissue$chr!="X","age"]),color="#035782",linetype="dashed")+
      geom_vline(xintercept = median(focal.tissue[focal.tissue$chr=="X","age"]),color="#9c0b1f",linetype="dashed")+
      guides(color=FALSE)+
      annotate("text",x=1300,y=0.00002,label=paste0("P=",format(wil$p.value,5)),size=5)
    ggsave(p,filename = file.path("/home/qians/MamDC/Result/Genorigin",Species,"Picture/Frequency.age.AX.pdf"),
           width = 5,height = 4)
    topptx(p,filename = file.path("/home/qians/MamDC/Result/Genorigin",Species,"Picture/Frequency.age.AX.pptx"),
           width = 5,height = 4)
  }
  
  
}