#Distribution of the organ in which maximum expression is observed for X-linked genes and autosome ones
## 1. CardosoKaessmann2019NatureRNAseq
{
  rm(list = ls())
  wd = "~/MamDC/Data/Seqdata/CardosoKaessmann2019NatureRNAseq/"
  library(dplyr)
  library(tidyr)
  library(ggplotify)
  library(eoffice)
  dis.tissue.all = c()
  for( Species in c("Human","Mouse")) {
    a <- fread(paste0(wd,Species,"/",Species,".allgene.fpkm.txt")) %>% as.data.frame()
    colnames(a) %<>% gsub("KidneyTestis_4w_Male","Kidney_4w_Male",.) %>%
      gsub("Forebrain","Brain",.) %>% gsub("Hindbrain","Cerebellum",.)
    b <- read.csv(file.path("~/MamDC/Data/Ref",Species,paste0("gene",Species,".bed")),header = F,sep = "\t") 
    b$V5 %<>% gsub(" ","",.)
    a[,c("chr","type")] <- b[match(a$V1,b$V4),c(1,5)]
    a %<>% filter(chr %in% c(1:100,"X") & type %in% "protein_coding")
    a[a$chr!="X","chr"]="A"
    
    a = a[,-c(1,ncol(a))]
    a[1:3,1:3]
    a$sample = colnames(a)[apply(a[,-ncol(a)], 1, which.max)]
    a$tissue = lapply(a$sample,function(x)(strsplit(x,"_")[[1]][1])) %>% unlist()
    #a$maxexpress = apply(a[,-c(ncol(a),ncol(a)-1,ncol(a)-2)], 1, max)
    dis.tissue = as.data.frame(t(table(a$chr,a$tissue)/as.numeric(table(a$chr))))
    dis.tissue$species = Species
    dis.tissue.all = rbind(dis.tissue.all,dis.tissue)
  }
  p = ggplot(dis.tissue.all,aes(Var2,Freq,fill=Var1))+geom_bar(stat = "identity",position = "fill")+
    theme_classic()+facet_grid(.~species)+ylab("Tissue distribution")+
    theme(axis.text.x = element_text(size=14),
          axis.text.y = element_text(size=14),
          axis.title.x = element_blank(),
          axis.title.y = element_text(size=16),
          strip.text.x = element_text(size = 12))+
    guides(fill=guide_legend(title=NULL))+
    scale_fill_manual(values = c("#3298c8","#33cafe","#c60202","#cb9900","#349a00","#cc3397","#ff6700"))
  ggsave(p, filename = "/home/qians/MamDC/Result/CardosoKaessmann2019NatureRNAseq/Human/Picture/001.HumanMouse.Organdistribution.pdf",
         device = "pdf",width = 4.2, height = 4)
  topptx(p,"/home/qians/MamDC/Result/CardosoKaessmann2019NatureRNAseq/Human/Picture/001.HumanMouse.Organdistribution.pptx",
         width = 4.2, height = 4)
}
    
## 2.1 WangKuster2019MSBRNAseq RNA-seq
{
  rm(list = ls())
  wd = "/home/qians/MamDC/Data/Seqdata/WangKuster2019MSBRNAseq"
  Species ="Human"
  a <- read.csv(file.path(wd,paste0(Species,".allgene.fpkm.txt")),header = T,sep = "\t",stringsAsFactors = FALSE)
  b <- read.csv(file.path("~/MamDC/Data/Ref",Species,paste0("gene",Species,".bed")),
                header = FALSE,sep = "\t",stringsAsFactors = F)
  b$V5 %<>% gsub(" ","",.)
  a[,c("chr","type")] <- b[match(row.names(a),b$V4),c(1,5)]
  a %<>% dplyr::filter(.,type == "protein_coding" & chr %in% c(1:22,"X")) %>% dplyr::select(.,-type)
  a[a$chr!="X","chr"]="A"
  a$sample = colnames(a)[apply(a[,-ncol(a)], 1, which.max)]
  a$tissue = lapply(a$sample,function(x)(strsplit(x,"_rep")[[1]][1])) %>% unlist()
  dis.tissue = as.data.frame(t(table(a$chr,a$tissue)/as.numeric(table(a$chr))))
  dis.tissue$Var1 %<>% as.character() %>% Hmisc::capitalize(.) %>% gsub("_"," ",.)
  dis.tissue$color = rep(c("grey","#DDD7C6"),nrow(dis.tissue)/2)
  dis.tissue[dis.tissue$Var1=="Cerebral cortex","color"]= "#3298c8"
  dis.tissue[dis.tissue$Var1=="Testis","color"]= "#ff6700"
  dis.tissue[dis.tissue$Var1=="Ovary","color"]= "#cc3397"
  dis.tissue[dis.tissue$Var1=="Esophagus","color"]= "black"
  p = ggplot(dis.tissue,aes(Var2,Freq,fill=Var1))+geom_bar(stat = "identity",position = "fill")+
    theme_classic()+ylab("Tissue distribution")+
    theme(axis.text.x = element_text(size=14),
          axis.text.y = element_text(size=14),
          axis.title.x = element_blank(),
          axis.title.y = element_text(size=16),
          strip.text.x = element_text(size = 12))+
    guides(fill=guide_legend(title=NULL))+
    scale_fill_manual(values = dis.tissue$color)
  topptx(p,"/home/qians/MamDC/Result/WangKuster2019MSBRNAseq/Human/Picture/001.Human.Organdistribution.RNAlevel.pptx",
         width = 7.2, height = 4)
  ggsave(p, filename = "/home/qians/MamDC/Result/WangKuster2019MSBRNAseq/Human/Picture/001.Human.Organdistribution.RNAlevel.pdf",
         device = "pdf",width = 6, height = 4)
  #0.00357 X Esophagus
  fisher.test(as.data.frame(matrix(c(222,841,3336,18853),nrow = 2)))
  for (tissue in unique(a$tissue)) {
    f = sum(a$chr=="A" & a$tissue=="esophagus")
  }
}

## 2.2 WangKuster2019MSBRNAseq protein
{
  rm(list = ls())
  wd = "/home/qians/MamDC/Data/Seqdata/WangKuster2019MSBRNAseq"
  a = fread(file.path(wd,"iBAQprotein.WangKuster2019MSBRNAseq.table1.csv")) %>% as.data.frame()
  Species="Human"
  b <- read.csv(file.path("~/MamDC/Data/Ref",Species,paste0("gene",Species,".bed")),
                header = FALSE,sep = "\t",stringsAsFactors = F)
  a$chr = b[match(a$`Gene ID`,b$V4),1]
  a = a[,-c(1:2,32:35)]
  a %<>% dplyr::filter(.,chr %in% c(1:22,"X"))
  a$sample = colnames(a)[apply(a[,-ncol(a)], 1, which.max)]
  a[a$chr!="X","chr"]="A"
  dis.tissue = as.data.frame(t(table(a$chr,a$sample)/as.numeric(table(a$chr))))
  dis.tissue$color = rep(c("grey","#DDD7C6"),nrow(dis.tissue)/2)
  dis.tissue[dis.tissue$Var1=="Brain","color"]= "#3298c8"
  dis.tissue[dis.tissue$Var1=="Testis","color"]= "#ff6700"
  dis.tissue[dis.tissue$Var1=="Ovary","color"]= "#cc3397"
  dis.tissue[dis.tissue$Var1=="Esophagus","color"]= "black"
  p = ggplot(dis.tissue,aes(Var2,Freq,fill=Var1))+geom_bar(stat = "identity",position = "fill")+
    theme_classic()+ylab("Tissue distribution")+
    theme(axis.text.x = element_text(size=14),
          axis.text.y = element_text(size=14),
          axis.title.x = element_blank(),
          axis.title.y = element_text(size=16),
          strip.text.x = element_text(size = 12))+
    guides(fill=guide_legend(title=NULL))+
    scale_fill_manual(values = dis.tissue$color)
  topptx(p,"/home/qians/MamDC/Result/WangKuster2019MSBRNAseq/Human/Picture/001.Human.Organdistribution.proteinlevel.pptx",
         width = 4.2, height = 4)
  ggsave(p, filename = "/home/qians/MamDC/Result/WangKuster2019MSBRNAseq/Human/Picture/001.Human.Organdistribution.proteinlevel.pdf",
         device = "pdf",width = 6, height = 4)
  # 0.00395 X Gallbladder
  # 0.00593 X Esophagus
}
  