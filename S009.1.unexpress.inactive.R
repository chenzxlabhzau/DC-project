{ #display merged results
  rm(list = ls())
  wd = "~/MamDC/Data/Seqdata/CardosoKaessmann2019NatureRNAseq/"
  library(dplyr)
  library(tidyr)
  library(ggplotify)
  library(eoffice)
  for( Species in c("Human","Mouse")) {
    a <- fread(paste0(wd,Species,"/",Species,".allgene.fpkm.txt")) %>% as.data.frame()
    colnames(a) %<>% gsub("KidneyTestis_4w_Male","Kidney_4w_Male",.) %>%
      gsub("Forebrain","Brain",.) %>% gsub("Hindbrain","Cerebellum",.)
    b <- read.csv(file.path("~/MamDC/Data/Ref",Species,paste0("gene",Species,".bed")),header = F,sep = "\t") 
    b$V5 %<>% gsub(" ","",.)
    a[,c("chr","type")] <- b[match(a$V1,b$V4),c(1,5)]
    a %<>% filter(chr %in% c(1:100,"X") & type %in% "protein_coding") #chr & coding
    a[a$chr!="X","chr"]="A"
    
    
    tapply(a$max<1, a$chr, sum)/table(a$chr)
    tapply(a$max, a$chr, median)
    a$max = apply(a[,-c(1,ncol(a),ncol(a)-1)],1,max)
    a$max =a$max/median(a[a$chr=="A","max"])
    tapply(a$max, a$chr, median)
    
    p = ggplot(a,aes(chr,log2(max),fill=chr))+geom_boxplot(notch = TRUE,outlier.colour = "white")+
      theme_classic()+ylab("Relative maximum expression level (log2)")+
      theme(axis.text.x = element_text(size=10),
            axis.text.y = element_text(size=10),
            axis.title.x = element_blank(),
            axis.title.y = element_text(size=12))+
      scale_fill_manual(values = c("#035782","#9c0b1f"))+
      coord_cartesian(ylim = c(-6,6))+
      geom_hline(yintercept = 0,color="black",linetype="dashed")+
      geom_hline(yintercept = -1,color="grey",linetype="dashed")+
      scale_y_continuous(breaks = c(-6,-3,-1,0,3,6),labels =  c(-6,-3,-1,0,3,6))+
      guides(fill=FALSE)
    ggsave(filename = paste0("/home/qians/MamDC/Result/CardosoKaessmann2019NatureRNAseq/",Species,"/Picture/S009.",Species,".unexpress.inactive.allstages.pdf"),
           p,device = "pdf", width = 4,height = 3.75)
                        