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
    row.names(a) = a$V1
    coldata = read.csv(paste0("~/Pseudo/Result/",Species,"/Savedata/coldata.csv"),row.names = 1,sep = ",")
    if (Species=="Human") {
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
      coldata <- coldata[order(coldata$V9),]
      coldata$name %<>% gsub("KidneyTestis_4w_Male","Kidney_4w_Male",.) %>%
        gsub("Forebrain","Brain",.) %>% gsub("Hindbrain","Cerebellum",.)
      coldata$rep = lapply(coldata$name,function(x)(strsplit(x,"_")[[1]][4])) %>% unlist()
      #coldata$name2 = paste(coldata$condition,coldata$stage2,coldata$sex,coldata$rep,sep = "_")
    }
    
    tissue = coldata[coldata$V9 %in% c(28,501,6485) & coldata$condition %in% c("Brain","Liver","Testis"),]
    express = a[,tissue$name]
    mbg = t(apply(express, 1, function(x)tapply(x,paste(tissue$condition,tissue$stage2,sep = "_"),mean))) %>% as.data.frame()
    b <- read.csv(file.path("~/MamDC/Data/Ref",Species,paste0("gene",Species,".bed")),header = F,sep = "\t") 
    b$V5 %<>% gsub(" ","",.)
    mbg[,c("chr","type")] <- b[match(row.names(mbg),b$V4),c(1,5)]
    mbg %<>% filter(chr %in% c(1:100,"X") & type %in% "protein_coding")
    mbg[mbg$chr!="X","chr"]="A"
    
    if (Species=="Human") {
      age = read.csv("~/MamDC/Data/Seqdata/Genorigin/Human/Homo_sapiens.csv",header = TRUE,sep = ",",stringsAsFactors = FALSE)
      age$gene_age %<>% gsub(">","",.) %>% as.numeric() 
      age$group = cut(as.numeric(age$gene_age),breaks = c(0,332,645,5000),labels = c("Young","Middle","Old"))
      mbg$group = age[match(row.names(mbg),age$ensembl_id),4]
    }
    
    for (T in c("Brain","Liver","Testis")) {
      focal.express = mbg %>% dplyr::select(.,c(starts_with(T),chr,group))
      focal.express = focal.express[apply(focal.express[,1:3], 1, min) > 1,]
      colnames(focal.express)[1:3] %<>% gsub(paste0(T,"_"),"",.) %>% 
        gsub("17ypb","Adult",.) %>% gsub("221dpb","Infant",.) %>% gsub("4wpc","Embryo",.)
      for (i in 1:3) { #normalize
        for (group in c("Young","Middle","Old")) {
          focal.express[focal.express$group == group,i] = focal.express[focal.express$group == group,i] / 
            median(focal.express[focal.express$group == group & focal.express$chr=="A",i])
        }
      }
      focal.express %<>% pivot_longer(.,cols=1:3) %>% as.data.frame()
      testP = c()
      for (i in c("Young","Middle","Old")) {
        for (j in c("Embryo","Infant","Adult")) {
          test = wilcox.test(focal.express[focal.express$group==i & focal.express$name==j & focal.express$chr=="A","value"],
                             focal.express[focal.express$group==i & focal.express$name==j & focal.express$chr=="X","value"])
          testP = c(testP,test$p.value)
        }
      }
      focal.express$group %<>% factor(.,levels = c("Young","Middle","Old"))
      focal.express$name %<>% factor(.,levels = c("Embryo","Infant","Adult"))
      p = ggplot(focal.express,aes(group,log2(value)))+geom_boxplot(aes(fill=chr),notch = TRUE,outlier.alpha = 0)+
        ylab(paste0("X:A ratio in ",T))+xlab("Gene age")+theme_bw()+
        theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(),
              axis.title.x=element_text(size=12),
              axis.text.x = element_text(size=10,angle = 45, hjust = 1, vjust = 1),
              axis.title.y=element_text(size=12),axis.text.y = element_text(size=10))+
        #guides(fill=guide_legend(title=NULL))+theme(legend.position=c(.97,.08))+
        guides(fill=FALSE)+
        theme(panel.border = element_blank())+theme(axis.line = element_line(size=0.5, colour = "black"))+
        scale_fill_manual(values=c("#035782","#9c0b1f"))+
        facet_grid(. ~ name )+theme(strip.text.x = element_text(size=12),strip.background = element_rect(fill="white"))+
        geom_hline(yintercept = 0,color="black",linetype="dashed")+
        geom_hline(yintercept = 1,color="grey",linetype="dashed")+
        geom_hline(yintercept = -1,color="grey",linetype="dashed")+
        coord_cartesian(ylim = c(-3,9))+
        scale_y_continuous(breaks = c(-1,0,1),labels = c(-1,0,1))+
        geom_text(aes(x,y,label =lab),
                  data = data.frame(x = rep(1:3,each=3),
                                    y = ((testP>=0.05)/50)+8,
                                    lab = cut(testP,breaks = c(1,0.05,0.01,0.001,0),
                                              labels = rev(c("n.s.","*","**","***"))),
                                    name = factor(rep(c("Embryo","Infant","Adult"),3),levels = c("Embryo","Infant","Adult"))),
                  vjust = 1, size=6)
      
    ggsave(p,filename = file.path("/home/qians/MamDC/Result/Genorigin",Species,"Picture",paste0("S003.XAratio.age.stage.",T,".pdf")),
           width = 5,height = 4)
    topptx(p,filename = file.path("/home/qians/MamDC/Result/Genorigin",Species,"Picture",paste0("S003.XAratio.age.stage.",T,".pptx")),
           width = 5,height = 4)
    }
  }
}


