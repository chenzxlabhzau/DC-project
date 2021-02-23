{ #Tissue specificity
  rm(list = ls())
  wd = "~/MamDC/Data/Seqdata/CardosoKaessmann2019NatureRNAseq/"
  library(dplyr)
  library(ggplotify)
  library(eoffice)
  for( Species in c("Human","Mouse")) {
    a <- fread(paste0(wd,Species,"/",Species,".allgene.fpkm.txt")) %>% as.data.frame()
    colnames(a) %<>% gsub("KidneyTestis_4w_Male","Kidney_4w_Male",.)
    b <- read.csv(file.path("~/MamDC/Data/Ref",Species,paste0("gene",Species,".bed")),header = F,sep = "\t") 
    b$V5 %<>% gsub(" ","",.)
    a[,c("chr","type")] <- b[match(a$V1,b$V4),c(1,5)]
    a %<>% filter(chr %in% c(1:100,"X") & type %in% "protein_coding")
    sample <- colnames(a)[2:(ncol(a)-2)] %>% as.data.frame()
    sample$tissue <- lapply(as.character(sample$.),function(x)(strsplit(x,"_")[[1]][1])) %>% unlist() %>%
      gsub("Forebrain","Brain",.) %>% gsub("Hindbrain","Cerebellum",.) %>% 
      gsub("Testis","Gonad",.) %>% gsub("Ovary","Gonad",.)
    
    sample$sex <- lapply(as.character(sample$.),function(x)(strsplit(x,"_")[[1]][3])) %>% unlist()
    sample$TS <- paste(sample$tissue,sample$sex,sep = "_")
    row.names(a) <- a$V1
    mbg <- t(apply(a[,2:(ncol(a)-2)],1,function(x){tapply(x,sample$TS,mean)})) %>% as.data.frame()
    mbg.sex.tau.all = c()
    for (Sex in c("Female","Male")) {
      mbg.sex <- dplyr::select(mbg,ends_with(Sex,ignore.case = FALSE))
      colnames(mbg.sex) %<>% lapply(., function(x)(strsplit(x,"_")[[1]][1])) %>% unlist()
      tau <- function(x){
        x$sum <- rowSums(x)
        x <- x[x$sum>0,-ncol(x)]
        x$max <- apply(x,1,max)
        x <- 1-x/x$max
        x <- x[,-ncol(x)]
        x$tau <- apply(x,1,sum)/(ncol(x)-1)
        return(as.data.frame(x))
      }
      mbg.sex.tau <- tau(mbg.sex)
      mbg.sex.tau$chr <- a[match(row.names(mbg.sex.tau),row.names(a)),"chr"]
      add.A = mbg.sex.tau[mbg.sex.tau$chr!="X",]
      add.A$chr="A"
      mbg.sex.tau = rbind(mbg.sex.tau,add.A)
      mbg.sex.tau$chr = factor(mbg.sex.tau$chr,levels = c(1:(length(unique(mbg.sex.tau$chr))-2),"A","X"))
      mbg.sex.tau$sex = Sex
      test <- wilcox.test(mbg.sex.tau[mbg.sex.tau$chr=="X" & mbg.sex.tau$sex==Sex,"tau"],
                           mbg.sex.tau[mbg.sex.tau$chr=="A" & mbg.sex.tau$sex==Sex,"tau"],alternative = "greater")
      if (test$p.value < 0.001) {
        p = ggplot(mbg.sex.tau,aes(chr,tau))+
          geom_boxplot(notch = T,fill=c(rep("#8db8d5",(length(unique(sort(mbg.sex.tau$chr)))-2)),"#035782","#9c0b1f"))+
          ylab(paste0("Tissue specificity in ",tolower(Species)," ",tolower(Sex)))+theme_bw()+
          theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(),
                axis.title.x=element_blank(),
                axis.text.x = element_text(size=16,angle = 45, hjust = 1, vjust = 1),
                axis.title.y=element_text(size=18),axis.text.y = element_text(size=16))+
          guides(fill=guide_legend(title=NULL))+theme(legend.position=c(.5,.9),legend.direction = "horizontal")+
          theme(panel.border = element_blank())+theme(axis.line = element_line(size=0.5, colour = "black"))+
          annotate("text", x=length(unique(mbg.sex.tau$chr)), y=0.05, label="***",size=6)
        ggsave(p, filename = paste0("/home/qians/MamDC/Result/CardosoKaessmann2019NatureRNAseq/",Species,"/Picture/",Species,".TissueSpecificity.",Sex,".pdf"),
               device = "pdf",width = 7, height = 7)
      }
      mbg.sex.tau.all = rbind(mbg.sex.tau.all,mbg.sex.tau)
    }
    #merge
      test1 <- wilcox.test(mbg.sex.tau.all[mbg.sex.tau.all$chr=="X" & mbg.sex.tau.all$sex=="Male","tau"],
                           mbg.sex.tau.all[mbg.sex.tau.all$chr=="A" & mbg.sex.tau.all$sex=="Male","tau"],alternative = "greater")
      test2 <- wilcox.test(mbg.sex.tau.all[mbg.sex.tau.all$chr=="X" & mbg.sex.tau.all$sex=="Female","tau"],
                           mbg.sex.tau.all[mbg.sex.tau.all$chr=="A" & mbg.sex.tau.all$sex=="Female","tau"],alternative = "greater")
      if (test1$p.value < 0.001 & test2$p.value < 0.001) {
        p = ggplot(mbg.sex.tau.all,aes(chr,tau))+
          geom_boxplot(notch = T,fill=rep(c(rep("#8db8d5",(length(unique(sort(mbg.sex.tau$chr)))-2)),"#035782","#9c0b1f"),2))+
          ylab(paste0("Tissue specificity in ",tolower(Species)))+theme_bw()+
          theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(),
                axis.title.x=element_blank(),
                axis.text.x = element_text(size=10,angle = 45, hjust = 1, vjust = 1),
                axis.title.y=element_text(size=12),axis.text.y = element_text(size=10),
                strip.background = element_rect(fill="white"))+
          guides(fill=guide_legend(title=NULL))+theme(legend.position=c(.5,.9),legend.direction = "horizontal")+
          theme(panel.border = element_blank())+theme(axis.line = element_line(size=0.5, colour = "black"))+
          annotate("text", x=length(unique(mbg.sex.tau$chr)), y=0.05, label="***",size=3)+
          facet_grid(~sex)
        topptx(p,filename = paste0("/home/qians/MamDC/Result/CardosoKaessmann2019NatureRNAseq/",Species,"/Picture/",Species,".TissueSpecificity.pptx"),
               width = 7, height = 7)
        ggsave(p,filename = paste0("/home/qians/MamDC/Result/CardosoKaessmann2019NatureRNAseq/",Species,"/Picture/",Species,".TissueSpecificity.pdf"),
               width = 7, height = 7)
    }
  }
}

{ #Stage specificity
  rm(list = ls())
  wd = "~/MamDC/Data/Seqdata/CardosoKaessmann2019NatureRNAseq/"
  library(dplyr)
  for( Species in c("Human","Mouse")) {
    a <- fread(paste0(wd,Species,"/",Species,".allgene.fpkm.txt")) %>% as.data.frame()
    colnames(a) %<>% gsub("KidneyTestis_4w_Male","Kidney_4w_Male",.) %>%
      gsub("Forebrain","Brain",.) %>% gsub("Hindbrain","Cerebellum",.) %>%
      gsub("Ovary","Gonad",.) %>% gsub("Testis","Gonad",.)
    b <- read.csv(file.path("~/MamDC/Data/Ref",Species,paste0("gene",Species,".bed")),header = F,sep = "\t") 
    b$V5 %<>% gsub(" ","",.)
    a[,c("chr","type")] <- b[match(a$V1,b$V4),c(1,5)]
    a %<>% filter(chr %in% c(1:100,"X") & type %in% "protein_coding")
    sample <- colnames(a)[2:(ncol(a)-2)] %>% as.data.frame()
    sample$tissue <- lapply(as.character(sample$.),function(x)(strsplit(x,"_")[[1]][1])) %>% unlist()
    sample$sex <- lapply(as.character(sample$.),function(x)(strsplit(x,"_")[[1]][3])) %>% unlist()
    row.names(a) <- a$V1
    c <- c()
    for (Tissue in unique(sort(sample$tissue))) {
      for (Sex in c("Female","Male")) {
        b <- dplyr::select(a,starts_with(Tissue)) %>% dplyr::select(.,contains(Sex,ignore.case = FALSE))
        tau <- function(x){
          x$sum <- rowSums(x)
          x <- x[x$sum>0,-ncol(x)]
          x$max <- apply(x,1,max)
          x <- 1-x/x$max
          x <- x[,-ncol(x)]
          x$tau <- apply(x,1,sum)/(ncol(x)-1)
          return(as.data.frame(x))
        }
        b <- tau(b)
        b$chr <- a[match(row.names(b),row.names(a)),"chr"]
        b <- b[b$chr!="Y",c("tau","chr")]
        b[b$chr!="X","chr"]="A"
        b$tissue <- Tissue
        b$sex <- Sex
        c <- rbind(c,b)
      }
    }
    c$sex <- factor(c$sex,levels = c("Male","Female"))
    
    testP = c()
    for (i in as.character(unique(c$sex))) {
      for (j in unique(c$tissue)) {
        test = wilcox.test(c[c$sex==i&c$tissue==j&c$chr=="A",1],c[c$sex==i&c$tissue==j&c$chr=="X",1],
                           alternative = "less")
        testP = c(testP,test$p.value)
      }
    }
    p = ggplot(c,aes(tissue,tau))+geom_boxplot(aes(fill=chr),notch = TRUE)+ylab(paste0("Stage specificity in ",tolower(Species)))+theme_bw()+
      theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(),
            axis.title.x=element_blank(),
            axis.text.x = element_text(size=10,angle = 45, hjust = 1, vjust = 1),
            axis.title.y=element_text(size=12),axis.text.y = element_text(size=10))+
      guides(fill=guide_legend(title=NULL))+theme(legend.position=c(.97,.08))+
      theme(panel.border = element_blank())+theme(axis.line = element_line(size=0.5, colour = "black"))+
      scale_fill_manual(values=c("#035782","#9c0b1f"))+
      facet_grid(. ~ sex )+theme(strip.text.x = element_text(size=12),strip.background = element_rect(fill="white"))+
      geom_text(aes(x,y,label =lab),
                data = data.frame(x = rep(1:6,2),
                                  y = ((testP>=0.05)/50)+1.03,
                                  lab = cut(testP,breaks = c(1,0.05,0.01,0.001,0),
                                            labels = rev(c("n.s.","*","**","***"))),
                                  sex = factor(rep(c("Female","Male"),each=6)),levels = c("Male","Female")),
                vjust = 1, size=4)
    #ggsave(paste0("/home/qians/MamDC/Result/CardosoKaessmann2019NatureRNAseq/",Species,"/Picture/",Species,".StageSpecificity.bothsex.pdf"),device = "pdf",p,width = 12, height = 7)
    topptx(p,filename = paste0("/home/qians/MamDC/Result/CardosoKaessmann2019NatureRNAseq/",Species,"/Picture/",Species,".StageSpecificity.bothsex.pptx"),
           width = 7, height = 7)
    }
}
