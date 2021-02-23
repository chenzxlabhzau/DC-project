{ # For human
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
  
  for (T in unique(coldata$condition)) {
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
    
    df = matrix(0,nrow = 10,ncol = 2) %>% as.data.frame()
    colnames(df) = c("A","X") #,"Direction"
    for (n in 1:5) {
      df[n,] = c(table(focal.express[focal.express$R> (n+4)*0.1 &focal.express$P<0.05,"chr"])[1],
                 table(focal.express[focal.express$R> (n+4)*0.1 &focal.express$P<0.05,"chr"])[2])
        #as.numeric(table(focal.express[focal.express$R>=(n+4)*0.1 &focal.express$P<0.05,"chr"]))
      df[n+5,] = c(table(focal.express[focal.express$R< -(n+4)*0.1 &focal.express$P<0.05,"chr"])[1],
                   table(focal.express[focal.express$R< -(n+4)*0.1 &focal.express$P<0.05,"chr"])[2])
    }
    df[is.na(df)]=0
    df$ratio = df$X/(df$A+df$X)
    df$cutoff = paste0("R>",rep(seq(0.5,0.9,by = 0.1),2))
    df$Direction = rep(c("Positive","Negative"),each=5) %>% factor(.,levels = c("Positive","Negative"))
    
    average = table(focal.express$chr)[2]/nrow(focal.express)*100
    p = ggplot(df,aes(cutoff,ratio*100,fill=Direction))+geom_bar(stat = "identity",position = "dodge")+
      theme_classic()+xlab("Threshold on correlation coefficient")+
      ylab(paste0("(%) X-linked genes in ",tolower(T)))+
      theme(axis.text.x = element_text(size=12),axis.title.x = element_text(size=14),
            axis.text.y = element_text(size=12),axis.title.y = element_text(size=14),
            legend.position=c(.35,.9),legend.direction = "horizontal",
            legend.text = element_text(size = 10))+
      guides(fill=guide_legend(title=NULL))+
      coord_cartesian(ylim = c(0,7))+scale_y_continuous(expand = c(0, 0))+
      geom_hline(yintercept = average,color="grey",linetype="dashed")
    if (T=="Liver") {
      p = ggplot(df[c(1:4,6:9),],aes(cutoff,ratio*100,fill=Direction))+geom_bar(stat = "identity",position = "dodge")+
        theme_classic()+xlab("Threshold on correlation coefficient")+
        ylab(paste0("(%) X-linked genes in ",tolower(T)))+
        theme(axis.text.x = element_text(size=12),axis.title.x = element_text(size=14),
              axis.text.y = element_text(size=12),axis.title.y = element_text(size=14),
              legend.position=c(.35,.9),legend.direction = "horizontal",
              legend.text = element_text(size = 10))+
        guides(fill=guide_legend(title=NULL))+
        coord_cartesian(ylim = c(0,7))+scale_y_continuous(expand = c(0, 0))+
        geom_hline(yintercept = average,color="grey",linetype="dashed")
    }
    if (T=="Testis") {
      p = ggplot(df,aes(cutoff,ratio*100,fill=Direction))+geom_bar(stat = "identity",position = "dodge")+
        theme_classic()+xlab("Threshold on correlation coefficient")+
        ylab(paste0("(%) X-linked genes in ",tolower(T)))+
        theme(axis.text.x = element_text(size=12),axis.title.x = element_text(size=14),
              axis.text.y = element_text(size=12),axis.title.y = element_text(size=14),
              legend.position=c(.35,.9),legend.direction = "horizontal",
              legend.text = element_text(size = 10))+
        guides(fill=guide_legend(title=NULL))+
        coord_cartesian(ylim = c(0,15))+scale_y_continuous(expand = c(0, 0))+
        geom_hline(yintercept = average,color="grey",linetype="dashed")
    }
      #scale_y_continuous(breaks = c(0,2,4,6),labels = c(0%,2%,4%,6%))
    ggsave(p,filename = paste0("/home/qians/MamDC/Result/CardosoKaessmann2019NatureRNAseq/",S,"/Picture/003.",S,
                             ".RatioXlinked.correlation.",T,".pdf"),device = "pdf",width = 5.5, height = 4.5)
    topptx(p,filename = paste0("/home/qians/MamDC/Result/CardosoKaessmann2019NatureRNAseq/",S,"/Picture/003.",S,
                               ".RatioXlinked.correlation.",T,".pptx"),width = 4.5, height = 3.5)
    
    dynamic = focal.express[abs(focal.express$R) > 0.8&focal.express$P<0.05,] #positive & negative
    library(org.Hs.eg.db)
    library(clusterProfiler)
    ego.positive <- enrichGO(gene = row.names(dynamic[dynamic$R > 0.8,]),
                             OrgDb = org.Hs.eg.db,
                             keyType = "ENSEMBL",
                             ont = "BP",
                             pAdjustMethod = "BH",
                             pvalueCutoff  = 0.01,
                             qvalueCutoff  = 0.05)
    write.csv(as.data.frame(ego.positive),
              file = paste0("/home/qians/MamDC/Result/CardosoKaessmann2019NatureRNAseq/",S,"/Picture/003.",S,
                                         ".Positivecorrelation.GObp.",T,".csv"))
    ego.negative <- enrichGO(gene = row.names(dynamic[dynamic$R < -0.8,]),
                             OrgDb = org.Hs.eg.db,
                             keyType = "ENSEMBL",
                             ont = "BP",
                             pAdjustMethod = "BH",
                             pvalueCutoff  = 0.01,
                             qvalueCutoff  = 0.05)
    write.csv(as.data.frame(ego.negative),
              file = paste0("/home/qians/MamDC/Result/CardosoKaessmann2019NatureRNAseq/",S,"/Picture/003.",S,
                            ".Negativecorrelation.GObp.",T,".csv")) #testis Rho< -0.7 as no GO result under rho< -0.8 
    
    dynamic = focal.express[focal.express$R > 0.8&focal.express$P<0.05,] #positive
    library(tidyverse)
    dynamic$gene = row.names(dynamic)
    dynamic = pivot_longer(dynamic,cols = 1:nrow(focal.tissue))
    dynamic$stage = coldata[match(dynamic$name,coldata$name),9]
    #dynamic %<>% dplyr::filter(.,chr=="X")
    dynamic$stage %<>% factor(.,levels = unique(.))
    p = ggplot(dynamic,aes(stage,log10(value+1),fill=stage))+geom_boxplot(outlier.colour = "white",notch = TRUE)+
      geom_smooth(se=TRUE, aes(group=1),color="#FDFC93")+
      theme_classic()+xlab("From young to old")+ylab("Expression level (log10)")+
      theme(axis.text.x = element_blank(),axis.title.x = element_text(size=16),
            axis.text.y = element_text(size=14),axis.title.y = element_text(size=16),axis.ticks.x = element_blank())+
      theme(axis.line = element_line(arrow = arrow(length = unit(0.4, 'cm'))))+
      guides(fill=FALSE)
    ggsave(p,filename = paste0("/home/qians/MamDC/Result/CardosoKaessmann2019NatureRNAseq/",S,"/Picture/003.",S,
                             ".Positivecorrelation.ExpressStage.",T,".pdf"),device = "pdf",width = 4.5, height = 5)
    topptx(p,filename = paste0("/home/qians/MamDC/Result/CardosoKaessmann2019NatureRNAseq/",S,"/Picture/003.",S,
                               ".Positivecorrelation.ExpressStage.",T,".pptx"),width = 4.5, height = 5)
    
    dynamic = focal.express[focal.express$R < -0.8&focal.express$P<0.05,] #negative
    library(tidyverse)
    dynamic$gene = row.names(dynamic)
    dynamic = pivot_longer(dynamic,cols = 1:nrow(focal.tissue))
    dynamic$stage = coldata[match(dynamic$name,coldata$name),9]
    #dynamic %<>% dplyr::filter(.,chr=="X")
    dynamic$stage %<>% factor(.,levels = unique(.))
    p = ggplot(dynamic,aes(stage,log10(value+1),fill=stage))+geom_boxplot(outlier.colour = "white",notch = TRUE)+
      geom_smooth(se=TRUE, aes(group=1),color="#FDFC93")+
      theme_classic()+xlab("From young to old")+ylab("Expression level (log10)")+
      theme(axis.text.x = element_blank(),axis.title.x = element_text(size=16),
            axis.text.y = element_text(size=14),axis.title.y = element_text(size=16),axis.ticks.x = element_blank())+
      theme(axis.line = element_line(arrow = arrow(length = unit(0.4, 'cm'))))+
      guides(fill=FALSE)
    ggsave(p,filename = paste0("/home/qians/MamDC/Result/CardosoKaessmann2019NatureRNAseq/",S,"/Picture/003.",S,
                             ".Negativecorrelation.ExpressStage.",T,".pdf"),device = "pdf",width = 4.5, height = 5)
    topptx(p,filename = paste0("/home/qians/MamDC/Result/CardosoKaessmann2019NatureRNAseq/",S,"/Picture/003.",S,
                               ".Negativecorrelation.ExpressStage.",T,".pptx"),width = 4.5, height = 5)
  }
}
  