{ #Sampling N genes from Auto (N=number of X-linked genes) 
  #Ratio = mean(X)/mean(sampledAuto)
  ########### Note: Hindbrain_11w_Male_rep1 & Cerebellum_11w_Male_rep1
  rm(list = ls())
  wd = "~/MamDC/Data/Seqdata/CardosoKaessmann2019NatureRNAseq/"
  library(dplyr)
  library(tidyr)
  for( Species in c("Human","Mouse")) {
    a <- fread(paste0(wd,Species,"/",Species,".allgene.fpkm.txt")) %>% as.data.frame()
    colnames(a) %<>% gsub("KidneyTestis_4w_Male","Kidney_4w_Male",.)
    b <- read.csv(file.path("~/MamDC/Data/Ref",Species,paste0("gene",Species,".bed")),header = F,sep = "\t") 
    b$V5 %<>% gsub(" ","",.)
    a[,c("chr","type")] <- b[match(a$V1,b$V4),c(1,5)]
    a %<>% filter(chr %in% c(1:100,"X") & type %in% "protein_coding") #chr & coding
    a[a$chr!="X","chr"]="A"
    
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
      coldata$name %<>% gsub("KidneyTestis_4w_Male","Kidney_4w_Male",.) 
    }
    
    row.names(a) = a$V1
    a = a[,coldata$name]
    a[1:3,1:3]
    if (all(colnames(a)==coldata$name)) { 
      if (Species=="Human") {
        mbg = apply(a, 1, function(x)tapply(x, paste(coldata$condition,coldata$stage2,coldata$sex,sep = "_"), mean)) %>% 
          t() %>% as.data.frame()
      }else {
        mbg = apply(a, 1, function(x)tapply(x, paste(coldata$condition,coldata$stage,coldata$sex,sep = "_"), mean)) %>% 
          t() %>% as.data.frame()
      }
    }
    mbg$chr = b[match(row.names(mbg),b$V4),1] #FPKM_geneinfo table
    
    #colnames(mbg) %<>% gsub("Testis","Gonad",.) %>% gsub("Ovary","Gonad",.) 
    mbg[mbg$chr!="X","chr"]="A"
    focal.x = mbg[mbg$chr=="X",-ncol(mbg)]
    focal.auto = mbg[mbg$chr=="A",-ncol(mbg)]
    
    cutoff = 1
    mean.x = t(do.call(rbind, lapply(1:ncol(focal.x), function(i){
      mean(focal.x[focal.x[,i]>cutoff,i])}))) %>% as.data.frame()
    mean.auto = t(do.call(rbind, lapply(1:ncol(focal.auto), function(i){
      mean(focal.auto[focal.auto[,i]>cutoff,i])}))) %>% as.data.frame()
    ratio = as.data.frame(t(mean.x/mean.auto))
    row.names(ratio) = colnames(mbg)[-ncol(mbg)]
    ratio$tissue = lapply(row.names(ratio),function(x)strsplit(x,split = "_",fixed = T)[[1]][1]) %>% unlist()
    ratio$stage = lapply(row.names(ratio),function(x)strsplit(x,split = "_",fixed = T)[[1]][2]) %>% unlist()
    ratio$sex = lapply(row.names(ratio),function(x)strsplit(x,split = "_",fixed = T)[[1]][3]) %>% unlist()
    ratio$time = coldata[match(ratio$stage,coldata$stage2),"V9"]
    ratio = ratio[order(ratio$time),]
    ratio$tissue2 = row.names(ratio) %>% factor(.,levels = .)
    ratio$stage %<>% factor(.,levels = unique(.))
    ggplot(ratio,aes(stage,V1))+geom_point(aes(stage,V1,group=tissue,color=tissue))+
      theme(axis.text.x = element_blank())+
      #facet_grid(.~sex)+
      theme_classic()+
      geom_smooth(aes(stage,V1,group=tissue,color=tissue),position = "identity",method = "lm")
    
    ggplot(ratio,aes(tissue2,V1))+geom_point(aes(group=tissue,color=tissue))+
      facet_grid(.~sex)+geom_smooth(aes(group=tissue,color=tissue))+
      theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))+
      theme_classic()
    
    get.x.auto.ratio <- function(T,cutoff,trim=0){
      focal.one<- dplyr::select(mbg,c(chr,starts_with(T)))
      e <- matrix(rep(0,(ncol(focal.one)-1)*4),ncol = 4) %>% as.data.frame()
      colnames(e) <- c("sample","mean","CI0.5","CI0.95")
      e$sample <- colnames(focal.one)[-1]
      #bootstrap
      for (i in 2:ncol(focal.one)) {
        c <- focal.one[,c(1,i)]
        colnames(c) <- c("chr","express")
        c <- c[c$express>cutoff,] #filter gene
        d <- list()
        for (j in 1:1000) {
          d[[j]] <- mean(c[c$chr=="X","express"],trim = trim)/mean(sample(c[c$chr=="A","express"],
                                                                          size = sum(c$chr=="X"), replace = FALSE),trim = trim)
          boot.mean <- unlist(lapply(d, mean))
        }
        quantile(boot.mean, probs = c(0.05, 0.95))
        e[i-1,2:4] = c(mean(c[c$chr=="X","express"],trim = trim)/mean(c[c$chr=="A","express"],trim = trim),
                       quantile(boot.mean, probs = c(0.05, 0.95))[1],
                       quantile(boot.mean, probs = c(0.05, 0.95))[2]) %>% as.numeric()
      }
      return(e)
    }
    
    for (T in c("Brain","Cerebellum","Heart","Kidney","Liver","Gonad")) {
      test = get.x.auto.ratio(T,1)
      test$stage = lapply(test$sample,function(x)strsplit(x,split = "_",fixed = T)[[1]][2]) %>% unlist()
      
      if (Species=="Human") {
        test$stage2 = coldata[match(test$stage,coldata$stage2),"V9"]
        test$sex = lapply(test$sample,function(x)strsplit(x,split = "_",fixed = T)[[1]][3]) %>% unlist()
        test = test[order(test$stage2),]
        test$stage %<>% factor(.,levels = unique(.))
      } else {
        test$sex = lapply(test$sample,function(x)strsplit(x,split = "_",fixed = T)[[1]][3]) %>% unlist()
        test$stage %<>% factor(.,levels = c("e10.5","e11.5","e12.5","e13.5","e14.5","e15.5","e16.5","e17.5","e18.5","0dpb","3dpb","2wpb","4wpb","9wpb"))
        test$stage %<>% droplevels(.)
      }
      ggplot(test,aes(stage,mean))+
        geom_point(size=3.0,position = position_dodge(0.5))+
        geom_errorbar(aes(x = stage, ymax=CI0.5, ymin=CI0.95,width =0.3),position = position_dodge(0.5))+
        xlab("Developmental stage")+ylab("Expression ratio")+theme_bw()+facet_grid(sex ~ .)+
        theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(),
              panel.background = element_blank(),axis.title.x=element_text(size=20),
              axis.title.y=element_text(size=20),axis.text.x = element_text(size=18,angle = 45, hjust = 1),
              axis.text.y = element_text(size=18),strip.text.y = element_text(size = 16, angle = 90))+
        geom_hline(aes(yintercept=1),color="gray",linetype="dashed")+
        ylim(c(0,max(test$CI0.95)))
      ggsave(filename = paste0("/home/qians/MamDC/Result/CardosoKaessmann2019NatureRNAseq/",Species,"/Picture/002.",Species,
                               ".RatioStage.SampleAuto.",T,".pdf"),device = "pdf",width = 12, height = 7)
    }
  }
}



