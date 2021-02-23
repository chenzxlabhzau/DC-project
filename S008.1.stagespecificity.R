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
    }else {
      for (i in 1:nrow(coldata)) {
        if (substr(coldata[i,2],1,1)=="e") {
          coldata[i,6] = substr(coldata[i,2],2,5) %>% as.numeric()
        }
        if (substr(coldata[i,2],2,4)== "dpb") {
          coldata[i,6] = substr(coldata[i,2],1,1) %>% as.numeric() +20 #20 day
        }
        if (substr(coldata[i,2],2,4)== "wpb") {
          coldata[i,6] = substr(coldata[i,2],1,1) %>% as.numeric() * 7 +20 #20 day
        }
      }
    }
    
    for (T in c("Brain","Cerebellum","Heart","Kidney","Liver","Testis","Ovary")) {
      focal.tissue <- dplyr::select(a,c(chr,starts_with(T)))
      ratio = matrix(0,ncol = 12,nrow = ncol(focal.tissue)-1) %>% as.data.frame()
      colnames(ratio) = c("All",(0:10)*0.1)
      row.names(ratio) = colnames(focal.tissue)[-1]
      
      focal.tissue.x = focal.tissue[focal.tissue$chr=="X",-1]
      focal.tissue.auto = focal.tissue[focal.tissue$chr=="A",-1]
      
      for (n in 1:12) {
        cutoff = (n-2)*0.1
        mean.x = t(do.call(rbind, lapply(1:ncol(focal.tissue.x), function(i){
          mean(focal.tissue.x[focal.tissue.x[,i]>cutoff,i])}))) %>% as.data.frame()
        mean.auto = t(do.call(rbind, lapply(1:ncol(focal.tissue.auto), function(i){
          mean(focal.tissue.auto[focal.tissue.auto[,i]>cutoff,i])}))) %>% as.data.frame()
        ratio[,n] = as.numeric(mean.x/mean.auto )
      }
      ratio$tissue = row.names(ratio)# 
      
      #changename = coldata[grepl("CS",coldata$stage),2:3] %>% dplyr::distinct(.,stage,stage2)
      changename = coldata %>% dplyr::distinct(.,stage,stage2) %>% dplyr::filter(.,stage != stage2)
      
      if (Species=="Human") {
        ratio$time = coldata[match(ratio$tissue,coldata$name),9]
        ratio[grepl("11w",ratio$tissue),"time"]=77 # Cerebellum_11w_Male_rep1.1
        ratio$tissue %<>% gsub("rep1.1","rep3",.) %>% gsub(paste0(T,"_"),"",.)
        for (i in 1:nrow(changename)) {
          ratio$tissue %<>% gsub(paste0(changename[i,1],"_"),paste0(changename[i,2],"_"),.)
        }
        
        ratio$group = lapply(ratio$tissue,function(x)strsplit(x,split = "_",fixed = TRUE)[[1]][1]) %>% unlist()
        df = apply(ratio[,1:12], 2, function(x)tapply(x,ratio$group,mean)) %>% as.data.frame()
        df$time = ratio[match(row.names(df),ratio$group),"time"]
        df = df[order(df$time),]
        df$time2 = row.names(df) %>% factor(.,levels = .)
        born = max(grep("wpc",df$time2))
        for (i in 1:12) {
          print(cor.test(df[,i],log10(df$time),method = "spearman"))
        }
        cor = cor.test(df[,12],df$time,method = "spearman")
        df %<>% pivot_longer(.,cols=1:12)
        cor.test(df$value,log10(df$time))
        df$name %<>% factor(.,levels = c("All",seq(0,1,0.1)))
        p = ggplot(df,aes(time2,value))+geom_line(aes(group = name,color=name))+theme_classic()+
          theme(axis.title.x= element_blank(),
                axis.text.x = element_text(size=10,angle = 45, hjust = 1, vjust = 1),
                axis.title.y=element_text(size=18),axis.text.y = element_text(size=16))+
          ylab(paste0("X:A ratio in ",tolower(T)))+coord_cartesian(ylim = c(0,1.8))+
          geom_hline(yintercept = 1,color="black",linetype="dashed")+
          geom_hline(yintercept = 0.5,color="black",linetype="dashed")+
          geom_vline(xintercept = born+1,color="grey")+
          scale_color_manual(values=c("black","#fc9272",rep("grey",9),"#A40C0E"),name="Cutoff")+
          annotate(geom="text", x=5, y=1.6, size=5,
                   label=paste0("rho=",round(cor$estimate,2),",\n P=",format(cor$p.value,2,scientific = TRUE)))
        #ggsave(filename = paste0("/home/qians/MamDC/Result/CardosoKaessmann2019NatureRNAseq/",Species,"/Picture/002.",Species,".ExpressRatio.FPKM.",T,".pdf"),device = "pdf",width = 12, height = 7)
        topptx(p,
               paste0("/home/qians/MamDC/Result/CardosoKaessmann2019NatureRNAseq/",Species,"/Picture/002.",Species,".ExpressRatio.FPKM.",T,".merged.pptx"),
               width = 7,height = 5)
        ggsave(filename = paste0("/home/qians/MamDC/Result/CardosoKaessmann2019NatureRNAseq/",Species,"/Picture/002.",Species,".ExpressRatio.FPKM.",T,".merged.pdf"),
               p,device = "pdf", width = 7,height = 5)
      }
      if (Species=="Mouse") {
        #ratio$tissue %<>% gsub(paste0(T,"_"),"",.) %>% factor(.,levels = coldata[coldata$condition==T,"name"])
        ratio$time = coldata[match(ratio$tissue,coldata$name),6]
        ratio$tissue %<>% gsub(paste0(T,"_"),"",.)
        ratio$group = lapply(ratio$tissue,function(x)strsplit(x,split = "_",fixed = TRUE)[[1]][1]) %>% unlist()
        df = apply(ratio[,1:12], 2, function(x)tapply(x,ratio$group,mean)) %>% as.data.frame()
        df$time = ratio[match(row.names(df),ratio$group),"time"]
        df = df[order(df$time),]
        embryo.length = sum(!grepl("pb",row.names(df)))
        df$time2 = row.names(df) %>% gsub("s","",.) %>% paste0(rep(c("e",""),c(embryo.length,nrow(df)-embryo.length)),.) %>% factor(.,levels = .)
        born = min(grep("pb",df$time2))
        for (i in 1:12) {
          print(cor.test(df[,i],log10(df$time),method = "spearman"))
        }
        cor = cor.test(df[,12],df$time,method = "spearman")
        df %<>% pivot_longer(.,cols=1:12)
        cor.test(df$value,log10(df$time))
        df$name %<>% factor(.,levels = c("All",seq(0,1,0.1)))
        p = ggplot(df,aes(time2,value))+geom_line(aes(group = name,color=name))+theme_classic()+
          theme(axis.title.x= element_blank(),
                axis.text.x = element_text(size=10,angle = 45, hjust = 1, vjust = 1),
                axis.title.y=element_text(size=18),axis.text.y = element_text(size=16))+
          ylab(paste0("X:A ratio in ",tolower(T)))+coord_cartesian(ylim = c(0,1.8))+
          geom_hline(yintercept = 1,color="black",linetype="dashed")+
          geom_hline(yintercept = 0.5,color="black",linetype="dashed")+
          geom_vline(xintercept = born,color="grey")+
          scale_color_manual(values=c("black","#fc9272",rep("grey",9),"#A40C0E"),name="Cutoff")+
          annotate(geom="text", x=4, y=1.6, size=5,
                   label=paste0("rho=",round(cor$estimate,2),",\n P=",format(cor$p.value,2,scientific = TRUE)))
        #ggsave(filename = paste0("/home/qians/MamDC/Result/CardosoKaessmann2019NatureRNAseq/",Species,"/Picture/002.",Species,".ExpressRatio.FPKM.",T,".pdf"),device = "pdf",width = 12, height = 7)
        topptx(p,
               paste0("/home/qians/MamDC/Result/CardosoKaessmann2019NatureRNAseq/",Species,"/Picture/002.",Species,".ExpressRatio.FPKM.",T,".merged.pptx"),
               width = 7,height = 5)
        ggsave(filename = paste0("/home/qians/MamDC/Result/CardosoKaessmann2019NatureRNAseq/",Species,"/Picture/002.",Species,".ExpressRatio.FPKM.",T,".merged.pdf"),
               p,device = "pdf", width = 7,height = 5)
      }
    }
  }
}
