{ #trim from 0% to 5%(both side)
  rm(list = ls())
  wd = "~/MamDC/Data/Seqdata/CardosoKaessmann2019NatureRNAseq/"
  library(dplyr)
  library(tidyr)
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
    }
    
    for (T in c("Brain","Cerebellum","Heart","Kidney","Liver","Testis","Ovary")) {
      focal.tissue <- dplyr::select(a,c(chr,starts_with(T)))
      ratio = matrix(0,ncol = 11,nrow = ncol(focal.tissue)-1) %>% as.data.frame()
      colnames(ratio) = paste0("trim",(0:10)*0.01)
      row.names(ratio) = colnames(focal.tissue)[-1]
      for (n in 1:11) {
        trim = (n-1)*0.01/2
        mean.x = apply(focal.tissue[focal.tissue$chr=="X",-1], 2, function(x)mean(x,trim = trim))
        mean.aotu = apply(focal.tissue[focal.tissue$chr=="A",-1], 2, function(x)mean(x,trim = trim))
        ratio[,n] = as.numeric(mean.x/mean.aotu )
      }
      ratio$tissue = row.names(ratio)
      
     if (Species=="Human") {
       ratio$time = coldata[match(ratio$tissue,coldata$name),9]
       ratio = ratio[order(ratio$time),]
       ratio$tissue %<>% factor(.,levels = .)
     }
      if (Species=="Mouse") {
        ratio$tissue %<>% factor(.,levels = coldata[coldata$condition==T,"name"])
      }
      
      ratio %<>% pivot_longer(.,cols=1:11)
      ggplot(ratio,aes(tissue,value))+geom_line(aes(group = name,color=name))+theme_classic()+
        theme(axis.title.x= element_blank(),
              axis.text.x = element_text(size=14,angle = 90, hjust = 1, vjust = 0.5),
              axis.title.y=element_text(size=18),axis.text.y = element_text(size=16))+
        ylab("X:A ratio")+coord_cartesian(ylim = c(0,1.5))+
        geom_hline(yintercept = 1,color="black",linetype="dashed")+
        scale_color_manual(values=c("#61439A",rep("grey",4),"#EE1F26",rep("grey",4),"#862461"),name="Cutoff")
      ggsave(filename = paste0("/home/qians/MamDC/Result/CardosoKaessmann2019NatureRNAseq/",Species,"/Picture/002.",Species,
                               ".ExpressRatio.trim.",T,".pdf"),device = "pdf",width = 11, height = 7)
    }
  }
}

  
      