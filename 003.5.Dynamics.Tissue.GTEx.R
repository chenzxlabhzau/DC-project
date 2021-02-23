#GTEx RNA-seq data for DC analysis
{ #generate X:A ratio
  rm(list = ls())
  S="Human"
  a <- read.csv("~/Pseudo/Data/Seqdata/GTEx/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct",header = F,sep = "\t",stringsAsFactors = F )
  colnames(a) <- a[3,]
  a <- a[-c(1:3),-c(which(colnames(a)=="Cells - Cultured fibroblasts"),
                    which(colnames(a)=="Cells - EBV-transformed lymphocytes"),
                    which(colnames(a)=="Whole Blood"))]
  a$Name <- lapply(a$Name,function(x)strsplit(x,".",fixed=T)[[1]][1]) %>% unlist()
  b <- read.csv(file.path("/home/qians/Pseudo/Data/Ref",S,paste0("gene",S,".bed")),
                header = FALSE,sep = "\t",stringsAsFactors = F)
  a[,c("chr","type")] = b[match(a$Name,b$V4),c(1,5)]
  a %<>% dplyr::filter(.,chr %in% c(1:100,"X") & type == " protein_coding") %>% dplyr::select(.,-type)
  a[a$chr!="X","chr"] ="A"
  a = a[,-(1:2)]
  apply(a[,-ncol(a)], 2, function(x)tapply(as.numeric(x), a$chr, mean))
  
  e <- matrix(0,nrow = ncol(a)-1,ncol = 4) %>% as.data.frame()
  colnames(e) <- c("sample","mean","CI5","CI95")
  e$sample <- colnames(a)[-ncol(a)]
  cutoff = 1
  for (T in e$sample) {
    b = a[,c("chr",T)]
    colnames(b)[2] = "expression"
    b$expression %<>% as.numeric()
    b$expression <- b$expression/mean(b[b$chr=="A","expression"]) #Allexpress/Aexpress
    d <- list()
    for (j in 1:1000) {
      d[[j]] <- sample(b[b$chr=="X","expression"],size = floor(nrow(b
                                                                    [b$chr=="X",])/2), replace = TRUE)
      boot.mean <- unlist(lapply(d, mean))
    }
    quantile(boot.mean, probs = c(0.05, 0.95))
    e[e$sample==T,2:4] <- c(mean(b[b$chr=="X","expression"]),quantile(boot.mean, probs = c(0.05, 0.95))[1],quantile(boot.mean, probs = c(0.05, 0.95))[2]) %>% as.numeric()
  }
  write.table(e,file = "~/MamDC/Result/GTEx/Savedata/Human_tissues_XAratio.txt",sep = "\t",row.names = FALSE, col.names = TRUE,quote = FALSE)
}

{
  rm(list = ls())
  library(ggplotify)
  library(eoffice)
  e = read.csv("~/MamDC/Result/GTEx/Savedata/Human_tissues_XAratio.txt",
               sep = "\t",header = TRUE,stringsAsFactors = F)
  #e$sample %<>% Hmisc::capitalize(.) %>% gsub("_"," ",.) %>% gsub("."," ",.,fixed = T)
  e = e[order(e$mean,decreasing = FALSE),]
  e$sample %<>% factor(.,levels = .)
  p = ggplot(e,aes(sample,mean))+geom_point(size=3.0)+
    geom_errorbar(aes(x = sample, ymax=CI5, ymin=CI95,width =0.3))+
    ylab("X:A ratio")+theme_classic()+
    theme(axis.title.y=element_text(size=18),axis.title.x=element_blank(),
          axis.text.y = element_text(size=16),axis.text.x = element_text(size=12,angle = 90,hjust = 1,vjust = 0.5))+
    #scale_y_continuous(breaks = c(0,0.5,1,2),labels = c(0,0.5,1,2))+
    geom_hline(aes(yintercept=1),color="black",linetype="dashed")+
    geom_hline(aes(yintercept=0.5),color="black",linetype="dashed")
  topptx(p,"~/MamDC/Result/GTEx/Picture/Expressionratio.alltissues.GTEx.pptx",
         width = 7, height = 7)
}


