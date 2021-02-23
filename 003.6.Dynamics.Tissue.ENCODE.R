#GTEx RNA-seq data for DC analysis
{#generate FPKM
  rm(list = ls())
  S="Mouse"
  a = fread("~/Pseudo/Data/Seqdata/ENCODE/counts.tsv") %>% as.data.frame()
  a$V1 <- lapply(a$V1,function(x)strsplit(x,".",fixed=T)[[1]][1]) %>% unlist()
  b <- read.csv(file.path("/home/qians/Pseudo/Data/Ref",S,paste0("gene",S,".bed")),
                header = FALSE,sep = "\t",stringsAsFactors = F)
  
  row.names(a) = a$V1
  a = a[,-1]
  #exon length
  length <- read.csv(file = file.path("~/Pseudo/Data/Ref",S,paste0(S,".nonreExonlength.bed")),
                     header = FALSE,sep = "\t",stringsAsFactors = FALSE)
  #nonredundant gene length
  gene.length <- tapply(length$V5, length$V4, sum) %>% as.data.frame() 
  
  a$length = gene.length[match(row.names(a),row.names(gene.length)),1]
  
  totalcounts <- colSums(a[,-ncol(a)])
  rpkm <- t(do.call(rbind, lapply(1:length(totalcounts), function(i){
    10^9*a[,i]/a$length/totalcounts[i]}))) %>% as.data.frame()
  row.names(rpkm) = row.names(a)
  colnames(rpkm) = colnames(a)[-ncol(a)]
  write.table(rpkm,file = "/home/qians/MamDC/Result/ENCODE/Savedata/Mouse.allgene.fpkm.txt",
              sep = "\t",quote = F,row.names = T,col.names = T)
}

{#ratio
  rm(list = ls())
  gc()
  rm(list = ls())
  S="Mouse"
  a = read.csv("/home/qians/MamDC/Result/ENCODE/Savedata/Mouse.allgene.fpkm.txt",
               header = T,sep = "\t",stringsAsFactors = F)
  b <- read.csv(file.path("/home/qians/Pseudo/Data/Ref",S,paste0("gene",S,".bed")),
                header = FALSE,sep = "\t",stringsAsFactors = F)
  a[,c("chr","type")] = b[match(row.names(a),b$V4),c(1,5)]
  a %<>% dplyr::filter(.,chr %in% c(1:100,"X") & type == " protein_coding") %>% dplyr::select(.,-type)
  a[a$chr!="X","chr"] ="A"
  info = fread("~/Pseudo/Data/Seqdata/ENCODE/sample.srr.link",header = F) %>% as.data.frame()
  name = colnames(a)[-ncol(a)] %>% as.data.frame()
  name$V2 = info[match(name$.,info$V2),1] 
  name$V2 %<>% gsub("CNS","CentralNervousSystem",.) %>%
    gsub("Bladder","UrinaryBladder",.) %>%  gsub("SubcFatPad","SubcutaneousAdiposeTissue",.) %>% gsub("GenitalFatPad","GenitalAdiposeTissue",.) %>%
    gsub("LgIntestine","LargeIntestine",.) %>% gsub("Adrenal","AdrenalGland",.) %>%
    gsub("SmIntestine","SmallIntestine",.)
  if (all(colnames(a)[-ncol(a)]==name$.)) {
    mbg = t(apply(a[,-ncol(a)], 1, function(x)tapply(x, name$V2, mean))) %>% as.data.frame()
  }
  if (all(row.names(a)==row.names(mbg))) {
    mbg$chr = a$chr
  }
  
  e <- matrix(0,nrow = ncol(mbg)-1,ncol = 4) %>% as.data.frame()
  colnames(e) <- c("sample","mean","CI5","CI95")
  e$sample <- colnames(mbg)[-ncol(mbg)]
  cutoff = 1
  for (T in e$sample) {
    b = mbg[,c("chr",T)]
    colnames(b)[2] = "expression"
    b$expression %<>% as.numeric()
    b$expression <- b$expression/mean(b[b$chr=="A","expression"]) #Allexpress/Aexpress
    d <- list()
    for (j in 1:1000) {
      d[[j]] <- sample(b[b$chr=="X","expression"],size = floor(nrow(b[b$chr=="X",])/2), replace = TRUE)
      boot.mean <- unlist(lapply(d, mean))
    }
    quantile(boot.mean, probs = c(0.05, 0.95))
    e[e$sample==T,2:4] <- c(mean(b[b$chr=="X","expression"]),quantile(boot.mean, probs = c(0.05, 0.95))[1],quantile(boot.mean, probs = c(0.05, 0.95))[2]) %>% as.numeric()
  }
  write.table(e,file = "~/MamDC/Result/ENCODE/Savedata/Mouse_tissues_XAratio.txt",sep = "\t",row.names = FALSE, col.names = TRUE,quote = FALSE)
  }

{#
  rm(list = ls())
  library(ggplotify)
  library(eoffice)
  e = read.csv("~/MamDC/Result/ENCODE/Savedata/Mouse_tissues_XAratio.txt",
               sep = "\t",header = TRUE,stringsAsFactors = F)
  e$sample %<>% gsub("GenitalAdiposeTissue","GAT",.) %>% 
    gsub("SubcutaneousAdiposeTissue","SCAT",.) %>% gsub("WholeBrain","Brain",.) 
    
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
  topptx(p,"~/MamDC/Result/ENCODE/Picture/Expressionratio.alltissues.GTEx.pptx",
         width = 7, height = 7)
}
