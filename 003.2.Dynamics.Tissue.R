{ # generate fpkm
  rm(list = ls())
  wd = "/home/qians/MamDC/Data/Seqdata/WangKuster2019MSBRNAseq"
  Species = "Human"
  a <- fread(file.path(wd,"counts.tsv")) %>% as.data.frame()
  a[1:3,1:3]
  a <- dplyr::filter(a,!grepl("_PAR_Y",V1))
  a$V1 <- lapply(a$V1,function(x)strsplit(x,".",fixed=T)[[1]][1]) %>% unlist()
  b <- read.csv(file.path("~/MamDC/Data/Ref",Species,paste0("gene",Species,".bed")),
                header = FALSE,sep = "\t",stringsAsFactors = F)
  a$type <- b[match(a$V1,b$V4),5]
  row.names(a) <- a$V1
  a = a[,-1]
    
  #exon length
  length <- read.csv(file = file.path("~/Pseudo/Data/Ref",Species,paste0(Species,".nonreExonlength.bed")),
                     header = FALSE,sep = "\t",stringsAsFactors = FALSE)
  #nonredundant gene length
  gene.length <- tapply(length$V5, length$V4, sum) %>% as.data.frame() 
  
  a$length = gene.length[match(row.names(a),row.names(gene.length)),1]
  a = na.omit(a)
  
  #allgene fpkm
  totalcounts <- colSums(a[,-c(ncol(a)-1,ncol(a))])
  rpkm <- t(do.call(rbind, lapply(1:length(totalcounts), function(i){
    10^9*a[,i]/a$length/totalcounts[i]}))) %>% as.data.frame()
  row.names(rpkm) = row.names(a)
  colnames(rpkm) = colnames(a)[1:(ncol(a)-2)]
  write.table(rpkm,file = file.path(wd,paste0(Species,".allgene.fpkm.txt")),
              sep = "\t",quote = F,row.names = T,col.names = T)

}

{ #generate ratio
  
  rm(list = ls())
  wd = "/home/qians/MamDC/Data/Seqdata/WangKuster2019MSBRNAseq"
  Species ="Human"
  a <- read.csv(file.path(wd,paste0(Species,".allgene.fpkm.txt")),header = T,sep = "\t",stringsAsFactors = FALSE)
  sample = colnames(a) %>% lapply(.,function(x)strsplit(x,split = "_rep")[[1]][1]) %>% unlist() #From Brain_s0dpb_Female_rep1 → Brain_s0dpb_Female
  mbg = t(apply(a, 1, function(x)tapply(x, sample, median))) %>% as.data.frame()
  
  b <- read.csv(file.path("~/MamDC/Data/Ref",Species,paste0("gene",Species,".bed")),
                header = FALSE,sep = "\t",stringsAsFactors = F)
  b$V5 %<>% gsub(" ","",.)
  mbg[,c("chr","type")] <- b[match(row.names(mbg),b$V4),c(1,5)]
  mbg <- dplyr::filter(mbg,type == "protein_coding" & chr %in% c(1:22,"X")) %>% dplyr::select(.,-type)
  dim(mbg)
  mbg[mbg$chr!="X","chr"]="A"
  
  mbg %<>% pivot_longer(cols=1:(ncol(mbg)-1), names_to= "samples", values_to = "expression")
  mbg[mbg$chr!="X","chr"]="A"
  
  e <- matrix(0,nrow = length(unique(sort(mbg$samples))),ncol = 4) %>% as.data.frame()
  colnames(e) <- c("sample","mean","CI5","CI95")
  e$sample <- unique(sort(mbg$samples))
  cutoff = 1
  for (T in unique(sort(mbg$samples))) {
    b <- dplyr::filter(mbg,samples %in% T & expression > cutoff) %>% as.data.frame()
    b$expression <- b$expression/mean(b[b$chr=="A","expression"]) #Allexpress/Aexpress
    d <- list()
    for (j in 1:1000) {
      d[[j]] <- sample(b[b$chr=="X","expression"],size = floor(nrow(b[b$chr=="X",])/2), replace = TRUE)
      boot.mean <- unlist(lapply(d, mean))
    }
    quantile(boot.mean, probs = c(0.05, 0.95))
    e[e$sample==T,2:4] <- c(mean(b[b$chr=="X","expression"]),quantile(boot.mean, probs = c(0.05, 0.95))[1],quantile(boot.mean, probs = c(0.05, 0.95))[2]) %>% as.numeric()
  }
  write.table(e,file = "~/MamDC/Result/WangKuster2019MSBRNAseq/Human/Savedata/Human_tissues_XAratio.txt",sep = "\t",row.names = FALSE, col.names = TRUE,quote = FALSE)
}

{ #X:A ratio & tissue specific genes
  #horizontal
  library(Hmisc)
  library(ggplotify)
  library(eoffice)
  e = read.csv("~/MamDC/Result/WangKuster2019MSBRNAseq/Human/Savedata/Human_tissues_XAratio.txt",
               sep = "\t",header = TRUE,stringsAsFactors = F)
  e$sample %<>% Hmisc::capitalize(.) %>% gsub("_"," ",.) %>% gsub("."," ",.,fixed = T)
  e = e[order(e$mean,decreasing = TRUE),]
  e$sample %<>% factor(.,levels = .)
  p = ggplot(e,aes(sample,mean))+geom_point(size=3.0)+geom_errorbar(aes(x = sample, ymax=CI5, ymin=CI95,width =0.3))+
    ylab("X:A ratio")+theme_bw()+theme(axis.title.x=element_text(size=18),axis.title.y=element_blank(),
                                       axis.text.x = element_text(size=16),axis.text.y = element_text(size=14))+
    scale_y_continuous(breaks = c(0,0.5,1,2),labels = c(0,0.5,1,2))+
    geom_hline(aes(yintercept=1),color="black",linetype="dashed")+
    geom_hline(aes(yintercept=0.5),color="black",linetype="dashed")+coord_flip()
  ggsave(p,filename = "~/MamDC/Result/WangKuster2019MSBRNAseq/Human/Picture/Expressionratio.alltissues.pdf",
         device = "pdf", width = 5.8, height = 7)
  topptx(p,filename = "~/MamDC/Result/WangKuster2019MSBRNAseq/Human/Picture/Expressionratio.alltissues.pptx",
         width = 5.8, height = 7)
}
  
{ #X:A ratio & tissue specific genes
  #vertical
  library(Hmisc)
  library(ggplotify)
  library(eoffice)
  e = read.csv("~/MamDC/Result/WangKuster2019MSBRNAseq/Human/Savedata/Human_tissues_XAratio.txt",
               sep = "\t",header = TRUE,stringsAsFactors = F)
  e$sample %<>% Hmisc::capitalize(.) %>% gsub("_"," ",.) %>% gsub("."," ",.,fixed = T)
  e = e[order(e$mean,decreasing = FALSE),]
  e$sample %<>% factor(.,levels = .)
  p= ggplot(e,aes(sample,mean))+geom_point(size=3.0)+geom_errorbar(aes(x = sample, ymax=CI5, ymin=CI95,width =0.3))+
    ylab("X:AA ratio")+theme_bw()+
    theme(axis.title.x=element_blank(),axis.title.y=element_text(size=18),
          axis.text.x = element_text(size=16,angle = 90,hjust = 1,vjust = 0.5),axis.text.y = element_text(size=14))+
    scale_y_continuous(breaks = c(0,0.5,1,2),labels = c(0,0.5,1,2))+
    geom_hline(aes(yintercept=1),color="black",linetype="dashed")+
    geom_hline(aes(yintercept=0.5),color="black",linetype="dashed")
  ggsave(p,filename = "~/MamDC/Result/WangKuster2019MSBRNAseq/Human/Picture/Expressionratio.alltissues.vertical.pdf",
         device = "pdf", width = 7, height = 5.8)
  topptx(p,filename = "~/MamDC/Result/WangKuster2019MSBRNAseq/Human/Picture/Expressionratio.alltissues.vertical.pptx",
         width = 7, height = 5.8)
}  

  ####### tissue specific gene
{
  rm(list = ls())
  wd = "/home/qians/MamDC/Data/Seqdata/WangKuster2019MSBRNAseq"
  Species ="Human"
  a <- read.csv(file.path(wd,paste0(Species,".allgene.fpkm.txt")),header = T,sep = "\t",stringsAsFactors = FALSE)
  sample = colnames(a) %>% lapply(.,function(x)strsplit(x,split = "_rep")[[1]][1]) %>% unlist() #From Brain_s0dpb_Female_rep1 → Brain_s0dpb_Female
  mbg = t(apply(a, 1, function(x)tapply(x, sample, median))) %>% as.data.frame()
  b <- read.csv(file.path("~/MamDC/Data/Ref",Species,paste0("gene",Species,".bed")),
                header = FALSE,sep = "\t",stringsAsFactors = F)
  b$V5 %<>% gsub(" ","",.)
  mbg[,c("chr","type")] <- b[match(row.names(mbg),b$V4),c(1,5)]
  mbg <- dplyr::filter(mbg,type == "protein_coding" & chr %in% c(1:22,"X")) %>% dplyr::select(.,-type)
  dim(mbg)
  mbg[mbg$chr!="X","chr"]="A"
  mbg %<>% dplyr::select(.,-"chr")
  
  tissue_gene <- as.data.frame(t(mbg))
  judge_matrix <- matrix(nrow = length(colnames(mbg)),ncol = length(rownames(mbg)))
  
  vec_quan <- c()
  for (tissue_num in seq(1,nrow(tissue_gene),1)) {
    quan=quantile(as.numeric(tissue_gene[tissue_num,]),0.6)
    vec_quan <- append(vec_quan,as.numeric(quan))
  }
  
  for (col_i in seq(1,ncol(tissue_gene),1)) {
    is_top5=rank(tissue_gene[,col_i])>= nrow(tissue_gene)-(5-1)
    for (row_j in seq(1,ncol(mbg),1)) {
      if (is_top5[row_j]) {
        quan=vec_quan[row_j]
        if (tissue_gene[row_j,col_i] > quan) {
          judge_matrix[row_j,col_i]=1
        } else {
          judge_matrix[row_j,col_i]=0
        }
      } else {
        judge_matrix[row_j,col_i]=0
      }
    }
  }
  
  judge_matrix <- as.data.frame(judge_matrix)
  colnames(judge_matrix) <- colnames(tissue_gene)
  rownames(judge_matrix) <- rownames(tissue_gene)
  #转置回来
  gene_ts <- as.data.frame(t(judge_matrix))
  gene_ts <- gene_ts[!(rowSums(gene_ts)==0),]
  
  gene_ts$chr <- b[match(rownames(gene_ts),b$V4),1]
  
  #c <- as.numeric(colSums(gene_ts[gene_ts$chr=="X",-ncol(gene_ts)]))/as.numeric(colSums(gene_ts[,-ncol(gene_ts)])) %>% as.data.frame()
  
  c <- matrix(1:(ncol(gene_ts)-1)) %>% as.data.frame
  c$V1 <- colSums(gene_ts[gene_ts$chr=="X",-ncol(gene_ts)]) #x linked
  c$V2 <- colSums(gene_ts[,-ncol(gene_ts)]) #all
  c$tissue <- colnames(gene_ts)[-ncol(gene_ts)]
  p = ggplot(c,aes(tissue,V1))+geom_point()+coord_flip()+theme_bw()+ylab("# X-linked tissue-specific genes")+
    theme(axis.title.x=element_text(size=18),axis.title.y=element_blank(),
          axis.text.x = element_text(size=16),axis.text.y = element_text(size=14))
  ggsave(p, filename = "~/MamDC/Result/WangKuster2019MSBRNAseq/Human/Picture/NumofXlinkedTSgenes.alltissues.pdf",
         device = "pdf", width = 6.5, height = 7)
  
  c$type <- "N"
  for (i in 1:nrow(c)) {
    if (c[i,"tissue"] %in% c("liver","pancreas","saliva.secreting_gland","skeletal_muscle_tissue")) {
      c[i,"type"]= "No DC"
    } 
    else {
      c[i,"type"]= "DC"
    }
  }
  
  c$tissue %<>% gsub("_"," ",.,fixed = TRUE) %>% capitalize() %>% gsub("."," ",.,fixed = T)
  
  e = read.csv("~/MamDC/Result/WangKuster2019MSBRNAseq/Human/Savedata/Human_tissues_XAratio.txt",
               sep = "\t",header = TRUE,stringsAsFactors = F)
  e$sample %<>% Hmisc::capitalize(.) %>% gsub("_"," ",.) %>% gsub("."," ",.,fixed = T)
  e = e[order(e$mean,decreasing = TRUE),]
  e$sample %<>% factor(.,levels = .)
  
  c$mean = e[match(c$tissue,e$sample),"mean"]
  cor = cor.test(c$mean,c$V1)
  p = ggplot(c,aes(mean,V1))+geom_point()+geom_smooth(method = "lm",se = F)+
    theme_classic()+xlab("X:A ratio")+ylab("# X-linked tissue-specific genes")+
    theme(axis.title.x=element_text(size=16),axis.title.y=element_text(size=16),
          axis.text.x = element_text(size=14),axis.text.y = element_text(size=14))+
    annotate(geom="text", x=0.4, y=150, size=4,
             label=paste0("R=",round(cor$estimate,2),", P value=",format(cor$p.value,scientific = TRUE,digits = 2)))
  ggsave(p,filename = "~/MamDC/Result/WangKuster2019MSBRNAseq/Human/Picture/NumofXlinkedTSgenes.XAratio.pdf",
         device = "pdf", width = 4.8, height = 4)
  topptx(p,filename = "~/MamDC/Result/WangKuster2019MSBRNAseq/Human/Picture/NumofXlinkedTSgenes.XAratio.pptx",
         width = 4.8, height = 4)
  
  c$propor = c$V1/c$V2*100 #proportion of Xlinked TS genes
  cor = cor.test(c$mean,c$propor)
  p = ggplot(c,aes(mean,propor))+geom_point()+geom_smooth(method = "lm",se = F)+
    theme_classic()+xlab("X:A ratio")+ylab("(%) X-linked tissue-specific genes")+
    theme(axis.title.x=element_text(size=16),axis.title.y=element_text(size=16),
          axis.text.x = element_text(size=14),axis.text.y = element_text(size=14))+
    annotate(geom="text", x=0.4, y=6, size=4,
             label=paste0("R=",round(cor$estimate,2),", P value=",format(cor$p.value,scientific = TRUE,digits = 2)))
  #ggsave(p,filename = "~/MamDC/Result/WangKuster2019MSBRNAseq/Human/Picture/NumofXlinkedTSgenes.XAratio.pdf",
   #      device = "pdf", width = 4.8, height = 4)
  topptx(p,filename = "~/MamDC/Result/WangKuster2019MSBRNAseq/Human/Picture/ProporofXlinkedTSgenes.XAratio.pptx",
         width = 4.8, height = 4)
  ggsave(p,filename = "~/MamDC/Result/WangKuster2019MSBRNAseq/Human/Picture/ProporofXlinkedTSgenes.XAratio.pdf",
         width = 4.8, height = 4)
  
  c$type <- factor(c$type,levels = c("No DC","DC"))
  wilcox.test(c[c$type=="DC",1],c[c$type!="DC",1],alternative = "greater")
  p = ggplot(c,aes(type,V1,fill=type))+geom_boxplot()+geom_point(position="jitter",pch=16,cex=1)+
    ylab("# X-linked tissue-specific genes")+theme_classic()+
    theme(axis.title.x=element_blank(),axis.title.y=element_text(size=16),
          axis.text.x = element_text(size=14),axis.text.y = element_text(size=14))+
    scale_fill_manual(values = c("#E39E3E","#C7533B"))+
    guides(fill=FALSE)+
    geom_segment(aes(x=1, y=159, xend=2, yend=159))+
    annotate("text", x=1.5, y=161, label="**",size=5)
  ggsave(p, filename = "~/MamDC/Result/WangKuster2019MSBRNAseq/Human/Picture/NumofXlinkedTSgenes.DCtype.pdf",
         device = "pdf", width = 4, height = 3.8)
  topptx(p,filename = "~/MamDC/Result/WangKuster2019MSBRNAseq/Human/Picture/NumofXlinkedTSgenes.DCtype.pptx",
         width = 4, height = 3.8)
}