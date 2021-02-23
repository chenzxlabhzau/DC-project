{ #number
  rm(list = ls())
  for (S in c("Human","Mouse")) {
    
    b <- read.csv(file.path("/home/qians/Pseudo/Data/Ref",S,paste0("gene",S,".bed")),
                  header = FALSE,sep = "\t",stringsAsFactors = F)
    c = read.csv(file.path("/home/qians/Pseudo/Result/GTRDforTF",S,paste0(S,".promoter.number.bed")),header = FALSE,sep = "\t",stringsAsFactors = F)
    c$type= b[match(c$V4,b$V4),5]
    c %<>% dplyr::filter(.,type == " protein_coding" & V1 %in% paste0("chr",c(1:100,"X")))
    c[c$V1!="chrX","V1"]="A"
    c[c$V1=="chrX","V1"]="X"
    if (S=="Human") {
      age = read.csv("~/MamDC/Data/Seqdata/Genorigin/Human/Homo_sapiens.csv",header = TRUE,sep = ",",stringsAsFactors = FALSE)
    }
    c = c[,c(1,ncol(c)-1)]
    colnames(c) = c("chr","number")
    wilcox.test(c[c$chr=="X",2],c[c$chr!="X",2])
    if (S=="Human") {
      p = ggplot(c,aes(chr,log10(number+1),fill=chr))+geom_boxplot(outlier.alpha = 0,notch = TRUE)+theme_bw()+
        ylab("TF-binding sites (log10)")+
        theme(axis.title.x=element_blank(),legend.position='none',
              axis.text.x = element_text(size=12),
              axis.title.y=element_text(size=14),axis.text.y = element_text(size=12))+
        scale_fill_manual(values = c("#035782","#9c0b1f"))+
        annotate("text", x=1.5, y=3.8, label="***",size=8)+
        coord_cartesian(ylim = c(0.5,4))
    }
    if (S=="Mouse") {
      p = ggplot(c,aes(chr,log10(number+1),fill=chr))+geom_boxplot(outlier.alpha = 0,notch = TRUE)+theme_bw()+
        ylab("TF-binding sites (log10)")+
        theme(axis.title.x=element_blank(),legend.position='none',
              axis.text.x = element_text(size=12),
              axis.title.y=element_text(size=14),axis.text.y = element_text(size=12))+
        scale_fill_manual(values = c("#035782","#9c0b1f"))+
        annotate("text", x=1.5, y=3.5, label="***",size=8)
    }
    topptx(p,filename = file.path("~/MamDC/Result/GTRDforTF",S,"Picture/TF.promoter.numberAX.pptx"),
           width = 3,height = 3)
  }
}

{ #type
  rm(list = ls())
  for (S in c("Human","Mouse")) {
    
    b <- read.csv(file.path("/home/qians/Pseudo/Data/Ref",S,paste0("gene",S,".bed")),
                  header = FALSE,sep = "\t",stringsAsFactors = F)
    c =fread(file.path("~/Pseudo/Result/GTRDforTF",S,paste0(S,".promoter.type3.bed")),header = FALSE) %>% as.data.frame()
    c[,c("type","chr")]= b[match(c$V2,b$V4),c(5,1)]
    
    c %<>% dplyr::filter(.,type == " protein_coding" & chr %in% c(1:100,"X"))
    c[c$chr!="X","chr"]="A"
    c[c$chr=="X","chr"]="X"
    wilcox.test(c[c$chr=="X",1],c[c$chr!="X",1])
    if (S=="Human") {
      p = ggplot(c,aes(chr,log10(V1+1),fill=chr))+geom_boxplot(outlier.alpha = 0,notch = TRUE)+theme_bw()+
      ylab("TF-binding diversity (log10)")+
      theme(axis.title.x=element_blank(),legend.position='none',
            axis.text.x = element_text(size=12),
            axis.title.y=element_text(size=14),axis.text.y = element_text(size=12))+
      scale_fill_manual(values = c("#035782","#9c0b1f"))+
      annotate("text", x=1.5, y=3.1, label="***",size=8)+
      coord_cartesian(ylim = c(0.5,3.2))
    }
    if (S=="Mouse") {
      p = ggplot(c,aes(chr,log10(V1+1),fill=chr))+geom_boxplot(outlier.alpha = 0,notch = TRUE)+theme_bw()+
        ylab("TF-binding diversity (log10)")+
        theme(axis.title.x=element_blank(),legend.position='none',
              axis.text.x = element_text(size=12),
              axis.title.y=element_text(size=14),axis.text.y = element_text(size=12))+
        scale_fill_manual(values = c("#035782","#9c0b1f"))+
        annotate("text", x=1.5, y=2.9, label="***",size=8)+
        coord_cartesian(ylim = c(0,3))
    }
    topptx(p,filename = file.path("~/MamDC/Result/GTRDforTF",S,"Picture/TF.promoter.typeAX.pptx"),
           width = 3,height = 3)
    
   
  }
}

{ #venn TF type
  rm(list = ls())
  a = fread("~/Pseudo/Result/GTRDforTF/Mouse/Mouse.promoter.type2.bed",header = FALSE) %>% as.data.frame()
  a$gene = lapply(a$V1,function(x)strsplit(as.character(x),"NAME",fixed=T)[[1]][1]) %>% unlist()
  a$tf = lapply(a$V1,function(x)strsplit(as.character(x),"NAME",fixed=T)[[1]][2]) %>% unlist()
  S="Mouse"
  b <- read.csv(file.path("/home/qians/Pseudo/Data/Ref",S,paste0("gene",S,".bed")),
                header = FALSE,sep = "\t",stringsAsFactors = F)
  a$chr = b[match(a$gene,b$V4),1]
  a1 = a %>% dplyr::filter(.,chr %in% c(1:100))
  a2 = a %>% dplyr::filter(.,chr == "X")
  t1 = table(a1$tf) %>% as.data.frame()
  t1 = t1[order(t1$Freq,decreasing = TRUE),]
  t2 = table(a2$tf) %>% as.data.frame()
  t2 = t2[order(t2$Freq,decreasing = TRUE),]
  write.table(t1,file = "~/Pseudo/Result/GTRDforTF/Mouse/Mouse.promoter.TFtype.A.txt",quote = FALSE,col.names = FALSE,row.names = FALSE)
  write.table(t2,file = "~/Pseudo/Result/GTRDforTF/Mouse/Mouse.promoter.TFtype.X.txt",quote = FALSE,col.names = FALSE,row.names = FALSE)
  
  if (FALSE) {
    library(VennDiagram)
    setwd(file.path("~/MamDC/Result/GTRDforTF",S,"Picture"))
    pdf(file.path("~/MamDC/Result/GTRDforTF",S,"Picture",paste0(S,".TFtype.venn.XA.pdf")))
    grid.draw(venn.diagram(list("X" = head(as.character(t2$Var1),nrow(t2)*0.2),
                                "Autosome" = head(as.character(t1$Var1),nrow(t1)*0.2)),
                           NULL,
                           imagetype ="svg",
                           fill = c("#d95f0d","#fc9272"),
                           cat.pos=c(0,3),
                           cat.cex=c(2,2),
                           cex = 2,
                           width = 4,height = 4))
    dev.off()
  }
}


############ correlate with age
{ #number
  rm(list = ls())
  for (S in c("Human","Mouse")) {
    b <- read.csv(file.path("/home/qians/Pseudo/Data/Ref",S,paste0("gene",S,".bed")),
                  header = FALSE,sep = "\t",stringsAsFactors = F)
    c = read.csv(file.path("/home/qians/Pseudo/Result/GTRDforTF",S,paste0(S,".promoter.number.bed")),header = FALSE,sep = "\t",stringsAsFactors = F)
    c$type= b[match(c$V4,b$V4),5]
    c %<>% dplyr::filter(.,type == " protein_coding" & V1 %in% paste0("chr",c(1:100,"X")))
    c[c$V1!="chrX","V1"]="A"
    c[c$V1=="chrX","V1"]="X"
    if (S=="Human") {
      age = read.csv("~/MamDC/Data/Seqdata/Genorigin/Human/Homo_sapiens.csv",header = TRUE,sep = ",",stringsAsFactors = FALSE)
    }
    if (S=="Mouse") {
      age = read.csv("~/MamDC/Data/Seqdata/Genorigin/Mouse/Mus_musculus.csv",header = TRUE,sep = ",",stringsAsFactors = FALSE)
      c$V5 = c$V6
    }
    age$gene_age %<>% gsub(">","",.) %>% as.numeric()
    c$age = age[match(c$V4,age$ensembl_id),2] 
    c$group = cut(c$age,breaks = c(0,332,645,5000),labels = c("Young","Middle","Old"))
    
    my_comparisons <- list(c("Young","Middle"),c("Middle","Old"),c("Young","Old"))
    p = ggplot(c,aes(group,log10(V5+1),color=group))+geom_boxplot(notch = TRUE,outlier.alpha = 0)+
      stat_compare_means(comparisons = my_comparisons)+theme_classic()+
      ylab("Number of TF-binding sites (log10)")+
      theme(axis.title.x=element_blank(),
            axis.text.x = element_text(size=10),
            axis.title.y=element_text(size=12),
            axis.text.y = element_text(size=10),
            strip.background  = element_blank(),
            strip.text = element_text(size=12))+
      guides(color=FALSE)+
      scale_color_manual(values = c("#D6B6BB","#CB7582","#9c0b1f"))+
      facet_wrap(.~V1)
    
    topptx(p,filename = file.path("~/MamDC/Result/GTRDforTF",S,"Picture/TF.promoter.numberAX.age.pptx"),
           width = 3,height = 3)
  }
}

{ #type
  rm(list = ls())
  for (S in c("Human","Mouse")) {
    
    b <- read.csv(file.path("/home/qians/Pseudo/Data/Ref",S,paste0("gene",S,".bed")),
                  header = FALSE,sep = "\t",stringsAsFactors = F)
    c =fread(file.path("~/Pseudo/Result/GTRDforTF",S,paste0(S,".promoter.type3.bed")),header = FALSE) %>% as.data.frame()
    c[,c("type","chr")]= b[match(c$V2,b$V4),c(5,1)]
    
    c %<>% dplyr::filter(.,type == " protein_coding" & chr %in% c(1:100,"X"))
    c[c$chr!="X","chr"]="A"
    c[c$chr=="X","chr"]="X"
    if (S=="Human") {
      age = read.csv("~/MamDC/Data/Seqdata/Genorigin/Human/Homo_sapiens.csv",header = TRUE,sep = ",",stringsAsFactors = FALSE)
    }
    if (S=="Mouse") {
      age = read.csv("~/MamDC/Data/Seqdata/Genorigin/Mouse/Mus_musculus.csv",header = TRUE,sep = ",",stringsAsFactors = FALSE)
      
    }
    age$gene_age %<>% gsub(">","",.) %>% as.numeric()
    c$age = age[match(c$V2,age$ensembl_id),2] 
    c$group = cut(c$age,breaks = c(0,332,645,5000),labels = c("Young","Middle","Old"))
    
    my_comparisons <- list(c("Young","Middle"),c("Middle","Old"),c("Young","Old"))
    p = ggplot(c,aes(group,log10(V1+1),color=group))+geom_boxplot(notch = TRUE,outlier.alpha = 0)+
      stat_compare_means(comparisons = my_comparisons)+theme_classic()+
      ylab("TF-binding Diversity (log10)")+
      theme(axis.title.x=element_blank(),
            axis.text.x = element_text(size=10),
            axis.title.y=element_text(size=12),
            axis.text.y = element_text(size=10),
            strip.background  = element_blank(),
            strip.text = element_text(size=12))+
      guides(color=FALSE)+
      scale_color_manual(values = c("#D6B6BB","#CB7582","#9c0b1f"))+
      facet_wrap(.~chr)
    topptx(p,filename = file.path("~/MamDC/Result/GTRDforTF",S,"Picture/TF.promoter.typeAX.age.pptx"),
           width = 3,height = 3)
  }
}
