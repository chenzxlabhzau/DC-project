#
rm(list = ls())
gc()
rm(list = ls())
#ref.species = "Platypus" #"Opossum" #"Chicken" 
#datatype = "ribo"
for (ref.species in c("Opossum","Chicken","Platypus")) {
  for (datatype in c("ribo","rna")) {
    focal.species = "Human"
    wd <- "~/MamDC/Data/Seqdata/WangKaessmann2020NatureRNARiboseq"
    base.object = ls()
    #focal species
    focal.dir = file.path(wd,tolower(focal.species))
    focal.files <- grep(paste0(datatype,"_[1-9].txt$"),list.files(focal.dir),value=TRUE)  
    focal.filePath <- sapply(focal.files, function(x){paste(focal.dir,x,sep='/')})   
    focal.data <- lapply(focal.filePath, function(x){ fread(x)})  
    a <- focal.data[1] %>% as.data.frame()
    a = a[,c(1,ncol(a))]
    colnames(a) <- c("Gene",strsplit(colnames(a)[1],split=".",fix=TRUE)[[1]][1]) # geneID, sample name
    for (i in 2:length(focal.data)) {
      b <- focal.data[i] %>% as.data.frame()
      b = b[,c(1,ncol(b))]
      colnames(b) <- c("Gene",strsplit(colnames(b)[1],split=".",fix=TRUE)[[1]][1])
      a <- merge(a,b,by="Gene")
    }
    row.names(a) = a$Gene
    a[is.na(a)]=0
    group = lapply(colnames(a)[-1], function(x)strsplit(x,split = "_",fixed = TRUE)[[1]][2]) %>% unlist()
    focal.mbg = t(apply(a[,-1], 1, function(x)tapply(x, group, mean))) %>%as.data.frame()
    with(focal.mbg,cor.test(brain,testis,method = "spearman"))
    rm(list = setdiff(ls(),c(base.object,"focal.mbg")))  
    
    #ref.species
    ref.dir = file.path(wd,tolower(ref.species))
    ref.files <- grep(paste0(datatype,"_[1-9].txt$"),list.files(ref.dir),value=TRUE)  
    ref.filePath <- sapply(ref.files, function(x){paste(ref.dir,x,sep='/')})   
    ref.data <- lapply(ref.filePath, function(x){ fread(x)})  
    a <- ref.data[1] %>% as.data.frame()
    a = a[,c(1,ncol(a))]
    colnames(a) <- c("Gene",strsplit(colnames(a)[1],split=".",fix=TRUE)[[1]][1]) # geneID, sample name
    for (i in 2:length(ref.data)) {
      b <- ref.data[i] %>% as.data.frame()
      b = b[,c(1,ncol(b))]
      colnames(b) <- c("Gene",strsplit(colnames(b)[1],split=".",fix=TRUE)[[1]][1])
      a <- merge(a,b,by="Gene")
    }
    row.names(a) = a$Gene
    a[is.na(a)]=0
    group = lapply(colnames(a)[-1], function(x)strsplit(x,split = "_",fixed = TRUE)[[1]][2]) %>% unlist()
    ref.mbg = t(apply(a[,-1], 1, function(x)tapply(x, group, mean))) %>%as.data.frame()
    with(ref.mbg,cor.test(brain,testis,method = "spearman"))
    rm(list = setdiff(ls(),grep("mbg|species|datatype",ls(),value = TRUE)))
    
    #build homo genes
    homo = read.csv(file.path("~/MamDC/Data/Ref",focal.species,
                              paste(focal.species,ref.species,"ortholog.csv",sep = ".")))
    colnames(homo) = c("gene1","chr1","gene2","chr2","homologytype")
    focal.gene = read.csv(file.path("/home/qians/MamDC/Data/Ref",focal.species,paste0(focal.species,".biomart.gene102.csv")))
    ref.gene = read.csv(file.path("/home/qians/MamDC/Data/Ref",ref.species,paste0(ref.species,".biomart.gene102.csv")))
    
    homo$type1 = focal.gene[match(homo$gene1,focal.gene$GeneID),"Genetype"]
    homo$type2 = ref.gene[match(homo$gene2,ref.gene$GeneID),"Genetype"]
    if (ref.species=="Chicken") {
      homo %<>% dplyr::filter(chr1 %in% c(1:100,"X") & 
                                chr2 %in% c(1:100) &
                                homologytype %in% "ortholog_one2one" &
                                type1 == "protein_coding" &
                                type2 == "protein_coding" )
      homo$chrtype = "Others"
      homo[homo$chr1=="X" & homo$chr2 %in% c(1,4),"chrtype"]="XX"
      homo[homo$chr1!="X","chrtype"]="AA"
    }
    if (ref.species=="Opossum") {
      homo %<>% dplyr::filter(chr1 %in% c(1:100,"X") & 
                                chr2 %in% c(1:100,"X") &
                                homologytype %in% "ortholog_one2one" &
                                type1 == "protein_coding" &
                                type2 == "protein_coding" )
      homo$chrtype = "Others"
      homo[homo$chr1=="X" & homo$chr2 %in% c(4,7,"X"),"chrtype"]="XX"
      homo[homo$chr1!="X","chrtype"]="AA"
    }
    if (ref.species=="Platypus") {
      homo %<>% dplyr::filter(chr1 %in% c(1:100,"X") & 
                                chr2 %in% c(1:100) &
                                homologytype %in% "ortholog_one2one" &
                                type1 == "protein_coding" &
                                type2 == "protein_coding" )
      
      homo$chrtype = "Others"
      homo[homo$chr1=="X" & homo$chr2 %in% c(6,15,18),"chrtype"]="XX" #415 
      homo[homo$chr1!="X","chrtype"]="AA" #10740
    }
    
    homo %<>% dplyr::select(.,c("gene1","gene2","chrtype")) %>% dplyr::filter(.,chrtype %in% c("AA","XX"))
    
    #calculate pearson correlation
    library(ggpointdensity)
    library(viridis)
    for (Tissue in colnames(focal.mbg)) {
      homo$t1 = focal.mbg[match(homo$gene1,row.names(focal.mbg)),Tissue]
      homo$t2 = ref.mbg[match(homo$gene2,row.names(ref.mbg)),Tissue]
      cor1 = homo %>% na.omit() %>% dplyr::filter(.,chrtype == "AA") %>% with(.,cor.test(t1,t2,method = "pearson"))
      cor2 = homo %>% na.omit() %>% dplyr::filter(.,chrtype == "XX") %>% with(.,cor.test(t1,t2,method = "pearson"))
      
      cor = homo %>% na.omit() %>% with(.,cor.test(t1,t2,method = "pearson"))
      #homo %>% na.omit() %>% dplyr::filter(.,t1>=1 & t2 >= 1)%>% with(.,cor.test(t1,t2,method = "pearson"))
      
      p = homo %>% na.omit() %>% ggplot(aes(log(t1+1),log(t2+1)))+geom_pointdensity() +
        scale_color_viridis()+theme_classic()+
        xlab(paste(focal.species,Tissue,sep = " "))+
        ylab(paste(ref.species,Tissue,sep = " "))+
        theme(axis.title.x = element_text(size = 14),
              axis.title.y = element_text(size = 14),
              axis.text.x = element_text(size = 12),
              axis.text.y = element_text(size = 12))+
        #annotate(geom="text", x= 3, y=homo %>% na.omit() %>% with(.,log(quantile(t2)[5])),
        #         size=5, label=paste0("R=",round(cor$estimate,2)))+
        ggtitle(paste0("R=",round(cor$estimate,2)))+
        #theme(plot.title = element_text(vjust = - 8,hjust = 0.5))+
        facet_wrap(.~chrtype)+
        geom_text(aes(x,y,label =lab),
                  data = data.frame(x = 1:2,
                                    y = 8.5,
                                    lab = paste0("R=",c(round(cor1$estimate,2),round(cor2$estimate,2))),
                                    chrtype = factor(c("AA","XX")),levels = c("AA","XX")),
                  vjust = 1, size=6)
      
      ggsave(paste0("~/MamDC/Result/WangKaessmann2020NatureRNARiboseq/",
                    focal.species,"/Picture/",focal.species,".",ref.species,".",Tissue,".",datatype,".correlation.pdf"),
             p,device = "pdf",width = 7,height = 4)
      topptx(p,filename = paste0("~/MamDC/Result/WangKaessmann2020NatureRNARiboseq/",
                                 focal.species,"/Picture/",focal.species,".",ref.species,".",Tissue,".",datatype,".correlation.pptx"),
             width = 7,height = 4)
    }
    
    # X:XX A:AA Ratio
    e = data.frame()
    for (Tissue in colnames(focal.mbg)) {
      homo$t1 = focal.mbg[match(homo$gene1,row.names(focal.mbg)),Tissue]
      homo$t2 = ref.mbg[match(homo$gene2,row.names(ref.mbg)),Tissue]
      
      #if (FALSE) { #He et al method
      #  homo.x = homo %>% na.omit() %>% dplyr::filter(.,chrtype == "XX")
      #  homo.x$tissue = Tissue
      #  homo.x[homo.x$t1<0.01,"t1"]=0.01
      #  homo.x[homo.x$t2<0.01,"t2"]=0.01
      #  homo.auto = homo %>% na.omit() %>% dplyr::filter(.,chrtype == "AA")
      #  homo.auto$tissue = Tissue
      #  homo.auto[homo.auto$t1<0.01,"t1"]=0.01
      #  homo.auto[homo.auto$t2<0.01,"t2"]=0.01
      #}
      
      homo.x = homo %>% na.omit() %>% dplyr::filter(.,t1>=1 & t2 >= 1 & chrtype == "XX")
      homo.x$tissue = Tissue
      homo.auto = homo %>% na.omit() %>% dplyr::filter(.,t1>=1 & t2 >= 1 & chrtype == "AA")
      homo.auto$tissue = Tissue
      #scale factor
      scale.factor = median(homo.auto$t1/homo.auto$t2)
      homo.auto$t1 = homo.auto$t1 / scale.factor
      homo.x$t1 = homo.x$t1 / scale.factor
      e = rbind(e,homo.x,homo.auto)
    }
    e$ratio = e$t1/e$t2
    
    wilcox.test(e[e$tissue=="testis" & e$chrtype=="AA","t1"],
                e[e$tissue=="testis" & e$chrtype=="XX","t1"],
                alternative = "greater")
    tapply(e$ratio, paste(e$tissue,e$chrtype,sep = "_"), median)
    
    p = ggplot(e,aes(tissue,log2(ratio),fill=chrtype))+geom_boxplot(notch = T,outlier.colour = "white")+
      theme_classic()+
      coord_cartesian(ylim = c(-4,3))+
      scale_y_continuous(breaks = c(-2,-1,0,1,2),labels = c(-2,-1,0,1,2))+
      scale_x_discrete(breaks = unique(e$tissue),labels = unique(e$tissue) %>%capitalize())+
      #xlab(paste(focal.species,Tissue,sep = " "))+
      ylab(paste0("Log2(",focal.species,":",ref.species,")"))+
      theme(axis.title.x = element_blank(),
            axis.title.y = element_text(size = 14),
            axis.text.x = element_text(size = 14,angle = 45,hjust = 1),
            axis.text.y = element_text(size = 12))+
      scale_fill_manual(values = c("#035782","#9c0b1f"))+
      geom_hline(yintercept = 0,color="black",linetype="dashed")+
      geom_hline(yintercept = -1,color="black",linetype="dashed")+
      guides(fill=guide_legend(title=NULL))+theme(legend.position=c(.2,.92),legend.direction = "horizontal")
    ggsave(paste0("/home/qians/MamDC/Result/WangKaessmann2020NatureRNARiboseq/",
                  focal.species,"/Picture/006.",focal.species,".",ref.species,".",datatype,".log2ratio.pdf"),
           p,device = "pdf",width = 4, height = 4)
    topptx(p,filename = paste0("/home/qians/MamDC/Result/WangKaessmann2020NatureRNARiboseq/",
                               focal.species,"/Picture/006.",focal.species,".",ref.species,".",datatype,".log2ratio.pptx"),
           width = 4, height = 4)
    
    if (FALSE) { #XX:AA
      for (i in seq(2,6,2)) {
        print(apply(e[,c("t1","t2")], 2, function(x)tapply(x, paste(e$tissue,e$chrtype,sep = "_"), mean))[i,]/
                apply(e[,c("t1","t2")], 2, function(x)tapply(x, paste(e$tissue,e$chrtype,sep = "_"), mean))[i-1,])
      }
      #human XX:AA
      tmp1 = matrix(0,nrow = length(unique(e$tissue)),ncol = 5) %>% as.data.frame()
      colnames(tmp1) <- c("sample","mean","CI0.5","CI0.95","species")
      tmp1$sample = unique(e$tissue)
      tmp1$species ="Human"
      for (i in 1:3) {
        c <- e[e$tissue==unique(e$tissue)[i],c("chrtype","t1")]
        colnames(c) <- c("chr","express")
        #c <- c[c$express>cutoff,] # has done in e
        c$express <- c$express/mean(c[c$chr=="AA","express"]) #Allexpress/Aexpress
        d <- list()
        for (j in 1:1000) {
          d[[j]] <- sample(c[c$chr=="XX","express"],size = floor(nrow(c[c$chr=="XX",])/2), replace = FALSE)
          boot.mean <- unlist(lapply(d, median)) #trim?
        }
        quantile(boot.mean, probs = c(0.05, 0.95))
        tmp1[i,2:4] = c(median(c[c$chr=="XX","express"]),
                        quantile(boot.mean, probs = c(0.05, 0.95))[1],
                        quantile(boot.mean, probs = c(0.05, 0.95))[2]) %>% as.numeric()
      }
      #chicken XX:AA
      tmp2 = matrix(0,nrow = length(unique(e$tissue)),ncol = 5) %>% as.data.frame()
      colnames(tmp2) <- c("sample","mean","CI0.5","CI0.95","species")
      tmp2$sample = unique(e$tissue)
      tmp2$species ="Chicken"
      for (i in 1:3) {
        c <- e[e$tissue==unique(e$tissue)[i],c("chrtype","t2")]
        colnames(c) <- c("chr","express")
        #c <- c[c$express>cutoff,] # has done in e
        c$express <- c$express/median(c[c$chr=="AA","express"]) #Allexpress/Aexpress
        d <- list()
        for (j in 1:1000) {
          d[[j]] <- sample(c[c$chr=="XX","express"],size = floor(nrow(c[c$chr=="XX",])/2), replace = FALSE)
          boot.mean <- unlist(lapply(d, median)) #trim?
        }
        quantile(boot.mean, probs = c(0.05, 0.95))
        tmp2[i,2:4] = c(median(c[c$chr=="XX","express"]),
                        quantile(boot.mean, probs = c(0.05, 0.95))[1],
                        quantile(boot.mean, probs = c(0.05, 0.95))[2]) %>% as.numeric()
      }
      tmp = rbind(tmp1,tmp2)
      tmp$sample %<>% Hmisc::capitalize(.)
      ggplot(tmp,aes(sample,mean,shape=species))+geom_point(position = position_dodge(0.5),size=3.0)+
        geom_errorbar(aes(x = sample, ymax=CI0.5, ymin=CI0.95,width =0.3),position = position_dodge(0.5))+
        ylab("X:A ratio")+theme_classic()+
        theme(axis.title.y=element_text(size=18),axis.title.x=element_blank(),
              axis.text.y = element_text(size=16),axis.text.x = element_text(size=12,angle = 45,hjust = 1,vjust = 1))+
        #scale_y_continuous(breaks = c(0,0.5,1,2),labels = c(0,0.5,1,2))+
        geom_hline(aes(yintercept=1),color="black",linetype="dashed")+
        geom_hline(aes(yintercept=0.5),color="black",linetype="dashed")+
        scale_y_continuous(expand = c(0, 0))+ coord_cartesian(ylim = c(0,max(tmp$CI0.95)*1.1))
    }
    }
  rm(list = ls())
}
