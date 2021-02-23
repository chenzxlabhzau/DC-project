rm(list = ls())
gc()
rm(list = ls())
#mouse data
{ #generate fpkm
  Species = "Mouse"
  a = fread(file.path("/home/qians/MamDC/Data/Seqdata/Ourlabdata",Species,paste0(tolower(Species),"_all_reference.count"))) %>% as.data.frame()
  row.names(a) <- a[,1]
  a = a[,-1]
  colnames(a) %<>%  gsub("mouse_male_heart_03","mouse_male_cerebellum_04",.) %>%  gsub("intestine","colon",.)
  #exon length
  length <- read.csv(file = file.path("~/Pseudo/Data/Ref",Species,paste0(Species,".nonreExonlength.bed")),
                     header = FALSE,sep = "\t",stringsAsFactors = FALSE)
  #nonredundant gene length
  gene.length <- tapply(length$V5, length$V4, sum) %>% as.data.frame() 
  
  a$length = gene.length[match(row.names(a),row.names(gene.length)),1]
  a = na.omit(a)
  totalcounts <- colSums(a[,-ncol(a)])
  rpkm <- t(do.call(rbind, lapply(1:length(totalcounts), function(i){
    10^9*a[,i]/a$length/totalcounts[i]}))) %>% as.data.frame()
  row.names(rpkm) = row.names(a)
  colnames(rpkm) = colnames(a)[1:(ncol(a)-1)]
  write.table(rpkm,file = file.path("~/MamDC/Data/Seqdata/Ourlabdata",Species,paste0(Species,".pseudogene.fpkm.txt")),
              sep = "\t",quote = F,row.names = T,col.names = T)
}

{ # normalize A median to 1, see X distribution in both sex (X, A balance)
  rm(list = ls())
  Species = "Mouse"
  a = fread(file.path("~/MamDC/Data/Seqdata/Ourlabdata",Species,paste0(Species,".pseudogene.fpkm.txt"))) %>% as.data.frame()
  b <- read.csv(file.path("/home/qians/Pseudo/Data/Ref",Species,paste0("gene",Species,".bed")),
                header = FALSE,sep = "\t",stringsAsFactors = F)
  a[,c("chr","type")] <- b[match(a$V1,b$V4),c(1,5)]
  a %<>% dplyr::filter(.,chr %in% c(1:100,"X") & type == " protein_coding")
  a[a$chr!="X","chr"]="A"
  row.names(a) = a$V1
  a = a[,-c(1,ncol(a))]
  #correlate = cor(a)
  #pheatmap(correlate)
  sex = lapply(colnames(a)[-ncol(a)], function(x)strsplit(x,split = "_",fixed = TRUE)[[1]][2]) %>% unlist()
  tissue = lapply(colnames(a)[-ncol(a)], function(x)strsplit(x,split = "_",fixed = TRUE)[[1]][3]) %>% unlist()
  mbg = t(apply(a[,-ncol(a)], 1, function(x)tapply(x, paste(sex,tissue,sep = "_"), mean))) %>% as.data.frame()
  mbg$chr = b[match(row.names(mbg),b$V4),1]
  
  cutoff = 1
  focal.tissue.auto = mbg[mbg$chr!="X",-ncol(mbg)]
  mean.auto = do.call(rbind, lapply(1:ncol(focal.tissue.auto), function(i){
    median(focal.tissue.auto[focal.tissue.auto[,i]>cutoff,i])})) %>% as.data.frame()
  row.names(mean.auto) = colnames(focal.tissue.auto)
  
  focal.tissue.x = mbg[mbg$chr=="X",-ncol(mbg)] %>% pivot_longer(.,cols=1:(ncol(mbg)-1))
  focal.tissue.x %<>% dplyr::filter(.,value > cutoff)
  focal.tissue.x$denominator = mean.auto[match(focal.tissue.x$name,row.names(mean.auto)),1]
  focal.tissue.x$normal.x = focal.tissue.x$value/focal.tissue.x$denominator
  focal.tissue.x$sex = lapply(focal.tissue.x$name,function(x)strsplit(x,split = "_",fixed = TRUE)[[1]][1]) %>% 
    unlist() %>% factor(.,levels = c("male","female"),labels = c("Male","Female"))
  focal.tissue.x$tissue = lapply(focal.tissue.x$name,function(x)strsplit(x,split = "_",fixed = TRUE)[[1]][2]) %>% 
    unlist() %>% Hmisc::capitalize() %>% factor(.,levels = c("Brain","Cerebellum","Heart","Colon","Gonad"))
  p = ggplot(focal.tissue.x,aes(log2(normal.x),group=sex))+geom_density(aes(color=sex))+facet_grid(.~tissue)+theme_classic()+
    scale_color_manual(values = c("#3596ca","#a24f51"))+
    geom_vline(xintercept = 0,linetype=2,color="grey")+
    geom_vline(data=filter(focal.tissue.x, tissue=="Gonad"), aes(xintercept=-1), linetype=2,colour="grey") +
    xlab("Log2(FPKM)") + ylab("Density")
  ggsave(p,filename = file.path("~/MamDC/Result/Ourdata",Species,"Picture/bothsex.medianA1.X.density.pdf"),
         width = 9,height = 4)
}

{ # (male, female balance)
  rm(list = ls())
  Species = "Mouse"
  a = fread(file.path("~/MamDC/Data/Seqdata/Ourlabdata",Species,paste0(Species,".pseudogene.fpkm.txt"))) %>% as.data.frame()
  b <- read.csv(file.path("/home/qians/Pseudo/Data/Ref",Species,paste0("gene",Species,".bed")),
                header = FALSE,sep = "\t",stringsAsFactors = F)
  a[,c("chr","type")] <- b[match(a$V1,b$V4),c(1,5)]
  a %<>% dplyr::filter(.,chr %in% c(1:100,"X") & type == " protein_coding")
  a[a$chr!="X","chr"]="A"
  row.names(a) = a$V1
  a = a[,-c(1,ncol(a))]
  #correlate = cor(a)
  #pheatmap(correlate)
  sex = lapply(colnames(a)[-ncol(a)], function(x)strsplit(x,split = "_",fixed = TRUE)[[1]][2]) %>% unlist()
  tissue = lapply(colnames(a)[-ncol(a)], function(x)strsplit(x,split = "_",fixed = TRUE)[[1]][3]) %>% unlist()
  mbg = t(apply(a[,-ncol(a)], 1, function(x)tapply(x, paste(sex,tissue,sep = "_"), mean))) %>% as.data.frame()
  mbg$chr = b[match(row.names(mbg),b$V4),1]
  #group = lapply(colnames(mbg)[-ncol(mbg)], function(x)strsplit(x,split = "_",fixed = TRUE)[[1]][2]) %>% unlist()
  
  c = c()
  cutoff = 1
  for (T in unique(tissue)) {
    focal.tissue = mbg %>% dplyr::select(.,c(paste0("female_",T),paste0("male_",T)))
    focal.tissue = focal.tissue[apply(focal.tissue, 1, function(x){min(x)>cutoff}),]
    focal.tissue$ratio = focal.tissue[,2] / focal.tissue[,1]
    focal.tissue$chr = b[match(row.names(focal.tissue),b$V4),1]
    focal.tissue[focal.tissue$chr!="X","chr"]="A"
    focal.tissue$ratio = focal.tissue$ratio/median(focal.tissue[focal.tissue$chr=="A","ratio"])
    focal.tissue$tissue = T
    focal.tissue = focal.tissue[,c("ratio","tissue","chr")]
    c = rbind(c,focal.tissue)
  }
  c$tissue %<>% Hmisc::capitalize() %>% factor(.,levels = c("Brain","Cerebellum","Heart","Colon","Gonad"))
  p = ggplot(c,aes(log2(ratio),group=chr))+geom_density(aes(color=chr))+facet_grid(.~tissue)+theme_classic()+
    scale_color_manual(values = c("#3596ca","#a24f51"))+
    geom_vline(xintercept = 0,linetype=2,color="grey")+
    xlab("Log2(FPKM)") + ylab("Density")+
    coord_cartesian(xlim = c(-2,2))
  ggsave(p,filename = file.path("~/MamDC/Result/Ourdata",Species,"Picture/male2female.medianA1.X.density.pdf"),
         width = 9,height = 4)
}