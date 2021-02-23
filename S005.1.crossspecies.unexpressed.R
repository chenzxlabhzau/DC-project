rm(list = ls());gc();rm(list = ls())
for (ref.species in c("Chicken","Platypus","Opossum")) {
  wd <- "~/MamDC/Data/Seqdata/WangKaessmann2020NatureRNARiboseq"
  datatype = "rna"
  focal.species = "Human"
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
  homo[,paste(focal.species,c("Brain","Liver","Testis"),sep = "_")] = focal.mbg[match(homo$gene1,row.names(focal.mbg)),c("brain","liver","testis")]
  homo[,paste(ref.species,c("Brain","Liver","Testis"),sep = "_")] = ref.mbg[match(homo$gene2,row.names(ref.mbg)),c("brain","liver","testis")]
  homo %<>% na.omit()
  e = matrix(0,nrow = 3,ncol = 5) %>% as.data.frame()
  colnames(e)  = c("Tissue",
                   paste(rep(c(focal.species,ref.species),each=2),
                         rep(c("X","A"),2),
                         sep = "_"))
  e$Tissue = c("Brain","Liver","Testis")
  for (cutoff in c(0.01,1)) {
    #cutoff=1
    for (Tissue in c("Brain","Liver","Testis")) {
      focal.express = homo %>% dplyr::select(.,c(chrtype,ends_with(Tissue)))
      colnames(focal.express)[2:3] = c("Human","Ref") 
      e[e$Tissue==Tissue,2:3] = as.numeric(tapply(focal.express$Human, focal.express$chrtype, function(x)sum(x<= cutoff))[2:1])
      e[e$Tissue==Tissue,4:5] = as.numeric(tapply(focal.express$Ref, focal.express$chrtype, function(x)sum(x<= cutoff))[2:1])
    }
    e$XX = sum(homo$chrtype=="XX")
    e$AA = sum(homo$chrtype=="AA")
    
    e$P1 = 1
    #Human X,A
    for (i in 1:nrow(e)) {
      fish = fisher.test(data.frame(matrix(as.numeric(e[i,c(2,3,6,7)]),nrow = 2)))
      e[i,"P1"] = fish$p.value
    }
    e$P2=1
    #Ref X,A
    for (i in 1:nrow(e)) {
      fish = fisher.test(data.frame(matrix(as.numeric(e[i,4:7]),nrow = 2)))
      e[i,"P2"] = fish$p.value
    }
    e$P3=1
    #Human Ref
    for (i in 1:nrow(e)) {
      fish = fisher.test(data.frame(matrix(as.numeric(e[i,2:5]),nrow = 2)))
      e[i,"P3"] = fish$p.value
    }
  
    e[,c(2,4)] = e[,c(2,4)] / sum(homo$chrtype=="XX") *100 # Note: Not X or A
    e[,c(3,5)] = e[,c(3,5)] / sum(homo$chrtype=="AA") *100
    
    elonger = e %>% pivot_longer(.,cols=2:5) %>% as.data.frame()
    elonger$Species = lapply(elonger$name, function(x)strsplit(x,split = "_",fixed = TRUE)[[1]][1]) %>% unlist()
    elonger$chrtype = lapply(elonger$name, function(x)strsplit(x,split = "_",fixed = TRUE)[[1]][2]) %>% unlist()
    
    elonger$name %<>% gsub("_"," ",.) %>% factor(.,levels = gsub("_"," ",colnames(e)[2:5]))
    
      p = ggplot(elonger,aes(Tissue,value,fill=name))+
        geom_bar(stat = "identity",position = "dodge")+
        theme_classic()+ylab("Frequency of unexpressed genes (%)")+
        theme(axis.text.x = element_text(size=12,angle = 45, hjust = 1, vjust = 1),
              axis.text.y = element_text(size=12),
              axis.title.x = element_blank(),
              axis.title.y = element_text(size=14))+
        scale_y_continuous(expand = c(0, 0))+ coord_cartesian(ylim = c(0,ceiling(max(elonger$value)*1.1)))+
        scale_fill_manual(values = c("#9c0b1f","#035782","#DC7684","#8db8d5"))+
        geom_text(aes(x,y,label =lab),inherit.aes=FALSE,
                  data = data.frame(x = 1:3,
                                    y = ceiling(max(elonger$value)*1.08),
                                    lab = cut(e$P1,breaks = c(1,0.05,0.01,0.001,0),
                                              labels = rev(c("n.s.","*","**","***"))),
                                    Tissue = factor(e$Tissue,levels = e$Tissue)), vjust = 1, size=7)+
        geom_text(aes(x,y,label =lab),inherit.aes=FALSE,
                  data = data.frame(x = 1:3,
                                    y = ceiling(max(elonger$value)*1.05),
                                    lab = cut(e$P2,breaks = c(1,0.05,0.01,0.001,0),
                                              labels = rev(c("n.s.","*","**","***"))),
                                    Tissue = factor(e$Tissue,levels = e$Tissue)), vjust = 0.7, size=5)+
        geom_text(aes(x,y,label =lab),inherit.aes=FALSE,
                  data = data.frame(x = 1:3,
                                    y = ceiling(max(elonger$value)*1.01),
                                    lab = cut(e$P3,breaks = c(1,0.05,0.01,0.001,0),
                                              labels = rev(c("n.s.","*","**","***"))),
                                    Tissue = factor(e$Tissue,levels = e$Tissue)), vjust = 1, size=7)+
        guides(fill=guide_legend(title=NULL))+theme(legend.position="bottom",legend.direction = "horizontal")
    
    ggsave(paste0("~/MamDC/Result/WangKaessmann2020NatureRNARiboseq/Human/Picture/S005.1.Freq.unexpressed.human",ref.species,"XA.FPKM",cutoff,".pdf"),
           p,device = "pdf",width = 6.5,height = 6)
    topptx(p,paste0("~/MamDC/Result/WangKaessmann2020NatureRNARiboseq/Human/Picture/S005.1.Freq.unexpressed.human",ref.species,"XA.FPKM",cutoff,".pptx"),
           width = 6.5,height = 6)
  }
  rm(list = ls())
}
