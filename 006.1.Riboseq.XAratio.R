#Using Ribo-seq data from "Transcriptome and translatome co-evolution in mammals"
rm(list = ls())
library(ggplotify)
library(eoffice)
for (Species in c("Human","Mouse")) {
  wd <- "~/MamDC/Data/Seqdata/WangKaessmann2020NatureRNARiboseq"   #current directory
  #Species = tolower(Species)
  directory = file.path(wd,tolower(Species))
  Files <- grep("ribo_[1-9].txt$",list.files(directory),value=TRUE) %>% 
    grep("brain|liver|testis",.,value=TRUE)
  filePath <- sapply(Files, function(x){paste(directory,x,sep='/')})   
  data <- lapply(filePath, function(x){ fread(x)})  
  a <- data[1] %>% as.data.frame()
  a = a[,c(1,ncol(a))]
  colnames(a) <- c("Gene",strsplit(colnames(a)[1],split=".",fix=TRUE)[[1]][1]) # geneID, sample name
  for (i in 2:length(data)) {
    b <- data[i] %>% as.data.frame()
    b = b[,c(1,ncol(b))]
    colnames(b) <- c("Gene",strsplit(colnames(b)[1],split=".",fix=TRUE)[[1]][1])
    a <- merge(a,b,by="Gene")
  }
  row.names(a) = a$Gene
  group = colnames(a)[-1] %>% lapply(., function(x)strsplit(x,split = "_",fixed = TRUE)[[1]][2]) %>% unlist()
  mbg = t(apply(a[,-1], 1, function(x)tapply(x, group, mean))) %>% as.data.frame()
  
  library(Hmisc)
  b= read.csv(file.path("~/MamDC/Data/Ref",capitalize(Species),paste0("gene",capitalize(Species),".bed")),
              header = FALSE,sep = "\t",stringsAsFactors = F)
  mbg$chr = b[match(row.names(mbg),b$V4),1]
  mbg %<>% dplyr::filter(.,chr %in% c(1:100,"X"))
  mbg[mbg$chr!="X","chr"]="A"
  mbg %<>% pivot_longer(.,cols=1:(ncol(.)-1)) %>% na.omit() %>% dplyr::filter(.,value > 1) %>% as.data.frame()
  mbg$name %<>% capitalize(.)
  
  library(Rmisc)
  #smy = ddply(mbg,.(chr,name),summarise, ci = log2(CI(value)))
  smy = ddply(mbg,.(chr,name),summarise, ci = log2(median(value)))
  #smy$interval = rep(c("upper","mean","lower"),nrow(smy)/3)
  smy$interval = "median"
  smy.wider = pivot_wider(smy,names_from = chr,values_from = ci) #correspond to facet
  #### fail example: pivot_wider(smy,names_from = name,values_from = mean)
  if (FALSE) { #area
    smy.wider %<>% dplyr::filter(.,interval %in% c("upper","lower"))
    ggplot(mbg,aes(log2(value)))+geom_density(aes(color=chr))+
      xlab("log2(FPKM)")+ylab("Frequency")+
      facet_wrap(.~name)+theme_classic()+
      scale_color_manual(values = c("#8db8d5","#9c0b1f")) +
      geom_vline(data = smy.wider,mapping = aes(xintercept = A),color="#8db8d5",linetype="dashed")+
      geom_vline(data = smy.wider,mapping = aes(xintercept = X),color="#9c0b1f",linetype="dashed")+
      geom_area(data=smy.wider,aes(x=A,y=0.2),fill="8db8d5")
  }
  
  #smy.wider %<>% dplyr::filter(.,interval %in% "median")
  testP = c()
  for (i in unique(mbg$name)) {
    test = wilcox.test(mbg[mbg$name==i&mbg$chr=="X",3],mbg[mbg$name==i&mbg$chr=="A",3],
                       alternative = "less")
    testP = c(testP,test$p.value)
  }
  #
  if (Species == "Human") {
    ymax = 0.25
  }
  if (Species == "Mouse") {
    ymax = 0.25
  }
  #threetissues = c("Brain","Liver","Testis")
  #mbg$name %<>% factor(.,levels = c(threetissues,setdiff(unique(mbg$name),threetissues)))
  p1 = ggplot(mbg,aes(log2(value)))+geom_density(aes(color=chr))+
    xlab("Log2(FPKM)")+ylab("Frequency")+
    facet_wrap(.~name)+theme_classic()+
    scale_color_manual(values = c("#8db8d5","#9c0b1f")) +
    theme(axis.title.x = element_text(size = 14),
          axis.title.y = element_text(size = 14),
          strip.text = element_text(size=12),
          legend.position=c(0.95,0.80))+
          #legend.position=c(0.23,0.9),legend.direction = "horizontal")+
    geom_vline(data = smy.wider,mapping = aes(xintercept = A),color="#8db8d5",linetype="dashed")+
    geom_vline(data = smy.wider,mapping = aes(xintercept = X),color="#9c0b1f",linetype="dashed")+
    scale_y_continuous(expand = c(0, 0))+coord_cartesian(ylim = c(0,ymax))+
    geom_text(aes(x,y,label =lab),
              data = data.frame(x = apply(smy.wider[,-(1:2)], 1,max)+3.2,
                                y = max(density(mbg$value)$y)*1.2,
                                lab = cut(testP,breaks = c(1,0.05,0.01,0.001,0),
                                          labels = rev(c("n.s.","*","**","***"))),
                                #sex = factor(c("brain","liver","testis")),levels = c("brain","liver","testis")),
                                name = factor(unique(mbg$name)),levels = unique(mbg$name)), #facet use name
              vjust = 1, size=5)
  #cumulative curve
  double = mbg[mbg$chr=="X",]
  double$value = double$value*2
  double$chr = "2X"
  mbg = rbind(mbg,double)
  mbg$chr %<>% factor(.,levels = c("A","X","2X"))
  
  p2 = ggplot(mbg,aes(log2(value)))+stat_ecdf(aes(color=chr,linetype=chr))+facet_wrap(.~name)+
    scale_color_manual(values = c("#8db8d5","#9c0b1f","#9c0b1f")) +theme_classic()+
    scale_linetype_manual(values=c("solid","solid","dotted"))+
    theme(axis.title.x = element_text(size = 14),
          axis.title.y = element_text(size = 14),
          strip.text = element_text(size=12),
          legend.position=c(0.95,0.25))+
    ylab("Cumulative frequency (%)")+
    guides(fill=guide_legend(title=NULL))
    

    
  if (Species=="Human") {
    ggsave(paste0("/home/qians/MamDC/Result/WangKaessmann2020NatureRNARiboseq/",Species,"/Picture/006.",Species,".Riboseq.density.pdf"),
           p1,device = "pdf",width = 7, height = 4)
    topptx(p1,paste0("/home/qians/MamDC/Result/WangKaessmann2020NatureRNARiboseq/",Species,"/Picture/006.",Species,".Riboseq.density.pptx"),
           width = 5, height = 3)
    topptx(p2,paste0("/home/qians/MamDC/Result/WangKaessmann2020NatureRNARiboseq/",Species,"/Picture/006.",Species,".Riboseq.cumulate.pptx"),
           width = 5, height = 3)
  } else{
    ggsave(paste0("/home/qians/MamDC/Result/WangKaessmann2020NatureRNARiboseq/",Species,"/Picture/006.",Species,".Riboseq.density.pdf"),
           p1,device = "pdf",width = 7, height = 6)
    topptx(p1,paste0("/home/qians/MamDC/Result/WangKaessmann2020NatureRNARiboseq/",Species,"/Picture/006.",Species,".Riboseq.density.pptx"),
           width = 5, height = 3)
    topptx(p2,paste0("/home/qians/MamDC/Result/WangKaessmann2020NatureRNARiboseq/",Species,"/Picture/006.",Species,".Riboseq.cumulate.pptx"),
           width = 5, height = 3)
  }
}