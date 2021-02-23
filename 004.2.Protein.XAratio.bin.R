#X:A ratio at protein level using expression bins
rm(list = ls())
wd = "/home/qians/MamDC/Data/Seqdata/WangKuster2019MSBRNAseq"
a = fread(file.path(wd,"iBAQprotein.WangKuster2019MSBRNAseq.table1.csv")) %>% as.data.frame()
Species="Human"
b <- read.csv(file.path("~/MamDC/Data/Ref",Species,paste0("gene",Species,".bed")),
              header = FALSE,sep = "\t",stringsAsFactors = F)
a$chr = b[match(a$`Gene ID`,b$V4),1]
a = a[,-c(1:2,32:35)]
a %<>% dplyr::filter(.,chr %in% c(1:22,"X"))
a[a$chr!="X","chr"]="A"
a[1:3,1:3]

#27 tissues
for (n in c(9,18,27)) {
  e = matrix(0,nrow = 25,ncol = 2*9+1) %>% as.data.frame()
  focal.tissue = a[,colnames(a)[c((n-8):n,ncol(a))]] %>% pivot_longer(.,cols=1:(ncol(.)-1))
  for (i in 1:18) {
    colnames(e)[i] = tapply(focal.tissue$value, paste(focal.tissue$chr,focal.tissue$name,sep = "_"), function(x)rev(quantile(x,seq(0.76,1,0.01))))[i] %>%
      names()
    e[,i] = tapply(focal.tissue$value, paste(focal.tissue$chr,focal.tissue$name,sep = "_"), function(x)rev(quantile(x,seq(0.76,1,0.01))))[[i]] %>%
      as.numeric()
  }
  if (all(lapply(colnames(e)[1:9], function(x)strsplit(x,split = "_",fixed = T)[[1]][2]) %>% unlist() == 
          lapply(colnames(e)[10:18], function(x)strsplit(x,split = "_",fixed = T)[[1]][2]) %>% unlist())) {
    e[,ncol(e)] = 1:25
    colnames(e)[ncol(e)]= "bins"
  }
  e %<>% pivot_longer(.,cols=1:(ncol(.)-1)) %>% as.data.frame()
  for (i in seq(1,50,by = 2)) {
    e[(9 * i  - 8):(9 * (i + 1)),4] = rep(e[((9 * (i + 1) - 8):(9 * (i + 1))), 3]/e[((9 * i  - 8):(9 * i)), 3],2)
  }
  
  e$chr = lapply(e$name,function(x)strsplit(x,split = "_",fixed = T)[[1]][1]) %>% unlist()
  e$tissue = lapply(e$name,function(x)strsplit(x,split = "_",fixed = T)[[1]][2]) %>% unlist()
  all(e$name== paste(e$chr,e$tissue,sep = "_"))
  
  scale.factor = 9
  p = ggplot()+
    geom_point(data = e,aes(bins,log10(value),color=chr),size=2,shape=1)+
    scale_color_manual(values = c("#8db8d5","#9c0b1f")) +
    geom_point(data = e,aes(bins,scale.factor*V4),shape=18,size=4)+
    #scale_shape_manual(values = 18) +
    scale_y_continuous(name = "Normalized iBAQ (log10)", 
                       sec.axis = sec_axis(trans = ~./scale.factor,
                                           name = "X:A ratio"))+
    facet_wrap(.~tissue)+theme_classic()+
    geom_hline(yintercept = 0.5*scale.factor,color="gray",linetype="dotted")+
    geom_hline(yintercept = 1.0*scale.factor,color="gray",linetype="dotted")+
    geom_smooth(data = e,aes(bins,scale.factor*V4))
  ggsave(filename = paste0("~/MamDC/Result/WangKuster2019MSBRNAseq/Human/Picture/004.iBAQXAratio.bins.",n,".pdf"),
         p,device = "pdf", width = 11,height = 6)
  topptx(p,filename = paste0("~/MamDC/Result/WangKuster2019MSBRNAseq/Human/Picture/004.iBAQXAratio.bins.",n,".pptx"),
                  width = 11,height = 6)
  
  if (n==1) {
    ggplot()+
      geom_point(data = e[e$tissue=="Brain",],aes(bins,log10(value),color=chr),size=2,shape=1)+
      scale_color_manual(values = c("#8db8d5","#9c0b1f")) +
      #geom_point(data = e[e$tissue=="Brain",],aes(bins,scale.factor*V4),shape=18,size=4)+
      #scale_shape_manual(values = 18) +
      scale_y_continuous(name = "Normalized iBAQ (log10)")+
      facet_wrap(.~tissue)+theme_classic()+xlab("Protein exprssion level bins")+
      geom_hline(yintercept = 0.5*scale.factor,color="gray",linetype="dotted")+
      geom_hline(yintercept = 1.0*scale.factor,color="gray",linetype="dotted")
    
    ggplot()+
      geom_point(data = e[e$tissue=="Brain",],aes(bins,log10(value),color=chr),size=2,shape=1)+
      scale_color_manual(values = c("#8db8d5","#9c0b1f")) +
      geom_point(data = e[e$tissue=="Brain",],aes(bins,scale.factor*V4),shape=18,size=4)+
      #scale_shape_manual(values = 18) +
      scale_y_continuous(name = "Normalized iBAQ (log10)", 
                         sec.axis = sec_axis(trans = ~./scale.factor,
                                             name = "X:A ratio"))+
      facet_wrap(.~tissue)+theme_classic()+xlab("Protein exprssion level bins")+
      geom_hline(yintercept = 0.5*scale.factor,color="gray",linetype="dotted")+
      geom_hline(yintercept = 1.0*scale.factor,color="gray",linetype="dotted")+
      geom_smooth(data = e[e$tissue=="Brain",],aes(bins,scale.factor*V4))
  }
}

#last 2 tissues
if (TRUE) {
  e = matrix(0,nrow = 25,ncol = 2*2+1) %>% as.data.frame()
  focal.tissue = a[,colnames(a)[28:30]] %>% pivot_longer(.,cols=1:(ncol(.)-1))
  for (i in 1:4) {
    colnames(e)[i] = tapply(focal.tissue$value, paste(focal.tissue$chr,focal.tissue$name,sep = "_"), function(x)rev(quantile(x,seq(0.76,1,0.01))))[i] %>%
      names()
    e[,i] = tapply(focal.tissue$value, paste(focal.tissue$chr,focal.tissue$name,sep = "_"), function(x)rev(quantile(x,seq(0.76,1,0.01))))[[i]] %>%
      as.numeric()
  }
  if (all(lapply(colnames(e)[1:2], function(x)strsplit(x,split = "_",fixed = T)[[1]][2]) %>% unlist() == 
          lapply(colnames(e)[3:4], function(x)strsplit(x,split = "_",fixed = T)[[1]][2]) %>% unlist())) {
    e[,ncol(e)] = 1:25
    colnames(e)[ncol(e)]= "bins"
  }
  e %<>% pivot_longer(.,cols=1:(ncol(.)-1)) %>% as.data.frame()
  for (i in seq(1,50,by = 2)) {
    e[(2 * i  - 1):(2 * (i + 1)),4] = rep(e[((2 * (i + 1) - 1):(2 * (i + 1))), 3]/e[((2 * i  - 1):(2 * i)), 3],2)
  }
  
  e$chr = lapply(e$name,function(x)strsplit(x,split = "_",fixed = T)[[1]][1]) %>% unlist()
  e$tissue = lapply(e$name,function(x)strsplit(x,split = "_",fixed = T)[[1]][2]) %>% unlist()
  all(e$name== paste(e$chr,e$tissue,sep = "_"))
  
  scale.factor = 9
  p = ggplot()+
    geom_point(data = e,aes(bins,log10(value),color=chr),size=2,shape=1)+
    scale_color_manual(values = c("#8db8d5","#9c0b1f")) +
    geom_point(data = e,aes(bins,scale.factor*V4),shape=18,size=4)+
    #scale_shape_manual(values = 18) +
    scale_y_continuous(name = "Normalized iBAQ (log10)", 
                       sec.axis = sec_axis(trans = ~./scale.factor,
                                           name = "X:A ratio"))+
    facet_wrap(.~tissue)+theme_classic()+
    geom_hline(yintercept = 0.5*scale.factor,color="gray",linetype="dotted")+
    geom_hline(yintercept = 1.0*scale.factor,color="gray",linetype="dotted")+
    geom_smooth(data = e,aes(bins,scale.factor*V4))
  ggsave(filename = paste0("~/MamDC/Result/WangKuster2019MSBRNAseq/Human/Picture/004.iBAQXAratio.bins.",29,".pdf"),
         p,device = "pdf", width = 8,height = 4)
  topptx(p,filename = paste0("~/MamDC/Result/WangKuster2019MSBRNAseq/Human/Picture/004.iBAQXAratio.bins.",29,".pptx"),
         width = 8,height = 4)
}

#Display brain liver and testis
if (TRUE) {
  display = c("brain","liver","testis")
  e = matrix(0,nrow = 25,ncol = 2*length(display)+1) %>% as.data.frame()
  focal.tissue = a[,c(capitalize(display),"chr")] %>% pivot_longer(.,cols=1:(ncol(.)-1))
  for (i in 1:(2*length(display))) {
    colnames(e)[i] = tapply(focal.tissue$value, paste(focal.tissue$chr,focal.tissue$name,sep = "_"), function(x)rev(quantile(x,seq(0.76,1,0.01))))[i] %>%
      names()
    e[,i] = tapply(focal.tissue$value, paste(focal.tissue$chr,focal.tissue$name,sep = "_"), function(x)rev(quantile(x,seq(0.76,1,0.01))))[[i]] %>%
      as.numeric()
  }
  if (all(lapply(colnames(e)[1:length(display)], function(x)strsplit(x,split = "_",fixed = T)[[1]][2]) %>% unlist() == 
          lapply(colnames(e)[(length(display)+1):(2*length(display))], function(x)strsplit(x,split = "_",fixed = T)[[1]][2]) %>% unlist())) {
    e[,ncol(e)] = 1:25
    colnames(e)[ncol(e)]= "bins"
  }
  e %<>% pivot_longer(.,cols=1:(ncol(.)-1)) %>% as.data.frame()
  for (i in seq(1,50,by = 2)) {
    e[(length(display) * i  - (length(display)-1)):(length(display) * (i + 1)),4] = 
      rep(e[((length(display) * (i + 1) - (length(display)-1)):(length(display) * (i + 1))), 3]/
            e[((length(display) * i  - (length(display)-1)):(length(display) * i)), 3],2)
  }
  
  e$chr = lapply(e$name,function(x)strsplit(x,split = "_",fixed = T)[[1]][1]) %>% unlist()
  e$tissue = lapply(e$name,function(x)strsplit(x,split = "_",fixed = T)[[1]][2]) %>% unlist()
  all(e$name== paste(e$chr,e$tissue,sep = "_"))
  
  scale.factor = 9
  p = ggplot()+
    geom_point(data = e,aes(bins,log10(value),color=chr),size=2,shape=1)+
    scale_color_manual(values = c("#8db8d5","#9c0b1f")) +
    xlab("Protein expression level bins")+
    geom_point(data = e,aes(bins,scale.factor*V4),shape=18,size=4)+
    #scale_shape_manual(values = 18) +
    scale_y_continuous(name = "Normalized iBAQ (log10)", 
                       sec.axis = sec_axis(trans = ~./scale.factor,
                                           name = "X:A ratio"))+
    facet_wrap(.~tissue)+theme_classic()+
    geom_hline(yintercept = 0.5*scale.factor,color="gray",linetype="dotted")+
    geom_hline(yintercept = 1.0*scale.factor,color="gray",linetype="dotted")+
    geom_smooth(data = e,aes(bins,scale.factor*V4))+
    theme(axis.title.x = element_text(size = 14),
          axis.title.y = element_text(size = 14),
          legend.position=c(0.15,0.9),legend.direction = "horizontal")+
    guides(fill=guide_legend(title=NULL))
  
  ggsave(filename = paste0("~/MamDC/Result/WangKuster2019MSBRNAseq/Human/Picture/004.iBAQXAratio.bins.",paste(display,collapse = ""),".pdf"),
         p,device = "pdf", width = 11,height = 6)
  topptx(p,filename = paste0("~/MamDC/Result/WangKuster2019MSBRNAseq/Human/Picture/004.iBAQXAratio.bins.",paste(display,collapse = ""),".pptx"),
         width = 8,height = 4)
}