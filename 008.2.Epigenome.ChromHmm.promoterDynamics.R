#(0~15) state dynamics across X/A
rm(list = ls())
gene = fread("/home/qians/MamDC/Result/Roadmap/Chromhmm/Savedata/all.mnemonics.promoter.bed",header = FALSE) %>% as.data.frame()
n = paste(gene$V5,gene$V9,sep = ".") %>% unique() %>% as.data.frame()
n$gene = lapply(n$., function(x)strsplit(x,".",fixed=T)[[1]][2]) %>% unlist()
freq = table(n$gene) %>% as.data.frame()
freq$chr = gene[match(freq$Var1,gene$V9),1]

freq %<>% dplyr::filter(.,chr != "chrY")
freq[freq$chr!="chrX","chr"]="A"
freq[freq$chr=="chrX","chr"]="X"
p1 = ggplot(freq,aes(Freq))+geom_bar(aes(y = ..prop..,group=chr,fill=chr),position = "dodge")+
  xlab("Number of ChromHMM states")+ylab("Proportion")+theme_bw()+ 
  theme(axis.title.x=element_text(size=14),
        axis.text.x = element_text(size=12),
        axis.title.y=element_text(size=14),
        axis.text.y = element_text(size=12),
        legend.position=c(.2,.8),legend.direction = "horizontal",
        legend.title = element_blank(),legend.background = element_rect(fill="transparent"))+
  #guides(fill=FALSE)+
  scale_fill_manual(values = c("#035782","#9c0b1f"))
topptx(p1,filename = "~/MamDC/Result/Roadmap/Chromhmm/Picture/Human.promoter.dynamics.pptx",
       width = 4.5,height = 4)
ggsave(p1,filename = "~/MamDC/Result/Roadmap/Chromhmm/Picture/Human.promoter.dynamics.pdf",
       width = 4.5,height = 4)
#freq$Freq2 = apply(freq[,-1], 1, function(x){(3^((x[2]=="X")+1-1)-2)*as.numeric(x[1])})

p2 = ggplot(freq,aes(chr,Freq,fill=chr))+geom_boxplot(notch = TRUE,outlier.colour = "white")+theme_classic()+
  xlab("Chromosome")+ylab("Number of epigenetic state")+
  theme(axis.title.x=element_text(size=14),
        axis.text.x = element_text(size=12),
        axis.title.y=element_text(size=14),
        axis.text.y = element_text(size=12))+
  scale_fill_manual(values = c("#035782","#9c0b1f"))+
  geom_segment(aes(x=1, y=15.2, xend=2, yend=15.2))+
  annotate("text", x=1.5, y=15.3, label="*",size=5)+
  guides(fill=FALSE)
ggsave(p2,filename = "~/MamDC/Result/Roadmap/Chromhmm/Picture/Human.promoter.dynamics.boxplot.pdf",
       width = 2.5,height = 4)

  ggplot(freq,aes(Freq))+geom_bar(aes(y = ..prop..,group=chr,fill=chr),position = "dodge")+
  xlab("Number of ChromHMM states")+ylab("Proportion")+theme_bw()+ 
  theme(axis.title.x=element_text(size=14),
        axis.text.x = element_text(size=12),
        axis.title.y=element_text(size=14),
        axis.text.y = element_text(size=12),
        legend.position=c(.2,.8),legend.direction = "horizontal",
        legend.title = element_blank(),legend.background = element_rect(fill="transparent"))