rm(list = ls())
wd = "/home/qians/MamDC/Data/Seqdata/WangKuster2019MSBRNAseq"
a = fread(file.path(wd,"iBAQprotein.WangKuster2019MSBRNAseq.table1.csv")) %>% as.data.frame()
Species="Human"
b <- read.csv(file.path("~/MamDC/Data/Ref",Species,paste0("gene",Species,".bed")),
              header = FALSE,sep = "\t",stringsAsFactors = F)
a$chr = b[match(a$`Gene ID`,b$V4),1]
a = a[,-c(1:2,32:35)]
a %<>% dplyr::filter(.,chr %in% c(1:22,"X"))

focal.tissue = a[,colnames(a)[c(1:9,ncol(a))]] %>% pivot_longer(.,cols=1:(ncol(.)-1))
focal.tissue$chr %<>% factor(.,levels = c(1:22,"X"))
ggplot(focal.tissue,aes(chr,log10(value+1)))+geom_boxplot(aes(fill=chr),outlier.colour = "white",notch = T)+
  facet_wrap(.~name)+theme_classic()+xlab("Chomosome")+ylab("Normalized iBAQ (log10)")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))+
  scale_fill_manual(values=c(rep("#8db8d5",length(unique(focal.tissue$chr))-1),"#9c0b1f"))+
  guides(fill=FALSE)
ggsave("~/MamDC/Result/WangKuster2019MSBRNAseq/Human/Picture/004.iBAQ.9.chr.pdf",device = "pdf",width = 11,height = 6)

focal.tissue = a[,colnames(a)[c(10:18,ncol(a))]] %>% pivot_longer(.,cols=1:(ncol(.)-1))
focal.tissue$chr %<>% factor(.,levels = c(1:22,"X"))
ggplot(focal.tissue,aes(chr,log10(value+1)))+geom_boxplot(aes(fill=chr),outlier.colour = "white",notch = T)+
  facet_wrap(.~name)+theme_classic()+xlab("Chomosome")+ylab("Normalized iBAQ (log10)")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))+
  scale_fill_manual(values=c(rep("#8db8d5",length(unique(focal.tissue$chr))-1),"#9c0b1f"))+
  guides(fill=FALSE)
ggsave("~/MamDC/Result/WangKuster2019MSBRNAseq/Human/Picture/004.iBAQ.18.chr.pdf",device = "pdf",width = 11,height = 6)

focal.tissue = a[,colnames(a)[c(19:27,ncol(a))]] %>% pivot_longer(.,cols=1:(ncol(.)-1))
focal.tissue$chr %<>% factor(.,levels = c(1:22,"X"))
ggplot(focal.tissue,aes(chr,log10(value+1)))+geom_boxplot(aes(fill=chr),outlier.colour = "white",notch = T)+
  facet_wrap(.~name)+theme_classic()+xlab("Chomosome")+ylab("Normalized iBAQ (log10)")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))+
  scale_fill_manual(values=c(rep("#8db8d5",length(unique(focal.tissue$chr))-1),"#9c0b1f"))+
  guides(fill=FALSE)
ggsave("~/MamDC/Result/WangKuster2019MSBRNAseq/Human/Picture/004.iBAQ.27.chr.pdf",device = "pdf",width = 11,height = 6)

focal.tissue = a[,colnames(a)[c(28:29,ncol(a))]] %>% pivot_longer(.,cols=1:(ncol(.)-1))
focal.tissue$chr %<>% factor(.,levels = c(1:22,"X"))
ggplot(focal.tissue,aes(chr,log10(value+1)))+geom_boxplot(aes(fill=chr),outlier.colour = "white",notch = T)+
  facet_wrap(.~name)+theme_classic()+xlab("Chomosome")+ylab("Normalized iBAQ (log10)")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))+
  scale_fill_manual(values=c(rep("#8db8d5",length(unique(focal.tissue$chr))-1),"#9c0b1f"))+
  guides(fill=FALSE)
ggsave("~/MamDC/Result/WangKuster2019MSBRNAseq/Human/Picture/004.iBAQ.29.chr.pdf",device = "pdf",width = 8,height = 4)


apply(a[,-30], 2, function(x)tapply(x, a$chr=="X", mean))[2,]/apply(a[,-30], 2, function(x)tapply(x, a$chr=="X", mean))[1,]

###########################
#Sample X
e <- matrix(0,nrow = ncol(a)-1,ncol = 4) %>% as.data.frame()
colnames(e) <- c("sample","mean","CI5","CI95")
e$sample <- colnames(a)[-ncol(a)]
for (T in colnames(a)[-ncol(a)]) {
  b <- a[,c("chr",T)]
  colnames(b)[2] = "expression"
  b[b$chr!="X","chr"]="A"
  b$expression <- b$expression/mean(b[b$chr=="A","expression"]) #Allexpress/Aexpress
  d <- list()
  for (j in 1:1000) {
    d[[j]] <- sample(b[b$chr=="X","expression"],size = floor(nrow(b[b$chr=="X",])/2), replace = FALSE)
    boot.mean <- unlist(lapply(d, mean))
  }
  quantile(boot.mean, probs = c(0.25, 0.75))
  e[e$sample==T,2:4] <- c(mean(b[b$chr=="X","expression"]),quantile(boot.mean, probs = c(0.25, 0.75))[1],quantile(boot.mean, probs = c(0.25, 0.75))[2]) %>% as.numeric()
}
e = e[order(e$mean,decreasing = TRUE),]
e$sample %<>% factor(.,levels = .)
ggplot(e,aes(sample,mean))+geom_point(size=3.0)+geom_errorbar(aes(x = sample, ymax=CI5, ymin=CI95,width =0.3))+
  ylab("X:A ratio")+theme_bw()+theme(axis.title.x=element_text(size=18),axis.title.y=element_blank(),
                                     axis.text.x = element_text(size=16),axis.text.y = element_text(size=14))+
  scale_y_continuous(breaks = c(0,0.5,1,2),labels = c(0,0.5,1,2))+
  geom_hline(aes(yintercept=1),color="black",linetype="dashed")+
  geom_hline(aes(yintercept=0.5),color="black",linetype="dashed")+coord_flip()
ggsave(filename = "~/MamDC/Result/WangKuster2019MSBRNAseq/Human/Picture/iBAQ.XAratio.alltissues.sampleX.pdf",
       device = "pdf", width = 5.8, height = 7)
###########################

###########################
#Sample Auto
e <- matrix(0,nrow = ncol(a)-1,ncol = 4) %>% as.data.frame()
colnames(e) <- c("sample","mean","CI5","CI95")
e$sample <- colnames(a)[-ncol(a)]
for (T in colnames(a)[-ncol(a)]) {
  b <- a[,c("chr",T)]
  colnames(b)[2] = "expression"
  b[b$chr!="X","chr"]="A"
  d <- list()
  for (j in 1:1000) {
    d[[j]] <-  mean(b[b$chr=="X","expression"])/mean(sample(b[b$chr=="A","expression"],
                                                                     size = sum(b$chr=="X"), replace = FALSE))
    boot.mean <- unlist(lapply(d, mean))
  }
  #quantile(boot.mean, probs = c(0.25, 0.75))
  e[e$sample==T,2:4] <- c(quantile(boot.mean, probs = c(0.25, 0.5 ,0.75))[2],
                          quantile(boot.mean, probs = c(0.25, 0.5 ,0.75))[1],
                          quantile(boot.mean, probs = c(0.25, 0.5 ,0.75))[3]) %>% as.numeric()
}
e = e[order(e$mean,decreasing = TRUE),]
e$sample %<>% factor(.,levels = .)
ggplot(e,aes(sample,mean))+geom_point(size=3.0)+geom_errorbar(aes(x = sample, ymax=CI5, ymin=CI95,width =0.3))+
  ylab("X:A ratio")+theme_bw()+theme(axis.title.x=element_text(size=18),axis.title.y=element_blank(),
                                     axis.text.x = element_text(size=16),axis.text.y = element_text(size=14))+
  scale_y_continuous(breaks = c(0,0.5,1,2),labels = c(0,0.5,1,2))+
  geom_hline(aes(yintercept=1),color="black",linetype="dashed")+
  geom_hline(aes(yintercept=0.5),color="black",linetype="dashed")+coord_flip()
ggsave(filename = "~/MamDC/Result/WangKuster2019MSBRNAseq/Human/Picture/iBAQ.XAratio.alltissues.sampleA.pdf",
       device = "pdf", width = 5.8, height = 7)
###########################

a[a$chr!="X","chr"]="A"
a %<>% pivot_longer(.,cols=1:(ncol(.)-1))
ggplot(a,aes(chr,log(value+1)))+geom_boxplot(aes(fill=chr),notch = TRUE,outlier.colour = "white")+facet_wrap(.~name)+
  theme_classic()+scale_fill_manual(values=c("#8db8d5","#9c0b1f"))+
  xlab("Chromosome")+ylab("X:A ratio")+
  guides(fill=FALSE)


for (Tissue in colnames(a)[-ncol(a)]) {
  focal.tissue = a[,c(Tissue,"chr")]
  colnames(focal.tissue)[1] = "Express"
  focal.tissue$chr %<>% factor(.,levels = c(1:22,"X"))
  wilcox.test(focal.tissue[focal.tissue$chr=="X",1],focal.tissue[focal.tissue$chr!="X",1])
  ggplot(focal.tissue,aes(chr,log10(Express+1)))+geom_boxplot(outlier.colour = "white",notch = T)
}
ggplot(focal.tissue,aes(chr,log10(value+1)))+geom_boxplot()+facet_grid(.~name)