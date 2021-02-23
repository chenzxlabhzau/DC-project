#age with state number
rm(list = ls())
library(ggpubr)
age = read.csv("~/MamDC/Data/Seqdata/Genorigin/Human/Homo_sapiens.csv",header = TRUE,sep = ",",stringsAsFactors = FALSE)
age$gene_age %<>% gsub(">","",.) %>% as.numeric()

#15 state
gene = fread("/home/qians/MamDC/Result/Roadmap/Chromhmm/Savedata/all.mnemonics.promoter.bed",header = FALSE) %>% as.data.frame()
n = paste(gene$V5,gene$V9,sep = ".") %>% unique() %>% as.data.frame()
n$gene = lapply(n$., function(x)strsplit(x,".",fixed=T)[[1]][2]) %>% unlist()
freq = table(n$gene) %>% as.data.frame()
age$s15 = freq[match(age$ensembl_id,freq$Var1),2]

#18 state
s18 = fread("/home/qians/MamDC/Result/Roadmap/Chromhmm/Savedata/18state/all_18_core_K27ac_hg38lift_mnemonics.promoter.bed",header = FALSE) %>% as.data.frame()
n = paste(s18$V4,s18$V8,sep = ".") %>% unique() %>% as.data.frame()
n$gene = lapply(n$., function(x)strsplit(x,".",fixed=T)[[1]][2]) %>% unlist()
freq = table(n$gene) %>% as.data.frame()
age$s18 = freq[match(age$ensembl_id,freq$Var1),2]

#50 state
s50 = fread("/home/qians/MamDC/Result/Roadmap/Chromhmm/Savedata/50state/all.50_segments.hg38.promoter.bed",header = FALSE) %>% as.data.frame()
n = paste(s50$V4,s50$V8,sep = ".") %>% unique() %>% as.data.frame()
n$gene = lapply(n$., function(x)strsplit(x,".",fixed=T)[[1]][2]) %>% unlist()
freq = table(n$gene) %>% as.data.frame()
age$s50 = freq[match(age$ensembl_id,freq$Var1),2]

age = na.omit(age)
cor.test(as.numeric(age$gene_age),age$s18)
age$chr = gene[match(age$ensembl_id,gene$V9),1]
age %<>% dplyr::filter(., chr %in% paste0("chr",c(1:100,"X")))
age[age$chr!="chrX","chr"]="A"
age[age$chr=="chrX","chr"]="X"
#divide group
tmp = age[age$chr=="X",]
tmp = tmp[order(tmp$gene_age),]
tmp$group = rep(c("Young","Middle","Old"),c(285,327,218))
tmp[285,"gene_age"]
tmp[612,"gene_age"]

age$group = cut(age$gene_age,breaks = c(0,332,645,5000),labels = c("Young","Middle","Old"))
sum(table(age$gene_age,paste(age$chr,age$group,sep = "_"))[,5])
age = age[,c(2,4:8)] %>% pivot_longer(.,cols=2:4)

for (state in c("s15","s18","s50")) {
  age$group %<>% factor(.,levels = c("Young","Middle","Old"))
  my_comparisons <- list(c("Young","Middle"),c("Middle","Old"),c("Young","Old"))
  p = age %>% dplyr::filter(.,name %in% state) %>% 
    ggplot(., aes(group,value,color=group))+geom_boxplot(notch = TRUE,outlier.alpha = 0)+
    stat_compare_means(comparisons = my_comparisons)+theme_classic()+
    ylab("Number of state")+
    theme(axis.title.x=element_blank(),
          axis.text.x = element_text(size=12),
          axis.title.y=element_text(size=14),
          axis.text.y = element_text(size=12),
          strip.background  = element_blank(),
          strip.text = element_text(size=10))+
    guides(color=FALSE)+
    scale_color_manual(values = c("#D6B6BB","#CB7582","#9c0b1f"))+
    facet_grid(.~chr)
  topptx(p,filename = paste0("/home/qians/MamDC/Result/Roadmap/Chromhmm/Picture/Number.State.XA.",state,".pptx"),
         width = 4,height = 3)
  ggsave(p,filename = paste0("/home/qians/MamDC/Result/Roadmap/Chromhmm/Picture/Number.State.XA.",state,".pdf"),
         width = 5,height = 4)
}

###discard: calculate mean dynamic state value and correlate with age
if (FALSE) {
  tmp = age[age$chr=="chrX",]
  tmp = tmp[order(tmp$gene_age),]
  tmp$group = rep(c("Young","Middle","Old"),c(285,327,218))
  tmp$group %<>% factor(.,levels = c("Young","Middle","Old"))
  my_comparisons <- list(c("Young","Middle"),c("Middle","Old"),c("Young","Old"))
  ggplot(tmp, aes(group,s50,color=group))+geom_boxplot(notch = TRUE,outlier.alpha = 0)+
    stat_compare_means(comparisons = my_comparisons)+theme_bw()+
    scale_color_manual(values = c("#D6B6BB","#CB7582","#9c0b1f"))+
    ylab("Number of state")+
    theme(axis.title.x=element_blank(),
          axis.text.x = element_text(size=12),
          axis.title.y=element_text(size=14),
          axis.text.y = element_text(size=12))+
    guides(color=FALSE)
  
  
  freq$chr = gene[match(freq$Var1,gene$V9),1]
  
  freq %<>% dplyr::filter(.,chr != "chrY")
  freq[freq$chr!="chrX","chr"]="A"
  freq[freq$chr=="chrX","chr"]="X"
  freq$age = age[match(freq$Var1,age$ensembl_id),2] %>% gsub(">","",.) #%>% as.numeric()
  freq15 = freq
  rm(freq)
  
  
  
  #age with state number
  age = read.csv("~/MamDC/Data/Seqdata/Genorigin/Human/Homo_sapiens.csv",header = TRUE,sep = ",",stringsAsFactors = FALSE)
  age$s15 = freq[match(age$ensembl_id,freq$Var1),2]
  age$gene_age %<>% gsub(">","",.) %>% factor(.,levels = c(1:5000))
  age$gene_age %<>% droplevels()
  age %<>% na.omit()
  ggplot(age,aes(gene_age,s15))+geom_boxplot()
  
  #18 state
  s18 = fread("/home/qians/MamDC/Result/Roadmap/Chromhmm/Savedata/18state/all_18_core_K27ac_hg38lift_mnemonics.promoter.bed",header = FALSE) %>% as.data.frame()
  n = paste(s18$V4,s18$V8,sep = ".") %>% unique() %>% as.data.frame()
  n$gene = lapply(n$., function(x)strsplit(x,".",fixed=T)[[1]][2]) %>% unlist()
  freq = table(n$gene) %>% as.data.frame()
  age$s18 = freq[match(age$ensembl_id,freq$Var1),2]
  
  #50 state
  s50 = fread("/home/qians/MamDC/Result/Roadmap/Chromhmm/Savedata/50state/all.50_segments.hg38.promoter.bed",header = FALSE) %>% as.data.frame()
  n = paste(s50$V4,s50$V8,sep = ".") %>% unique() %>% as.data.frame()
  n$gene = lapply(n$., function(x)strsplit(x,".",fixed=T)[[1]][2]) %>% unlist()
  freq = table(n$gene) %>% as.data.frame()
  age$s50 = freq[match(age$ensembl_id,freq$Var1),2]
  
  age$chr = gene[match(age$ensembl_id,gene$V9),1]
  tmp = age
  age %<>% dplyr::filter(.,chr == "chrX")
  
  t15 = tapply(age$s15, age$gene_age, function(x)mean(x)) %>% as.data.frame()
  t15$age = row.names(t15)
  t18 = tapply(age$s18, age$gene_age, function(x)mean(x)) %>% as.data.frame()
  t18$age = row.names(t18)
  t50 = tapply(age$s50, age$gene_age, function(x)mean(x,na.rm=TRUE)) %>% as.data.frame()
  t50$age = row.names(t50)
  
  t15$t18 = t18[match(t15$age,t18$age),1]
  t15$t50 = t50[match(t15$age,t50$age),1]
  t15 = t15[,c(2,1,3,4)]
  colnames(t15) = c("age","t15","t18","t50")
  cor.test(as.numeric(t15$age),t15$t18,method = "spearman")
  
  t15 %<>% na.omit()
  t15$t15 = t15$t15 / 15
  t15$t18 = t15$t18 / 18
  t15$t50 = t15$t50 /50
  t15$age %<>% factor(.,levels = 1:5000) %>% droplevels()
  t15 %<>% pivot_longer(.,cols=2:4)
  ggplot(t15,aes(age,value,color=name,group=name))+geom_point()+geom_smooth(method = "lm")
  #0.69 0.72 0.69
  
  cor.test(as.numeric(t15$age),t15$t50,method = "spearman")
  
  
  freq %<>% dplyr::filter(.,chr != "chrY")
  freq[freq$chr!="chrX","chr"]="A"
  freq[freq$chr=="chrX","chr"]="X"
  
  cor.test(as.numeric(age$gene_age),age$s50)
  
  freq$age = age[match(freq$Var1,age$ensembl_id),2] %>% gsub(">","",.) #%>% as.numeric()
  freq %<>% dplyr::filter(.,chr == "X")
  cor.test(freq$Freq,as.numeric(freq$age))
  
  age.freq = tapply(freq$Freq, freq$age, mean) %>% as.data.frame()
  age.freq$age = as.numeric(row.names(age.freq))
  cor.test(log10(age.freq$age),age.freq$.)
  ggplot(age.freq,aes(log10(age),.))+geom_point()+geom_smooth(method = "lm")
  
  freq$age %<>% factor(.,levels = as.character(sort(as.numeric(unique(freq$age)))))
  freq$age 
  age$num = freq[match(age$ensembl_id)]
}
