#Data from "Ageing hallmarks exhibit organ-specific temporal signatures"
## Stage
rm(list = ls())
a <- fread("~/MamDC/Data/Seqdata/NicholasTony2020NatureRNAseq/GSE132040_190214_A00111_0269_AHH3J3DSXX_190214_A00111_0270_BHHMFWDSXX.csv") %>% 
  as.data.frame()
a = a[1:(nrow(a)-5),]
row.names(a) = a$gene
cpm = apply(a[,-1], 2, function(x)x/sum(x)*10^6) %>% as.data.frame()
cpm[1:3,1:3]
colnames(cpm) %<>% gsub(".gencode.vM19","",.)
info = read.csv("~/MamDC/Data/Seqdata/NicholasTony2020NatureRNAseq/GSE132040_MACA_Bulk_metadata.csv")
info$characteristics..age %<>% as.numeric()
table(colnames(cpm) %in% info$Sample.name)
b = read.csv("~/MamDC/Data/Seqdata/NicholasTony2020NatureRNAseq/Mouse.gencode.gene.bed",
             header = FALSE,sep = "\t",stringsAsFactors = F)
cpm$type = b[match(row.names(cpm),b$V6),5]
cpm %<>% dplyr::filter(.,type=="protein_coding")

info$source.name %<>% lapply(.,function(x)strsplit(x,split = "_",fixed = T)[[1]][1]) %>% unlist() %>%
  gsub("Limb","Limb Muscle",.) %>% gsub("Small","Small Intestine",.)
info %<>% dplyr::filter(.,!grepl("NA",source.name))
info$sample = paste(info$source.name,info$characteristics..sex,info$characteristics..age,sep = "_") # tissue + sex + stage

df = matrix(0,nrow = length(unique(info$source.name)),ncol = 5) %>% as.data.frame()
colnames(df) = c(paste(rep(c("X","A"),2),rep(c("Positive","Negative"),each=2),sep = "_"),"Meanratio") #X_Up A_Up X_Down A_Down
row.names(df) = unique(info$source.name)
for (Tissue in unique(info$source.name)) {
  focal.tissue = dplyr::filter(info,source.name == Tissue)
  focal.express = cpm[,focal.tissue$Sample.name]
  if (all(focal.tissue$Sample.name==colnames(focal.express))) { # check order
    cor.r = t(do.call(rbind, lapply(1:nrow(focal.express), function(i){
      cor.test(as.numeric(focal.express[i,]), focal.tissue$characteristics..age)$estimate})))
    cor.p = t(do.call(rbind, lapply(1:nrow(focal.express), function(i){
      cor.test(as.numeric(focal.express[i,]), focal.tissue$characteristics..age)$p.value})))
    focal.express$R = as.numeric(cor.r)
    focal.express$P = as.numeric(cor.p)
    focal.express$chr = b[match(row.names(focal.express),b$V6),1]
    focal.express %<>% dplyr::filter(.,chr %in% paste0("chr",c(1:100,"X")))
    focal.express[focal.express$chr!="chrX","chr"]="chrA"
    focal.express[is.na(focal.express)]=0
    
    df[Tissue,"Meanratio"] = sum(focal.express$chr=="chrX")/nrow(focal.express)
    df[Tissue,"X_Positive"] = sum(focal.express$R> 0.4 & focal.express$P<0.05 & focal.express$chr=="chrX")
    df[Tissue,"A_Positive"] = sum(focal.express$R> 0.4 & focal.express$P<0.05 & focal.express$chr=="chrA")
    df[Tissue,"X_Negative"] = sum(focal.express$R< -0.4 & focal.express$P<0.05 & focal.express$chr=="chrX")
    df[Tissue,"A_Negative"] = sum(focal.express$R< -0.4 & focal.express$P<0.05 & focal.express$chr=="chrA")
  }
}
df$Pratio = df$X_Positive/(df$X_Positive+df$A_Positive)
df$Nratio = df$X_Negative/(df$X_Negative+df$A_Negative)
df$tissue = row.names(df)
df = df[sort(df$tissue),]

testP = c()
for (i in 1:nrow(df)) {
  testP[i] = fisher.test(as.data.frame(matrix(as.numeric(df[i,1:4]),nrow = 2)))$p.value
}
df %<>% dplyr::select(.,c(ends_with("ratio"),tissue)) %>% 
  pivot_longer(.,cols=1:3) %>% dplyr::filter(.,name != "Meanratio") %>% as.data.frame()
df$name %<>% factor(.,levels = c("Pratio","Nratio"))
ggplot(df,aes(tissue,value*100,fill=name))+geom_bar(stat = "identity",position = "dodge")+theme_classic()+
  theme(axis.text.x = element_text(size = 12,angle = 45,hjust = 1,vjust = 1),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 12),
        axis.title.y = element_text(size = 14))+
  xlab("Tissue")+ylab("(%) X-linked genes")+
  scale_y_continuous(expand = c(0, 0))+coord_cartesian(ylim = c(0,8))+
  geom_hline(yintercept = 4.4,color="grey",linetype = "dashed")+
  geom_text(aes(x,y,label =lab),inherit.aes=FALSE,
            data = data.frame(x = 1:length(unique(df$tissue)),
                              y = 7.5,
                              lab = cut(testP,breaks = c(1,0.05,0.01,0.001,0),
                                        labels = rev(c("","*","**","***")))),
            vjust = 1, size=7)
ggsave(filename = paste0("~/MamDC/Result/NicholasTony2020NatureRNAseq/Mouse/Picture/003.Mouse.RatioXlinked.Tissue.pdf")
       ,device = "pdf",width = 5.5, height = 4.5)
