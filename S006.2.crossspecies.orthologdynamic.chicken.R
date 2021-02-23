rm(list = ls())
gc()
rm(list = ls())
#ref.species = "Platypus" #"Opossum" #"Chicken" 
#for (ref.species in c("Opossum","Chicken","Platypus")) {
ref.species = "Chicken"
focal.species = "Human"
wd = "~/MamDC/Data/Seqdata/CardosoKaessmann2019NatureRNAseq"
base.object = ls()
#focal species
focal.express = fread(file.path(wd,focal.species,paste0(focal.species,".allgene.fpkm.txt"))) %>% as.data.frame()
#colnames(focal.express) %<>% gsub("KidneyTestis_4w_Male","Kidney_4w_Male",.)
b <- read.csv(file.path("~/MamDC/Data/Ref",focal.species,paste0("gene",focal.species,".bed")),header = F,sep = "\t") 
b$V5 %<>% gsub(" ","",.)
focal.express[,c("chr","type")] <- b[match(focal.express$V1,b$V4),c(1,5)]
focal.express %<>% filter(chr %in% c(1:100,"X") & type %in% "protein_coding") %>% select(.,-type)
row.names(focal.express) = focal.express$V1

s1 = c("12wpc","13wpc","19wpc","20wpc","2ypb","4ypb",
       "25ypb","28ypb","29ypb","32ypb","39ypb")
stage = data.frame(s1 = s1,
                   s2 = c("12wpc","13wpc","19wpc","20wpc","toddler","toddler",
                          "youngadult","youngadult","youngadult","youngadult","youngmiddle"))
coldata = read.csv(paste0("~/Pseudo/Result/",focal.species,"/Savedata/coldata.csv"),row.names = 1,sep = ",") 
coldata %<>% dplyr::filter(stage2 %in% s1)
coldata$s1 = coldata$stage2
for (i in 1:nrow(coldata)) {
  coldata[i,"s1"] %<>% gsub("32ypb", "youngadult",.) %>%  gsub("2ypb", "toddler",.,useBytes = TRUE) %>% 
    gsub("4ypb", "toddler",.) %>%
    gsub("25ypb", "youngadult",.) %>% gsub("28ypb", "youngadult",.) %>% 
    gsub("29ypb", "youngadult",.) %>% 
    gsub("39ypb", "youngmiddle",.) #%>%
   #gsub("46ypb", "oldmiddle",.) %>%
   #gsub("50ypb", "oldmiddle",.) %>%
   #gsub("53ypb", "oldmiddle",.) %>%
   #gsub("54ypb", "oldmiddle",.) %>%
   #gsub("55ypb", "senior",.) %>%
   #gsub("58ypb", "senior",.)
}
focal.express = focal.express[,coldata$name]
focal.mbg = apply(focal.express, 1, function(x)tapply(x, paste(coldata$condition,coldata$s1,sep = "_"), mean)) %>%
  t() %>% as.data.frame()

correlate = data.frame(focal = c("12wpc","13wpc","19wpc","20wpc","toddler",
                                 "youngadult","youngmiddle"),
                       ref=c("10","12","12","14","0dph","10wph","Adult"))

ref.express = fread(file.path("/home/qians/MamDC/Data/Seqdata/CardosoKaessmann2019NatureRNAseq",ref.species,paste0(ref.species,".allgene.fpkm.txt"))) %>% 
  as.data.frame()
row.names(ref.express) = ref.express$V1
ref.express = ref.express[,-1]
#colnames(ref.express) %<>% gsub("adult","Adult",.)
ref.coldata = colnames(ref.express) %>% as.data.frame()
ref.coldata$stage = ref.coldata$. %>% lapply(.,function(x)strsplit(x,split = "_",fixed = TRUE)[[1]][2]) %>% unlist()
ref.coldata %<>% dplyr::filter(.,stage %in% c(correlate$ref,"adult"))
ref.coldata$tissue = ref.coldata$. %>% lapply(.,function(x)strsplit(x,split = "_",fixed = TRUE)[[1]][1]) %>% unlist()
ref.express = ref.express[,ref.coldata$.]  
all(colnames(ref.express) == ref.coldata$.)
ref.mbg = apply(ref.express, 1, function(x)tapply(x, paste(ref.coldata$tissue,capitalize(ref.coldata$stage),sep = "_"), mean)) %>%
  t() %>% as.data.frame()

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

homo %<>% dplyr::select(.,c("gene1","gene2","chrtype")) %>% dplyr::filter(.,chrtype %in% c("AA","XX"))
homo[,4:(ncol(focal.mbg)+4-1)] = focal.mbg[match(homo$gene1,row.names(focal.mbg)),]
homo[,(ncol(homo)+1):(ncol(homo)+ncol(ref.mbg))] = ref.mbg[match(homo$gene2,row.names(ref.mbg)),]
#homo = homo[apply(homo[,-c(1:3)], 1, min)>1,]
row.names(homo) = homo$gene1

Tissue = "Brain" ###################################################################
express = data.frame()
testP = c()
for (i in 1:nrow(correlate)) {
  one = homo %>% dplyr::select(.,c(chrtype,starts_with(Tissue))) %>% na.omit()
  #one = one[apply(one[,-1], 1, min)>1,]
  #table(one[apply(one[,-1], 1, min) >1,1])
  onestage = one %>% dplyr::select(chrtype,matches(paste0(c(paste0(as.character(correlate[i,])[1],"$"),
                                                            paste0(as.character(correlate[i,])[2],"$")),collapse = "|"),ignore.case = FALSE))
  onestage = onestage[apply(onestage[,-1], 1, min) >1,]
  colnames(onestage)[2:3] = c("focal","ref")
  tapply(onestage$focal/onestage$ref, onestage$chrtype, median)
  median(onestage[onestage$chrtype=="AA","focal"]/onestage[onestage$chrtype=="AA","ref"])
  onestage$focal = onestage$focal / unname(tapply(onestage$focal/onestage$ref, onestage$chrtype, median)["AA"])
  onestage$ratio = onestage$focal / onestage$ref
  test =wilcox.test(onestage[onestage$chrtype=="AA","ratio"],onestage[onestage$chrtype=="XX","ratio"])
  testP = c(testP,test$p.value)
  onestage$stage = correlate[i,1]
  express = rbind(express,onestage)
}
express$stage %<>% factor(.,levels = c("12wpc","13wpc","19wpc","20wpc",
                                       "toddler","youngadult","youngmiddle"))
cor= cor.test(as.numeric(with(express[express$chrtype=="XX",],tapply(ratio,stage,median))),1:nrow(correlate),method = "spearman")
p1 = ggplot(express,aes(stage,log2(ratio),fill=chrtype))+geom_boxplot(notch = TRUE,outlier.alpha = 0)+
  theme_classic()+coord_cartesian(ylim = c(-4,4))+
  xlab("Development stage") + ylab(paste0("Expression ratio in ",tolower(Tissue)))+
  theme(axis.text.x = element_text(size=10,angle = 45, hjust = 1, vjust = 1),
        axis.title.x =  element_text(size=12),
        axis.text.y = element_text(size=10),
        axis.title.y =  element_text(size=12))+
  scale_fill_manual(values = c("#035782","#9c0b1f")) +
  geom_hline(yintercept = 0,color="black",linetype="dashed")+
  geom_hline(yintercept = -1,color="grey",linetype="dashed")+
  scale_y_continuous(breaks = c(-4,-2,-1,0,1,2,4),labels =  c(-4,-2,-1,0,1,2,4))+
  annotate(geom="text", x=2.5, y=3.7, size=5,
           label=paste0("rho=",round(cor$estimate,2),",\n P=",format(cor$p.value,2)))
ggsave(paste0("~/MamDC/Result/WangKaessmann2020NatureRNARiboseq/Human/Picture/S006.1.Orthologratio.dynamics.human.",ref.species,".",Tissue,".pdf"),
       p1,device = "pdf",width = 7,height = 6)
topptx(p1,
       paste0("~/MamDC/Result/WangKaessmann2020NatureRNARiboseq/Human/Picture/S006.1.Orthologratio.dynamics.human.",ref.species,".",Tissue,".pptx"),
       width = 7,height = 6)

Tissue = "Liver" ###################################################################
express = data.frame()
for (i in c(1:nrow(correlate))) {
  one = homo %>% dplyr::select(.,c(chrtype,starts_with(Tissue))) %>% na.omit()
  #one = one[apply(one[,-1], 1, min)>1,]
  #table(one[apply(one[,-1], 1, min) >1,1])
  onestage = one %>% dplyr::select(chrtype,matches(paste0(c(paste0(as.character(correlate[i,])[1],"$"),
                                                            paste0(as.character(correlate[i,])[2],"$")),collapse = "|"),ignore.case = FALSE))
  onestage = onestage[apply(onestage[,-1], 1, min) >1,]
  colnames(onestage)[2:3] = c("focal","ref")
  tapply(onestage$focal/onestage$ref, onestage$chrtype, median)
  median(onestage[onestage$chrtype=="AA","focal"]/onestage[onestage$chrtype=="AA","ref"])
  onestage$focal = onestage$focal / unname(tapply(onestage$focal/onestage$ref, onestage$chrtype, median)["AA"])
  onestage$ratio = onestage$focal / onestage$ref
  onestage$stage = correlate[i,1]
  express = rbind(express,onestage)
}
express$stage %<>% factor(.,levels = c("12wpc","13wpc","19wpc","20wpc",
                                       "toddler","youngadult","youngmiddle"))
cor= cor.test(as.numeric(with(express[express$chrtype=="XX",],tapply(ratio,stage,median))),1:nrow(correlate),method = "spearman")
p2 = ggplot(express,aes(stage,log2(ratio),fill=chrtype))+geom_boxplot(notch = TRUE,outlier.alpha = 0)+
  theme_classic()+coord_cartesian(ylim = c(-4,4))+
  xlab("Development stage") +  ylab(paste0("Expression ratio in ",tolower(Tissue)))+
  theme(axis.text.x = element_text(size=10,angle = 45, hjust = 1, vjust = 1),
        axis.title.x =  element_text(size=12),
        axis.text.y = element_text(size=10),
        axis.title.y =  element_text(size=12))+
  scale_fill_manual(values = c("#035782","#9c0b1f")) +
  geom_hline(yintercept = 0,color="black",linetype="dashed")+
  geom_hline(yintercept = -1,color="grey",linetype="dashed")+
  scale_y_continuous(breaks = c(-4,-2,-1,0,1,2,4),labels =  c(-4,-2,-1,0,1,2,4))+
  annotate(geom="text", x=2.5, y=3.7, size=5,
           label=paste0("rho=",round(cor$estimate,2),",\n P=",format(cor$p.value,2)))
ggsave(paste0("~/MamDC/Result/WangKaessmann2020NatureRNARiboseq/Human/Picture/S006.1.Orthologratio.dynamics.human.",ref.species,".",Tissue,".pdf"),
       p2,device = "pdf",width = 7,height = 6)
topptx(p2,
       paste0("~/MamDC/Result/WangKaessmann2020NatureRNARiboseq/Human/Picture/S006.1.Orthologratio.dynamics.human.",ref.species,".",Tissue,".pptx"),
       width = 7,height = 6)

Tissue = "Testis" ###################################################################
express = data.frame()
for (i in c(1:3,6,7)) {
  one = homo %>% dplyr::select(.,c(chrtype,starts_with(Tissue))) %>% na.omit()
  #one = one[apply(one[,-1], 1, min)>1,]
  #table(one[apply(one[,-1], 1, min) >1,1])
  onestage = one %>% dplyr::select(chrtype,matches(paste0(c(paste0(as.character(correlate[i,])[1],"$"),
                                                            paste0(as.character(correlate[i,])[2],"$")),collapse = "|"),ignore.case = FALSE))
  onestage = onestage[apply(onestage[,-1], 1, min) >1,]
  colnames(onestage)[2:3] = c("focal","ref")
  tapply(onestage$focal/onestage$ref, onestage$chrtype, median)
  median(onestage[onestage$chrtype=="AA","focal"]/onestage[onestage$chrtype=="AA","ref"])
  onestage$focal = onestage$focal / unname(tapply(onestage$focal/onestage$ref, onestage$chrtype, median)["AA"])
  onestage$ratio = onestage$focal / onestage$ref
  onestage$stage = correlate[i,1]
  express = rbind(express,onestage)
}
express$stage %<>% factor(.,levels = c("12wpc","13wpc","19wpc","youngadult","youngmiddle"))
cor= cor.test(as.numeric(with(express[express$chrtype=="XX",],tapply(ratio,stage,median))),1:5,method = "spearman")
p3 = ggplot(express,aes(stage,log2(ratio),fill=chrtype))+geom_boxplot(notch = TRUE,outlier.alpha = 0)+
  theme_classic()+coord_cartesian(ylim = c(-4,4))+
  xlab("Development stage") +  ylab(paste0("Expression ratio in ",tolower(Tissue)))+
  theme(axis.text.x = element_text(size=10,angle = 45, hjust = 1, vjust = 1),
        axis.title.x =  element_text(size=12),
        axis.text.y = element_text(size=10),
        axis.title.y =  element_text(size=12))+
  scale_fill_manual(values = c("#035782","#9c0b1f")) +
  geom_hline(yintercept = 0,color="black",linetype="dashed")+
  geom_hline(yintercept = -1,color="grey",linetype="dashed")+
  scale_y_continuous(breaks = c(-4,-2,-1,0,1,2,4),labels =  c(-4,-2,-1,0,1,2,4))+
  annotate(geom="text", x=2.5, y=3.7, size=5,
           label=paste0("rho=",round(cor$estimate,2),",\n P=",format(cor$p.value,2)))
ggsave(paste0("~/MamDC/Result/WangKaessmann2020NatureRNARiboseq/Human/Picture/S006.1.Orthologratio.dynamics.human.",ref.species,".",Tissue,".pdf"),
       p3,device = "pdf",width = 7,height = 6)
topptx(p3,
       paste0("~/MamDC/Result/WangKaessmann2020NatureRNARiboseq/Human/Picture/S006.1.Orthologratio.dynamics.human.",ref.species,".",Tissue,".pptx"),
       width = 7,height = 6)