rm(list = ls())
gc()

## Raw FPKM
rm(list = ls())
a = read.csv("~/MamDC/Data/Seqdata/LinHe2012PNASRNAseq/LinHe.HumanChicken.raw.csv",header = TRUE)
a$chrtype = "Other"
a[a$Human_Chr=="X" & a$Chicken_Chr %in% c(1,4),"chrtype"]="XX"
a[a$Human_Chr!="X","chrtype"]="AA"
colnames(a) %<>% gsub("Male[0-9]","Male",.)

e = matrix(0,nrow = length(colnames(a)[3:12]),ncol = 5) %>% as.data.frame()
colnames(e)  = c("Tissue",
                 paste(rep(c("Human","Chicken"),each=2),
                       rep(c("X","A"),2),
                       sep = "_"))
e$Tissue = gsub("Human_","",colnames(a)[3:12])
for (cutoff in c(0.01,1)) {
  for (Tissue in gsub("Human_","",colnames(a)[3:12])) {
    focal.express = a %>% dplyr::select(.,c(chrtype,ends_with(Tissue)))
    colnames(focal.express)[2:3] = c("Human","Chicken") 
    e[e$Tissue==Tissue,2:3] = as.numeric(tapply(focal.express$Human, focal.express$chrtype, function(x)sum(x<= cutoff))[2:1])
    e[e$Tissue==Tissue,4:5] = as.numeric(tapply(focal.express$Chicken, focal.express$chrtype, function(x)sum(x<= cutoff))[2:1])
  }
  e$XX = sum(a$chrtype=="XX")
  e$AA = sum(a$chrtype=="AA")
  
  e$P1 = 1
  #Human X,A
  for (i in 1:nrow(e)) {
    fish = fisher.test(data.frame(matrix(as.numeric(e[i,c(2,3,6,7)]),nrow = 2)))
    e[i,"P1"] = fish$p.value
  }
  
  e$P2=1
  #Chicken X,A
  for (i in 1:nrow(e)) {
    fish = fisher.test(data.frame(matrix(as.numeric(e[i,4:7]),nrow = 2)))
    e[i,"P2"] = fish$p.value
  }
  
  e$P3=1
  #Human Chicken
  for (i in 1:nrow(e)) {
    fish = fisher.test(data.frame(matrix(as.numeric(e[i,2:5]),nrow = 2)))
    e[i,"P3"] = fish$p.value
  }
  
  e[,c(2,4)] = e[,c(2,4)] / sum(a$chrtype=="XX") *100 # Note: Not X or A
  e[,c(3,5)] = e[,c(3,5)] / sum(a$chrtype=="AA") *100
  
  elonger = e %>% pivot_longer(.,cols=2:5) %>% as.data.frame()
  elonger$Species = lapply(elonger$name, function(x)strsplit(x,split = "_",fixed = TRUE)[[1]][1]) %>% unlist()
  elonger$chrtype = lapply(elonger$name, function(x)strsplit(x,split = "_",fixed = TRUE)[[1]][2]) %>% unlist()
  
  elonger$name %<>% gsub("_"," ",.) %>% factor(.,levels = gsub("_"," ",colnames(e)[2:5]))
  elonger$Tissue %<>% gsub("_Female"," (F)",.) %>% gsub("_Male"," (M)",.)
  
  if (cutoff=1) {
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
                data = data.frame(x = 1:10,
                                  y = ceiling(max(elonger$value)*1.08),
                                  lab = cut(e$P1,breaks = c(1,0.05,0.01,0.001,0),
                                            labels = rev(c("n.s.","*","**","***"))),
                                  Tissue = factor(e$Tissue,levels = e$Tissue)), vjust = 1, size=7)+
      geom_text(aes(x,y,label =lab),inherit.aes=FALSE,
                data = data.frame(x = 1:10,
                                  y = ceiling(max(elonger$value)*1.05),
                                  lab = cut(e$P2,breaks = c(1,0.05,0.01,0.001,0),
                                            labels = rev(c("n.s.","*","**","***"))),
                                  Tissue = factor(e$Tissue,levels = e$Tissue)), vjust = 0.7, size=5)+
      geom_text(aes(x,y,label =lab),inherit.aes=FALSE,
                data = data.frame(x = 1:10,
                                  y = ceiling(max(elonger$value)*1.01),
                                  lab = cut(e$P3,breaks = c(1,0.05,0.01,0.001,0),
                                            labels = rev(c("n.s.","*","**","***"))),
                                  Tissue = factor(e$Tissue,levels = e$Tissue)), vjust = 1, size=7)+
      guides(fill=guide_legend(title=NULL))+theme(legend.position="bottom",legend.direction = "horizontal")
  }
  if (cutoff=0.01) {
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
                data = data.frame(x = 1:10,
                                  y = ceiling(max(elonger$value)*1.08),
                                  lab = cut(e$P1,breaks = c(1,0.05,0.01,0.001,0),
                                            labels = rev(c("n.s.","*","**","***"))),
                                  Tissue = factor(e$Tissue,levels = e$Tissue)), vjust = 1, size=6)+
      geom_text(aes(x,y,label =lab),inherit.aes=FALSE,
                data = data.frame(x = 1:10,
                                  y = ceiling(max(elonger$value)*1.00),
                                  lab = cut(e$P2,breaks = c(1,0.05,0.01,0.001,0),
                                            labels = rev(c("n.s.","*","**","***"))),
                                  Tissue = factor(e$Tissue,levels = e$Tissue)), vjust = 0.7, size=6)+
      geom_text(aes(x,y,label =lab),inherit.aes=FALSE,
                data = data.frame(x = 1:10,
                                  y = ceiling(max(elonger$value)*0.92),
                                  lab = cut(e$P3,breaks = c(1,0.05,0.01,0.001,0),
                                            labels = rev(c("n.s.","*","**","***"))),
                                  Tissue = factor(e$Tissue,levels = e$Tissue)), vjust = 1, size=6)+
      guides(fill=guide_legend(title=NULL))+theme(legend.position="bottom",legend.direction = "horizontal")
  }
  ggsave(paste0("~/MamDC/Result/LinHe2012PNASRNAseq/Picture/006.3.Freq.unexpressed.humanchickenXA.FPKM",cutoff,".pdf"),
         p,device = "pdf",width = 6.5,height = 6)
  topptx(p,paste0("~/MamDC/Result/LinHe2012PNASRNAseq/Picture/006.3.Freq.unexpressed.humanchickenXA.FPKM",cutoff,".pptx"),
         width = 6.5,height = 6)
}