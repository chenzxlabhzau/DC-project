#state 1 to 2
rm(list = ls())
gene = fread("/home/qians/MamDC/Result/Roadmap/Chromhmm/Savedata/all.mnemonics.promoter.bed",header = FALSE) %>% as.data.frame()
n = table(gene$V5,gene$V9) %>% as.data.frame()
number = tapply(n$Freq, n$Var2, sum) %>% as.data.frame()
n$all = number[match(n$Var2,row.names(number)),1]
n$prop = n$Freq/n$all
n$number = lapply(n$Var1,function(x)strsplit(as.character(x),"_",fixed = TRUE)[[1]][1]) %>% unlist() %>% as.numeric()



for (state1 in 1:15) {  
  df = matrix(1,nrow = length(unique(n$Var2)),ncol = 15) %>% as.data.frame()
  row.names(df) = unique(n$Var2) %>% as.character()
  colnames(df) = paste(state1,1:15,sep = "_")
  for (i in row.names(df)) {
    #id = rownames(df)[i]
    if (n[n$Var2==i & n$number==1,"Freq"]<=1) {
      df[i,1]=0
    }else {
      df[i,1] = (n[n$Var2==i & n$number==1,"Freq"]-1) / (n[n$Var2==i & n$number==1,"all"]-1)
    }
  }
  for (state2 in 2:15) {
    df[i,state2] = n[n$Var2==i & n$number==state2,"Freq"] / (n[n$Var2==i & n$number==state2,"all"]-1)
  }
  write
}

lapply(1:nrow(df), function(i){
  id = df[i,1];
  if (n[n$Var2==id & n$number==1,"Freq"]==0) {
    df[i,-1] =0
  }else {
    for (m in 2:15) {
      df[i,m] = n[n$Var2==id & n$number== m,3]/(n[n$Var2==id & n$number== m,4]-n[n$Var2==id & n$number==1,"Freq"])
    }
    
  }
    
})

rpkm <- t(do.call(rbind, lapply(1:length(totalcounts), function(i){
  10^9*pseu[,i]/pseu$length/totalcounts[i]}))) %>% as.data.frame()