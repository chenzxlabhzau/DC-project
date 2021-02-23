rm(list = ls())
wd = "~/MamDC/Result/CardosoKaessmann2019NatureRNAseq/Human/Picture/"
for (cor in c("Positive","Negative")) {
  for (Tissue in c("Brain","Liver")) {
    a = read.csv(paste0(wd,"003.Human.",cor,"correlation.GObp.",Tissue,".csv"))
    go_id = a$ID
    mat = GO_similarity(go_id)
    pdf(file = paste0(wd,"003.Human.",cor,"correlation.GObp.",Tissue,".pdf"),width = 10,height = 7)
    df = simplifyGO(mat, column_title = paste(nrow(a),"GO terms clustered on",Tissue,cor,"genesets",sep = " "))
    dev.off()
  }
}