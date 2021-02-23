#cd ~/MamDC/Data/Seqdata/RobinDavid2015CurrentBioRNAseq
#wget ftp://ftp.ensembl.org/pub/release-65/gtf/mus_musculus/Mus_musculus.NCBIM37.65.gtf.gz
#awk -F ";" '{print $1}' Mus_musculus.NCBIM37.65.gtf | awk -F "\t" '{print $1"NF}'| sed 's/gene_id//g' | sed 's/"//g' | sed 's/ //g' | sort | uniq > chr.gene.txt

rm(list = ls())
Num = "S007.1."

a = read.csv("~/MamDC/Data/Seqdata/RobinDavid2015CurrentBioRNAseq/GSE64960_8-10week_gonad_counts.csv")
b = read.csv("~/MamDC/Data/Seqdata/RobinDavid2015CurrentBioRNAseq/chr.gene.txt",sep = "\t",header = FALSE,stringsAsFactors = FALSE)
a[,c("chr","type")] = b[match(a$X,b$V3),1:2]
row.names(a) = a$X
a %<>% dplyr::filter(., chr %in% c(1:100,"X") & type == "protein_coding") %>% dplyr::select(!c(X,type))
a[a$chr!="X","chr"]="A"
cpm = apply(a[,-ncol(a)], 2, function(x)x/sum(x)* 10^6) %>% as.data.frame()
#cpm$chr = a$chr
group = rep(c("DMEf","WTf","WTm"),each=2)
mbg = apply(cpm,1,function(x)tapply(x, group, mean)) %>% t() %>% as.data.frame()

mbg$chr = a[match(row.names(mbg),row.names(a)),"chr"]
mbg = mbg[apply(mbg[,-ncol(mbg)], 1, function(x)min(x) >0),]

mbg[,1] = mbg[,1]/median(mbg[mbg$chr=="A",1])
mbg[,2] = mbg[,2]/median(mbg[mbg$chr=="A",2])
mbg[,3] = mbg[,3]/median(mbg[mbg$chr=="A",3])

mbg %<>% pivot_longer(.,cols=1:3)
#mbg %<>% dplyr::filter(.,value > 10) %>% as.data.frame()
#mbg[mbg$name=="DMEf",3] = mbg[mbg$name=="DMEf",3]/ median(mbg[mbg$chr=="A"&mbg$name=="DMEf",3])
#mbg[mbg$name=="WTf",3] = mbg[mbg$name=="WTf",3] / median(mbg[mbg$chr=="A"&mbg$name=="WTf",3])
#mbg[mbg$name=="WTm",3] = mbg[mbg$name=="WTm",3] / median(mbg[mbg$chr=="A"&mbg$name=="WTm",3])
ggplot(mbg,aes(chr,log2(value+1),fill=name))+geom_boxplot(notch = T)

