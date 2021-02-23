# body plot
rm(list = ls())
library(gganatogram)
library(dplyr)
library(viridis)
library(gridExtra)
library(ggplotify)
library(eoffice)
# for human male
if (TRUE) {
  prim = hgMale_key
  prim$value = 100
  wang = read.csv("~/MamDC/Result/WangKuster2019MSBRNAseq/Human/Savedata/Human_tissues_XAratio.txt",
                  sep = "\t",header = TRUE,stringsAsFactors = F)
  gtex = read.csv("~/MamDC/Result/GTEx/Savedata/Human_tissues_XAratio.txt",
                  sep = "\t",header = TRUE,stringsAsFactors = F)
  if (TRUE) {
    prim[prim$organ=="adipose_tissue","value"] = mean(c(wang[wang$sample=="adipose_tissue",2], gtex[gtex$sample=="Adipose - Subcutaneous",2]))
    prim[prim$organ=="adrenal_gland","value"] = mean(c(wang[wang$sample=="adrenal_gland",2], gtex[gtex$sample=="Adrenal Gland",2]))
    prim[prim$organ=="amygdala","value"] = gtex[gtex$sample=="Brain - Amygdala",2]
    prim[prim$organ=="aorta","value"] = gtex[gtex$sample=="Artery - Aorta",2]
    prim[prim$organ=="appendix","value"] = wang[wang$sample=="vermiform_appendix",2]
    prim[prim$organ=="atrial_appendage","value"] = gtex[gtex$sample=="Heart - Atrial Appendage",2]
    prim[prim$organ=="bone","value"] = wang[wang$sample=="bone_marrow",2]
    prim[prim$organ=="bone_marrow","value"] = wang[wang$sample=="bone_marrow",2]
    prim[prim$organ=="brain","value"] = mean(c(wang[wang$sample=="cerebral_cortex",2],gtex[gtex$sample=="Brain - Cortex",2]))
    prim[prim$organ=="breast","value"] = gtex[gtex$sample=="Breast - Mammary Tissue",2]
    prim[prim$organ=="cerebellar_hemisphere","value"] = gtex[gtex$sample=="Brain - Cerebellar Hemisphere",2]
    prim[prim$organ=="cerebellum","value"] = gtex[gtex$sample=="Brain - Cerebellum",2]
    prim[prim$organ=="colon","value"] = mean(c(wang[wang$sample=="colon",2],gtex[gtex$sample=="Colon - Sigmoid",2],gtex[gtex$sample=="Colon - Transverse",2]))
    prim[prim$organ=="coronary_artery","value"] = gtex[gtex$sample=="Artery - Coronary",2]
    prim[prim$organ=="duodenum","value"] = wang[wang$sample=="duodenum",2]
    prim[prim$organ=="esophagus","value"] = mean(c(wang[wang$sample=="esophagus",2], gtex[gtex$sample=="Esophagus - Mucosa",2], gtex[gtex$sample=="Esophagus - Mucosa",2]))
    prim[prim$organ=="frontal_cortex","value"] = gtex[gtex$sample=="Brain - Frontal Cortex (BA9)",2]
    prim[prim$organ=="gall_bladder","value"] = mean(c(wang[wang$sample=="gall_bladder",2], gtex[gtex$sample=="Bladder",2]))
    prim[prim$organ=="gastroesophageal_junction","value"] = gtex[gtex$sample=="Esophagus - Gastroesophageal Junction",2]
    prim[prim$organ=="heart","value"] = mean(c(wang[wang$sample=="heart",2], gtex[gtex$sample=="Heart - Atrial Appendage",2], gtex[gtex$sample=="Heart - Left Ventricle",2]))
    prim[prim$organ=="hippocampus","value"] = gtex[gtex$sample=="Brain - Hippocampus",2]
    prim[prim$organ=="ileum","value"] = gtex[gtex$sample=="Small Intestine - Terminal Ileum",2]
    prim[prim$organ=="kidney","value"] = mean(c(wang[wang$sample=="kidney",2], gtex[gtex$sample=="Kidney - Cortex",2], gtex[gtex$sample=="Kidney - Medulla",2]))
    prim[prim$organ=="left_ventricle","value"] = gtex[gtex$sample=="Heart - Left Ventricle",2]
    prim[prim$organ=="liver","value"] = mean(c(wang[wang$sample=="liver",2], gtex[gtex$sample=="Liver",2]))
    prim[prim$organ=="lung","value"] = mean(c(wang[wang$sample=="lung",2], gtex[gtex$sample=="Lung",2]))
    prim[prim$organ=="nerve","value"] = gtex[gtex$sample=="Nerve - Tibial",2]
    prim[prim$organ=="pancreas","value"] = mean(c(wang[wang$sample=="pancreas",2], gtex[gtex$sample=="Pancreas",2]))
    prim[prim$organ=="pituitary_gland","value"] = gtex[gtex$sample=="Pituitary",2]
    prim[prim$organ=="prostate","value"] =  mean(c(wang[wang$sample=="prostate_gland",2], gtex[gtex$sample=="Prostate",2]))
    prim[prim$organ=="rectum","value"] = wang[wang$sample=="rectum",2]
    prim[prim$organ=="salivary_gland","value"] = mean(c(wang[wang$sample=="saliva.secreting_gland",2],gtex[gtex$sample=="Minor Salivary Gland",2]))
    prim[prim$organ=="skeletal_muscle","value"] = mean(c(wang[wang$sample=="skeletal_muscle_tissue",2],gtex[gtex$sample=="Muscle - Skeletal",2]))
    prim[prim$organ=="skin","value"] = mean(c(wang[wang$sample=="zone_of_skin",2],gtex[gtex$sample=="Skin - Not Sun Exposed (Suprapubic)",2],gtex[gtex$sample=="Skin - Sun Exposed (Lower leg)",2]))
    prim[prim$organ=="small_intestine","value"] = mean(c(wang[wang$sample=="small_intestine",2], gtex[gtex$sample=="Small Intestine - Terminal Ileum",2]))
    prim[prim$organ=="smooth_muscle","value"] = wang[wang$sample=="smooth_muscle_tissue",2]
    prim[prim$organ=="spinal_cord","value"] = gtex[gtex$sample=="Brain - Spinal cord (cervical c-1)",2]
    prim[prim$organ=="spleen","value"] = mean(c(wang[wang$sample=="spleen",2],gtex[gtex$sample=="Spleen",2]))
    prim[prim$organ=="stomach","value"] = mean(c(wang[wang$sample=="stomach",2],gtex[gtex$sample=="Stomach",2]))  
    prim[prim$organ=="testis","value"] = mean(c(wang[wang$sample=="testis",2],gtex[gtex$sample=="Testis",2])) 
    prim[prim$organ=="thyroid_gland","value"] = mean(c(wang[wang$sample=="thyroid_gland",2],gtex[gtex$sample=="Thyroid",2])) 
    prim[prim$organ=="tonsil","value"] = wang[wang$sample=="tonsil",2]
    prim[prim$organ=="urinary_bladder","value"] = wang[wang$sample=="urinary_bladder",2]
  }
  
  except = c("bronchus","caecum","cartilage","diaphragm","epididymis","left_atrium",
             "leukocyte","lymph_node","mitral_valve","nasal_pharynx","nose","parotid_gland",
             "penis","pleura","prefrontal_cortex","pulmonary_valve","renal_cortex",
             "seminal_vesicle","submandibular_gland","temporal_lobe","throat","tongue",
             "trachea","tricuspid_valve","vas_deferens")
  prim %<>% dplyr::filter(., !organ %in%  except)
  #ref: https://gotellilab.github.io/GotelliLabMeetingHacks/NickGotelli/ViridisColorPalette.html
  p1 = gganatogram(data=prim, fillOutline='white', organism='human', sex='male', fill="value") + 
    theme_void() + scale_fill_viridis_c(option="magma")
  #+guides(fill=guide_legend(title="X:A ratio"))
  #ggtitle("Dosage compensation state across 43 tissues")
  ggsave(p1,filename = "~/MamDC/Result/GTEx/Picture/Human.bodyplot.male.pdf",
         width = 4,height = 5)
  topptx(p1,"~/MamDC/Result/GTEx/Picture/Human.bodyplot.male.pptx",
         width = 4,height = 5)
}

# for human female
if (TRUE) {
  prim2 = hgFemale_key
  prim2$value = 100
  prim2$value = c(prim$value,head(prim$value,nrow(prim2)-nrow(prim)))
  
  prim2[prim2$organ=="ovary","value"] = mean(c(wang[wang$sample=="ovary",2],gtex[gtex$sample=="Ovary",2]))
  prim2[prim2$organ=="endometrium","value"] = wang[wang$sample=="endometrium",2]
  prim2[prim2$organ=="vagina","value"] = gtex[gtex$sample=="vagina",2]
  prim2[prim2$organ=="fallopian_tube","value"] = mean(c(wang[wang$sample=="fallopian_tube",2],gtex[gtex$sample=="Fallopian Tube",2]))
  prim2[prim2$organ=="uterus","value"] = gtex[gtex$sample=="Uterus",2]
  prim2[prim2$organ=="uterine_cervix","value"] = gtex[gtex$sample=="Cervix - Endocervix",2]
  prim2[prim2$organ=="ectocervix ","value"] = gtex[gtex$sample=="Cervix - Ectocervix",2]
  
  prim2 = prim2[c(which.max(prim2$value),which.min(prim2$value),which(prim2$organ %in% c("ovary","endometrium","vagina",
                                                                                         "fallopian_tube","uterus",
                                                                                         "uterine_cervix","ectocervix"))),]
  #ref: https://gotellilab.github.io/GotelliLabMeetingHacks/NickGotelli/ViridisColorPalette.html
  p3 = gganatogram(data=prim2, fillOutline='white', organism='human', sex='female', fill="value") + 
    theme_void() + scale_fill_viridis_c(option="magma")
  #+guides(fill=guide_legend(title="X:A ratio"))
  #ggtitle("Dosage compensation state across 43 tissues")
  ggsave(p3,filename = "~/MamDC/Result/GTEx/Picture/Human.bodyplot.female.pdf",
         width = 4,height = 5)
  topptx(p3,"~/MamDC/Result/GTEx/Picture/Human.bodyplot.male.pptx",
         width = 4,height = 5)
}

#for mouse male
rm(list = setdiff(ls(),"p1"))
if (TRUE) {
  prim = mmMale_key
  prim$value = 100
  nicho = read.csv("~/MamDC/Result/NicholasTony2020NatureRNAseq/Mouse/Savedata/Mouse_tissues_XAratio.txt",
                   sep = "\t",header = TRUE,stringsAsFactors = F)
  nicho = tapply(nicho$ratio, nicho$tissue, mean) %>% as.data.frame()
  nicho$sample = row.names(nicho)
  encode= read.csv("~/MamDC/Result/ENCODE/Savedata/Mouse_tissues_XAratio.txt",
                   sep = "\t",header = TRUE,stringsAsFactors = F)
  if (TRUE) {
    prim[prim$organ=="adrenal_gland","value"] = encode[encode$sample=="AdrenalGland",2]
    prim[prim$organ=="bone_marrow","value"] = nicho[nicho$sample=="Marrow",1]
    prim[prim$organ=="brain","value"] = mean(c(encode[encode$sample=="WholeBrain",2],nicho[nicho$sample=="Brain",1]))
    prim[prim$organ=="brown_adipose_tissue","value"] = nicho[nicho$sample=="BAT",1]
    prim[prim$organ=="colon","value"] = encode[encode$sample=="Colon",2]
    prim[prim$organ=="duodenum","value"] = encode[encode$sample=="Duodenum",2]
    prim[prim$organ=="femur","value"] = nicho[nicho$sample=="Bone",1]
    prim[prim$organ=="heart","value"] = mean(c(encode[encode$sample=="Heart",2],nicho[nicho$sample=="Heart",1]))
    prim[prim$organ=="hindlimb","value"] = encode[encode$sample=="Limb",2]
    prim[prim$organ=="kidney","value"] = mean(c(encode[encode$sample=="Kidney",2],nicho[nicho$sample=="Kidney",1]))
    prim[prim$organ=="liver","value"] = mean(c(encode[encode$sample=="Liver",2],nicho[nicho$sample=="Liver",1]))
    prim[prim$organ=="lung","value"] = mean(c(encode[encode$sample=="Lung",2],nicho[nicho$sample=="Lung",1]))
    prim[prim$organ=="pancreas","value"] = nicho[nicho$sample=="Pancreas",1]
    prim[prim$organ=="skin","value"] = nicho[nicho$sample=="Skin",1]
    prim[prim$organ=="small_intestine","value"] = mean(c(encode[encode$sample=="SmallIntestine",2],nicho[nicho$sample=="Small Intestine",1]))
    prim[prim$organ=="spleen","value"] = mean(c(encode[encode$sample=="Spleen",2],nicho[nicho$sample=="Spleen",1]))
    prim[prim$organ=="stomach","value"] = encode[encode$sample=="Stomach",2]
    prim[prim$organ=="testis","value"] = encode[encode$sample=="Testis",2]
    prim[prim$organ=="thymus","value"] = encode[encode$sample=="Thymus",2]
    prim[prim$organ=="urinary_bladder","value"] = encode[encode$sample=="UrinaryBladder",2]
  }
  except = c("aorta","blood_vessel","caecum","cartilage","circulatory_system","diaphragm","epididymis","esophagus","eye",
             "gall_bladder","ileum","intestinal_mucosa","jejunum","lymph_node","penis","peripheral_nervous_system","prostate_gland",
             "quadriceps_femoris","sciatic_nerve","seminal_vesicle","skeletal_muscle","spinal_cord","trachea","trigeminal_nerve","vas_deferens")
  prim %<>% dplyr::filter(., !organ %in%  except)
  p2 = gganatogram(data=prim, fillOutline='white', organism='mouse', sex='male', fill="value") + 
    theme_void() + scale_fill_viridis_c(option="magma")
  #+guides(fill=guide_legend(title="X:A ratio"))
  #ggtitle("Dosage compensation state across 43 tissues")
  topptx(p2,"~/MamDC/Result/ENCODE/Picture/Mouse.bodyplot.male.pptx",
         width = 4,height = 5)
  ggsave(p2,filename = "~/MamDC/Result/ENCODE/Picture/Mouse.bodyplot.male.pdf",
         width = 4,height = 5)
  p = p1+p2+plot_layout(ncol = 2, widths =  c(3, 2),heights = c(2,1))
  topptx(p,"~/MamDC/Result/GTEx/Picture/HumanMouse.bodyplot.pptx",
         width = 4,height = 5)
  ggsave("~/MamDC/Result/GTEx/Picture/HumanMouse.bodyplot.pdf",
         p,width = 6,height = 5)
}

#for mouse female
if (TRUE) {
  prim2 = mmFemale_key
  prim2$value = 100
  prim2$value = c(rep(prim$value,2),head(prim$value,nrow(prim2)-nrow(prim)*2))
  
  prim2[prim2$organ=="mammary_gland","value"] = encode[encode$sample=="MammaryGland",2]
  prim2[prim2$organ=="reproductive_system","value"] = encode[encode$sample=="Ovary",2]
  prim2 = prim2[c(which.max(prim2$value),which.min(prim2$value),which(prim2$organ %in% c("mammary_gland",
                                                                                         "reproductive_system"))),]
  p4 = gganatogram(data=prim2, fillOutline='white', organism='mouse', sex='female', fill="value") + 
    theme_void() + scale_fill_viridis_c(option="magma")
  #+guides(fill=guide_legend(title="X:A ratio"))
  #ggtitle("Dosage compensation state across 43 tissues")
  topptx(p4,"~/MamDC/Result/ENCODE/Picture/Mouse.bodyplot.female.pptx",
         width = 4,height = 5)
  ggsave(p4,filename = "~/MamDC/Result/ENCODE/Picture/Mouse.bodyplot.female.pdf",
         width = 4,height = 5)
}
  