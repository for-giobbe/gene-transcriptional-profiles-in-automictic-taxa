library(ggfortify)
library(data.table)
library(ggplot2)
library(gridExtra)

myCol <- c("#f9d47c","#00b4d7","#71ceea","#f78b88")

###################### GONADS ############################################################################################

BAT <- read.csv("../input_files/BAT_RSEM.TMM.EXPR.matrix.3sp.OGonly.csv")
BRO <- read.csv("../input_files/BRO_RSEM.TMM.EXPR.matrix.3sp.OGonly.csv")
BGM <- read.csv("../input_files/BGM_RSEM_mf.TMM.EXPR.matrix.3sp.OGonly.csv")

tmp <- merge(BAT,BRO,by="transcript")
df <- merge(tmp,BGM,by="transcript")

BRO <- read.table(file = "../input_files/new_DE_BGM_VS_BRO_single.lst", sep = " ", header= TRUE)
BAT <- read.table(file = "../input_files/new_DE_BGM_VS_BAT_single.lst", sep = " ", header= TRUE)

BGM <- (subset(BRO, BGM_F_padj<0.05 & BGM_F_logFC>1))$OG
BRO <- (subset(BRO, BRO_padj<0.05 & BRO_logFC>1))$OG
BAT <- (subset(BAT, BAT_padj<0.05 & BAT_logFC>1))$OG

tot <- c(BGM,BRO,BAT)

de <- df[df$transcript %in% tot,]

de_gonads <- de[,c(2,4,6,9,11,12,14,16,18,20,22,24,26,28,30,32,34,36,38,40,42,44,46,48)]
de_gonads <- as.data.frame(t(de_gonads))
rownames(de_gonads)

pca_res <- prcomp(de_gonads)

de_gonads$condition <- c("atticus","atticus","atticus","atticus","atticus","atticus",
                         "rossius","rossius","rossius","rossius","rossius","rossius",
                         "grandii female","grandii female","grandii female","grandii female","grandii female","grandii female",
                         "grandii male","grandii male","grandii male","grandii male","grandii male","grandii male")

gnds <- autoplot(pca_res, data =  de_gonads, colour ='condition', frame = TRUE, frame.type = 'norm') + 
  theme_light() + coord_fixed() + scale_color_manual(values = myCol) + theme(legend.position="none") + xlim(1,-1) + ylim(1,-1)

###################### LEGS ############################################################################################

BAT <- read.csv("../input_files/BAT_RSEM.TMM.EXPR.matrix.3sp.OGonly.csv")
BRO <- read.csv("../input_files/BRO_RSEM.TMM.EXPR.matrix.3sp.OGonly.csv")
BGM <- read.csv("../input_files/BGM_RSEM_mf.TMM.EXPR.matrix.3sp.OGonly.csv")

tmp <- merge(BAT,BRO,by="transcript")
df <- merge(tmp,BGM,by="transcript")

BRO <- read.table(file = "../input_files/new_DE_BGM_VS_BRO_single.lst", sep = " ", header= TRUE)
BAT <- read.table(file = "../input_files/new_DE_BGM_VS_BAT_single.lst", sep = " ", header= TRUE)

BGM <- (subset(BRO, BGM_F_padj<0.05 & BGM_F_logFC<1))$OG
BRO <- (subset(BRO, BRO_padj<0.05 & BRO_logFC<1))$OG
BAT <- (subset(BAT, BAT_padj<0.05 & BAT_logFC<1))$OG

tot <- c(BGM,BRO,BAT)

de <- df[df$transcript %in% tot,]

colnames(de)

de_legs <- de[,c(3,5,7,8,10,13,15,17,19,21,23,25,27,29,31,33,35,37,39,41,43,45,47,49)]
de_legs <- as.data.frame(t(de_legs))
rownames(de_legs)

pca_res <- prcomp(de_legs)

de_legs$condition <- c("atticus","atticus","atticus","atticus","atticus","atticus",
                         "rossius","rossius","rossius","rossius","rossius","rossius",
                         "grandii female","grandii female","grandii female","grandii female","grandii female","grandii female",
                         "grandii male","grandii male","grandii male","grandii male","grandii male","grandii male")

legs <- autoplot(pca_res, data =  de_legs, colour ='condition', frame = TRUE, frame.type = 'norm') + 
  theme_light() + coord_fixed() + scale_color_manual(values = myCol) + theme(legend.position="none") + xlim(1,-1) + ylim(1,-1)

###################### GONADSALL ############################################################################################

BAT <- read.csv("../input_files/BAT_RSEM.TMM.EXPR.matrix.3sp.OGonly.csv")
BRO <- read.csv("../input_files/BRO_RSEM.TMM.EXPR.matrix.3sp.OGonly.csv")
BGM <- read.csv("../input_files/BGM_RSEM_mf.TMM.EXPR.matrix.3sp.OGonly.csv")

tmp <- merge(BAT,BRO,by="transcript")
df <- merge(tmp,BGM,by="transcript")

de_gonads <- df[,c(2,4,6,9,11,12,14,16,18,20,22,24,26,28,30,32,34,36,38,40,42,44,46,48)]
de_gonads <- as.data.frame(t(de_gonads))
rownames(de_gonads)

pca_res <- prcomp(de_gonads)

de_gonads$condition <- c("atticus","atticus","atticus","atticus","atticus","atticus",
                         "rossius","rossius","rossius","rossius","rossius","rossius",
                         "grandii female","grandii female","grandii female","grandii female","grandii female","grandii female",
                         "grandii male","grandii male","grandii male","grandii male","grandii male","grandii male")

gndsall <- autoplot(pca_res, data =  de_gonads, colour ='condition', frame = TRUE, frame.type = 'norm') + 
  theme_light() + coord_fixed() + scale_color_manual(values = myCol) + theme(legend.position="none") + xlim(1.1,-1.1) + ylim(1.1,-1.1)


###################### LEGSALL ############################################################################################

BAT <- read.csv("../input_files/BAT_RSEM.TMM.EXPR.matrix.3sp.OGonly.csv")
BRO <- read.csv("../input_files/BRO_RSEM.TMM.EXPR.matrix.3sp.OGonly.csv")
BGM <- read.csv("../input_files/BGM_RSEM_mf.TMM.EXPR.matrix.3sp.OGonly.csv")

tmp <- merge(BAT,BRO,by="transcript")
df <- merge(tmp,BGM,by="transcript")

de_legs <- df[,c(3,5,7,8,10,13,15,17,19,21,23,25,27,29,31,33,35,37,39,41,43,45,47,49)]
de_legs <- as.data.frame(t(de_legs))
rownames(de_legs)

pca_res <- prcomp(de_legs)

de_legs$condition <- c("atticus","atticus","atticus","atticus","atticus","atticus",
                       "rossius","rossius","rossius","rossius","rossius","rossius",
                       "grandii female","grandii female","grandii female","grandii female","grandii female","grandii female",
                       "grandii male","grandii male","grandii male","grandii male","grandii male","grandii male")

legsall <- autoplot(pca_res, data =  de_legs, colour ='condition', frame = TRUE, frame.type = 'norm') + 
  theme_light() + coord_fixed() + scale_color_manual(values = myCol) + theme(legend.position="none") + xlim(1.1,-1.1) + ylim(1.1,-1.1)

###################### plot ############################################################################################

grid.arrange(gndsall,gnds,legsall,legs)

