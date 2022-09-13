library(ggfortify)
library(data.table)
library(ggplot2)
library(gridExtra)

myCol <- c("#f9d47c","#71ceea","blue","#f78b88")

###################### GONADS ############################################################################################

BAT <- read.csv("../intermediate_files//BAT_RSEM.TMM.EXPR.matrix.3sp.OGonly.csv")
BRO <- read.csv("../intermediate_files/BRO_RSEM.TMM.EXPR.matrix.3sp.OGonly.csv")
BGM <- read.csv("../intermediate_files/BGM_RSEM_mf.TMM.EXPR.matrix.3sp.OGonly.csv")

tmp <- merge(BAT,BRO,by="transcript")
df <- merge(tmp,BGM,by="transcript")

BRO <- read.table(file = "../intermediate_files/new_DE_BGM_VS_BRO_single.lst", sep = " ", header= TRUE)
BAT <- read.table(file = "../intermediate_files/new_DE_BGM_VS_BAT_single.lst", sep = " ", header= TRUE)

BGM <- (subset(BRO, BGM_F_padj<0.05 & BGM_F_logFC>1))$OG
BRO <- (subset(BRO, BRO_padj<0.05 & BRO_logFC>1))$OG
BAT <- (subset(BAT, BAT_padj<0.05 & BAT_logFC>1))$OG

tot <- c(BGM,BRO,BAT)

de <- df[df$transcript %in% tot,]

de_gonads <- de[,c(2,4,6,9,11,12,14,16,18,20,22,24,26,28,30,32,34,36,38,40,42,44,46,48)]
de_gonads <- as.data.frame(t(de_gonads))
rownames(de_gonads)

pca_res <- prcomp(de_gonads)

summ <- summary(pca_res)

de_gonads$condition <- c("atticus","atticus","atticus","atticus","atticus","atticus",
                         "rossius","rossius","rossius","rossius","rossius","rossius",
                         "grandii female","grandii female","grandii female","grandii female","grandii female","grandii female",
                         "grandii male","grandii male","grandii male","grandii male","grandii male","grandii male")

df <- cbind(pca_res$x[,1:2], de_gonads$condition) %>% as.data.frame()
df$V3 <- as.factor(df$V3)
df$PC1 <- as.numeric(df$PC1)
df$PC2 <- as.numeric(df$PC2)

gnds <- ggplot(df, aes(PC1, PC2, colour = V3)) +
  geom_point(size = 3) + scale_color_manual(values = myCol) + coord_fixed() +
  stat_ellipse(geom = "polygon", aes(fill = after_scale(alpha(colour, 0))), size = 0.75, type = "norm", linetype = 2) + theme_bw() + 
  theme(legend.position = "none") +
  labs(x=paste("PC1 (",round(summ$importance[2,1],2)*100,"%)", sep="")) + labs(y=paste("PC2 (",round(summ$importance[2,2],2)*100,"%)", sep="")) +
  theme(aspect.ratio=1) + theme(axis.ticks.x = element_blank(), axis.text.x = element_blank(),axis.ticks.y = element_blank(), axis.text.y = element_blank())


###################### LEGS ############################################################################################

BAT <- read.csv("../intermediate_files/BAT_RSEM.TMM.EXPR.matrix.3sp.OGonly.csv")
BRO <- read.csv("../intermediate_files/BRO_RSEM.TMM.EXPR.matrix.3sp.OGonly.csv")
BGM <- read.csv("../intermediate_files/BGM_RSEM_mf.TMM.EXPR.matrix.3sp.OGonly.csv")

tmp <- merge(BAT,BRO,by="transcript")
df <- merge(tmp,BGM,by="transcript")

BRO <- read.table(file = "../intermediate_files/new_DE_BGM_VS_BRO_single.lst", sep = " ", header= TRUE)
BAT <- read.table(file = "../intermediate_files/new_DE_BGM_VS_BAT_single.lst", sep = " ", header= TRUE)

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

summ <- summary(pca_res)

de_legs$condition <- c("atticus","atticus","atticus","atticus","atticus","atticus",
                         "rossius","rossius","rossius","rossius","rossius","rossius",
                         "grandii female","grandii female","grandii female","grandii female","grandii female","grandii female",
                         "grandii male","grandii male","grandii male","grandii male","grandii male","grandii male")

df <- cbind(pca_res$x[,1:2], de_legs$condition) %>% as.data.frame()
df$V3 <- as.factor(df$V3)
df$PC1 <- as.numeric(df$PC1)
df$PC2 <- as.numeric(df$PC2)

legs <- ggplot(df, aes(PC1, PC2, colour = V3)) +
  geom_point(size = 3) + scale_color_manual(values = myCol) + coord_fixed() +
  stat_ellipse(geom = "polygon", aes(fill = after_scale(alpha(colour, 0))), size = 0.75, type = "norm", linetype = 2) + theme_bw() + 
  theme(legend.position = "none") +
  labs(x=paste("PC1 (",round(summ$importance[2,1],2)*100,"%)", sep="")) + labs(y=paste("PC2 (",round(summ$importance[2,2],2)*100,"%)", sep="")) +
  theme(aspect.ratio=1) + theme(axis.ticks.x = element_blank(), axis.text.x = element_blank(),axis.ticks.y = element_blank(), axis.text.y = element_blank())

###################### GONADSALL ############################################################################################

BAT <- read.csv("../intermediate_files/BAT_RSEM.TMM.EXPR.matrix.3sp.OGonly.csv")
BRO <- read.csv("../intermediate_files/BRO_RSEM.TMM.EXPR.matrix.3sp.OGonly.csv")
BGM <- read.csv("../intermediate_files/BGM_RSEM_mf.TMM.EXPR.matrix.3sp.OGonly.csv")

tmp <- merge(BAT,BRO,by="transcript")
df <- merge(tmp,BGM,by="transcript")

de_gonads <- df[,c(2,4,6,9,11,12,14,16,18,20,22,24,26,28,30,32,34,36,38,40,42,44,46,48)]
de_gonads <- as.data.frame(t(de_gonads))
rownames(de_gonads)

pca_res <- prcomp(de_gonads)

summ <- summary(pca_res)

de_gonads$condition <- c("atticus","atticus","atticus","atticus","atticus","atticus",
                         "rossius","rossius","rossius","rossius","rossius","rossius",
                         "grandii female","grandii female","grandii female","grandii female","grandii female","grandii female",
                         "grandii male","grandii male","grandii male","grandii male","grandii male","grandii male")

df <- cbind(pca_res$x[,1:2], de_gonads$condition) %>% as.data.frame()
df$V3 <- as.factor(df$V3)
df$PC1 <- as.numeric(df$PC1)
df$PC2 <- as.numeric(df$PC2)

gndsall <- ggplot(df, aes(PC1, PC2, colour = V3)) +
  geom_point(size = 3) + scale_color_manual(values = myCol) + coord_fixed() +
  stat_ellipse(geom = "polygon", aes(fill = after_scale(alpha(colour, 0))), size = 0.75, type = "norm", linetype = 2) + theme_bw() + 
  theme(legend.position = "none") +
  labs(x=paste("PC1 (",round(summ$importance[2,1],2)*100,"%)", sep="")) + labs(y=paste("PC2 (",round(summ$importance[2,2],2)*100,"%)", sep="")) +
  theme(aspect.ratio=1) + theme(axis.ticks.x = element_blank(), axis.text.x = element_blank(),axis.ticks.y = element_blank(), axis.text.y = element_blank())

###################### LEGSALL ############################################################################################

BAT <- read.csv("../intermediate_files/BAT_RSEM.TMM.EXPR.matrix.3sp.OGonly.csv")
BRO <- read.csv("../intermediate_files/BRO_RSEM.TMM.EXPR.matrix.3sp.OGonly.csv")
BGM <- read.csv("../intermediate_files/BGM_RSEM_mf.TMM.EXPR.matrix.3sp.OGonly.csv")

tmp <- merge(BAT,BRO,by="transcript")
df <- merge(tmp,BGM,by="transcript")

de_legs <- df[,c(3,5,7,8,10,13,15,17,19,21,23,25,27,29,31,33,35,37,39,41,43,45,47,49)]
de_legs <- as.data.frame(t(de_legs))
rownames(de_legs)

pca_res <- prcomp(de_legs)

summ <- summary(pca_res)

de_legs$condition <- c("atticus","atticus","atticus","atticus","atticus","atticus",
                       "rossius","rossius","rossius","rossius","rossius","rossius",
                       "grandii female","grandii female","grandii female","grandii female","grandii female","grandii female",
                       "grandii male","grandii male","grandii male","grandii male","grandii male","grandii male")

df <- cbind(pca_res$x[,1:2], de_legs$condition) %>% as.data.frame()
df$V3 <- as.factor(df$V3)
df$PC1 <- as.numeric(df$PC1)
df$PC2 <- as.numeric(df$PC2)

legsall <- ggplot(df, aes(PC1, PC2, colour = V3)) +
  geom_point(size = 3) + scale_color_manual(values = myCol) + coord_fixed() +
  stat_ellipse(geom = "polygon", aes(fill = after_scale(alpha(colour, 0))), size = 0.75, type = "norm", linetype = 2) + theme_bw() + 
  theme(legend.position = "none") +
  labs(x=paste("PC1 (",round(summ$importance[2,1],2)*100,"%)", sep="")) + labs(y=paste("PC2 (",round(summ$importance[2,2],2)*100,"%)", sep="")) +
  theme(aspect.ratio=1) + theme(axis.ticks.x = element_blank(), axis.text.x = element_blank(),axis.ticks.y = element_blank(), axis.text.y = element_blank())

###################### plot ############################################################################################

grid.arrange(gndsall,gnds,legsall,legs)
