library(topGO)
library(tidyr)
library(tidyverse)
require(gridExtra)
 
geneID2GO <- readMappings(file = "../intermediate_files/def_GO_universe_nonredundant")
geneUniverse <- names(geneID2GO)

############################################################################################################      male_upregulated_sexual_asexual_upregulated ################## 

genesOfInterest <- read.table("../intermediate_files/upregulated_asex_&_upregulated_BGM_M.txt",header=FALSE)
genesOfInterest <- as.character(genesOfInterest$V1)
geneList <- factor(as.integer(geneUniverse %in% genesOfInterest))
names(geneList) <- geneUniverse

myGOdata <- new("topGOdata", description="My project", ontology="BP", allGenes=geneList,  annot = annFUN.gene2GO, gene2GO = geneID2GO, nodeSize = 1)
sp4_bgm_male_upregulated_convergent_shifts_weight_fisher <- runTest(myGOdata, algorithm="weight", statistic="fisher")
sp4_bgm_male_upregulated_convergent_shifts_elim_fisher <- runTest(myGOdata, algorithm="elim", statistic="fisher")
sp4_bgm_male_upregulated_convergent_shifts_allRes <- GenTable(myGOdata, 
                                                       Elim_Fisher = sp4_bgm_male_upregulated_convergent_shifts_elim_fisher, 
                                                       Weight_Fisher = sp4_bgm_male_upregulated_convergent_shifts_weight_fisher, 
                                                       orderBy = "Elim_Fisher", ranksOf = "Weight_Fisher", topNodes = 500)
sp4_bgm_male_upregulated_convergent_shifts_allRes <- subset(sp4_bgm_male_upregulated_convergent_shifts_allRes, Elim_Fisher<0.05)
write.table(sp4_bgm_male_upregulated_convergent_shifts_allRes, 'sex_male_upreg_asex_upreg_BP.csv'  , append = F, sep=',', row.names = F, quote = FALSE)

myGOdata <- new("topGOdata", description="My project", ontology="MF", allGenes=geneList,  annot = annFUN.gene2GO, gene2GO = geneID2GO, nodeSize = 10)
sp4_bgm_male_upregulated_convergent_shifts_weight_fisher <- runTest(myGOdata, algorithm="weight", statistic="fisher")
sp4_bgm_male_upregulated_convergent_shifts_elim_fisher <- runTest(myGOdata, algorithm="parentchild", statistic="fisher")
sp4_bgm_male_upregulated_convergent_shifts_allRes <- GenTable(myGOdata, 
                                                              Elim_Fisher = sp4_bgm_male_upregulated_convergent_shifts_elim_fisher, 
                                                              Weight_Fisher = sp4_bgm_male_upregulated_convergent_shifts_weight_fisher, 
                                                              orderBy = "Elim_Fisher", ranksOf = "Weight_Fisher", topNodes = 500)
sp4_bgm_male_upregulated_convergent_shifts_allRes <- subset(sp4_bgm_male_upregulated_convergent_shifts_allRes, Elim_Fisher<0.05)
write.table(sp4_bgm_male_upregulated_convergent_shifts_allRes, 'sex_male_upreg_asex_upreg_MF.csv', append = F, sep=',', row.names = F, quote = FALSE)

############################################################################################################      male_upregulated_sexual_asexual_upregulated_BRO ################## 

genesOfInterest <- read.table("../intermediate_files/upregulated_BRO_&_upregulated_BGM_M.txt",header=FALSE)
genesOfInterest <- as.character(genesOfInterest$V1)
geneList <- factor(as.integer(geneUniverse %in% genesOfInterest))
names(geneList) <- geneUniverse

myGOdata <- new("topGOdata", description="My project", ontology="BP", allGenes=geneList,  annot = annFUN.gene2GO, gene2GO = geneID2GO, nodeSize = 10)
sp4_bgm_male_upregulated_convergent_shifts_weight_fisher <- runTest(myGOdata, algorithm="weight", statistic="fisher")
sp4_bgm_male_upregulated_convergent_shifts_elim_fisher <- runTest(myGOdata, algorithm="parentchild", statistic="fisher")
sp4_bgm_male_upregulated_convergent_shifts_allRes <- GenTable(myGOdata, 
                                                              Elim_Fisher = sp4_bgm_male_upregulated_convergent_shifts_elim_fisher,
                                                              Weight_Fisher = sp4_bgm_male_upregulated_convergent_shifts_weight_fisher, 
                                                              orderBy = "Elim_Fisher", ranksOf = "Weight_Fisher", topNodes = 500)
sp4_bgm_male_upregulated_convergent_shifts_allRes <- subset(sp4_bgm_male_upregulated_convergent_shifts_allRes, Elim_Fisher<0.05)
write.table(sp4_bgm_male_upregulated_convergent_shifts_allRes, 'sex_male_upreg_BRO_upreg_BP.csv', append = F, sep=',', row.names = F, quote = FALSE)

myGOdata <- new("topGOdata", description="My project", ontology="MF", allGenes=geneList,  annot = annFUN.gene2GO, gene2GO = geneID2GO, nodeSize = 10)
sp4_bgm_male_upregulated_convergent_shifts_weight_fisher <- runTest(myGOdata, algorithm="weight", statistic="fisher")
sp4_bgm_male_upregulated_convergent_shifts_elim_fisher <- runTest(myGOdata, algorithm="parentchild", statistic="fisher")
sp4_bgm_male_upregulated_convergent_shifts_allRes <- GenTable(myGOdata, 
                                                              Elim_Fisher = sp4_bgm_male_upregulated_convergent_shifts_elim_fisher,
                                                              Weight_Fisher = sp4_bgm_male_upregulated_convergent_shifts_weight_fisher, 
                                                              orderBy = "Elim_Fisher", ranksOf = "Weight_Fisher", topNodes = 500)
sp4_bgm_male_upregulated_convergent_shifts_allRes <- subset(sp4_bgm_male_upregulated_convergent_shifts_allRes, Elim_Fisher<0.05)
write.table(sp4_bgm_male_upregulated_convergent_shifts_allRes, 'sex_male_upreg_BRO_upreg_MF.csv', append = F, sep=',', row.names = F, quote = FALSE)

############################################################################################################      male_upregulated_sexual_asexual_upregulated_BAT ################## 

genesOfInterest <- read.table("../intermediate_files/upregulated_BAT_&_upregulated_BGM_M.txt",header=FALSE)
genesOfInterest <- as.character(genesOfInterest$V1)
geneList <- factor(as.integer(geneUniverse %in% genesOfInterest))
names(geneList) <- geneUniverse

myGOdata <- new("topGOdata", description="My project", ontology="BP", allGenes=geneList,  annot = annFUN.gene2GO, gene2GO = geneID2GO, nodeSize = 10)
sp4_bgm_male_upregulated_convergent_shifts_weight_fisher <- runTest(myGOdata, algorithm="weight", statistic="fisher")
sp4_bgm_male_upregulated_convergent_shifts_elim_fisher <- runTest(myGOdata, algorithm="parentchild", statistic="fisher")
sp4_bgm_male_upregulated_convergent_shifts_allRes <- GenTable(myGOdata, 
                                                              Elim_Fisher = sp4_bgm_male_upregulated_convergent_shifts_elim_fisher, 
                                                              Weight_Fisher = sp4_bgm_male_upregulated_convergent_shifts_weight_fisher,
                                                              orderBy = "Elim_Fisher", ranksOf = "Weight_Fisher", topNodes = 500)
sp4_bgm_male_upregulated_convergent_shifts_allRes <- subset(sp4_bgm_male_upregulated_convergent_shifts_allRes, Elim_Fisher<0.05)
write.table(sp4_bgm_male_upregulated_convergent_shifts_allRes, 'sex_male_upreg_BAT_upreg_BP.csv', append = F, sep=',', row.names = F, quote = FALSE)

myGOdata <- new("topGOdata", description="My project", ontology="MF", allGenes=geneList,  annot = annFUN.gene2GO, gene2GO = geneID2GO, nodeSize = 10)
sp4_bgm_male_upregulated_convergent_shifts_weight_fisher <- runTest(myGOdata, algorithm="weight", statistic="fisher")
sp4_bgm_male_upregulated_convergent_shifts_elim_fisher <- runTest(myGOdata, algorithm="parentchild", statistic="fisher")
sp4_bgm_male_upregulated_convergent_shifts_allRes <- GenTable(myGOdata, 
                                                              Elim_Fisher = sp4_bgm_male_upregulated_convergent_shifts_elim_fisher, 
                                                              Weight_Fisher = sp4_bgm_male_upregulated_convergent_shifts_weight_fisher, 
                                                              orderBy = "Elim_Fisher", ranksOf = "Weight_Fisher", topNodes = 500)
sp4_bgm_male_upregulated_convergent_shifts_allRes <- subset(sp4_bgm_male_upregulated_convergent_shifts_allRes, Elim_Fisher<0.05)
write.table(sp4_bgm_male_upregulated_convergent_shifts_allRes, 'sex_male_upreg_BAT_upreg_MF.csv', append = F, sep=',', row.names = F, quote = FALSE)

############################################################################################################      not_upregulated_sexual_asexual_upregulated ################## 

genesOfInterest <- read.table("../intermediate_files/upregulated_asex_&_not_BGM_B.txt",header=FALSE)
genesOfInterest <- as.character(genesOfInterest$V1)
geneList <- factor(as.integer(geneUniverse %in% genesOfInterest))
names(geneList) <- geneUniverse

myGOdata <- new("topGOdata", description="My project", ontology="BP", allGenes=geneList,  annot = annFUN.gene2GO, gene2GO = geneID2GO, nodeSize = 10)
not_upregulated_sexual_asexual_upregulated_weight_fisher <- runTest(myGOdata, algorithm="weight", statistic="fisher")
not_upregulated_sexual_asexual_upregulated_elim_fisher <- runTest(myGOdata, algorithm="parentchild", statistic="fisher")
not_upregulated_sexual_asexual_upregulated_allRes <- GenTable(myGOdata, 
                                                              Elim_Fisher = not_upregulated_sexual_asexual_upregulated_elim_fisher,
                                                              Weight_Fisher = not_upregulated_sexual_asexual_upregulated_weight_fisher, 
                                                              orderBy = "Elim_Fisher", ranksOf = "Weight_Fisher", topNodes = 100)
not_upregulated_sexual_asexual_upregulated_allRes <- subset(not_upregulated_sexual_asexual_upregulated_allRes, Elim_Fisher<0.05)
write.table(not_upregulated_sexual_asexual_upregulated_allRes, 'sex_not_upreg_asex_upreg_BP.csv', append = F, sep=',', row.names = F, quote = FALSE)

myGOdata <- new("topGOdata", description="My project", ontology="MF", allGenes=geneList,  annot = annFUN.gene2GO, gene2GO = geneID2GO, nodeSize = 10)
not_upregulated_sexual_asexual_upregulated_weight_fisher <- runTest(myGOdata, algorithm="weight", statistic="fisher")
not_upregulated_sexual_asexual_upregulated_elim_fisher <- runTest(myGOdata, algorithm="parentchild", statistic="fisher")
not_upregulated_sexual_asexual_upregulated_allRes <- GenTable(myGOdata, 
                                                              Elim_Fisher = not_upregulated_sexual_asexual_upregulated_elim_fisher, 
                                                              Weight_Fisher = not_upregulated_sexual_asexual_upregulated_weight_fisher, 
                                                              orderBy = "Elim_Fisher", ranksOf = "Weight_Fisher", topNodes = 100)
not_upregulated_sexual_asexual_upregulated_allRes <- subset(not_upregulated_sexual_asexual_upregulated_allRes, Elim_Fisher<0.05)
write.table(not_upregulated_sexual_asexual_upregulated_allRes, 'sex_not_upreg_asex_upreg_MF.csv', append = F, sep=',', row.names = F, quote = FALSE)

############################################################################################################      not_upregulated_sexual_asexual_upregulated_BRO ################## 

genesOfInterest <- read.table("../intermediate_files/upregulated_BRO_&_not_BGM_B.txt",header=FALSE)
genesOfInterest <- as.character(genesOfInterest$V1)
geneList <- factor(as.integer(geneUniverse %in% genesOfInterest))
names(geneList) <- geneUniverse

myGOdata <- new("topGOdata", description="My project", ontology="BP", allGenes=geneList,  annot = annFUN.gene2GO, gene2GO = geneID2GO, nodeSize = 10)
not_upregulated_sexual_asexual_upregulated_weight_fisher <- runTest(myGOdata, algorithm="weight", statistic="fisher")
not_upregulated_sexual_asexual_upregulated_elim_fisher <- runTest(myGOdata, algorithm="parentchild", statistic="fisher")
not_upregulated_sexual_asexual_upregulated_allRes <- GenTable(myGOdata, 
                                                              Elim_Fisher = not_upregulated_sexual_asexual_upregulated_elim_fisher, 
                                                              Weight_Fisher = not_upregulated_sexual_asexual_upregulated_weight_fisher, 
                                                              orderBy = "Elim_Fisher", ranksOf = "Weight_Fisher", topNodes = 100)
not_upregulated_sexual_asexual_upregulated_allRes <- subset(not_upregulated_sexual_asexual_upregulated_allRes, Elim_Fisher<0.05)
write.table(not_upregulated_sexual_asexual_upregulated_allRes, 'sex_not_upreg_asex_upreg_BRO_BP.csv', append = F, sep=',', row.names = F, quote = FALSE)

myGOdata <- new("topGOdata", description="My project", ontology="MF", allGenes=geneList,  annot = annFUN.gene2GO, gene2GO = geneID2GO, nodeSize = 10)
not_upregulated_sexual_asexual_upregulated_weight_fisher <- runTest(myGOdata, algorithm="weight", statistic="fisher")
not_upregulated_sexual_asexual_upregulated_elim_fisher <- runTest(myGOdata, algorithm="parentchild", statistic="fisher")
not_upregulated_sexual_asexual_upregulated_allRes <- GenTable(myGOdata, 
                                                              Elim_Fisher = not_upregulated_sexual_asexual_upregulated_elim_fisher, 
                                                              Weight_Fisher = not_upregulated_sexual_asexual_upregulated_weight_fisher, 
                                                              orderBy = "Elim_Fisher", ranksOf = "Weight_Fisher", topNodes = 100)
not_upregulated_sexual_asexual_upregulated_allRes <- subset(not_upregulated_sexual_asexual_upregulated_allRes, Elim_Fisher<0.05)
write.table(not_upregulated_sexual_asexual_upregulated_allRes, 'sex_not_upreg_asex_upreg_BRO_MF.csv', append = F, sep=',', row.names = F, quote = FALSE)

############################################################################################################      not_upregulated_sexual_asexual_upregulated_BAT ################## 

genesOfInterest <- read.table("../intermediate_files/upregulated_BAT_&_not_BGM_B.txt",header=FALSE)
genesOfInterest <- as.character(genesOfInterest$V1)
geneList <- factor(as.integer(geneUniverse %in% genesOfInterest))
names(geneList) <- geneUniverse

myGOdata <- new("topGOdata", description="My project", ontology="BP", allGenes=geneList,  annot = annFUN.gene2GO, gene2GO = geneID2GO, nodeSize = 10)
not_upregulated_sexual_asexual_upregulated_weight_fisher <- runTest(myGOdata, algorithm="weight", statistic="fisher")
not_upregulated_sexual_asexual_upregulated_elim_fisher <- runTest(myGOdata, algorithm="parentchild", statistic="fisher")
not_upregulated_sexual_asexual_upregulated_allRes <- GenTable(myGOdata, 
                                                              Elim_Fisher = not_upregulated_sexual_asexual_upregulated_elim_fisher, 
                                                              Weight_Fisher = not_upregulated_sexual_asexual_upregulated_weight_fisher, 
                                                              orderBy = "Elim_Fisher", ranksOf = "Weight_Fisher", topNodes = 100)
not_upregulated_sexual_asexual_upregulated_allRes <- subset(not_upregulated_sexual_asexual_upregulated_allRes, Elim_Fisher<0.05)
write.table(not_upregulated_sexual_asexual_upregulated_allRes, 'sex_not_upreg_asex_upreg_BAT_BP.csv'  , append = F, sep=',', row.names = F, quote = FALSE)

myGOdata <- new("topGOdata", description="My project", ontology="MF", allGenes=geneList,  annot = annFUN.gene2GO, gene2GO = geneID2GO, nodeSize = 10)
not_upregulated_sexual_asexual_upregulated_weight_fisher <- runTest(myGOdata, algorithm="weight", statistic="fisher")
not_upregulated_sexual_asexual_upregulated_elim_fisher <- runTest(myGOdata, algorithm="parentchild", statistic="fisher")
not_upregulated_sexual_asexual_upregulated_allRes <- GenTable(myGOdata, 
                                                              Elim_Fisher = not_upregulated_sexual_asexual_upregulated_elim_fisher, 
                                                              Weight_Fisher = not_upregulated_sexual_asexual_upregulated_weight_fisher, 
                                                              orderBy = "Elim_Fisher", ranksOf = "Weight_Fisher", topNodes = 100)
not_upregulated_sexual_asexual_upregulated_allRes <- subset(not_upregulated_sexual_asexual_upregulated_allRes, Elim_Fisher<0.05)
write.table(not_upregulated_sexual_asexual_upregulated_allRes, 'sex_not_upreg_asex_upreg_BAT_MF.csv'  , append = F, sep=',', row.names = F, quote = FALSE)

############################################################################################################      gonad_upregulated_sexual_asexual_notregulated ################## 

genesOfInterest <- read.table("../intermediate_files/notregulated_asex_&_upregulated_BGM_B.txt",header=FALSE)
genesOfInterest <- as.character(genesOfInterest$V1)
geneList <- factor(as.integer(geneUniverse %in% genesOfInterest))
names(geneList) <- geneUniverse

myGOdata <- new("topGOdata", description="My project", ontology="BP", allGenes=geneList,  annot = annFUN.gene2GO, gene2GO = geneID2GO, nodeSize = 10)
sp4_bgm_not_upregulated_convergent_shifts_weight_fisher <- runTest(myGOdata, algorithm="weight", statistic="fisher")
sp4_bgm_not_upregulated_convergent_shifts_elim_fisher <- runTest(myGOdata, algorithm="parentchild", statistic="fisher")
sp4_bgm_not_upregulated_convergent_shifts_allRes <- GenTable(myGOdata, 
                                                             Elim_Fisher = sp4_bgm_not_upregulated_convergent_shifts_elim_fisher, 
                                                             Weight_Fisher = sp4_bgm_not_upregulated_convergent_shifts_weight_fisher, 
                                                             orderBy = "Elim_Fisher", ranksOf = "Weight_Fisher", topNodes = 100)
sp4_bgm_not_upregulated_convergent_shifts_allRes <- subset(sp4_bgm_not_upregulated_convergent_shifts_allRes, Elim_Fisher<0.05)
write.table(sp4_bgm_not_upregulated_convergent_shifts_allRes, 'sex_gonad_upreg_asex_not_BP.csv', append = F, sep=',', row.names = F, quote = FALSE)

myGOdata <- new("topGOdata", description="My project", ontology="MF", allGenes=geneList,  annot = annFUN.gene2GO, gene2GO = geneID2GO, nodeSize = 10)
sp4_bgm_not_upregulated_convergent_shifts_weight_fisher <- runTest(myGOdata, algorithm="weight", statistic="fisher")
sp4_bgm_not_upregulated_convergent_shifts_elim_fisher <- runTest(myGOdata, algorithm="parentchild", statistic="fisher")
sp4_bgm_not_upregulated_convergent_shifts_allRes <- GenTable(myGOdata, 
                                                             Elim_Fisher = sp4_bgm_not_upregulated_convergent_shifts_elim_fisher, 
                                                             Weight_Fisher = sp4_bgm_not_upregulated_convergent_shifts_weight_fisher, 
                                                             orderBy = "Elim_Fisher", ranksOf = "Weight_Fisher", topNodes = 100)

sp4_bgm_not_upregulated_convergent_shifts_allRes <- subset(sp4_bgm_not_upregulated_convergent_shifts_allRes, Elim_Fisher<0.05)
write.table(sp4_bgm_not_upregulated_convergent_shifts_allRes, 'sex_gonad_upreg_asex_not_MF.csv', append = F, sep=',', row.names = F, quote = FALSE)

############################################################################################################      gonad_upregulated_sexual_asexual_notregulated_BRO ################## 

genesOfInterest <- read.table("../intermediate_files/notregulated_BRO_&_upregulated_BGM_B.txt",header=FALSE)
genesOfInterest <- as.character(genesOfInterest$V1)
geneList <- factor(as.integer(geneUniverse %in% genesOfInterest))
names(geneList) <- geneUniverse

myGOdata <- new("topGOdata", description="My project", ontology="BP", allGenes=geneList,  annot = annFUN.gene2GO, gene2GO = geneID2GO, nodeSize = 10)
sp4_bgm_not_upregulated_convergent_shifts_weight_fisher <- runTest(myGOdata, algorithm="weight", statistic="fisher")
sp4_bgm_not_upregulated_convergent_shifts_elim_fisher <- runTest(myGOdata, algorithm="parentchild", statistic="fisher")
sp4_bgm_not_upregulated_convergent_shifts_allRes <- GenTable(myGOdata, 
                                                             Elim_Fisher = sp4_bgm_not_upregulated_convergent_shifts_elim_fisher, 
                                                             Weight_Fisher = sp4_bgm_not_upregulated_convergent_shifts_weight_fisher, 
                                                             orderBy = "Elim_Fisher", ranksOf = "Weight_Fisher", topNodes = 100)
sp4_bgm_not_upregulated_convergent_shifts_allRes <- subset(sp4_bgm_not_upregulated_convergent_shifts_allRes, Elim_Fisher<0.05)
write.table(sp4_bgm_not_upregulated_convergent_shifts_allRes, 'sex_gonad_upreg_asex_not_BRO_BP.csv', append = F, sep=',', row.names = F, quote = FALSE)

myGOdata <- new("topGOdata", description="My project", ontology="MF", allGenes=geneList,  annot = annFUN.gene2GO, gene2GO = geneID2GO, nodeSize = 10)
sp4_bgm_not_upregulated_convergent_shifts_weight_fisher <- runTest(myGOdata, algorithm="weight", statistic="fisher")
sp4_bgm_not_upregulated_convergent_shifts_elim_fisher <- runTest(myGOdata, algorithm="parentchild", statistic="fisher")
sp4_bgm_not_upregulated_convergent_shifts_allRes <- GenTable(myGOdata, 
                                                             Elim_Fisher = sp4_bgm_not_upregulated_convergent_shifts_elim_fisher, 
                                                             Weight_Fisher = sp4_bgm_not_upregulated_convergent_shifts_weight_fisher, 
                                                             orderBy = "Elim_Fisher", ranksOf = "Weight_Fisher", topNodes = 100)

sp4_bgm_not_upregulated_convergent_shifts_allRes <- subset(sp4_bgm_not_upregulated_convergent_shifts_allRes, Elim_Fisher<0.05)
write.table(sp4_bgm_not_upregulated_convergent_shifts_allRes, 'sex_gonad_upreg_asex_not_BRO_MF.csv', append = F, sep=',', row.names = F, quote = FALSE)

############################################################################################################      gonad_upregulated_sexual_asexual_notregulated_BAT ################## 

genesOfInterest <- read.table("../intermediate_files/notregulated_BAT_&_upregulated_BGM_B.txt",header=FALSE)
genesOfInterest <- as.character(genesOfInterest$V1)
geneList <- factor(as.integer(geneUniverse %in% genesOfInterest))
names(geneList) <- geneUniverse

myGOdata <- new("topGOdata", description="My project", ontology="BP", allGenes=geneList,  annot = annFUN.gene2GO, gene2GO = geneID2GO, nodeSize = 10)
sp4_bgm_not_upregulated_convergent_shifts_weight_fisher <- runTest(myGOdata, algorithm="weight", statistic="fisher")
sp4_bgm_not_upregulated_convergent_shifts_elim_fisher <- runTest(myGOdata, algorithm="parentchild", statistic="fisher")
sp4_bgm_not_upregulated_convergent_shifts_allRes <- GenTable(myGOdata, 
                                                             Elim_Fisher = sp4_bgm_not_upregulated_convergent_shifts_elim_fisher, 
                                                             Weight_Fisher = sp4_bgm_not_upregulated_convergent_shifts_weight_fisher, 
                                                             orderBy = "Elim_Fisher", ranksOf = "Weight_Fisher", topNodes = 100)
sp4_bgm_not_upregulated_convergent_shifts_allRes <- subset(sp4_bgm_not_upregulated_convergent_shifts_allRes, Elim_Fisher<0.05)
write.table(sp4_bgm_not_upregulated_convergent_shifts_allRes, 'sex_gonad_upreg_asex_not_BAT_BP.csv', append = F, sep=',', row.names = F, quote = FALSE)

myGOdata <- new("topGOdata", description="My project", ontology="MF", allGenes=geneList,  annot = annFUN.gene2GO, gene2GO = geneID2GO, nodeSize = 10)
sp4_bgm_not_upregulated_convergent_shifts_weight_fisher <- runTest(myGOdata, algorithm="weight", statistic="fisher")
sp4_bgm_not_upregulated_convergent_shifts_elim_fisher <- runTest(myGOdata, algorithm="parentchild", statistic="fisher")
sp4_bgm_not_upregulated_convergent_shifts_allRes <- GenTable(myGOdata, 
                                                             Elim_Fisher = sp4_bgm_not_upregulated_convergent_shifts_elim_fisher, 
                                                             Weight_Fisher = sp4_bgm_not_upregulated_convergent_shifts_weight_fisher, 
                                                             orderBy = "Elim_Fisher", ranksOf = "Weight_Fisher", topNodes = 100)

sp4_bgm_not_upregulated_convergent_shifts_allRes <- subset(sp4_bgm_not_upregulated_convergent_shifts_allRes, Elim_Fisher<0.05)
write.table(sp4_bgm_not_upregulated_convergent_shifts_allRes, 'sex_gonad_upreg_asex_not_BAT_MF.csv', append = F, sep=',', row.names = F, quote = FALSE)

############################################################################################################      female_upregulated_sexual_asexual_notregulated ################## 

genesOfInterest <- read.table("../intermediate_files/notregulated_asex_&_upregulated_BGM_F.txt",header=FALSE)
genesOfInterest <- as.character(genesOfInterest$V1)
geneList <- factor(as.integer(geneUniverse %in% genesOfInterest))
names(geneList) <- geneUniverse

myGOdata <- new("topGOdata", description="My project", ontology="BP", allGenes=geneList,  annot = annFUN.gene2GO, gene2GO = geneID2GO, nodeSize = 10)
myGOdata
sp4_bgm_not_upregulated_convergent_shifts_weight_fisher <- runTest(myGOdata, algorithm="weight", statistic="fisher")
sp4_bgm_not_upregulated_convergent_shifts_elim_fisher <- runTest(myGOdata, algorithm="parentchild", statistic="fisher")
sp4_bgm_not_upregulated_convergent_shifts_allRes <- GenTable(myGOdata, 
                                                             Elim_Fisher = sp4_bgm_not_upregulated_convergent_shifts_elim_fisher, 
                                                             Weight_Fisher = sp4_bgm_not_upregulated_convergent_shifts_weight_fisher, 
                                                             orderBy = "Elim_Fisher", ranksOf = "Weight_Fisher", topNodes = 100)
sp4_bgm_not_upregulated_convergent_shifts_allRes <- subset(sp4_bgm_not_upregulated_convergent_shifts_allRes, Elim_Fisher<0.05)
write.table(sp4_bgm_not_upregulated_convergent_shifts_allRes, 'sex_female_upreg_asex_not_BP.csv', append = F, sep=',', row.names = F, quote = FALSE)

myGOdata <- new("topGOdata", description="My project", ontology="MF", allGenes=geneList,  annot = annFUN.gene2GO, gene2GO = geneID2GO, nodeSize = 10)
myGOdata
sp4_bgm_not_upregulated_convergent_shifts_weight_fisher <- runTest(myGOdata, algorithm="weight", statistic="fisher")
sp4_bgm_not_upregulated_convergent_shifts_elim_fisher <- runTest(myGOdata, algorithm="parentchild", statistic="fisher")
sp4_bgm_not_upregulated_convergent_shifts_allRes <- GenTable(myGOdata, 
                                                             Elim_Fisher = sp4_bgm_not_upregulated_convergent_shifts_elim_fisher, 
                                                             Weight_Fisher = sp4_bgm_not_upregulated_convergent_shifts_weight_fisher, 
                                                             orderBy = "Elim_Fisher", ranksOf = "Weight_Fisher", topNodes = 100)
sp4_bgm_not_upregulated_convergent_shifts_allRes <- subset(sp4_bgm_not_upregulated_convergent_shifts_allRes, Elim_Fisher<0.05)
write.table(sp4_bgm_not_upregulated_convergent_shifts_allRes, 'sex_female_upreg_asex_not_MF.csv', append = F, sep=',', row.names = F, quote = FALSE)

############################################################################################################      female_upregulated_sexual_asexual_notregulated_BRO ################## 

genesOfInterest <- read.table("../intermediate_files/notregulated_BRO_&_upregulated_BGM_F.txt",header=FALSE)
genesOfInterest <- as.character(genesOfInterest$V1)
geneList <- factor(as.integer(geneUniverse %in% genesOfInterest))
names(geneList) <- geneUniverse

myGOdata <- new("topGOdata", description="My project", ontology="BP", allGenes=geneList,  annot = annFUN.gene2GO, gene2GO = geneID2GO, nodeSize = 10)
myGOdata
sp4_bgm_not_upregulated_convergent_shifts_weight_fisher <- runTest(myGOdata, algorithm="weight", statistic="fisher")
sp4_bgm_not_upregulated_convergent_shifts_elim_fisher <- runTest(myGOdata, algorithm="parentchild", statistic="fisher")
sp4_bgm_not_upregulated_convergent_shifts_allRes <- GenTable(myGOdata, 
                                                             Elim_Fisher = sp4_bgm_not_upregulated_convergent_shifts_elim_fisher, 
                                                             Weight_Fisher = sp4_bgm_not_upregulated_convergent_shifts_weight_fisher, 
                                                             orderBy = "Elim_Fisher", ranksOf = "Weight_Fisher", topNodes = 100)
sp4_bgm_not_upregulated_convergent_shifts_allRes <- subset(sp4_bgm_not_upregulated_convergent_shifts_allRes, Elim_Fisher<0.05)
write.table(sp4_bgm_not_upregulated_convergent_shifts_allRes, 'sex_female_upreg_asex_not_BRO_BP.csv', append = F, sep=',', row.names = F, quote = FALSE)

myGOdata <- new("topGOdata", description="My project", ontology="MF", allGenes=geneList,  annot = annFUN.gene2GO, gene2GO = geneID2GO, nodeSize = 10)
myGOdata
sp4_bgm_not_upregulated_convergent_shifts_weight_fisher <- runTest(myGOdata, algorithm="weight", statistic="fisher")
sp4_bgm_not_upregulated_convergent_shifts_elim_fisher <- runTest(myGOdata, algorithm="parentchild", statistic="fisher")
sp4_bgm_not_upregulated_convergent_shifts_allRes <- GenTable(myGOdata, 
                                                             Elim_Fisher = sp4_bgm_not_upregulated_convergent_shifts_elim_fisher, 
                                                             Weight_Fisher = sp4_bgm_not_upregulated_convergent_shifts_weight_fisher, 
                                                             orderBy = "Elim_Fisher", ranksOf = "Weight_Fisher", topNodes = 100)
sp4_bgm_not_upregulated_convergent_shifts_allRes <- subset(sp4_bgm_not_upregulated_convergent_shifts_allRes, Elim_Fisher<0.05)
write.table(sp4_bgm_not_upregulated_convergent_shifts_allRes, 'sex_female_upreg_asex_not_BRO_MF.csv', append = F, sep=',', row.names = F, quote = FALSE)

############################################################################################################      female_upregulated_sexual_asexual_notregulated_BAT ################## 

genesOfInterest <- read.table("../intermediate_files/notregulated_BAT_&_upregulated_BGM_F.txt",header=FALSE)
genesOfInterest <- as.character(genesOfInterest$V1)
geneList <- factor(as.integer(geneUniverse %in% genesOfInterest))
names(geneList) <- geneUniverse

myGOdata <- new("topGOdata", description="My project", ontology="BP", allGenes=geneList,  annot = annFUN.gene2GO, gene2GO = geneID2GO, nodeSize = 10)
myGOdata
sp4_bgm_not_upregulated_convergent_shifts_weight_fisher <- runTest(myGOdata, algorithm="weight", statistic="fisher")
sp4_bgm_not_upregulated_convergent_shifts_elim_fisher <- runTest(myGOdata, algorithm="parentchild", statistic="fisher")
sp4_bgm_not_upregulated_convergent_shifts_allRes <- GenTable(myGOdata, 
                                                             Elim_Fisher = sp4_bgm_not_upregulated_convergent_shifts_elim_fisher, 
                                                             Weight_Fisher = sp4_bgm_not_upregulated_convergent_shifts_weight_fisher, 
                                                             orderBy = "Elim_Fisher", ranksOf = "Weight_Fisher", topNodes = 100)
sp4_bgm_not_upregulated_convergent_shifts_allRes <- subset(sp4_bgm_not_upregulated_convergent_shifts_allRes, Elim_Fisher<0.05)
write.table(sp4_bgm_not_upregulated_convergent_shifts_allRes, 'sex_female_upreg_asex_not_BAT_BP.csv', append = F, sep=',', row.names = F, quote = FALSE)

myGOdata <- new("topGOdata", description="My project", ontology="MF", allGenes=geneList,  annot = annFUN.gene2GO, gene2GO = geneID2GO, nodeSize = 10)
myGOdata
sp4_bgm_not_upregulated_convergent_shifts_weight_fisher <- runTest(myGOdata, algorithm="weight", statistic="fisher")
sp4_bgm_not_upregulated_convergent_shifts_elim_fisher <- runTest(myGOdata, algorithm="parentchild", statistic="fisher")
sp4_bgm_not_upregulated_convergent_shifts_allRes <- GenTable(myGOdata, 
                                                             Elim_Fisher = sp4_bgm_not_upregulated_convergent_shifts_elim_fisher, 
                                                             Weight_Fisher = sp4_bgm_not_upregulated_convergent_shifts_weight_fisher, 
                                                             orderBy = "Elim_Fisher", ranksOf = "Weight_Fisher", topNodes = 100)
sp4_bgm_not_upregulated_convergent_shifts_allRes <- subset(sp4_bgm_not_upregulated_convergent_shifts_allRes, Elim_Fisher<0.05)
write.table(sp4_bgm_not_upregulated_convergent_shifts_allRes, 'sex_female_upreg_asex_not_BAT_MF.csv', append = F, sep=',', row.names = F, quote = FALSE)

############################################################################################################      positive_selection BRO ################## 

genesOfInterest <- read.table("../intermediate_files/BRO_branchsite_genes.lst",header=FALSE)
BRO <-read.table(file = "../intermediate_files/BRO_RSEM.TMM.EXPR.matrix.4sp.reformat.og_gonad_only", sep = " ", header= TRUE)
BRO$BRO_mean <- rowMeans(BRO[, 2:7])
BRO_exp <- (subset(BRO, BRO_mean > 10)$OG)
genesOfInterest <- intersect(BRO_exp,genesOfInterest$V1)

myGOdata <- new("topGOdata", description="My project", ontology="BP", allGenes=geneList,  annot = annFUN.gene2GO, gene2GO = geneID2GO, nodeSize = 10)
positive_selection_BRO_weight_fisher <- runTest(myGOdata, algorithm="weight", statistic="fisher")
positive_selection_BRO_elim_fisher <- runTest(myGOdata, algorithm="parentchild", statistic="fisher")
positive_selection_BRO_allRes <- GenTable(myGOdata, 
                                          Elim_Fisher = positive_selection_BRO_elim_fisher, 
                                          Weight_Fisher = positive_selection_BRO_weight_fisher, 
                                          orderBy = "Elim_Fisher", ranksOf = "Weight_Fisher", topNodes = 100)

positive_selection_BRO_allRes <- subset(positive_selection_BRO_allRes, Elim_Fisher<0.05)
write.table(positive_selection_BRO_allRes, 'positive_selection_BRO_BP.csv', append = F, sep=',', row.names = F, quote = FALSE)

myGOdata <- new("topGOdata", description="My project", ontology="MF", allGenes=geneList,  annot = annFUN.gene2GO, gene2GO = geneID2GO, nodeSize = 10)
positive_selection_BRO_weight_fisher <- runTest(myGOdata, algorithm="weight", statistic="fisher")
positive_selection_BRO_elim_fisher <- runTest(myGOdata, algorithm="parentchild", statistic="fisher")
positive_selection_BRO_allRes <- GenTable(myGOdata, 
                                          Elim_Fisher = positive_selection_BRO_elim_fisher, 
                                          Weight_Fisher = positive_selection_BRO_weight_fisher, 
                                          orderBy = "Elim_Fisher", ranksOf = "Weight_Fisher", topNodes = 100)

positive_selection_BRO_allRes <- subset(positive_selection_BRO_allRes, Elim_Fisher<0.05)
write.table(positive_selection_BRO_allRes, 'positive_selection_BRO_MF.csv', append = F, sep=',', row.names = F, quote = FALSE)

############################################################################################################      positive_selection_BAT ################## 

genesOfInterest <- read.table("../intermediate_files/BAT_branchsite_genes.lst",header=FALSE)
BAT <-read.table(file = "../intermediate_files/BAT_RSEM.TMM.EXPR.matrix.4sp.reformat.og_gonad_only", sep = " ", header= TRUE)
BAT$BAT_mean <- rowMeans(BAT[, 2:7])
BAT_exp <- (subset(BAT, BAT_mean > 10)$OG)
genesOfInterest <- intersect(BAT_exp,genesOfInterest$V1)

geneList <- factor(as.integer(geneUniverse %in% genesOfInterest))
names(geneList) <- geneUniverse

myGOdata <- new("topGOdata", description="My project", ontology="BP", allGenes=geneList,  annot = annFUN.gene2GO, gene2GO = geneID2GO, nodeSize = 10)
positive_selection_BAT_weight_fisher <- runTest(myGOdata, algorithm="weight", statistic="fisher")
positive_selection_BAT_elim_fisher <- runTest(myGOdata, algorithm="parentchild", statistic="fisher")
positive_selection_BAT_allRes <- GenTable(myGOdata, 
                                                             Elim_Fisher = positive_selection_BAT_elim_fisher, 
                                                             Weight_Fisher = positive_selection_BAT_weight_fisher, 
                                                             orderBy = "Elim_Fisher", ranksOf = "Weight_Fisher", topNodes = 100)

positive_selection_BAT_allRes <- subset(positive_selection_BAT_allRes, Elim_Fisher<0.05)
write.table(positive_selection_BAT_allRes, 'positive_selection_BAT_BP.csv', append = F, sep=',', row.names = F, quote = FALSE)

myGOdata <- new("topGOdata", description="My project", ontology="MF", allGenes=geneList,  annot = annFUN.gene2GO, gene2GO = geneID2GO, nodeSize = 10)
positive_selection_BAT_weight_fisher <- runTest(myGOdata, algorithm="weight", statistic="fisher")
positive_selection_BAT_elim_fisher <- runTest(myGOdata, algorithm="parentchild", statistic="fisher")
positive_selection_BAT_allRes <- GenTable(myGOdata, 
                                          Elim_Fisher = positive_selection_BAT_elim_fisher, 
                                          Weight_Fisher = positive_selection_BAT_weight_fisher, 
                                          orderBy = "Elim_Fisher", ranksOf = "Weight_Fisher", topNodes = 100)

positive_selection_BAT_allRes <- subset(positive_selection_BAT_allRes, Elim_Fisher<0.05)
write.table(positive_selection_BAT_allRes, 'positive_selection_BAT_MF.csv', append = F, sep=',', row.names = F, quote = FALSE)

############################################################################################################      positive_selection CONVERGENT ################## 

BAT <-read.table(file = "../intermediate_files/BAT_RSEM.TMM.EXPR.matrix.4sp.reformat.og_gonad_only", sep = " ", header= TRUE)
BRO <-read.table(file = "../intermediate_files/BRO_RSEM.TMM.EXPR.matrix.4sp.reformat.og_gonad_only", sep = " ", header= TRUE)
BGM <-read.table(file = "../intermediate_files/BGM_RSEM_mf.TMM.EXPR.matrix.4sp.reformat.ogonly", sep = " ", header= TRUE)
BGM_G = subset(BGM, select = -c(BGM_F_L_RNA10,BGM_F_L_RNA12,BGM_F_L_RNA14,BGM_F_L_RNA15,BGM_F_L_RNA17,BGM_F_L_RNA19,BGM_M_L_RNA04,BGM_M_L_RNA06,BGM_M_L_RNA08,BGM_M_L_RNA20,BGM_M_L_RNA21,BGM_M_L_RNA22) )
BGM_F = subset(BGM_G, select = -c(BGM_M_G_RNA04,BGM_M_G_RNA06,BGM_M_G_RNA08,BGM_M_G_RNA20,BGM_M_G_RNA21,BGM_M_G_RNA22))
BGM_M = subset(BGM_G, select = -c(BGM_F_G_RNA10,BGM_F_G_RNA12,BGM_F_G_RNA14,BGM_F_G_RNA15,BGM_F_G_RNA17,BGM_F_G_RNA19))

curated_positive_selection <- read.table("../intermediate_files/curated_positive_selection_gblocked.txt",header=TRUE, sep =",")

curated_positive_selection_BAT <- merge(curated_positive_selection, BAT, by="OG")
curated_positive_selection_BAT$BAT_mean <- rowMeans(curated_positive_selection_BAT[, 4:9])
curated_positive_selection_BAT$BAT_sd <- apply(curated_positive_selection_BAT[,4:9],1,sd)

curated_positive_selection_BRO <- merge(curated_positive_selection, BRO, by="OG")
curated_positive_selection_BRO$BRO_mean <- rowMeans(curated_positive_selection_BRO[, 4:9])
curated_positive_selection_BRO$BRO_sd <- apply(curated_positive_selection_BRO[,4:9],1,sd)

curated_positive_selection_BGM_F <- merge(curated_positive_selection, BGM_F, by="OG")
curated_positive_selection_BGM_F$BGM_F_mean <- rowMeans(curated_positive_selection_BGM_F[, 4:9])
curated_positive_selection_BGM_F$BGM_F_sd <- apply(curated_positive_selection_BGM_F[,4:9],1,sd)

curated_positive_selection_BGM_M <- merge(curated_positive_selection, BGM_M, by="OG")
curated_positive_selection_BGM_M$BGM_M_mean <- rowMeans(curated_positive_selection_BGM_M[, 4:9])
curated_positive_selection_BGM_M$BGM_M_sd <- apply(curated_positive_selection_BGM_M[,4:9],1,sd)

#

mf <-read.table(file = "../intermediate_files/mf_tab_reduced.txt", sep = " ", header= TRUE)
limits <- aes(ymax = mean + sd, ymin = mean - sd)

ptmm <- ggplot(mf, aes(x=OG, y=mean, width=0.7, fill=sp)) + geom_col(position = "dodge") +
  geom_errorbar(limits, position = position_dodge(0.7), width = 0.25, size = 0.6) + coord_flip() + 
  theme_minimal() +   theme(axis.text.y = element_text(angle = 20, hjust = 1, vjust=0.3, size=8)) + scale_fill_manual(values = c("#ffc037", "#f55a4f", "#BFE1B0","#39A96B") ) +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black", size = 0.6)) + 
  theme(aspect.ratio=1) + theme(legend.position = "none") + theme(axis.title.y=element_blank()) 

psites <- ggplot(curated_positive_selection, aes(x=OG, y=sites, width=0.7, fill=species)) + geom_col(position = "dodge") + coord_flip() + 
          theme_minimal() +   theme(axis.text.y = element_text(angle = 20, hjust = 1, vjust=0.3, size=8)) + scale_fill_manual(values = c("#ffc037", "#f55a4f") ) +
          theme(panel.border = element_blank(), panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), axis.line = element_line(colour = "black", size = 0.6)) + 
          theme(aspect.ratio=1) + theme(legend.position = "none") + theme(axis.title.y=element_blank()) 

grid.arrange(psites, ptmm, ncol = 2)


