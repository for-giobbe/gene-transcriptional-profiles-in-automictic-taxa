require(WGCNA)
require(data.table)
require(tidyverse)
require(grid)
require(gridExtra)
require(ggplot2)
require(dplyr)    
require(tidyr)
require(dplyr)
require(networkD3)
require(VennDiagram)
require(GeneOverlap)


##########################################################################################################################  final BRO easy sankey ############################

tot <- read.table(file = "../inputs/new_DE_BGM_VS_BRO_single.lst", sep = " ", header= TRUE)

dim(tot)

#   non DE in both                                        DE in both                                            DE in M                                               DE in F                                               
dim(subset(tot, BGM_M_padj>0.05 & BGM_F_padj>0.05)) + dim(subset(tot, BGM_M_padj<0.05 & BGM_F_padj<0.05)) + dim(subset(tot, BGM_M_padj<0.05 & BGM_F_padj>0.05)) + dim(subset(tot, BGM_M_padj>0.05 & BGM_F_padj<0.05))
#   non DE in both                                        DN in both                                                                                DN in M                                                               DN in F                                                                UP in M                                                                UP in F                                                                UP in both
dim(subset(tot, BGM_M_padj>0.05 & BGM_F_padj>0.05)) + dim(subset(tot, BGM_M_padj<0.05 & BGM_M_logFC<1 & BGM_F_padj<0.05 & BGM_F_logFC<1)) + dim(subset(tot, BGM_M_padj<0.05 & BGM_M_logFC<1 & BGM_F_padj>0.05)) + dim(subset(tot, BGM_M_padj>0.05 & BGM_F_padj<0.05 & BGM_F_logFC<1 )) + dim(subset(tot, BGM_F_padj>0.05  & BGM_M_padj<0.05 & BGM_M_logFC>1)) + dim(subset(tot, BGM_M_padj>0.05  & BGM_F_padj<0.05 & BGM_F_logFC>1)) + dim(subset(tot, BGM_M_padj<0.05 & BGM_M_logFC>1 & BGM_F_padj<0.05 & BGM_F_logFC>1)) + dim(subset(tot, BGM_M_padj<0.05 & BGM_M_logFC>1 & BGM_F_padj<0.05 & BGM_F_logFC<1)) + dim(subset(tot, BGM_M_padj<0.05 & BGM_M_logFC<1 & BGM_F_padj<0.05 & BGM_F_logFC>1))


#     non DE in both                    
dim(subset(tot, BGM_M_padj>0.05 & BGM_F_padj>0.05))
#     DN in both 
dim(subset(tot, BGM_M_padj<0.05 & BGM_M_logFC<1 & BGM_F_padj<0.05 & BGM_F_logFC<1))
#     DN in M
dim(subset(tot, BGM_M_padj<0.05 & BGM_M_logFC<1 & BGM_F_padj>0.05))
#     DN in F
dim(subset(tot, BGM_M_padj>0.05 & BGM_F_padj<0.05 & BGM_F_logFC<1 ))
#     UP in M
dim(subset(tot, BGM_F_padj>0.05  & BGM_M_padj<0.05 & BGM_M_logFC>1))
#     UP in F
dim(subset(tot, BGM_M_padj>0.05  & BGM_F_padj<0.05 & BGM_F_logFC>1))
#     UP in both
dim(subset(tot, BGM_M_padj<0.05 & BGM_M_logFC>1 & BGM_F_padj<0.05 & BGM_F_logFC>1))
#     UP in M & DN in F
dim(subset(tot, BGM_M_padj<0.05 & BGM_M_logFC>1 & BGM_F_padj<0.05 & BGM_F_logFC<1))
#     DN in M & UP in F
dim(subset(tot, BGM_M_padj<0.05 & BGM_M_logFC<1 & BGM_F_padj<0.05 & BGM_F_logFC>1))





#     non_DE_in_both                    
NDB_NR <- nrow(subset(tot, BGM_M_padj>0.05 & BGM_F_padj>0.05 & BRO_padj>0.05))
NDB_UR <- nrow(subset(tot, BGM_M_padj>0.05 & BGM_F_padj>0.05 & BRO_padj<0.05 & BRO_logFC>1))
NDB_DR <- nrow(subset(tot, BGM_M_padj>0.05 & BGM_F_padj>0.05 & BRO_padj<0.05 & BRO_logFC<1))

#     DN_in_both 
DNB_NR <- nrow(subset(tot, BGM_M_padj<0.05 & BGM_M_logFC<1 & BGM_F_padj<0.05 & BGM_F_logFC<1 & BRO_padj>0.05))
DNB_UR <- nrow(subset(tot, BGM_M_padj<0.05 & BGM_M_logFC<1 & BGM_F_padj<0.05 & BGM_F_logFC<1 & BRO_padj<0.05 & BRO_logFC>1))
DNB_DR <- nrow(subset(tot, BGM_M_padj<0.05 & BGM_M_logFC<1 & BGM_F_padj<0.05 & BGM_F_logFC<1 & BRO_padj<0.05 & BRO_logFC<1))

#     DN_in_M
DNM_NR <- nrow(subset(tot, BGM_M_padj<0.05 & BGM_M_logFC<1 & BGM_F_padj>0.05 & BRO_padj>0.05))
DNM_UR <- nrow(subset(tot, BGM_M_padj<0.05 & BGM_M_logFC<1 & BGM_F_padj>0.05 & BRO_padj<0.05 & BRO_logFC>1))
DNM_DR <- nrow(subset(tot, BGM_M_padj<0.05 & BGM_M_logFC<1 & BGM_F_padj>0.05 & BRO_padj<0.05 & BRO_logFC<1))

#     DN_in_F
DNF_NR <- nrow(subset(tot, BGM_M_padj>0.05 & BGM_F_padj<0.05 & BGM_F_logFC<1 & BRO_padj>0.05))
DNF_UR <- nrow(subset(tot, BGM_M_padj>0.05 & BGM_F_padj<0.05 & BGM_F_logFC<1 & BRO_padj<0.05 & BRO_logFC>1))
DNF_DR <- nrow(subset(tot, BGM_M_padj>0.05 & BGM_F_padj<0.05 & BGM_F_logFC<1 & BRO_padj<0.05 & BRO_logFC<1))

#     UP_in_M
UPM_NR <- nrow(subset(tot, BGM_F_padj>0.05  & BGM_M_padj<0.05 & BGM_M_logFC>1 & BRO_padj>0.05))
UPM_UR <- nrow(subset(tot, BGM_F_padj>0.05  & BGM_M_padj<0.05 & BGM_M_logFC>1 & BRO_padj<0.05 & BRO_logFC>1))
UPM_DR <- nrow(subset(tot, BGM_F_padj>0.05  & BGM_M_padj<0.05 & BGM_M_logFC>1 & BRO_padj<0.05 & BRO_logFC<1))

#     UP_in_F
UPF_NR <- nrow(subset(tot, BGM_M_padj>0.05  & BGM_F_padj<0.05 & BGM_F_logFC>1 & BRO_padj>0.05))
UPF_UR <- nrow(subset(tot, BGM_M_padj>0.05  & BGM_F_padj<0.05 & BGM_F_logFC>1 & BRO_padj<0.05 & BRO_logFC>1))
UPF_DR <- nrow(subset(tot, BGM_M_padj>0.05  & BGM_F_padj<0.05 & BGM_F_logFC>1 & BRO_padj<0.05 & BRO_logFC<1))

#     UP_in_both
UPB_NR <- nrow(subset(tot, BGM_M_padj<0.05 & BGM_M_logFC>1 & BGM_F_padj<0.05 & BGM_F_logFC>1 & BRO_padj>0.05))
UPB_UR <- nrow(subset(tot, BGM_M_padj<0.05 & BGM_M_logFC>1 & BGM_F_padj<0.05 & BGM_F_logFC>1 & BRO_padj<0.05 & BRO_logFC>1))
UPB_DR <- nrow(subset(tot, BGM_M_padj<0.05 & BGM_M_logFC>1 & BGM_F_padj<0.05 & BGM_F_logFC>1 & BRO_padj<0.05 & BRO_logFC<1))

#     UP_in_M_&_DN_in_F
UMDF_NR <- nrow(subset(tot, BGM_M_padj<0.05 & BGM_M_logFC>1 & BGM_F_padj<0.05 & BGM_F_logFC<1 & BRO_padj>0.05))
UMDF_UR <- nrow(subset(tot, BGM_M_padj<0.05 & BGM_M_logFC>1 & BGM_F_padj<0.05 & BGM_F_logFC<1 & BRO_padj<0.05 & BRO_logFC>1))
UMDF_DR <- nrow(subset(tot, BGM_M_padj<0.05 & BGM_M_logFC>1 & BGM_F_padj<0.05 & BGM_F_logFC<1 & BRO_padj<0.05 & BRO_logFC<1))

#     DN_in_M_&_UP_in_F
DMUP_NR <- nrow(subset(tot, BGM_M_padj<0.05 & BGM_M_logFC<1 & BGM_F_padj<0.05 & BGM_F_logFC>1 & BRO_padj>0.05))
DMUP_UR <- nrow(subset(tot, BGM_M_padj<0.05 & BGM_M_logFC<1 & BGM_F_padj<0.05 & BGM_F_logFC>1 & BRO_padj<0.05 & BRO_logFC>1))
DMUP_DR <- nrow(subset(tot, BGM_M_padj<0.05 & BGM_M_logFC<1 & BGM_F_padj<0.05 & BGM_F_logFC>1 & BRO_padj<0.05 & BRO_logFC<1))






UPR_UPF <- UPF_UR + DMUP_UR
UPR_UPM <- UPM_UR + UMDF_UR
UPR_UPB <- UPB_UR
UPR_NOG <- NDB_UR + DNB_UR + DNM_UR + DNF_UR

NOR_UGF <- UPF_NR + UPF_DR + DMUP_NR + DMUP_DR
NOR_UGM <- UPM_NR + UPM_DR + UMDF_NR + UMDF_DR
NOR_UGB <- UPB_NR + UPB_DR
NOR_NOG <- NDB_NR + NDB_DR + DNB_NR + DNB_DR + DNM_NR + DNM_DR + DNF_NR + DNF_DR


links_BRO <- ""


links_BRO <-  data.frame(
    source=c("BGM_B_UP",         "BGM_B_UP",        "BGM_F_UP",        "BGM_F_UP",        "BGM_M_UP",        "BGM_M_UP",        "BGM_B_NO",        "BGM_B_NO"), 
    target=c("BRO_UP",           "BRO_NO",          "BRO_UP",          "BRO_NO",          "BRO_UP",          "BRO_NO",          "BRO_UP",          "BRO_NO"), 
  value=c(    UPR_UPB,            NOR_UGB,           UPR_UPF,           NOR_UGF,           UPR_UPM,           NOR_UGM,           UPR_NOG,           NOR_NOG)
)

# From these flows we need to create a node data frame: it lists every entities involved in the flow
nodes_BRO <- data.frame(
  name=c(as.character(links_BRO$source), 
         as.character(links_BRO$target)) %>% unique()
)

# With networkD3, connection must be provided using id, not using real name like in the links dataframe.. So we need to reformat it.
links_BRO$IDsource <- match(links_BRO$source, nodes_BRO$name)-1 
links_BRO$IDtarget <- match(links_BRO$target, nodes_BRO$name)-1

# Add a 'group' column to each connection:
links_BRO$group <- as.factor(c("a","b","a","b","a","b","a","b"))

# Add a 'group' column to each node. Here I decide to put all of them in the same group to make them grey
nodes_BRO$group <- as.factor(c("1"))

# Give a color for each group:
my_color <- 'd3.scaleOrdinal() .domain(["a", "b", "1"]) .range([ "#f55a4f", "lightgray", "white"])'

# composition of the genes upregulated in B. rossius gonads
BRO_SAN <-    sankeyNetwork(Links = links_BRO, Nodes = nodes_BRO,
                            Source = "IDsource", Target = "IDtarget",
                            Value = "value", NodeID = "name", 
                            sinksRight=TRUE, fontSize = 0, fontFamily = "Helvetica",
                            colourScale=my_color, LinkGroup="group", NodeGroup="group", nodePadding = 150, nodeWidth = 30,
                            iterations = 0)

##########################################################################################################################  final BAT easy sankey ############################

tot <- read.table(file = "new_DE_BGM_VS_BAT_single.lst", sep = " ", header= TRUE)


 na.omit(tot[3],tot[9])
       
#   non DE in both                                        DE in both                                            DE in M                                               DE in F                                               
dim(subset(tot, BGM_M_padj>0.05 & BGM_F_padj>0.05)) + dim(subset(tot, BGM_M_padj<0.05 & BGM_F_padj<0.05)) + dim(subset(tot, BGM_M_padj<0.05 & BGM_F_padj>0.05)) + dim(subset(tot, BGM_M_padj>0.05 & BGM_F_padj<0.05))
#     non DE in both                                        DN in both                                                                            DN in M                                                               DN in F                                                                UP in M                                                                UP in F                                                                UP in both
dim(subset(tot, BGM_M_padj>0.05 & BGM_F_padj>0.05)) + dim(subset(tot, BGM_M_padj<0.05 & BGM_M_logFC<1 & BGM_F_padj<0.05 & BGM_F_logFC<1)) + dim(subset(tot, BGM_M_padj<0.05 & BGM_M_logFC<1 & BGM_F_padj>0.05)) + dim(subset(tot, BGM_M_padj>0.05 & BGM_F_padj<0.05 & BGM_F_logFC<1 )) + dim(subset(tot, BGM_F_padj>0.05  & BGM_M_padj<0.05 & BGM_M_logFC>1)) + dim(subset(tot, BGM_M_padj>0.05  & BGM_F_padj<0.05 & BGM_F_logFC>1)) + dim(subset(tot, BGM_M_padj<0.05 & BGM_M_logFC>1 & BGM_F_padj<0.05 & BGM_F_logFC>1)) + dim(subset(tot, BGM_M_padj<0.05 & BGM_M_logFC>1 & BGM_F_padj<0.05 & BGM_F_logFC<1)) + dim(subset(tot, BGM_M_padj<0.05 & BGM_M_logFC<1 & BGM_F_padj<0.05 & BGM_F_logFC>1))





#     non DE in both                    
dim(subset(tot, BGM_M_padj>0.05 & BGM_F_padj>0.05))
#     DN in both 
dim(subset(tot, BGM_M_padj<0.05 & BGM_M_logFC<1 & BGM_F_padj<0.05 & BGM_F_logFC<1))
#     DN in M
dim(subset(tot, BGM_M_padj<0.05 & BGM_M_logFC<1 & BGM_F_padj>0.05))
#     DN in F
dim(subset(tot, BGM_M_padj>0.05 & BGM_F_padj<0.05 & BGM_F_logFC<1 ))
#     UP in M
dim(subset(tot, BGM_F_padj>0.05  & BGM_M_padj<0.05 & BGM_M_logFC>1))
#     UP in F
dim(subset(tot, BGM_M_padj>0.05  & BGM_F_padj<0.05 & BGM_F_logFC>1))
#     UP in both
dim(subset(tot, BGM_M_padj<0.05 & BGM_M_logFC>1 & BGM_F_padj<0.05 & BGM_F_logFC>1))
#     UP in M & DN in F
dim(subset(tot, BGM_M_padj<0.05 & BGM_M_logFC>1 & BGM_F_padj<0.05 & BGM_F_logFC<1))
#     DN in M & UP in F
dim(subset(tot, BGM_M_padj<0.05 & BGM_M_logFC<1 & BGM_F_padj<0.05 & BGM_F_logFC>1))





#     non_DE_in_both                    
NDB_NR <- nrow(subset(tot, BGM_M_padj>0.05 & BGM_F_padj>0.05 & BAT_padj>0.05))
NDB_UR <- nrow(subset(tot, BGM_M_padj>0.05 & BGM_F_padj>0.05 & BAT_padj<0.05 & BAT_logFC>1))
NDB_DR <- nrow(subset(tot, BGM_M_padj>0.05 & BGM_F_padj>0.05 & BAT_padj<0.05 & BAT_logFC<1))

#     DN_in_both 
DNB_NR <- nrow(subset(tot, BGM_M_padj<0.05 & BGM_M_logFC<1 & BGM_F_padj<0.05 & BGM_F_logFC<1 & BAT_padj>0.05))
DNB_UR <- nrow(subset(tot, BGM_M_padj<0.05 & BGM_M_logFC<1 & BGM_F_padj<0.05 & BGM_F_logFC<1 & BAT_padj<0.05 & BAT_logFC>1))
DNB_DR <- nrow(subset(tot, BGM_M_padj<0.05 & BGM_M_logFC<1 & BGM_F_padj<0.05 & BGM_F_logFC<1 & BAT_padj<0.05 & BAT_logFC<1))

#     DN_in_M
DNM_NR <- nrow(subset(tot, BGM_M_padj<0.05 & BGM_M_logFC<1 & BGM_F_padj>0.05 & BAT_padj>0.05))
DNM_UR <- nrow(subset(tot, BGM_M_padj<0.05 & BGM_M_logFC<1 & BGM_F_padj>0.05 & BAT_padj<0.05 & BAT_logFC>1))
DNM_DR <- nrow(subset(tot, BGM_M_padj<0.05 & BGM_M_logFC<1 & BGM_F_padj>0.05 & BAT_padj<0.05 & BAT_logFC<1))

#     DN_in_F
DNF_NR <- nrow(subset(tot, BGM_M_padj>0.05 & BGM_F_padj<0.05 & BGM_F_logFC<1 & BAT_padj>0.05))
DNF_UR <- nrow(subset(tot, BGM_M_padj>0.05 & BGM_F_padj<0.05 & BGM_F_logFC<1 & BAT_padj<0.05 & BAT_logFC>1))
DNF_DR <- nrow(subset(tot, BGM_M_padj>0.05 & BGM_F_padj<0.05 & BGM_F_logFC<1 & BAT_padj<0.05 & BAT_logFC<1))

#     UP_in_M
UPM_NR <- nrow(subset(tot, BGM_F_padj>0.05  & BGM_M_padj<0.05 & BGM_M_logFC>1 & BAT_padj>0.05))
UPM_UR <- nrow(subset(tot, BGM_F_padj>0.05  & BGM_M_padj<0.05 & BGM_M_logFC>1 & BAT_padj<0.05 & BAT_logFC>1))
UPM_DR <- nrow(subset(tot, BGM_F_padj>0.05  & BGM_M_padj<0.05 & BGM_M_logFC>1 & BAT_padj<0.05 & BAT_logFC<1))

#     UP_in_F
UPF_NR <- nrow(subset(tot, BGM_M_padj>0.05  & BGM_F_padj<0.05 & BGM_F_logFC>1 & BAT_padj>0.05))
UPF_UR <- nrow(subset(tot, BGM_M_padj>0.05  & BGM_F_padj<0.05 & BGM_F_logFC>1 & BAT_padj<0.05 & BAT_logFC>=1))
UPF_DR <- nrow(subset(tot, BGM_M_padj>0.05  & BGM_F_padj<0.05 & BGM_F_logFC>1 & BAT_padj<0.05 & BAT_logFC<1))

#     UP_in_both
UPB_NR <- nrow(subset(tot, BGM_M_padj<0.05 & BGM_M_logFC>1 & BGM_F_padj<0.05 & BGM_F_logFC>1 & BAT_padj>0.05))
UPB_UR <- nrow(subset(tot, BGM_M_padj<0.05 & BGM_M_logFC>1 & BGM_F_padj<0.05 & BGM_F_logFC>1 & BAT_padj<0.05 & BAT_logFC>1))
UPB_DR <- nrow(subset(tot, BGM_M_padj<0.05 & BGM_M_logFC>1 & BGM_F_padj<0.05 & BGM_F_logFC>1 & BAT_padj<0.05 & BAT_logFC<1))

#     UP_in_M_&_DN_in_F
UMDF_NR <- nrow(subset(tot, BGM_M_padj<0.05 & BGM_M_logFC>1 & BGM_F_padj<0.05 & BGM_F_logFC<1 & BAT_padj>0.05))
UMDF_UR <- nrow(subset(tot, BGM_M_padj<0.05 & BGM_M_logFC>1 & BGM_F_padj<0.05 & BGM_F_logFC<1 & BAT_padj<0.05 & BAT_logFC>1))
UMDF_DR <- nrow(subset(tot, BGM_M_padj<0.05 & BGM_M_logFC>1 & BGM_F_padj<0.05 & BGM_F_logFC<1 & BAT_padj<0.05 & BAT_logFC<1))

#     DN_in_M_&_UP_in_F
DMUP_NR <- nrow(subset(tot, BGM_M_padj<0.05 & BGM_M_logFC<1 & BGM_F_padj<0.05 & BGM_F_logFC>1 & BAT_padj>0.05))
DMUP_UR <- nrow(subset(tot, BGM_M_padj<0.05 & BGM_M_logFC<1 & BGM_F_padj<0.05 & BGM_F_logFC>1 & BAT_padj<0.05 & BAT_logFC>1))
DMUP_DR <- nrow(subset(tot, BGM_M_padj<0.05 & BGM_M_logFC<1 & BGM_F_padj<0.05 & BGM_F_logFC>1 & BAT_padj<0.05 & BAT_logFC<1))






UPR_UPF <- UPF_UR + DMUP_UR
UPR_UPM <- UPM_UR + UMDF_UR
UPR_UPB <- UPB_UR
UPR_NOG <- NDB_UR + DNB_UR + DNM_UR + DNF_UR

NOR_UGF <- UPF_NR + UPF_DR + DMUP_NR + DMUP_DR
NOR_UGM <- UPM_NR + UPM_DR + UMDF_NR + UMDF_DR
NOR_UGB <- UPB_NR + UPB_DR
NOR_NOG <- NDB_NR + NDB_DR + DNB_NR + DNB_DR + DNM_NR + DNM_DR + DNF_NR + DNF_DR


links_BAT <- ""


links_BAT <-  data.frame(
  source=c("BGM_B_UP",         "BGM_B_UP",        "BGM_F_UP",        "BGM_F_UP",        "BGM_M_UP",        "BGM_M_UP",        "BGM_B_NO",        "BGM_B_NO"), 
  target=c("BAT_UP",           "BAT_NO",          "BAT_UP",          "BAT_NO",          "BAT_UP",          "BAT_NO",          "BAT_UP",          "BAT_NO"), 
  value=c(    UPR_UPB,            NOR_UGB,           UPR_UPF,           NOR_UGF,           UPR_UPM,           NOR_UGM,           UPR_NOG,           NOR_NOG)
)

# From these flows we need to create a node data frame: it lists every entities involved in the flow
nodes_BAT <- data.frame(
  name=c(as.character(links_BAT$source), 
         as.character(links_BAT$target)) %>% unique()
)

# With networkD3, connection must be provided using id, not using real name like in the links dataframe.. So we need to reformat it.
links_BAT$IDsource <- match(links_BAT$source, nodes_BAT$name)-1 
links_BAT$IDtarget <- match(links_BAT$target, nodes_BAT$name)-1

# Add a 'group' column to each connection:
links_BAT$group <- as.factor(c("a","b","a","b","a","b","a","b"))

# Add a 'group' column to each node. Here I decide to put all of them in the same group to make them grey
nodes_BAT$group <- as.factor(c("1"))

# Give a color for each group:
my_color <- 'd3.scaleOrdinal() .domain(["a", "b", "1"]) .range([ "#ffc037", "lightgray", "white"])'

# composition of the genes upregulated in B. rossius gonads
BAT_SAN <-    sankeyNetwork(Links = links_BAT, Nodes = nodes_BAT,
                            Source = "IDsource", Target = "IDtarget",
                            Value = "value", NodeID = "name", 
                            sinksRight=TRUE, fontSize = 0, fontFamily = "Helvetica",
                            colourScale=my_color, LinkGroup="group", NodeGroup="group", nodePadding = 150, nodeWidth = 30,
                            iterations = 0)

##########################################################################################################################  final phylostratigraphy ############################


bat_totale_exp = read.csv("../input_files/BAT_RSEM.counts.matrix.BAT_GND_vs_BAT_LEG.DESeq2.DE_results", sep ="\t", header = TRUE)
bat_totale_gonad_specific <- subset(bat_totale_exp, padj<0.005 & log2FoldChange > 1)
write(bat_totale_exp[,1], file = "../input_files/bat_totale.txt")
write(bat_totale_gonad_specific[,1], file = "bat_gonad_specific.txt")

bro_totale_exp = read.csv("../input_files/BRO_RSEM.counts.matrix.BRO_GND_vs_BRO_LEG.DESeq2.DE_results", sep ="\t", header = TRUE)
bro_totale_gonad_specific <- subset(bro_totale_exp, padj<0.005 & log2FoldChange > 1)
write(bro_totale_exp[,1], file = "../input_files/bro_totale.txt")
write(bro_totale_gonad_specific[,1], file = "bro_gonad_specific.txt")

bgm_f_totale_exp = read.csv("../input_files/BGM_RSEM_f_only.counts.matrix.BGM_GND_vs_BGM_LEG.DESeq2.DE_results", sep ="\t", header = TRUE)
bgm_m_totale_exp = read.csv("../input_files/BGM_RSEM_m_only.counts.matrix.BGM_GND_vs_BGM_LEG.DESeq2.DE_results", sep ="\t", header = TRUE)

DE_bgm_B <- merge(bgm_f_totale_exp, bgm_m_totale_exp, by="transcript")
write(DE_bgm_B[,1], file = "../input_files/bgm_totale.txt")

bgm_B_totale_gonad_specic <- subset(DE_bgm_B, padj.x<0.005 & padj.y<0.005 & log2FoldChange.x > 1 & log2FoldChange.y > 1)
write(bgm_B_totale_gonad_specic[,1], file = "../input_files/bgm_B_gonad_specific.txt")

bgm_f_totale_gonad_specic <- subset(DE_bgm_B, padj.x<0.005 & log2FoldChange.x > 1 & padj.y>0.005 & log2FoldChange.y < 1)
write(bgm_f_totale_gonad_specic[,1], file = "../input_files/bgm_F_gonad_specific.txt")

bgm_m_totale_gonad_specic <- subset(DE_bgm_B, padj.y<0.005 & log2FoldChange.y > 1 & padj.x>0.005 & log2FoldChange.x < 1)
write(bgm_m_totale_gonad_specic[,1], file = "../input_files/bgm_M_gonad_specific.txt")


length(bgm_B_totale_gonad_specic[,1]) 
length(bgm_m_totale_gonad_specic[,1]) 
length(bgm_f_totale_gonad_specic[,1])
length(bat_totale_gonad_specific[,1])
length(bro_totale_gonad_specific[,1])


##########################################################################################################################  final phylostratigraphy ############################

# create a dataset
phylostrata = read.csv("../inputs/phylostrata.csv", sep =";", header = TRUE)
phylostrata <- as.data.table(phylostrata)
phylostrata_long <- melt(phylostrata, id.vars=c("type"))


# Grouped
age <- ggplot(phylostrata_long, aes(fill=type, y=variable, x=value)) + geom_bar(position="fill", stat="identity", width=0.9) + 
  theme(axis.text=element_text(size=10), axis.title=element_text(size=12,face="bold"), title=element_text(size=10,face="bold")) + 
  scale_fill_manual(values=c("#137177", "#137177", "#39A96B","#39A96B","#BFE1B0")) + theme(legend.position = "none") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_rect(colour = "black", size=1)) + 
  labs(title = "genes phylostratigraphy")+ theme_classic() +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
 
##########################################################################################################################  the new overlaps #######

tot <- read.table(file = "new_DE_BGM_VS_BAT_single.lst", sep = " ", header= TRUE)

#     non_DE_in_both                    
NDB_NR <- subset(tot, BGM_M_padj>0.05 & BGM_F_padj>0.05 & BAT_padj>0.05)
NDB_UR <- subset(tot, BGM_M_padj>0.05 & BGM_F_padj>0.05 & BAT_padj<0.05 & BAT_logFC>1)
NDB_DR <- subset(tot, BGM_M_padj>0.05 & BGM_F_padj>0.05 & BAT_padj<0.05 & BAT_logFC<1)
#     DN_in_both 
DNB_NR <- subset(tot, BGM_M_padj<0.05 & BGM_M_logFC<1 & BGM_F_padj<0.05 & BGM_F_logFC<1 & BAT_padj>0.05)
DNB_UR <- subset(tot, BGM_M_padj<0.05 & BGM_M_logFC<1 & BGM_F_padj<0.05 & BGM_F_logFC<1 & BAT_padj<0.05 & BAT_logFC>1)
DNB_DR <- subset(tot, BGM_M_padj<0.05 & BGM_M_logFC<1 & BGM_F_padj<0.05 & BGM_F_logFC<1 & BAT_padj<0.05 & BAT_logFC<1)
#     DN_in_M
DNM_NR <- subset(tot, BGM_M_padj<0.05 & BGM_M_logFC<1 & BGM_F_padj>0.05 & BAT_padj>0.05)
DNM_UR <- subset(tot, BGM_M_padj<0.05 & BGM_M_logFC<1 & BGM_F_padj>0.05 & BAT_padj<0.05 & BAT_logFC>1)
DNM_DR <- subset(tot, BGM_M_padj<0.05 & BGM_M_logFC<1 & BGM_F_padj>0.05 & BAT_padj<0.05 & BAT_logFC<1)
#     DN_in_F
DNF_NR <- subset(tot, BGM_M_padj>0.05 & BGM_F_padj<0.05 & BGM_F_logFC<1 & BAT_padj>0.05)
DNF_UR <- subset(tot, BGM_M_padj>0.05 & BGM_F_padj<0.05 & BGM_F_logFC<1 & BAT_padj<0.05 & BAT_logFC>1)
DNF_DR <- subset(tot, BGM_M_padj>0.05 & BGM_F_padj<0.05 & BGM_F_logFC<1 & BAT_padj<0.05 & BAT_logFC<1)
#     UP_in_M
UPM_NR <- subset(tot, BGM_F_padj>0.05  & BGM_M_padj<0.05 & BGM_M_logFC>1 & BAT_padj>0.05)
UPM_UR <- subset(tot, BGM_F_padj>0.05  & BGM_M_padj<0.05 & BGM_M_logFC>1 & BAT_padj<0.05 & BAT_logFC>1)
UPM_DR <- subset(tot, BGM_F_padj>0.05  & BGM_M_padj<0.05 & BGM_M_logFC>1 & BAT_padj<0.05 & BAT_logFC<1)
#     UP_in_F
UPF_NR <- subset(tot, BGM_M_padj>0.05  & BGM_F_padj<0.05 & BGM_F_logFC>1 & BAT_padj>0.05)
UPF_UR <- subset(tot, BGM_M_padj>0.05  & BGM_F_padj<0.05 & BGM_F_logFC>1 & BAT_padj<0.05 & BAT_logFC>=1)
UPF_DR <- subset(tot, BGM_M_padj>0.05  & BGM_F_padj<0.05 & BGM_F_logFC>1 & BAT_padj<0.05 & BAT_logFC<1)
#     UP_in_both
UPB_NR <- subset(tot, BGM_M_padj<0.05 & BGM_M_logFC>1 & BGM_F_padj<0.05 & BGM_F_logFC>1 & BAT_padj>0.05)
UPB_UR <- subset(tot, BGM_M_padj<0.05 & BGM_M_logFC>1 & BGM_F_padj<0.05 & BGM_F_logFC>1 & BAT_padj<0.05 & BAT_logFC>1)
UPB_DR <- subset(tot, BGM_M_padj<0.05 & BGM_M_logFC>1 & BGM_F_padj<0.05 & BGM_F_logFC>1 & BAT_padj<0.05 & BAT_logFC<1)
#     UP_in_M_&_DN_in_F
UMDF_NR <- subset(tot, BGM_M_padj<0.05 & BGM_M_logFC>1 & BGM_F_padj<0.05 & BGM_F_logFC<1 & BAT_padj>0.05)
UMDF_UR <- subset(tot, BGM_M_padj<0.05 & BGM_M_logFC>1 & BGM_F_padj<0.05 & BGM_F_logFC<1 & BAT_padj<0.05 & BAT_logFC>1)
UMDF_DR <- subset(tot, BGM_M_padj<0.05 & BGM_M_logFC>1 & BGM_F_padj<0.05 & BGM_F_logFC<1 & BAT_padj<0.05 & BAT_logFC<1)
#     DN_in_M_&_UP_in_F
DMUP_NR <- subset(tot, BGM_M_padj<0.05 & BGM_M_logFC<1 & BGM_F_padj<0.05 & BGM_F_logFC>1 & BAT_padj>0.05)
DMUP_UR <- subset(tot, BGM_M_padj<0.05 & BGM_M_logFC<1 & BGM_F_padj<0.05 & BGM_F_logFC>1 & BAT_padj<0.05 & BAT_logFC>1)
DMUP_DR <- subset(tot, BGM_M_padj<0.05 & BGM_M_logFC<1 & BGM_F_padj<0.05 & BGM_F_logFC>1 & BAT_padj<0.05 & BAT_logFC<1)

UPR_UPF_BAT <- c(UPF_UR$OG,DMUP_UR$OG)
UPR_UPM_BAT <- c(UPM_UR$OG,UMDF_UR$OG)
UPR_UPB_BAT <- UPB_UR$OG
UPR_NOG_BAT <- c(NDB_UR$OG,DNB_UR$OG,DNM_UR$OG,DNF_UR$OG)

NOR_UGF_BAT <- c(UPF_NR$OG,UPF_DR$OG,DMUP_NR$OG,DMUP_DR$OG)
NOR_UGM_BAT <- c(UPM_NR$OG,UPM_DR$OG,UMDF_NR$OG,UMDF_DR$OG)
NOR_UGB_BAT <- c(UPB_NR$OG,UPB_DR$OG)
NOR_NOG_BAT <- c(NDB_NR$OG,NDB_DR$OG,DNB_NR$OG,DNB_DR$OG,DNM_NR$OG,DNM_DR$OG,DNF_NR$OG,DNF_DR$OG)

tot <- read.table(file = "new_DE_BGM_VS_BRO_single.lst", sep = " ", header= TRUE)

#     non_DE_in_both                    
NDB_NR <- subset(tot, BGM_M_padj>0.05 & BGM_F_padj>0.05 & BRO_padj>0.05)
NDB_UR <- subset(tot, BGM_M_padj>0.05 & BGM_F_padj>0.05 & BRO_padj<0.05 & BRO_logFC>1)
NDB_DR <- subset(tot, BGM_M_padj>0.05 & BGM_F_padj>0.05 & BRO_padj<0.05 & BRO_logFC<1)
#     DN_in_both 
DNB_NR <- subset(tot, BGM_M_padj<0.05 & BGM_M_logFC<1 & BGM_F_padj<0.05 & BGM_F_logFC<1 & BRO_padj>0.05)
DNB_UR <- subset(tot, BGM_M_padj<0.05 & BGM_M_logFC<1 & BGM_F_padj<0.05 & BGM_F_logFC<1 & BRO_padj<0.05 & BRO_logFC>1)
DNB_DR <- subset(tot, BGM_M_padj<0.05 & BGM_M_logFC<1 & BGM_F_padj<0.05 & BGM_F_logFC<1 & BRO_padj<0.05 & BRO_logFC<1)
#     DN_in_M
DNM_NR <- subset(tot, BGM_M_padj<0.05 & BGM_M_logFC<1 & BGM_F_padj>0.05 & BRO_padj>0.05)
DNM_UR <- subset(tot, BGM_M_padj<0.05 & BGM_M_logFC<1 & BGM_F_padj>0.05 & BRO_padj<0.05 & BRO_logFC>1)
DNM_DR <- subset(tot, BGM_M_padj<0.05 & BGM_M_logFC<1 & BGM_F_padj>0.05 & BRO_padj<0.05 & BRO_logFC<1)
#     DN_in_F
DNF_NR <- subset(tot, BGM_M_padj>0.05 & BGM_F_padj<0.05 & BGM_F_logFC<1 & BRO_padj>0.05)
DNF_UR <- subset(tot, BGM_M_padj>0.05 & BGM_F_padj<0.05 & BGM_F_logFC<1 & BRO_padj<0.05 & BRO_logFC>1)
DNF_DR <- subset(tot, BGM_M_padj>0.05 & BGM_F_padj<0.05 & BGM_F_logFC<1 & BRO_padj<0.05 & BRO_logFC<1)
#     UP_in_M
UPM_NR <- subset(tot, BGM_F_padj>0.05  & BGM_M_padj<0.05 & BGM_M_logFC>1 & BRO_padj>0.05)
UPM_UR <- subset(tot, BGM_F_padj>0.05  & BGM_M_padj<0.05 & BGM_M_logFC>1 & BRO_padj<0.05 & BRO_logFC>1)
UPM_DR <- subset(tot, BGM_F_padj>0.05  & BGM_M_padj<0.05 & BGM_M_logFC>1 & BRO_padj<0.05 & BRO_logFC<1)
#     UP_in_F
UPF_NR <- subset(tot, BGM_M_padj>0.05  & BGM_F_padj<0.05 & BGM_F_logFC>1 & BRO_padj>0.05)
UPF_UR <- subset(tot, BGM_M_padj>0.05  & BGM_F_padj<0.05 & BGM_F_logFC>1 & BRO_padj<0.05 & BRO_logFC>=1)
UPF_DR <- subset(tot, BGM_M_padj>0.05  & BGM_F_padj<0.05 & BGM_F_logFC>1 & BRO_padj<0.05 & BRO_logFC<1)
#     UP_in_both
UPB_NR <- subset(tot, BGM_M_padj<0.05 & BGM_M_logFC>1 & BGM_F_padj<0.05 & BGM_F_logFC>1 & BRO_padj>0.05)
UPB_UR <- subset(tot, BGM_M_padj<0.05 & BGM_M_logFC>1 & BGM_F_padj<0.05 & BGM_F_logFC>1 & BRO_padj<0.05 & BRO_logFC>1)
UPB_DR <- subset(tot, BGM_M_padj<0.05 & BGM_M_logFC>1 & BGM_F_padj<0.05 & BGM_F_logFC>1 & BRO_padj<0.05 & BRO_logFC<1)
#     UP_in_M_&_DN_in_F
UMDF_NR <- subset(tot, BGM_M_padj<0.05 & BGM_M_logFC>1 & BGM_F_padj<0.05 & BGM_F_logFC<1 & BRO_padj>0.05)
UMDF_UR <- subset(tot, BGM_M_padj<0.05 & BGM_M_logFC>1 & BGM_F_padj<0.05 & BGM_F_logFC<1 & BRO_padj<0.05 & BRO_logFC>1)
UMDF_DR <- subset(tot, BGM_M_padj<0.05 & BGM_M_logFC>1 & BGM_F_padj<0.05 & BGM_F_logFC<1 & BRO_padj<0.05 & BRO_logFC<1)
#     DN_in_M_&_UP_in_F
DMUP_NR <- subset(tot, BGM_M_padj<0.05 & BGM_M_logFC<1 & BGM_F_padj<0.05 & BGM_F_logFC>1 & BRO_padj>0.05)
DMUP_UR <- subset(tot, BGM_M_padj<0.05 & BGM_M_logFC<1 & BGM_F_padj<0.05 & BGM_F_logFC>1 & BRO_padj<0.05 & BRO_logFC>1)
DMUP_DR <- subset(tot, BGM_M_padj<0.05 & BGM_M_logFC<1 & BGM_F_padj<0.05 & BGM_F_logFC>1 & BRO_padj<0.05 & BRO_logFC<1)

UPR_UPF_BRO <- c(UPF_UR$OG,DMUP_UR$OG)
UPR_UPM_BRO <- c(UPM_UR$OG,UMDF_UR$OG)
UPR_UPB_BRO <- UPB_UR$OG
UPR_NOG_BRO <- c(NDB_UR$OG,DNB_UR$OG,DNM_UR$OG,DNF_UR$OG)

NOR_UGF_BRO <- c(UPF_NR$OG,UPF_DR$OG,DMUP_NR$OG,DMUP_DR$OG)
NOR_UGM_BRO <- c(UPM_NR$OG,UPM_DR$OG,UMDF_NR$OG,UMDF_DR$OG)
NOR_UGB_BRO <- c(UPB_NR$OG,UPB_DR$OG)
NOR_NOG_BRO <- c(NDB_NR$OG,NDB_DR$OG,DNB_NR$OG,DNB_DR$OG,DNM_NR$OG,DNM_DR$OG,DNF_NR$OG,DNF_DR$OG)

###

myCol <- c("#f55a4f", "#ffc037")

venn.diagram(
  x = list(UPR_UPF_BRO, UPR_UPF_BAT),
  category.names = c("BRO" , "BAT "),
  filename = 'UPR_asex_UPF_sex.png',
  output=TRUE,
  imagetype="png",
  height = 2000, 
  width = 2000, 
  resolution = 300,
  compression = "lzw",
  # Circles
  lwd = 2,
  lty = 'blank',
  fill = myCol,
)
go.obj <- newGeneOverlap(UPR_UPF_BAT, 
                         UPR_UPF_BRO, 
                         15972)
go.obj <- testGeneOverlap(go.obj)
print(go.obj) 

venn.diagram(
  x = list(UPR_UPM_BRO, UPR_UPM_BAT),
  category.names = c("BRO" , "BAT "),
  filename = 'UPR_asex_UPM_sex.png',
  output=TRUE,
  imagetype="png" ,
  height = 2000, 
  width = 2000, 
  resolution = 300,
  compression = "lzw",
  # Circles
  lwd = 2,
  lty = 'blank',
  fill = myCol,
)
go.obj <- newGeneOverlap(UPR_UPM_BRO, 
                         UPR_UPM_BAT, 
                         15972)
go.obj <- testGeneOverlap(go.obj)
print(go.obj) 

venn.diagram(
  x = list(UPR_UPB_BRO, UPR_UPB_BAT),
  category.names = c("BRO" , "BAT "),
  filename = 'UPR_asex_UPB_sex.png',
  output=TRUE,
  imagetype="png" ,
  height = 2000, 
  width = 2000, 
  resolution = 300,
  compression = "lzw",
  # Circles
  lwd = 2,
  lty = 'blank',
  fill = myCol,
)
go.obj <- newGeneOverlap(UPR_UPB_BRO, 
                         UPR_UPB_BAT, 
                         15972)
go.obj <- testGeneOverlap(go.obj)
print(go.obj) 

venn.diagram(
  x = list(UPR_NOG_BRO, UPR_NOG_BAT),
  category.names = c("BRO" , "BAT "),
  filename = 'UPR_asex_NOR_sex.png',
  output=TRUE,
  imagetype="png" ,
  height = 2000, 
  width = 2000, 
  resolution = 300,
  compression = "lzw",
  # Circles
  lwd = 2,
  lty = 'blank',
  fill = myCol,
)
go.obj <- newGeneOverlap(UPR_NOG_BRO, 
                         UPR_NOG_BAT, 
                         15972)
go.obj <- testGeneOverlap(go.obj)
print(go.obj) 

#

venn.diagram(
  x = list(NOR_UGF_BRO, NOR_UGF_BAT),
  category.names = c("BRO" , "BAT "),
  filename = 'NOR_asex_UPF_sex.png',
  output=TRUE,
  imagetype="png" ,
  height = 2000, 
  width = 2000, 
  resolution = 300,
  compression = "lzw",
  # Circles
  lwd = 2,
  lty = 'blank',
  fill = myCol,
)
go.obj <- newGeneOverlap(NOR_UGF_BRO, 
                         NOR_UGF_BAT, 
                         15972)
go.obj <- testGeneOverlap(go.obj)
print(go.obj) 

venn.diagram(
  x = list(NOR_UGM_BRO, NOR_UGM_BAT),
  category.names = c("BRO" , "BAT "),
  filename = 'NOR_asex_UPM_sex.png',
  output=TRUE,
  imagetype="png" ,
  height = 2000, 
  width = 2000, 
  resolution = 300,
  compression = "lzw",
  # Circles
  lwd = 2,
  lty = 'blank',
  fill = myCol,
)
go.obj <- newGeneOverlap(NOR_UGM_BRO, 
                         NOR_UGM_BAT, 
                         15972)
go.obj <- testGeneOverlap(go.obj)
print(go.obj) 

venn.diagram(
  x = list(NOR_UGB_BRO, NOR_UGB_BAT),
  category.names = c("BRO" , "BAT "),
  filename = 'NOR_asex_UPB_sex.png',
  output=TRUE,
  imagetype="png" ,
  height = 2000, 
  width = 2000, 
  resolution = 300,
  compression = "lzw",
  # Circles
  lwd = 2,
  lty = 'blank',
  fill = myCol,
)
go.obj <- newGeneOverlap(NOR_UGB_BRO, 
                         NOR_UGB_BAT, 
                         15972)
go.obj <- testGeneOverlap(go.obj)
print(go.obj) 

venn.diagram(
  x = list(NOR_NOG_BRO, NOR_NOG_BAT),
  category.names = c("BRO" , "BAT "),
  filename = 'NOR_asex_NOR_sex.png',
  output=TRUE,
  imagetype="png" ,
  height = 2000, 
  width = 2000, 
  resolution = 300,
  compression = "lzw",
  # Circles
  lwd = 2,
  lty = 'blank',
  fill = myCol,
)
go.obj <- newGeneOverlap(NOR_NOG_BRO, 
                         NOR_NOG_BAT, 
                         15972)
go.obj <- testGeneOverlap(go.obj)
print(go.obj) 








##########################################################################################################################  the new gene list ######

UPR_UPF_BRO 
UPR_UPM_BRO 
UPR_UPB_BRO 
UPR_NOG_BRO 
NOR_UGF_BRO 
NOR_UGM_BRO 
NOR_UGB_BRO 
NOR_NOG_BRO 

UPR_UPF_BGM
UPR_UPM_BGM
UPR_UPB_BGM
UPR_NOG_BGM
NOR_UGF_BGM
NOR_UGM_BGM
NOR_UGB_BGM
NOR_NOG_BGM

write(intersect(UPR_NOG_BRO,UPR_NOG_BAT), file = "enrichment/upregulated_asex_&_not_BGM_B.txt")
write(intersect(UPR_UPM_BRO,UPR_UPM_BAT), file = "enrichment/upregulated_asex_&_upregulated_BGM_M.txt")
write(intersect(NOR_UGF_BRO,NOR_UGF_BAT), file = "enrichment/notregulated_asex_&_upregulated_BGM_F.txt")
write(intersect(NOR_UGB_BRO,NOR_UGB_BAT), file = "enrichment/notregulated_asex_&_upregulated_BGM_B.txt")

write(setdiff(as.vector(UPR_NOG_BRO),as.vector(intersect(UPR_NOG_BRO,UPR_NOG_BAT))), file = "enrichment/upregulated_BRO_&_not_BGM_B.txt")
write(setdiff(as.vector(UPR_NOG_BAT),as.vector(intersect(UPR_NOG_BRO,UPR_NOG_BAT))), file = "enrichment/upregulated_BAT_&_not_BGM_B.txt")

write(setdiff(as.vector(UPR_UPM_BRO),as.vector(intersect(UPR_UPM_BRO,UPR_UPM_BAT))), file = "enrichment/upregulated_BRO_&_upregulated_BGM_M.txt")
write(setdiff(as.vector(UPR_UPM_BAT),as.vector(intersect(UPR_UPM_BRO,UPR_UPM_BAT))), file = "enrichment/upregulated_BAT_&_upregulated_BGM_M.txt")

write(as.vector(setdiff(as.vector(NOR_UGF_BRO),as.vector(intersect(NOR_UGF_BRO,NOR_UGF_BAT)))), file = "enrichment/notregulated_BRO_&_upregulated_BGM_F.txt")
write(as.vector(setdiff(as.vector(NOR_UGF_BAT),as.vector(intersect(NOR_UGF_BRO,NOR_UGF_BAT)))), file = "enrichment/notregulated_BAT_&_upregulated_BGM_F.txt")

write(as.vector(setdiff(as.vector(NOR_UGB_BRO),as.vector(intersect(NOR_UGB_BRO,NOR_UGB_BAT)))), file = "enrichment/notregulated_BRO_&_upregulated_BGM_B.txt")
write(as.vector(setdiff(as.vector(NOR_UGB_BAT),as.vector(intersect(NOR_UGB_BRO,NOR_UGB_BAT)))), file = "enrichment/notregulated_BAT_&_upregulated_BGM_B.txt")

