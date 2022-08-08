require(tidyverse)
require(GeneOverlap)

df <- read.table("../intermediate_files/phylostrata.csv",sep=";",header=T)

df <- df %>%bind_rows(summarise_all(., ~if(is.numeric(.)) sum(.) else "Total"))

for(var in 2:9) {
  df[7,var] <- (df[5,as.numeric(var)]/df[6,as.numeric(var)]*100)
}

df[7,1] <- "proportion"

bat <- prop.test(x = c(26, 527), n = c(1169, 13703), alternative = "l")
bro <- prop.test(x = c(54, 773), n = c(1511, 14162), alternative = "l")
bgmf <- prop.test(x = c(4, 651), n = c(262, 1366), alternative = "l")
bgmm <- prop.test(x = c(60, 651), n = c(1615,1366), alternative = "l")
bgmb <- prop.test(x = c(3, 651), n = c(487, 1366), alternative = "l")
