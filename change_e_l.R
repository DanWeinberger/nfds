library(readxl)
cog_annotate <- read_excel("./Data/NIHMS74007-supplement-Supplementary_Dataset_1.xlsx")


d1 <- read.csv('./Data/data_MassUS.csv')

d1.cogs <- d1[, grep('CLS', names(d1))]

e_l_year <- aggregate(d1.cogs, by=list('year'=d1$Year), FUN=mean)

e_l_2001 <- as.vector(t(e_l_year[1,-1]))
e_l_2007 <- as.vector(t(e_l_year[3,-1]))

plot(e_l_2001,e_l_2007)
abline(a=0, b=1)
