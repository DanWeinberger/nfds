library(readxl)
cog_annotate <- read_excel("./Data/NIHMS74007-supplement-Supplementary_Dataset_1.xlsx")
#cog_annotate2 <- readRDS('./Data/cog_annotations.rds')

d1 <- read.csv('./Data/data_MassUS.csv')

d1.cogs <- d1[, grep('CLS', names(d1))]

e_l_year <- aggregate(d1.cogs, by=list('year'=d1$Year), FUN=mean)

cog_labels <- names(e_l_year)[-1]

e_l_2001 <- as.vector(t(e_l_year[1,-1]))
e_l_2007 <- as.vector(t(e_l_year[3,-1]))
e_l_ratio <- e_l_2007/e_l_2001
log_e_l_ratio <- log(e_l_ratio)

log_e_l_ratio <- cbind.data.frame(cog_labels,e_l_2001,e_l_2007,e_l_ratio, log_e_l_ratio)
log_e_l_ratio$cog_labels <- as.character(log_e_l_ratio$cog_labels)

#Seems to be a mismatch in labels between COG labels
log_e_l_ratio <- merge(log_e_l_ratio, cog_annotate, all=T, by.x='cog_labels', by.y='COG')

#We want locis that are ~1 (stable pre- and post-); ~0 on log scale
hist(log_e_l_ratio)


mod1 <- 


plot(e_l_2001,e_l_2007)
abline(a=0, b=1)


