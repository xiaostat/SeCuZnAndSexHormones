rm(list = ls())
library(foreign)
library(data.table)
library(fBasics)
library(nortest)
library(Hmisc)
library(questionr)
library(car)

setwd("/Users/xiao/Zero/RongLiu/TraceMetalSexSteroidHormone/data")
#setwd("E:/RongLiu/TraceMetalSexSteroidHormone")
cx.h <-  list.files(pattern = "*H.XPT")
cx.i <-  list.files(pattern = "*I.XPT")
num.h <- length(cx.h)
num.i <- length(cx.i)
index.h <- list()
index.i <- list()
#a <- list.files()
#file.remove(a[grep(".csv",a)])
for (i in 1:num.h) {
  index.h[[i]] <- data.frame(read.xport(cx.h[i]))
  #write.csv(index[[i]],paste(cx[i], ".csv", sep = ""))
}
for (j in 1:num.i) {
  index.i[[j]] <- data.frame(read.xport(cx.i[j]))
  #write.csv(index[[i]],paste(cx[i], ".csv", sep = ""))
}

slctvar.h <- Reduce(function(x, y) 
  merge(x, y, by=intersect(names(x),names(y)), all=TRUE), index.h, accumulate =F)
slctvar.i <- Reduce(function(x, y) 
  merge(x, y, by=intersect(names(x),names(y)), all=TRUE), index.i, accumulate =F)
slctvar <- merge(slctvar.h,slctvar.i,
                 by = intersect(names(slctvar.h),names(slctvar.i)),all = T)

#write.csv(slctvar,"tracemetal.csv")

vrb.anlys <- c('SEQN','RIDAGEYR','RIAGENDR','RIDRETH3','RIDEXMON','INDFMPIR',#'DMDEDUC2','SDDSRVYR','WTINT2YR','WTMEC2YR',#DEMO
               'LBDCOTLC','LBXCOT',#'LBXHCT','LBDHCTLC',,#COT
               'BMXBMI',#BMX
               'WTSA2YR','LBXSCU','LBXSSE','LBXSZN',#CUSEZN
               'PHDSESN',#FASTQX
               'LBXTST','LBXEST','LBXSHBG','LBDSHGLC','LBDTSTLC','LBDESTLC'#TST
               )
SE.all <- slctvar
SE.all <- SE.all[,colnames(SE.all) %in% vrb.anlys]
dim(SE.all)#26502    20
SE.all <- SE.all[which(SE.all$RIDAGEYR >= 6 & SE.all$RIDAGEYR <=19),]
dim(SE.all)#5451   20
# write.csv(SE,"tracemetal.csv")

#======= Remove NAs of primary variables in samples
if (nrow(SE.all[!complete.cases(SE.all),]) > 0) {
  SE <- SE.all[complete.cases(SE.all),]#na.omit(SE.all)
}
dim(SE)#1097   20
#======= Remove NAs of primary variables in samples

#======= Variables control
SE$chdado <- ifelse(SE$RIDAGEYR >= 6 & SE$RIDAGEYR <= 11,"A",'B')
SE$chdado <- factor(SE$chdado,labels = c('Children','Adolescent'))
SE$RIAGENDR<-factor(SE$RIAGENDR,labels=c('Male','Female'))
SE$MEC4YR <- 1/2 * SE$WTSA2YR
SE$group <- ifelse(SE$RIAGENDR == 'Male' & SE$chdado == 'Children',1,
                   ifelse(SE$RIAGENDR == 'Male' & SE$chdado == 'Adolescent',2,
                          ifelse(SE$RIAGENDR == 'Female' & SE$chdado == 'Children',3,
                                 ifelse(SE$RIAGENDR == 'Female' & SE$chdado == 'Adolescent',4,NA))))
SE$group <- factor(SE$group,labels = c('MaleChild','MaleAdole','FemaleChild','FemaleAdole'))

library(dplyr)
SE.mcma <- SE %>% mutate(group = as.character(group)) %>% 
  filter(group == 'MaleChild' | group == 'MaleAdole')
SE.mcfc <- SE %>% mutate(group = as.character(group)) %>% 
  filter(group == 'MaleChild' | group == 'FemaleChild')
SE.mcfa <- SE %>% mutate(group = as.character(group)) %>% 
  filter(group == 'MaleChild' | group == 'FemaleAdole')
SE.mafc <- SE %>% mutate(group = as.character(group)) %>% 
  filter(group == 'MaleAdole' | group == 'FemaleChild')
SE.mafa <- SE %>% mutate(group = as.character(group)) %>% 
  filter(group == 'MaleAdole' | group == 'FemaleAdole')
SE.fcfa <- SE %>% mutate(group = as.character(group)) %>% 
  filter(group == 'FemaleChild' | group == 'FemaleAdole')

SE$RIDRETH3 <- ifelse(SE$RIDRETH3 == 1,'C',
                      ifelse(SE$RIDRETH3 == 3,'A',
                             ifelse(SE$RIDRETH3 == 4,'B','D')))
SE$RIDRETH3 <- factor(SE$RIDRETH3,labels = c('Non-Hispanic White',
                                             'Non-Hispanic Black',
                                             'Mexican American',
                                             'Other Race'))

SE$pir <- ifelse(SE$INDFMPIR <= 1,'A','B')
SE$pir <- factor(SE$pir,labels = c('pir<=1','pir>1'))#ordered = T

SE$LBDCOTLC <- ifelse(SE$LBDCOTLC == 0,'B','A')
SE$LBDCOTLC <- factor(SE$LBDCOTLC,labels = c('Below','AtAbove'))#,ordered = T

SE$bmi <- ifelse(SE$BMXBMI < 25,'A',
                 ifelse(SE$BMXBMI <30,'B','C'))
SE$bmi <- factor(SE$bmi,labels = c('NormalWeight','Overweight','Obese'))#,ordered = T

SE$RIDEXMON <- factor(SE$RIDEXMON,labels = c('Nov2Apr','May2Oct'))

SE$PHDSESN <- factor(SE$PHDSESN,labels = c('Morning','Afternoon','Evening'))#,ordered = T

length(which(SE$LBXSCU <= 2.5/sqrt(2)))
length(which(SE$LBXSSE <= 2.9/sqrt(2)))
length(which(SE$LBXSZN <= 4.5/sqrt(2)))
length(which(SE$LBXTST <= 0.75/sqrt(2)))/length(SE$LBXTST)
length(which(SE$LBXEST <= 2.994/sqrt(2)))/length(SE$LBXEST)
length(which(SE$LBXSHBG <= 0.8/sqrt(2)))

length(which(SE[SE$group == 'MaleChild','LBXEST'] <= 2.994/sqrt(2)))/length(SE[SE$group == 'MaleChild','LBXEST'])
length(which(SE[SE$group == 'MaleAdole','LBXEST'] <= 2.994/sqrt(2)))/length(SE[SE$group == 'MaleChild','LBXEST'])
length(which(SE[SE$group == 'FemaleChild','LBXEST'] <= 2.994/sqrt(2)))/length(SE[SE$group == 'MaleChild','LBXEST'])
length(which(SE[SE$group == 'FemaleAdole','LBXEST'] <= 2.994/sqrt(2)))/length(SE[SE$group == 'MaleChild','LBXEST'])

#======= Variables control

#======= Pre-analysis
tab.vb <- c('group','RIDRETH3','pir','LBDCOTLC','bmi','RIDEXMON','PHDSESN')
a <- lapply(SE[,tab.vb],table)
lapply(a,addmargins)
rate.tab <- function(x){round(prop.table(x)*100,2)}
lapply(a,rate.tab)
b <- lapply(SE[,tab.vb],FUN = questionr::wtd.table, weights = SE$MEC4YR)
lapply(b,addmargins)
lapply(b,rate.tab)

lapply(SE[,tab.vb],FUN = function(vb,group){table(vb,group)},group = SE$group)
lapply(SE[,tab.vb],FUN = function(vb,group){rate.tab(table(vb,group))},group = SE$group)
lapply(SE[,tab.vb],FUN = function(vb,group){chisq.test(table(vb,group))},group = SE$group)

lapply(SE[,tab.vb],
       FUN = function(vb,group,weights){questionr::wtd.table(vb,group,weights)},
       group = SE$group,weights=SE$MEC4YR)
lapply(SE[,tab.vb],
       FUN = function(vb,group,weights){rate.tab(questionr::wtd.table(vb,group,weights))},
       group = SE$group,
       weights = SE$MEC4YR)
lapply(SE[,tab.vb],
       FUN = function(vb,group,weights){chisq.test(questionr::wtd.table(vb,group,weights))},
       group = SE$group,
       weights=SE$MEC4YR)

a <- chisq.test(questionr::wtd.table(SE.mcma[,tab.vb[2]],SE.mcma$group,weights=SE.mcma$MEC4YR))
a$p.value * 6
b <- chisq.test(questionr::wtd.table(SE.mcfc[,tab.vb[2]],SE.mcfc$group,weights=SE.mcfc$MEC4YR))
b$p.value * 6
c <- chisq.test(questionr::wtd.table(SE.mcfa[,tab.vb[2]],SE.mcfa$group,weights=SE.mcfa$MEC4YR))
c$p.value * 6
d <- chisq.test(questionr::wtd.table(SE.mafc[,tab.vb[2]],SE.mafc$group,weights=SE.mafc$MEC4YR))
d$p.value * 6
e <- chisq.test(questionr::wtd.table(SE.mafa[,tab.vb[2]],SE.mafa$group,weights=SE.mafa$MEC4YR))
e$p.value * 6
f <- chisq.test(questionr::wtd.table(SE.fcfa[,tab.vb[2]],SE.fcfa$group,weights=SE.fcfa$MEC4YR))
f$p.value * 6


shapiro.test(rep(SE.MaCh$RIDAGEYR,times=SE.MaCh$MEC4YR))
nortest::lillie.test(rep(SE.MaCh$RIDAGEYR,times=SE.MaCh$MEC4YR))
bartlett.test(rep(SE$RIDAGEYR,times=SE$MEC4YR) ~ rep(SE$group,times = SE$MEC4YR))#p-value < 2.2e-16
leveneTest(rep(SE$RIDAGEYR,times=SE$MEC4YR) ~ rep(SE$group,times = SE$MEC4YR))
bartlett.test(SE$RIDAGEYR ~ SE$group)#p-value < 2.2e-16
kruskal.test(SE$RIDAGEYR ~ SE$group)#p-value < 2.2e-16
kruskal.test(rep(SE$RIDAGEYR,times=SE$MEC4YR) ~ rep(SE$group,times = SE$MEC4YR))#p-value < 2.2e-16
anova(lm(SE$RIDAGEYR ~ SE$group,weights = SE$MEC4YR))#p-value < 2.2e-16
aov.year <- aov(SE$RIDAGEYR ~ SE$group,weights = SE$MEC4YR)#p-value < 2.2e-16
summary(aov.year)
TukeyHSD(aov.year)# Post Hoc.
#Welch's anova
oneway.test(rep(SE$RIDAGEYR,times=SE$MEC4YR) ~ rep(SE$group,times = SE$MEC4YR),var.equal = F)


vb <- c('RIDAGEYR','LBXSSE','LBXSCU','LBXSZN','LBXTST','LBXEST','LBXSHBG')
apply(SE[,vb],2,FUN = lillie.test)
skew <- data.frame(apply(SE[,vb],2,FUN = skewness))
colnames(skew) <- 'Total'
cat('Skewness of un-weighted continuous vb\n')
print(skew)
#kurtosis
null <- data.frame()
for (i in 1:length(vb)) {
  a <- t(data.frame(tapply(SE[,vb[i]],SE$group,FUN = skewness)))
  rownames(a) <- vb[i]
  null <- rbind(null,a)
}
#Skewness of un-weighted continuous vb
cbind(skew,null)
gc()
wtd.vb <- apply(SE[,vb],2,FUN = rep,times=SE$MEC4YR)# weight each vb
wtd.group <- rep(SE$group,times=SE$MEC4YR)
apply(wtd.vb,2,FUN = lillie.test)#weighted Kolmogorov-Smirnov normality test


wtd.skew <- data.frame(apply(wtd.vb,2,FUN = skewness))
colnames(wtd.skew) <- 'Total'
wtd.skew
null <- data.frame()
for (i in 1:length(vb)) {
  a <- t(data.frame(tapply(wtd.vb[,vb[i]],wtd.group,FUN = skewness)))
  rownames(a) <- vb[i]
  null <- rbind(null,a)
}
#Skewness of weighted continuous vb
cbind(wtd.skew,null)


IQR <- function(x){
  quantile(x,probs=0.75) - quantile(x,probs=0.25)
}

vb.median <- c('LBXTST','LBXEST','LBXSHBG')
md.iqr <- data.frame(apply(SE[,vb.median],2,FUN=median),apply(SE[,vb.median],2,FUN=IQR))
colnames(md.iqr) <- c('median','IQR')
print(md.iqr)
#          median     IQR
# LBXTST   19.00 112.950
# LBXEST   14.80  32.583
# LBXSHBG  55.05  58.770

null.md <- data.frame()
null.iqr <- data.frame()
for (i in 1:length(vb.median)) {
  a <- t(data.frame(tapply(SE[,vb.median[i]],SE$group,FUN = median)))
  b <- t(data.frame(tapply(SE[,vb.median[i]],SE$group,FUN = IQR)))
  rownames(a) <- vb.median[i]
  rownames(b) <- vb.median[i]
  null.md <- rbind(null.md,a)
  null.iqr <- rbind(null.iqr,b)
}
cbind(md.iqr[,'median'],null.md)
cbind(md.iqr[,'IQR'],null.iqr)

# weighted median and IQR
wtd.IQR <- function(x,weights){
  wtd.quantile(x,weights,probs=0.75) - wtd.quantile(x,weights,probs=0.25)
} 
wtd.md.iqr <- data.frame(apply(SE[,vb.median],2,FUN = wtd.quantile, weights = SE$MEC4YR, probs = 0.5),
                     apply(SE[,vb.median],2,FUN = wtd.IQR, weights = SE$MEC4YR))
colnames(wtd.md.iqr) <- c('wtd.median','wtd.iqr')
print(wtd.md.iqr)


null.wtd.md <- data.frame()
null.wtd.iqr <- data.frame()
for (i in 1:length(vb.median)) {
  a <- t(data.frame(tapply(rep(SE[,vb.median[i]],times=SE$MEC4YR),
                           rep(SE$group,times=SE$MEC4YR),FUN = median)))
  b <- t(data.frame(tapply(rep(SE[,vb.median[i]],times=SE$MEC4YR),
                           rep(SE$group,times=SE$MEC4YR),FUN = IQR)))
  rownames(a) <- vb.median[i]
  rownames(b) <- vb.median[i]
  null.wtd.md <- rbind(null.wtd.md,a)
  null.wtd.iqr <- rbind(null.wtd.iqr,b)
}
cbind(wtd.md.iqr[,'wtd.median'],null.wtd.md)
cbind(wtd.md.iqr[,'wtd.iqr'],null.wtd.iqr)


x <- rep(SE$LBXSSE,times=SE$MEC4YR)
y <- rep(SE$group,times=SE$MEC4YR)
leveneTest(x ~ y)
leveneTest(SE$LBXSSE ~ SE$group)#0.01465
bartlett.test(SE$LBXSSE ~ SE$group)#0.03454
kruskal.test(x ~ y)
kruskal.test(SE$LBXSSE~SE$group)
oneway.test(SE$LBXSSE ~ SE$group,var.equal = F) # Welch' ANOVA
oneway.test(rep(SE$LBXSSE,times=SE$MEC4YR) ~ rep(SE$group,times = SE$MEC4YR),var.equal = F)

vb.mean <- c('RIDAGEYR','LBXSSE','LBXSCU','LBXSZN')
mean.sd <- data.frame(apply(SE[,vb.mean],2,FUN=mean),apply(SE[,vb.mean],2,FUN=sd))
colnames(mean.sd) <- c('mean','sd')
print(mean.sd)
#              mean        sd
# RIDAGEYR  12.38013  3.878926
# LBXSSE   121.59325 14.290774
# LBXSCU   114.14148 26.133967
# LBXSZN    81.75205 14.424303

null.mean <- data.frame()
null.sd <- data.frame()
for (i in 1:length(vb.mean)) {
  a <- t(data.frame(tapply(SE[,vb.mean[i]],SE$group,FUN = mean)))
  b <- t(data.frame(tapply(SE[,vb.mean[i]],SE$group,FUN = sd)))
  rownames(a) <- vb.mean[i]
  rownames(b) <- vb.mean[i]
  null.mean <- rbind(null.mean,a)
  null.sd <- rbind(null.sd,b)
}
cbind(mean.sd[,'mean'],null.mean)
cbind(mean.sd[,'sd'],null.sd)


wtd.mean.sd <- data.frame(apply(SE[,vb.mean],2,FUN = wtd.mean, weights = SE$MEC4YR),
                          sqrt(apply(SE[,vb.mean],2,FUN = Hmisc:::wtd.var, weights = SE$MEC4YR)))
colnames(wtd.mean.sd) <- c('wtd.mean','wtd.sd')
print(wtd.mean.sd)


null.wtd.mean <- data.frame()
null.wtd.sd <- data.frame()
for (i in 1:length(vb.mean)) {
  a <- t(data.frame(tapply(rep(SE[,vb.mean[i]],times=SE$MEC4YR),
                           rep(SE$group,times=SE$MEC4YR),FUN = mean)))
  b <- t(data.frame(tapply(rep(SE[,vb.mean[i]],times=SE$MEC4YR),
                           rep(SE$group,times=SE$MEC4YR),FUN = sd)))
  rownames(a) <- vb.mean[i]
  rownames(b) <- vb.mean[i]
  null.wtd.mean <- rbind(null.wtd.mean,a)
  null.wtd.sd <- rbind(null.wtd.sd,b)
}
cbind(wtd.mean.sd[,'wtd.mean'],null.wtd.mean)
cbind(wtd.mean.sd[,'wtd.sd'],null.wtd.sd)


x <- rep(SE$LBXSCU,times=SE$MEC4YR)
y <- rep(SE$group,times=SE$MEC4YR)
leveneTest(x ~ y)
leveneTest(SE$LBXSCU ~ SE$group)#0.005278 **
bartlett.test(x ~ y)
kruskal.test(x ~ y)
anova(lm(rep(SE$LBXSCU,times=SE$MEC4YR) ~ rep(SE$group,times=SE$MEC4YR)))
leveneTest(SE$LBXSCU~SE$group)
bartlett.test(SE$LBXSCU~SE$group)
kruskal.test(SE$LBXSCU~SE$group)
anova(lm(SE$LBXSCU ~ SE$group))
oneway.test(SE$LBXSCU ~ SE$group,var.equal = F)
oneway.test(rep(SE$LBXSCU,times=SE$MEC4YR) ~ rep(SE$group,times = SE$MEC4YR),var.equal = F)


u <- rep(SE$LBXSZN,times=SE$MEC4YR)
v <- rep(SE$group,times=SE$MEC4YR)
leveneTest(u ~ v)
kruskal.test(u ~ v)
anova(lm(SE$LBXSZN ~ SE$group,weights = SE$MEC4YR))
leveneTest(SE$LBXSZN~SE$group)#0.1007
bartlett.test(SE$LBXSZN~SE$group)#0.07185
kruskal.test(SE$LBXSZN~SE$group)
anova(lm(SE$LBXSZN ~ SE$group,weights = SE$MEC4YR))
anova(lm(SE$LBXSZN ~ SE$group))
oneway.test(SE$LBXSZN ~ SE$group,var.equal = F)
oneway.test(rep(SE$LBXSZN,times=SE$MEC4YR) ~ rep(SE$group,times = SE$MEC4YR),var.equal = F)


anova(lm(SE$LBXTST ~ SE$group,weights = SE$MEC4YR))
kruskal.test(rep(SE$LBXTST,times=SE$MEC4YR) ~ rep(SE$group,times=SE$MEC4YR))
leveneTest(SE$LBXTST~SE$group)
bartlett.test(SE$LBXTST~SE$group)
kruskal.test(SE$LBXTST~SE$group)
oneway.test(SE$LBXTST ~ SE$group,var.equal = F)
oneway.test(rep(SE$LBXTST,times=SE$MEC4YR) ~ rep(SE$group,times = SE$MEC4YR),var.equal = F)


anova(lm(SE$LBXEST ~ SE$group,weights = SE$MEC4YR))
kruskal.test(rep(SE$LBXEST,times=SE$MEC4YR) ~ rep(SE$group,times=SE$MEC4YR))
leveneTest(SE$LBXEST~SE$group)
bartlett.test(SE$LBXEST~SE$group)
kruskal.test(SE$LBXEST~SE$group)
oneway.test(SE$LBXEST ~ SE$group,var.equal = F)
oneway.test(rep(SE$LBXEST,times=SE$MEC4YR) ~ rep(SE$group,times = SE$MEC4YR),var.equal = F)

anova(lm(SE$LBXSHBG ~ SE$group,weights = SE$MEC4YR))
kruskal.test(rep(SE$LBXSHBG,times=SE$MEC4YR) ~ rep(SE$group,times=SE$MEC4YR))
leveneTest(SE$LBXSHBG~SE$group)
bartlett.test(SE$LBXSHBG~SE$group)
kruskal.test(SE$LBXSHBG~SE$group)
oneway.test(SE$LBXSHBG ~ SE$group,var.equal = F)
oneway.test(rep(SE$LBXSHBG,times=SE$MEC4YR) ~ rep(SE$group,times = SE$MEC4YR),var.equal = F)


#=== post hoc
#vb.median <- c('LBXTST','LBXEST','LBXSHBG')
#vb.mean <- c('RIDAGEYR','LBXSSE','LBXSCU','LBXSZN')
library(multcomp)
fit <- aov(RIDAGEYR ~ group, data = SE)#weights=MEC4YR
summary(fit)
library(agricolae)
fit_snk <- SNK.test(fit,'group',group=FALSE)
fit_snk$comparison

fit <- aov(LBXSCU ~ group, data = SE)#weights=MEC4YR
summary(fit)
fit_snk <- SNK.test(fit,'group',group=FALSE)
fit_snk$comparison
# difference pvalue signif.  
# FemaleAdole - FemaleChild  -4.025169 0.0536  
# FemaleAdole - MaleAdole    16.916075 0.0000     *** 
# FemaleAdole - MaleChild    -9.890127 0.0000     *** 
# FemaleChild - MaleAdole    20.941244 0.0000     *** 
# FemaleChild - MaleChild    -5.864958 0.0050      ** 
# MaleAdole - MaleChild     -26.806202 0.0000     *** 

fit <- aov(LBXSSE ~ group, data = SE)#weights=MEC4YR
summary(fit)
fit_snk <- SNK.test(fit,'group',group=FALSE)
fit_snk$comparison
# difference pvalue signif.        LCL         UCL
# FemaleAdole - FemaleChild   7.727989 0.0000     *** 
# FemaleAdole - MaleAdole    -2.182112 0.0660       .  
# FemaleAdole - MaleChild     5.419393 0.0000     ***  
# FemaleChild - MaleAdole    -9.910101 0.0000     *** 
# FemaleChild - MaleChild    -2.308596 0.0518       .  
# MaleAdole - MaleChild       7.601505 0.0000     ***  
fit <- aov(LBXSZN ~ group, data = SE)#weights=MEC4YR
summary(fit)
fit_snk <- SNK.test(fit,'group',group=FALSE)
fit_snk$comparison
# difference pvalue signif.       LCL       UCL
# FemaleAdole - FemaleChild -1.01043095 0.6915         -3.906797  1.885935
# FemaleAdole - MaleAdole   -4.46920959 0.0017      ** -7.644591 -1.293828
# FemaleAdole - MaleChild   -0.04918067 0.9682         -2.470698  2.372336
# FemaleChild - MaleAdole   -3.45877864 0.0052      ** -5.880296 -1.037262
# FemaleChild - MaleChild    0.96125029 0.4362         -1.460267  3.382767
# MaleAdole - MaleChild      4.42002892 0.0010     ***  1.523663  7.316395

fit <- aov(LBXTST ~ group, data = SE)#weights=MEC4YR
summary(fit)
fit_snk <- SNK.test(fit,'group',group=FALSE)
fit_snk$comparison
# difference pvalue signif.         LCL        UCL
# FemaleAdole - FemaleChild   17.465131 0.1492           -4.504452   39.43471
# FemaleAdole - MaleAdole   -367.795090 0.0000     *** -386.162834 -349.42735
# FemaleAdole - MaleChild      9.317308 0.3198           -9.050436   27.68505
# FemaleChild - MaleAdole   -385.260221 0.0000     *** -409.346199 -361.17424
# FemaleChild - MaleChild     -8.147823 0.3843          -26.515567   10.21992
# MaleAdole - MaleChild      377.112398 0.0000     ***  355.142816  399.08198

fit <- aov(LBXEST ~ group, data = SE)#weights=MEC4YR
summary(fit)
fit_snk <- SNK.test(fit,'group',group=FALSE)
fit_snk$comparison
# difference pvalue signif.        LCL       UCL
# FemaleAdole - FemaleChild  70.103983 0.0000     ***  60.683066 79.524901
# FemaleAdole - MaleAdole    65.124016 0.0000     ***  57.247626 73.000406
# FemaleAdole - MaleChild    82.105083 0.0000     ***  71.776621 92.433545
# FemaleChild - MaleAdole    -4.979967 0.2150         -12.856357  2.896423
# FemaleChild - MaleChild    12.001099 0.0029      **   4.124710 19.877489
# MaleAdole - MaleChild      16.981067 0.0001     ***   7.560149 26.401984

fit <- aov(LBXSHBG ~ group, data = SE)#weights=MEC4YR
summary(fit)
fit_snk <- SNK.test(fit,'group',group=FALSE)
fit_snk$comparison
# difference pvalue signif.       LCL        UCL
# FemaleAdole - FemaleChild  -23.89110      0     *** -31.08018 -16.702030
# FemaleAdole - MaleAdole     18.10324      0     ***  10.91416  25.292313
# FemaleAdole - MaleChild    -40.17388      0     *** -48.77270 -31.575059
# FemaleChild - MaleAdole     41.99434      0     ***  33.39552  50.593166
# FemaleChild - MaleChild    -16.28278      0     *** -23.47185  -9.093703
# MaleAdole - MaleChild      -58.27712      0     *** -67.70429 -48.849948


#========= Weighted linear regression 
#reg.vb <- 'LBXTST'#'LBXEST','LBXSHBG'
Data   <- SE[SE$group=='FemaleAdole',]#MaleChild  MaleAdole FemaleChild FemaleAdole
reg.formula <- as.formula(LBXSHBG ~ RIDRETH3 + INDFMPIR + LBXCOT + BMXBMI 
                          + RIDEXMON + PHDSESN + LBXSSE + LBXSCU + LBXSZN)
reg.form.sqrt <- update(reg.formula,sqrt(.) ~ .)
reg.form.log  <- update(reg.formula,log(.) ~ .)


reg.origin <- lm(reg.formula, weights = MEC4YR, data = Data)
# Remove outliers
outlierResult <- car::outlierTest(reg.origin)# 离群点, 针对Y
while(outlierResult$bonf.p[1] < 0.05){
  exclusionRows <- names(outlierResult[[1]])
  inclusionRows <- !(rownames(Data) %in% exclusionRows)
  Data <- Data[inclusionRows,]
  reg.last <- lm(formula = reg.formula, weights = MEC4YR, data = Data)
  outlierResult <- car::outlierTest(reg.last)
}
# Remove outliers
# Remove High Leverage Points
leverage.con <- names(which(hatvalues(reg.last) > 
                              3 * length(coefficients(reg.last))/length(fitted(reg.last))))
while(length(leverage.con) > 0){
  exclusionRows <- leverage.con
  inclusionRows <- !(rownames(Data) %in% exclusionRows)
  Data <- Data[inclusionRows,]
  reg.last <- lm(reg.formula, weights = MEC4YR, data = Data)
  leverage.con <- names(which(hatvalues(reg.last) > 
                                3 * length(coefficients(reg.last))/length(fitted(reg.last))))
}
  
influential <- which(cooks.distance(reg.last) 
                    > 4/(nrow(Data)-length(reg.last$coefficients)-2))#4*mean(reg.last,na.rm=T)
while(length(influential) > 0){#
  exclusionRows <- rownames(Data)[influential]
  inclusionRows <- !(rownames(Data) %in% exclusionRows)
  Data <- Data[inclusionRows,]
  reg.last <- lm(reg.formula, weights = MEC4YR, data = Data)
  influential <- which(cooks.distance(reg.last) 
                        > 4/(nrow(Data)-length(reg.last$coefficients)-2))#4*mean(reg.last,na.rm=T)
    
}
summary(power <- car::powerTransform(reg.last))#变量变换使得残差正态
if (power$roundlam < 0.25) {
  reg.tsf <- lm(formula = reg.form.log, weights = MEC4YR,data=Data)
  cat('Box-Cox transformation is: log.')
  reg.final.result <- reg.tsf
} else if (power$roundlam < 0.75) {
  reg.tsf <- lm(formula = reg.form.sqrt, weights = MEC4YR,data=Data)
  cat('Box-Cox transformation is: square root.')
  reg.final.result <- reg.tsf
} else if (power$roundlam < 1.5) {
  reg.tsf <- reg.origin
  cat('Box-Cox transformation is: origin.')
  reg.final.result <- reg.origin
}

summary(reg.origin)
summary(reg.last)
summary(reg.tsf)
summary(reg.final.result)
car::vif(reg.last)
car::crPlots(reg.last)# 线性诊断
car::ncvTest(reg.last)# 方差齐性检验
car::ncvTest(reg.tsf)# 方差齐性检验
car::spreadLevelPlot(reg.last)# 推荐变换使得方差齐同
summary(unitest <- gvlma::gvlma(reg.last,alphalevel = 0.005))# 综合验证
summary(unitest <- gvlma::gvlma(reg.tsf,alphalevel = 0.005))# 综合验证
durbinWatsonTest(reg.last)#检验残差独立性
confint(reg.tsf)
boxTidwell(LBXTST ~ I(INDFMPIR+1) + LBXCOT + BMXBMI + LBXSSE + LBXSCU + LBXSZN, 
           ~ RIDRETH3 + RIDEXMON + PHDSESN, data=Data)

result <- as.data.frame(summary(reg.last)$coefficients[11:13,c('Estimate','Pr(>|t|)')])
result <- cbind(result,as.data.frame(confint(reg.last)[11:13,]))
colnames(result) <- c('Estimate', 'P value','Lower', 'Upper')
result
dim(Data)
sum(Data$MEC4YR)
