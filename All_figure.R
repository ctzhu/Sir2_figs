
library(survival)
setwd('/Users/user/Work/Sir2_figs')

pdf("All_figures.pdf", width=10,height=7.5)

df <- read.csv('Sir2_dup_del_stv.csv')
S1 <- survfit(Surv(time)~Food+Geno, conf.type="none", 
              data=df[(df$Food=='V')&((df$Geno=='22')|(df$Geno=='OO')),])
S2 <- survfit(Surv(time)~Food+Geno, conf.type="none",
              data=df[(df$Food=='IV')&((df$Geno=='22')|(df$Geno=='OO')),])
par(mfrow=c(2,1))
par(mar=c(2,2,2,2))
CLS <- c(1,2,3,4)
plot(S1, main="Starvation Survivorship after 1Y3S (high sugar) diet feeding", lty=1, lwd=1, col=CLS, xlim=c(0,135))
plot(S2, main="Starvation Survivorship after 3Y1S (high yeast) diet feeding", lty=1, lwd=1, col=CLS, xlim=c(0,135))
text(125, 0.05, 'Hours')
legend(70, 1., c('Duplication', 'Deletion'), 
       bty='n', lty=1, lwd=1, col=CLS)
mdl <- coxph(Surv(time)~Food*Geno, data=droplevels(df[((df$Geno=='22')|(df$Geno=='OO')),]))
S <- summary(mdl)
C <- printCoefmat(S$coefficient, digits=3)


library(survival)
setwd('/Users/user/Work/Sir2_figs')
df <- read.csv('Sir2_LS_Rep1.csv')
S1 <- survfit(Surv(time)~Food+Genotype, conf.type="none", data=df[df$Sex=='F',])
S2 <- survfit(Surv(time)~Food+Genotype, conf.type="none", data=df[df$Sex=='M',])
par(mfrow=c(2,1))
par(mar=c(2,2,2,2))
CLS <- c(1,2,3,4)
LTY <- rep(c(1,2), each=4)
plot(S1, main="Female", lty=LTY, lwd=1, col=CLS, xlim=c(0,135))
plot(S2, main="Male", lty=LTY, lwd=1, col=CLS, xlim=c(0,135))
text(125, 0.05, 'Days')
legend(90, 1.05, names(S1$strata), 
       bty='n', 
       lty=LTY, lwd=1, col=CLS)
mdl <- coxph(Surv(time)~Food*Genotype*Sex, data=df)
S <- summary(mdl)
C <- printCoefmat(S$coefficient, digits=3)

library(survival)
setwd('/Users/user/Work/Sir2_figs')
df <- read.csv('Sir2_LS_Rep2.csv')
S1 <- survfit(Surv(time)~Food+Genotype, conf.type="none", data=df[df$Sex=='F',])
S2 <- survfit(Surv(time)~Food+Genotype, conf.type="none", data=df[df$Sex=='M',])
par(mfrow=c(2,1))
par(mar=c(2,2,2,2))
CLS <- c(3,2,4,1)
LTY <- rep(c(1,2), each=4)
plot(S1, main="Female", lty=LTY, lwd=1, col=CLS, xlim=c(0,135))
plot(S2, main="Male", lty=LTY, lwd=1, col=CLS, xlim=c(0,135))
text(125, 0.05, 'Hours')
legend(90, 1.05, names(S1$strata), 
       bty='n', lty=LTY, lwd=1, col=CLS)
mdl <- coxph(Surv(time)~Food*Genotype*Sex, data=df)
S <- summary(mdl)
C <- printCoefmat(S$coefficient, digits=3)

library(survival)
setwd('/Users/user/Work/Sir2_figs')
df <- read.csv('Sir2_LS_Rep3_deldup.csv')
S1 <- survfit(Surv(time)~Food+Genotype, conf.type="none", data=df[df$Sex=='F',])
S2 <- survfit(Surv(time)~Food+Genotype, conf.type="none", data=df[df$Sex=='M',])
par(mfrow=c(2,1))
par(mar=c(2,2,2,2))
CLS <- c(1,2,3,4)
LTY <- rep(c(1,2), each=4)
plot(S1, main="Female", lty=LTY, lwd=1, col=CLS, xlim=c(0,135))
plot(S2, main="Male", lty=LTY, lwd=1, col=CLS, xlim=c(0,135))
text(125, 0.05, 'Days')
legend(90, 1.05, names(S1$strata), 
       bty='n', 
       lty=LTY, lwd=1, col=CLS)
mdl <- coxph(Surv(time)~Food*Genotype*Sex, data=df)
S <- summary(mdl)
C <- printCoefmat(S$coefficient, digits=3)

library(survival)
setwd('/Users/user/Work/Sir2_figs')
df <- read.csv('Sir2_LS_Rep3_overw1118.csv')
S1 <- survfit(Surv(time)~Food+Genotype, conf.type="none", data=df[df$Sex=='F',])
S2 <- survfit(Surv(time)~Food+Genotype, conf.type="none", data=df[df$Sex=='M',])
par(mfrow=c(2,1))
par(mar=c(2,2,2,2))
CLS <- c(3,2,4,1)
LTY <- rep(c(1,2), each=4)
plot(S1, main="Female", lty=LTY, lwd=1, col=CLS, xlim=c(0,135))
plot(S2, main="Male", lty=LTY, lwd=1, col=CLS, xlim=c(0,135))
text(125, 0.05, 'Hours')
legend(90, 1.05, names(S1$strata), 
       bty='n', lty=LTY, lwd=1, col=CLS)
mdl <- coxph(Surv(time)~Food*Genotype*Sex, data=df)
S <- summary(mdl)
C <- printCoefmat(S$coefficient, digits=3)

library(survival)
setwd('/Users/user/Work/Sir2_figs')
df <- read.csv('Sir2_dup_del_stv_mito.csv')
S1 <- survfit(Surv(Time1)~Sir2+mito, conf.type="none", data=df[(df$Food=='Y')&(df$Mother=='w1118'),])
S2 <- survfit(Surv(Time1)~Sir2+mito, conf.type="none", data=df[(df$Food=='S')&(df$Mother=='w1118'),])
S3 <- survfit(Surv(Time1)~Sir2+mito, conf.type="none", data=df[(df$Food=='Y')&(df$Mother=='mito'),])
S4 <- survfit(Surv(Time1)~Sir2+mito, conf.type="none", data=df[(df$Food=='S')&(df$Mother=='mito'),])
par(mfrow=c(2,2))
par(mar=c(2,2,2,2))
CLS <- c(1,2,3,4)
LTY <- rep(c(2,1), each=4)
plot(S1, main="w1118 mother, 3Y1S (high yeast) diet feeding", 
     lty=LTY, lwd=1, col=CLS, xlim=c(0,135))
plot(S2, main="w1118 mother, 1Y3S (high sugar) diet feeding",
     lty=LTY, lwd=1, col=CLS, xlim=c(0,135))
plot(S3, main="mito mother, 3Y1S (high yeast) diet feeding", 
     lty=LTY, lwd=1, col=CLS, xlim=c(0,135))
plot(S4, main="mito mother, 1Y3S (high sugar) diet feeding",
     lty=LTY, lwd=1, col=CLS, xlim=c(0,135))
text(125, 0.05, 'Hours')
legend(0, 0.55, names(S1$strata), 
       bty='n', 
       lty=LTY, lwd=1, col=CLS)


setwd('/Users/user/Work/Sir2_figs')
df   <- read.csv('Sir2_dup_del_eclosion_mito.csv')
lm1  <- lm(Counts~(Mito*Geno*Sex) %in% Block + Block, data=df)
aov1 <- summary(aov(lm1))
boxplot(Counts~Mito*Geno, data=df, ylab='No.offspring (5 mating pairs)', xlab='Mitotype-Genotype', 
        cex.lab=1.5, cex.axis=1)

setwd('/Users/user/Work/Sir2_figs')
df   <- read.csv('Sir2_dup_del_Biochem_VIV.CSV')
par(mfrow=c(3,3))
umug=expression(paste(mu, 'M/', mu, 'g'))
boxplot(fly_per~food*genotype, data=df, main='Fresh Weight per fly', ylab='mg/fly', cex.lab=1.5)
boxplot(BCA~food*genotype, data=df, main='Soluble Protein(BCA assay)', ylab=expression(paste(mu, 'g/L')), cex.lab=1.5)
boxplot(NADH_BCA~food*genotype, data=df, main='NADH/BCA', ylab=umug, cex.lab=1.5)
boxplot(NAD_BCA~food*genotype, data=df, main='NAD+/BCA', ylab=umug, cex.lab=1.5)
boxplot(NADPH_BCA~food*genotype, data=df, main='NADPH/BCA', ylab=umug, cex.lab=1.5)
boxplot(NADP_BCA~food*genotype, data=df, main='NADP+/BCA', ylab=umug, cex.lab=1.5)
boxplot(TAG_BCA~food*genotype, data=df, main='TAG/BCA', ylab=umug, cex.lab=1.5)
boxplot(Glc_BCA~food*genotype, data=df, main='Glucose/BCA', ylab=umug, cex.lab=1.5)
boxplot(Gly_BCA~food*genotype, data=df, main='Glycogen/BCA', ylab=umug, cex.lab=1.5)


dev.off()
