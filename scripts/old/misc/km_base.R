source('load.R')
library(survival)
pdf(file='~/Desktop/KM.pdf',width=15,height=10)
par(mfrow=c(2,2), mar=c(4,4,4,4))

# IPI Clonality
ipi_pre <- ipi[ipi$type=='PRE',]
ipi_diverse <- ipi_pre$Clonality<0.1
ipi_s <- Surv(ipi_pre$os)

ipi_fit <- survfit(ipi_s~ipi_diverse)
plot(ipi_fit,
     lwd=2,
     col="blue2",
     lty=c(5,1),
     xlab='Months After Diagnosis',
     ylab='Percent Alive',
     xlim=c(0,16),
     main='anti-CTLA4'
     )
legend('topright',c('Clonal','Diverse'), col='blue2',lty=c(5,1),bty="n",lwd=2)
text(8,0.3,"p=0.013")
coxph(ipi_s~ipi_diverse)


#IPI EXP CL
ipi_post <- ipi[ipi$type=='POST3',]
ipi_exp_cl <- ipi_post[,'Number of Expanded Clones vs PRE']>100
ipi_post_s <- Surv(ipi_post$os.m)
ipi_post_fit <- survfit(ipi_post_s ~ ipi_exp_cl)

plot(ipi_post_fit,
     lwd=2,
     col="blue2",
     lty=c(5,1),
     xlab='Months After Diagnosis',
     ylab='Percent Alive',
     xlim=c(0,16),
     main='anti-CTLA4'
     )
legend('topright',c('Low','High'), col='blue2',lty=c(5,1),bty="n", lwd=2)
coxph(ipi_post_s~ipi_exp_cl)
text(9,0.5,"p=0.007")

#PD1 Clonality
pd1_pre <- pd1[pd1$type=='PRE',]
pd1_diverse <- pd1_pre$Clonality<0.1
pd1_s <- Surv(pd1_pre$os,!pd1_pre$alive)

pd1_fit <- survfit(pd1_s~pd1_diverse)
plot(pd1_fit,
     lwd=2,
     col='darkgreen',
     lty=c(5,1),
     xlab='Months After Diagnosis',
     ylab='Percent Alive',
     mark.time=T,
     main='anti-PD1',
     xlim=c(0,16)
     )
legend('topright',c('Clonal','Diverse'), col='darkgreen',lty=c(5,1),bty="n", lwd=2)
coxph(pd1_s~pd1_diverse)


#PD1 Exp Cl
pd1_post <- pd1[pd1$type=='POST3',]
pd1_exp_cl <- pd1_post[,'Number of Expanded Clones vs PRE']>100
pd1_post_s <- Surv(pd1_post$os, !pd1_post$alive)
pd1_post_fit <- survfit(pd1_post_s ~ pd1_exp_cl)

plot(pd1_post_fit,
     lwd=2,
     col='darkgreen',
     lty=c(5,1),
     xlab='Months After Diagnosis',
     ylab='Percent Alive',
     mark.time=T,
     xlim=c(0,16),
     main='anti-PD1'
     )
legend('topright',c('Low','High'), col='darkgreen',lty=c(5,1),bty="n", lwd=2)
coxph(pd1_post_s~pd1_exp_cl)


dev.off()

# Now using pre-treatment ALC

pdf(file='~/Desktop/KM_alc.pdf',width=7.5,height=10)
par(mfrow=c(2,1), mar=c(4,4,4,4))

ipi_pre <- ipi[ipi$type=='PRE',]
ipi_alc <- ipi_pre$ALC>1000
ipi_s <- Surv(ipi_pre$os)

ipi_alc_fit <- survfit(ipi_s~ipi_alc)
plot(ipi_alc_fit,
     lwd=2,
     col="blue2",
     lty=c(5,1),
     xlab='Months After Diagnosis',
     ylab='Percent Alive',
     xlim=c(0,16),
     main='anti-CTLA4'
     )
legend('topright',c('Low ALC','High ALC'), col='blue2',lty=c(5,1),bty="n",lwd=2)
coxph(ipi_s~ipi_alc)



pd1_pre <- pd1[pd1$type=='PRE',]
pd1_alc <- pd1_pre$ALC>1000
pd1_s <- Surv(pd1_pre$os,!pd1_pre$alive)

pd1_alc_fit <- survfit(pd1_s~pd1_alc)
plot(pd1_alc_fit,
     lwd=2,
     col='darkgreen',
     lty=c(5,1),
     xlab='Months After Diagnosis',
     ylab='Percent Alive',
     mark.time=T,
     main='anti-PD1',
     xlim=c(0,16)
     )
legend('topright',c('Low ALC','High ALC'), col='darkgreen',lty=c(5,1),bty="n", lwd=2)
coxph(pd1_s~pd1_alc)

dev.off()


# by arm

pdf(file='~/Desktop/KM_arm.pdf',width=7.5,height=10)
par(mfrow=c(2,1), mar=c(4,4,4,4))

ipi_pre <- ipi[ipi$type=='PRE',]
ipi_s <- Surv(ipi_pre$os)

ipi_arm_fit <- survfit(ipi_s~ipi_pre$arm)
plot(ipi_arm_fit,
     lwd=2,
     col="blue2",
     lty=c(5,1),
     xlab='Months After Diagnosis',
     ylab='Percent Alive',
     xlim=c(0,16),
     main='anti-CTLA4'
     )
legend('topright',c('anti-CTLA4','anti-CTLA4 + GVAX'), col='blue2',lty=c(5,1),bty="n",lwd=2)
coxph(ipi_s~ipi_alc)



pd1_pre <- pd1[pd1$type=='PRE',]
pd1_s <- Surv(pd1_pre$os,!pd1_pre$alive)

pd1_arm_fit <- survfit(pd1_s~pd1_pre$arm)
plot(pd1_arm_fit,
     lwd=2,
     col='darkgreen',
     lty=c(1,5),
     xlab='Months After Diagnosis',
     ylab='Percent Alive',
     mark.time=T,
     main='anti-PD1',
     xlim=c(0,16)
     )
legend('topright',c('GVAX + LM', 'anti-PD1 + GVAX + LM'), col='darkgreen',lty=c(5,1),bty="n", lwd=2)
coxph(pd1_s~pd1_alc)

dev.off()


# trying to make a model with both (just for ipi)
tmp <- ipi_pre[,c('patient', 'Clonality', 'os.m')]
tmp2 <- ipi_post[,c('patient','Number of Expanded Clones vs PRE')]
both <- merge(tmp, tmp2, by='patient', all=T)
names(both)[4] <- 'exp_cl'
both$both <- both$Clonality<0.1 & both$exp_cl>100
both$s <- Surv(both$os.m)
w <- which(is.na(both$both))
both$both[w] <- 'Unknown'
both_fit <- survfit(both$s~both$both)

#plot 

plot(both_fit, 
     lwd=2,
     col="blue2",
     lty=c(5,1,3),
     xlab='Months After Diagnosis',
     ylab='Percent Alive',
     xlim=c(0,16),
     main='anti-CTLA4'
) 
legend('topright',c('Both','Not Both', 'Unknown'), col='blue2',lty=c(1,5,3),bty="n",lwd=2)
coxph(both$s~both$both)

