
##------------------------------------------------------------------------##
## This script reproduces the analyses done to Performance variables      ##
## using Full-subset analysis for each coral group on Central reefs       ##
## Model selection, variable importance scores and variance contribution  ##
## for the selected predictors in the Best Models                         ##
##------------------------------------------------------------------------##


library(dplyr)
library(RCurl)
require(car)
library(corrplot)
require(doBy)
require(gplots)
require(RColorBrewer)
require(gridExtra)
library(ggplot2)
library(dplyr)
library(gridExtra)
library(grid)
library(mgcv)
library(ggpubr)
library(FSSgam)
library(tidyr)
library(magrittr)

##----------------------------------------------
#Load data-set for analyses of Performance variables
##----------------------------------------------
load("CI_performance.Rdata")

##----------------------------------------------
#Central Performance variables
##----------------------------------------------

##Transform response for modeling with  beta-distribution
## transform fraction with 0 and 1 to <1 and >0
transform01 <- function(x) {
  (x * (length(x) - 1) + 0.5) / (length(x))
}


Centdat_scores_ndist%<>% mutate(
    ACBX_score=transform01(ACBX_score),
    ACTO_score=transform01(ACTO_score),
    CBRN_score=transform01(CBRN_score),
    MSE_score=transform01(MSE_score),
    Total_score=transform01(Total_score))%>%
    drop_na()


##-----------------------------------------------------------------
##Full Subset analyses run on HPC separately for each coral group-
##-----------------------------------------------------------------


Cent.preds=c("lagDHW","lagCOTS","lagCycl","TCI5_COTS","TCI5_cycl","TCI5","log.sal","sqrt.MUD","sqrt.DINe","Chla","mPAR","log.DINriv")
lnvars=c("TCI5_dhw","TCI3_dhw","lagCIacute")
factor.vars= c("zone","lagCOTSp")


ModelHC=gam(Total_score~s(lagDHW,bs='cr',k=3)+ s(time,bs='cr',k=5)+s(LATITUDE, k=5)+ 
              s(FULLREEF_ID,bs='re')+s(REEF_SITE_NO,bs='re'),family=betar(link="logit"),data=Centdat_scores_ndist)


CTotal_model.set=generate.model.set(use.dat=Centdat_scores_ndist,
                                   max.predictors=5, 
                                   test.fit=ModelHC, k=3, 
                                   pred.vars.cont=Cent.preds,
                                   linear.vars=lnvars,
                                   pred.vars.fact = factor.vars,
                                   factor.smooth.interactions=F,
                                   cyclic.vars=NA,
                                   smooth.smooth.interactions=F,
                                   cov.cutoff=0.28,                
                                   null.terms="s(time,bs='cr',k=5)+s(LATITUDE,k=5)+s(FULLREEF_ID,bs='re')+s(REEF_SITE_NO,bs='re')")

CTotal.list=fit.model.set(CTotal_model.set)## FITS mod
##extract Variable Importance Scores
TOTAL.var.imp=CTotal.list$variable.importance$aic$variable.weights.raw
# Model selection table and best models
mod.table.CTotal=CTotal.list$mod.data.out
mod.table.CTotal=mod.table.CTotal[order(mod.table.CTotal$AICc),]
CTotal.less.2AICc=mod.table.CTotal[which(mod.table.CTotal$delta.AICc<2),]


##-------------------
## Acropora tabular
##-------------------
Model.ACTO=gam(ACTO_score~s(lagDHW,bs='cr',k=3)+s(time, bs = "cr", k = 5) +
                 s(LATITUDE, k = 5)+ s(FULLREEF_ID,bs = "re")+ 
                 s(REEF_SITE_NO, bs = "re"),family=betar(link="logit"),data=Centdat_scores_ndist)


CActo_model.set=generate.model.set(use.dat=Centdat_scores_ndist,
                                         max.predictors=5, 
                                         test.fit=Model.ACTO, k=3, 
                                         linear.vars=lnvars,
                                         pred.vars.cont=Cent.preds,
                                         pred.vars.fact =factor.vars,
                                         factor.smooth.interactions= F,
                                         cyclic.vars=NA,
                                         smooth.smooth.interactions=F,
                                         cov.cutoff=0.28,              
                                         null.terms="s(time,bs='cr',k=5)+s(LATITUDE,k=5)+s(FULLREEF_ID,bs ='re')+s(REEF_SITE_NO,bs='re')")

CActo.list=fit.model.set(CActo_model.set)## FITS mod
##extract Variable Importance Scores
ACTO.var.imp=CActo.list$variable.importance$aic$variable.weights.raw
# Model selection table and best models
mod.table.CActo=CActo.list$mod.data.out
mod.table.CActo=mod.table.CActo[order(mod.table.CActo$AICc),]
CActo.less.2AICc=mod.table.CActo[which(mod.table.CActo$delta.AICc<2),]

##-------------------
## Acropora branching
##-------------------

Cent.preds=c("lagDHW","lagCOTS","lagCycl","lagCIacute","TCI5_COTS","TCI5_cycl","TCI3_dhw",
             "TCI5_dhw","TCI5","log.sal","sqrt.MUD","sqrt.DINe","Chla","mPAR","log.DINriv")


Model.ACBX=gam(ACBX_score~s(lagDHW,bs='cr',k=3)+s(time, bs = "cr", k = 5) +
                 s(LATITUDE, k = 5)+
                 s(FULLREEF_ID,bs = "re")+ 
                 s(REEF_SITE_NO, bs = "re"),family=betar(link="logit"),data=Centdat_scores_ndist)

CAcbr_model.set=generate.model.set(use.dat=Centdat_scores_ndist,
                                         max.predictors=5, 
                                         test.fit=Model.ACBX, k=3, 
                                         pred.vars.cont=Cent.preds,
                                         pred.vars.fact =factor.vars,
                                         factor.smooth.interactions= F,
                                         cyclic.vars=NA,
                                         smooth.smooth.interactions=F,
                                         cov.cutoff=0.28,              
                                         null.terms="s(time,bs='cr',k=5)+s(LATITUDE,k=5)+s(FULLREEF_ID,bs ='re')+s(REEF_SITE_NO,bs='re')")


CAcbr.list=fit.model.set(CAcbr_model.set)## FITS mod
##extract Variable Importance Scores
ACBR.var.imp=CAcbr.list$variable.importance$aic$variable.weights.raw
# Model selection table and best models
mod.table.CAcbr=CAcbr.list$mod.data.out
mod.table.CAcbr=mod.table.CAcbr[order(mod.table.CAcbr$AICc),]
CAcbr.less.2AICc=mod.table.CAcbr[which(mod.table.CAcbr$delta.AICc<2),]

##-------------------
## Other branching corals
##-------------------

Model.CBRN=gam(CBRN_score~s(lagDHW,bs='cr',k=3)+s(time, bs = "cr", k = 5) +
                 s(LATITUDE, k = 5)+
                 s(FULLREEF_ID,bs = "re")+ 
                 s(REEF_SITE_NO, bs = "re"),family=betar(link="logit"),data=Centdat_scores_ndist)


Central_modelCBRN.set=generate.model.set(use.dat=Centdat_scores_ndist,
                                         max.predictors=5, 
                                         test.fit=Model.CBRN, k=3, 
                                         pred.vars.cont=Cent.preds,
                                         pred.vars.fact =factor.vars,
                                         factor.smooth.interactions= F,
                                         cyclic.vars=NA,
                                         smooth.smooth.interactions=F,
                                         cov.cutoff=0.28,              
                                         null.terms="s(time,bs='cr',k=5)+s(LATITUDE,k=5)+s(FULLREEF_ID,bs ='re')+s(REEF_SITE_NO,bs='re')")

Ccbrn.list=fit.model.set(Ccbrn_model.set)## FITS mod
##extract Variable Importance Scores
CBRN.var.imp=Ccbrn.list$variable.importance$aic$variable.weights.raw
# Model selection table and best models
mod.table.Ccbrn=Ccbrn.list$mod.data.out
mod.table.Ccbrn=mod.table.Ccbrn[order(mod.table.Ccbrn$AICc),]
Ccbrn.less.2AICc=mod.table.Ccbrn[which(mod.table.Ccbrn$delta.AICc<2),]

##----------------------------
## Massive-sub-massive-encrusting
##----------------------------

Model.MSE=gam(MSE_score~s(lagDHW,bs='cr',k=3)+s(time, bs = "cr", k = 5) +
                s(LATITUDE, k = 5)+
                s(FULLREEF_ID,bs = "re")+ 
                s(REEF_SITE_NO, bs = "re"),family=betar(link="logit"),data=Centdat_scores_ndist)


CMse_model.set=generate.model.set(use.dat=Centdat_scores_ndist,
                                        max.predictors=5, 
                                        test.fit=Model.MSE, k=3, 
                                        pred.vars.cont=Cent.preds,
                                        pred.vars.fact =factor.vars,
                                        factor.smooth.interactions= F,
                                        cyclic.vars=NA,
                                        smooth.smooth.interactions=F,
                                        cov.cutoff=0.28,              
                                        null.terms="s(time,bs='cr',k=5)+s(LATITUDE,k=5)+s(FULLREEF_ID,bs ='re')+s(REEF_SITE_NO,bs='re')")

CMse.list=fit.model.set(CMse_model.set)## FITS mod
##extract Variable Importance Scores
MSE.var.imp=CMse.list$variable.importance$aic$variable.weights.raw
# Model selection table and best models
mod.table.CMse=CMse.list$mod.data.out
mod.table.CMse=mod.table.CMse[order(mod.table.CMse$AICc),]
CMse.less.2AICc=mod.table.CMse[which(mod.table.CMse$delta.AICc<2),]


#-----------------------------------
##Gather data from Models for Variable Importance Scores
#-----------------------------------

ALLcorals.vis<- do.call("cbind",mget(ls(pattern = "var.imp*")))
CPF_vis=as.data.frame(ALLcorals.vis)
names(CPF_vis)<- substr(colnames(CPF_vis), 1, 4)
CPF_vis$predictor <- rownames(CPF_vis)

CPF_vis$predictor=dplyr::recode(CPF_vis$predictor,
                                      TCI5_COTS="5yr-outbreaks",
                                      TCI5_cycl="5yr-storms",
                                      TCI3_dhw="3yr-bleaching risk",
                                      TCI5_dhw="5yr-bleaching risk",
                                      sqrt.DINe="Environ DIN",
                                      mPAR="PAR",
                                      log.DINriv="River DIN",
                                      log.sal="Salinity Index",
                                      sqrt.MUD="Sediments",
                                      TCI5="5yr-acute events",
                                      lagCOTSp="lagOutbreaks",
                                      lagCIacute="lagAcute")

rownames(CPF_vis)=CPF_vis$predictor
##Order columns and remove predictor column
CPF_vis=CPF_vis[c(2,1,3,4,5)]
##Convert to matrix
CPF_vis <- data.matrix(CPF_vis)
##Transpose for heatmap
CPF_vis<-t(CPF_vis)

ordpred=c("lagCOTS","lagCycl","lagDHW","lagOutbreaks","lagAcute","3yr-bleaching risk","5yr-outbreaks",
          "5yr-storms","5yr-bleaching risk","5yr-acute events","Chla","Environ DIN","PAR","River DIN","Salinity Index","Sediments","zone")

##----------------------------------------------------
##Heatmap of Variable Importance Scores: Figure S6.A
##----------------------------------------------------

dev.new()
ggsave(file="VarImp_CPF.jpeg",width = 9, height = 6, units = c("in"),dpi = 300,
       heatmap.2(CPF_vis[,ordpred],notecex=0.3,  dendrogram ="none",
                 col=colorRampPalette(c("white","#D7FADE", "#A1FC9C", "#17D136", "#2B9C09"))(15),
                 main="Central performance\n (1992-2017)",
                 trace="none",key.title = "",keysize=0.5,
                 notecol="black",key=T,
                 sepcolor = "black",margins=c(11,11), lhei=c(2,6),lwid=c(2,6),cexRow = 2,cexCol = 2,
                 Rowv=FALSE,Colv=FALSE))
dev.off()

#-----------------------------------
##Gather data for  Best Models
#-----------------------------------

CPF<- do.call("rbind",mget(ls(pattern = "2AICc*")))
CPFtable=CPF[c(1,3,9,11,5,7)]
CPFtable$Group <- substr(rownames(CPFtable), 1, 3)

##Table S3.Central
CPFtable%<>%as.data.frame(row.names = NULL)%>%
  rename(predictors=modname)%>%
  mutate(Model="CPF")%>%
  arrange(match(Group, c("Total","ACT", "ACB", "CBR","MSE", desc(delta.AICc))))%>%
  dplyr::select(8,7,everything())%>%
  glimpse()
##Select Only the Most parsimonious Models
CPFtable%<>% slice (18,2:3,10,13:15,16)


#-----------------------------------------------------------------------
##Re-fit Best mod and Estimate Variance Contribution to select variables 
#-----------------------------------------------------------------------

CPFtable$predictors[CPFtable$Group=="Total"]
#[1] lagCycl+TCI5

mod=gam(Total_score~ s(lagCycl, k = 3, bs = "cr")+ s(TCI5, k = 3, bs = "cr")+
          s(time, bs = "cr",k = 5)+ s(FULLREEF_ID, bs = "re")+ s(REEF_SITE_NO,bs = "re") + 
          s(LATITUDE, k = 5), family = betar(link = "logit"), 
        data=Centdat_scores_ndist,method="REML")

### remove each predictor at a time
mp_cyc=update(mod, ~.-s(lagCycl, k = 3, bs = "cr"),sp=mod$sp[-1])
mp_ci5=update(mod, ~.-s(TCI5, k = 3, bs = "cr"),sp=mod$sp[-2])

listmods=setNames(lapply(ls(pattern="\\mp_"), get), paste(ls(pattern="\\mp_")))
listmods=listmods[c(1:5)]
names(listmods) <- c('mcyc', 'mtc5')


#function to estimate Explained variance per term 
devexp<-function(x,y){
  ((summary(y)$dev.expl)-(summary(x)$dev.expl))/(summary(y)$dev.expl)
}


pred_dev= as.data.frame(matrix(nrow = 2, ncol = 2))
colnames(pred_dev)=c("pred","devex")
for (i in 1:2){
  dev=devexp(listmods[[i]],mod)*100
  pred_dev[i, 1] <- names(listmods)[i]
  pred_dev[i, 2] <- dev
}
pred_dev[3,1]="total"
pred_dev[3,2]=pred_dev%>%summarize_if(is.numeric, sum, na.rm=TRUE)
pred_dev[]

rm(list=setdiff(ls(), c("CPFtable","Centdat_scores_ndist","mod","devexp")))


##--------------------------------------------------------------
##Get Model predictions and Plot effects in order of Importance
##---------------------------------------------------------------

## pred Cycl
tdcycl <- expand.grid(lagCycl=seq(min(Centdat_scores_ndist$lagCycl),max(Centdat_scores_ndist$lagCycl),length.out = 50),
                      TCI5=median(mod$model$TCI5),
                      time=mean(mod$model$time),
                      REEF_SITE_NO=(mod$model$REEF_SITE_NO),
                      FULLREEF_ID=(mod$model$FULLREEF_ID),
                      LATITUDE=mean(mod$model$LATITUDE))%>%
  distinct()%>%
  glimpse()

fits <- predict.gam(mod, newdata=tdcycl, type='response', se.fit=T, exclude= c('s(REEF_SITE_NO)','s(FULLREEF_ID)'))

pcyl = tdcycl%>%data.frame(fits)%>%
  group_by(lagCycl)%>%
  summarise(response=mean(fit), lower=mean(fit-(1.96*se.fit)),
            upper=mean(fit+(1.96*se.fit)))%>%
  ungroup()

#####
tdtc5 <- expand.grid(TCI5=seq(min(Centdat_scores_ndist$TCI5),max(Centdat_scores_ndist$TCI5),length.out = 50),
                     lagCycl=median(mod$model$lagCycl),
                     time=mean(mod$model$time),
                     REEF_SITE_NO=(mod$model$REEF_SITE_NO),
                     FULLREEF_ID=(mod$model$FULLREEF_ID),
                     LATITUDE=mean(mod$model$LATITUDE))%>%
  distinct()%>%
  glimpse()

fits <- predict.gam(mod, newdata=tdtc5 , type='response', se.fit=T, exclude= c('s(REEF_SITE_NO)','s(FULLREEF_ID)'))

pci5 = tdtc5 %>%data.frame(fits)%>%
  group_by(TCI5)%>% 
  summarise(response=mean(fit), lower=mean(fit-(1.96*se.fit)),
            upper=mean(fit+(1.96*se.fit)))%>%
  ungroup()

##Backtransform response for plotting
backtransform.est <- function(x, n) {
  y <- (x * n - 0.5) / (n - 1)
  return(y)
}


listmods=setNames(lapply(ls(pattern="\\p"), get), paste(ls(pattern="\\p")))## 
names(listmods) <- c('pcyl', 'pti5')

for (i in names(listmods)){
  listmods[[i]]=mutate(listmods[[i]], response=(backtransform.est(response,nrow(Centdat_scores_ndist))),
                       lower=(backtransform.est(lower,nrow(Centdat_scores_ndist))),
                       upper=(backtransform.est(upper,nrow(Centdat_scores_ndist))))
}

# Plot effects in order of Importance
png(file="total_cpf_pred.png",
    width=900, height=250,res =100)
par(mfrow=c(1,4))

pcyl=listmods[[1]]
plot(pcyl$lagCycl^2, pcyl$response, type="n", lwd=3, ylim=(c(0,1)),
     main="",ylab="Performance",
     xlab = "Lagged-Cycl",cex.axis=2,cex.lab=2)
X=c(0,14)
Y1 <- c(0,0.41)
rect(X[1], Y1[1], X[2], Y1[2], border = "white", col = adjustcolor("red",alpha=0.2))
polygon(c(pcyl$lagCycl^2, rev(pcyl$lagCycl^2)), 
        c(pcyl$lower,rev(pcyl$upper)), col="darkorange",
        border=NA)
lines(pcyl$lagCycl^2, pcyl$response,  lwd=1)
quantile(Centdat_scores_ndist$lagCycl^2,probs=c(0.50,0.95))
abline(v=0, col="black", lty=2)
abline(v=5, col="black", lty=2)


pci5=listmods[[2]]
plot(pci5$TCI5, pci5$response, type="n", lwd=3, ylim=(c(0,1)),
     main="",ylab="Performance",
     xlab = "Lagged-acute events",cex.axis=2,cex.lab=2)
X=c(0,6.1)
Y1 <- c(0,0.41)
rect(X[1], Y1[1], X[2], Y1[2], border = "white", col = adjustcolor("red",alpha=0.2))
polygon(c(pci5$TCI5, rev(pci5$TCI5)), 
        c(pci5$lower,rev(pci5$upper)), col="grey",
        border=NA)
lines(pci5$TCI5, pci5$response,  lwd=1)
quantile(Centdat_scores_ndist$TCI5,probs=c(0.50,0.95))
abline(v=1, col="black", lty=2)
abline(v=4, col="black", lty=2)
dev.off()


rm(list=setdiff(ls(), c("CPFtable","Centdat_scores_ndist","backtransform.est","devexp")))

##----------------------------------
###ACROPORA TABULAR-Best Mod
##----------------------------------
CPFtable$predictors[CPFtable$Group=="ACT"]
#[1] lagCycl+log.sal+TCI5_COTS+TCI5_cycl

mod=gam(ACTO_score ~ s(lagCycl, k = 3, bs = "cr") + s(log.sal, k = 3, bs = "cr") + s(TCI5_COTS, k = 3, bs = "cr") + 
           s(TCI5_cycl, k = 3, bs = "cr") + s(time, bs = "cr", k = 5) + s(FULLREEF_ID, bs = "re") + s(REEF_SITE_NO,bs = "re")+ 
          s(LATITUDE, k = 5),family = betar(link = "logit"), data=Centdat_scores_ndist,method="REML")


mp_cyc=update(mod, ~.-s(lagCycl, k = 3, bs = "cr"),sp=mod$sp[-1])
mp_sal=update(mod, ~.-s(log.sal, k = 3, bs = "cr"),sp=mod$sp[-2])
mp_cots=update(mod, ~.-s(TCI5_COTS, k = 3, bs = "cr"),sp=mod$sp[-3])
mp_cyc5=update(mod, ~.-s(TCI5_cycl, k = 3, bs = "cr"),sp=mod$sp[-4])


## Contribution to Explained deviance
listmods=setNames(lapply(ls(pattern="\\mp_"), get), paste(ls(pattern="\\mp_")))## list of all model fits with names
listmods=listmods[c(1:5)]
names(listmods) <- c('mcots', 'mcyl', 'mcyc5', 'msal')

pred_dev= as.data.frame(matrix(nrow = 4, ncol = 2))
colnames(pred_dev)=c("pred","devex")
for (i in 1:4){
  dev=devexp(listmods[[i]],mod)*100
  pred_dev[i, 1] <- names(listmods)[i]
  pred_dev[i, 2] <- dev
}
pred_dev[5,1]="total"
pred_dev[5,2]=pred_dev%>%summarize_if(is.numeric, sum, na.rm=TRUE)

pred_dev[]


## Predict effects
tdcyc <- expand.grid(lagCycl=seq(min(Centdat_scores_ndist$lagCycl),max(Centdat_scores_ndist$lagCycl),length.out = 50),
                     log.sal=median(mod$model$log.sal),
                     TCI5_COTS=median(mod$model$TCI5_COTS),
                     TCI5_cycl=median(mod$model$TCI5_cycl),
                     time=median(mod$model$time),
                     REEF_SITE_NO=(mod$model$REEF_SITE_NO),
                     FULLREEF_ID=(mod$model$FULLREEF_ID),
                     LATITUDE=mean(mod$model$LATITUDE))%>%
  distinct()%>%
  glimpse()

fits <- predict.gam(mod, newdata=tdcyc, type='response', se.fit=T, exclude= c("s(REEF_SITE_NO)","s(FULLREEF_ID)"))

pcyl = tdcyc %>%data.frame(fits)%>%
  group_by(lagCycl)%>%
  summarise(response=mean(fit), lower=mean(fit-(1.96*se.fit)),
            upper=mean(fit+(1.96*se.fit)))%>%
  ungroup()


tdcyc2 <- expand.grid(TCI5_cycl=seq(min(Centdat_scores_ndist$TCI5_cycl),max(Centdat_scores_ndist$TCI5_cycl),length.out = 50),
                      log.sal=median(mod$model$log.sal),
                      TCI5_COTS=median(mod$model$TCI5_COTS),
                      lagCycl=median(mod$model$lagCycl),
                      time=median(mod$model$time),
                      REEF_SITE_NO=(mod$model$REEF_SITE_NO),
                      FULLREEF_ID=(mod$model$FULLREEF_ID),
                      LATITUDE=mean(mod$model$LATITUDE))%>%
  distinct()%>%
  glimpse()

fits <- predict.gam(mod, newdata=tdcyc2, type='response', se.fit=T, exclude= c("s(REEF_SITE_NO)","s(FULLREEF_ID)"))

pcyc5 = tdcyc2 %>%data.frame(fits)%>%
  group_by(TCI5_cycl)%>% 
  summarise(response=mean(fit), lower=mean(fit-(1.96*se.fit)),
            upper=mean(fit+(1.96*se.fit)))%>%
  ungroup()


tdcots <- expand.grid(TCI5_COTS=seq(min(Centdat_scores_ndist$TCI5_COTS),max(Centdat_scores_ndist$TCI5_COTS),length.out = 50),
                      log.sal=median(mod$model$log.sal),
                      TCI5_cycl=median(mod$model$TCI5_cycl),
                      lagCycl=median(mod$model$lagCycl),
                      time=median(mod$model$time),
                      REEF_SITE_NO=(mod$model$REEF_SITE_NO),
                      FULLREEF_ID=(mod$model$FULLREEF_ID),
                      LATITUDE=mean(mod$model$LATITUDE))%>%
  distinct()%>%
  glimpse()

fits <- predict.gam(mod, newdata=tdcots , type='response', se.fit=T, exclude= c("s(REEF_SITE_NO)","s(FULLREEF_ID)"))

pcots = tdcots  %>%data.frame(fits)%>%
  group_by(TCI5_COTS)%>% 
  summarise(response=mean(fit), lower=mean(fit-(1.96*se.fit)),
            upper=mean(fit+(1.96*se.fit)))%>%
  ungroup()



tdsal <- expand.grid(log.sal=seq(min(Centdat_scores_ndist$log.sal),max(Centdat_scores_ndist$log.sal),length.out = 50),
                     TCI5_COTS=median(mod$model$TCI5_COTS),
                     TCI5_cycl=median(mod$model$TCI5_cycl),
                     lagCycl=median(mod$model$lagCycl),
                     time=median(mod$model$time),
                     REEF_SITE_NO=(mod$model$REEF_SITE_NO),
                     FULLREEF_ID=(mod$model$FULLREEF_ID),
                     LATITUDE=mean(mod$model$LATITUDE))%>%
  distinct()%>%
  glimpse()

fits <- predict.gam(mod, newdata=tdsal , type='response', se.fit=T, exclude= c("s(REEF_SITE_NO)","s(FULLREEF_ID)"))

psal = tdsal %>%data.frame(fits)%>%
  group_by(log.sal)%>% 
  summarise(response=mean(fit), lower=mean(fit-(1.96*se.fit)),
            upper=mean(fit+(1.96*se.fit)))%>%
  ungroup()


listmods=setNames(lapply(ls(pattern="\\p"), get), paste(ls(pattern="\\p")))## 
names(listmods) <- c('pcots','pcyl','pcyl5','psal')

for (i in names(listmods)){
  listmods[[i]]=mutate(listmods[[i]], response=(backtransform.est(response,nrow(Centdat_scores_ndist))),
                       lower=(backtransform.est(lower,nrow(Centdat_scores_ndist))),
                       upper=(backtransform.est(upper,nrow(Centdat_scores_ndist))))
}

# Plot effects in order of Importance
png(file="acto_cpf_pred.png",
    width=900, height=250,res =100)
par(mfrow=c(1,4))


pcots=listmods[[1]]
plot(pcots$TCI5_COTS, pcots$response, type="n", lwd=3, ylim=(c(0,1)),
     main="",ylab="Performance",xlab = "5yr-COTS events",cex.axis=2,cex.lab=2)
X <- c(0,5)
Y1 <- c(0,0.41)
Y2 <- c(0.4,0.6)
rect(X[1], Y1[1], X[2], Y1[2], border = "white", col = adjustcolor("red",alpha=0.2))
polygon(c(pcots$TCI5_COTS, rev(pcots$TCI5_COTS)), 
        c(pcots$lower,rev(pcots$upper)), col="magenta",
        border=NA)
lines(pcots$TCI5_COTS, pcots$response,  lwd=1)
quantile(Centdat_scores_ndist$TCI5_COTS,probs=c(0.50,0.95))
abline(v=0, col="black", lty=2)
abline(v=2, col="black", lty=2)


pcyl=listmods[[2]]
plot(pcyl$lagCycl^2, pcyl$response, type="n", lwd=3, ylim=(c(0,1)),
     main="",ylab="Performance",
     xlab = "LagCyclone exposure",cex.axis=2,cex.lab=2)
polygon(c(pcyl$lagCycl^2, rev(pcyl$lagCycl^2)), 
        c(pcyl$lower,rev(pcyl$upper)), col="darkorange",
        border=NA)
lines(pcyl$lagCycl^2, pcyl$response,  lwd=1)
X=c(0,14)
rect(X[1], Y1[1], X[2], Y1[2], border = "white", col = adjustcolor("red",alpha=0.2))
quantile(Centdat_scores_ndist$lagCycl^2,probs=c(0.50,0.95))
abline(v=0, col="black", lty=2)
abline(v=5, col="black", lty=2)

pcyc5=listmods[[3]]
plot(pcyc5$TCI5_cycl, pcyc5$response, type="n", lwd=3, ylim=(c(0,1)),
     main="",ylab="Performance",
     xlab = "Lag 5yr-storms",cex.axis=2,cex.lab=2)
X=c(0,4)
rect(X[1], Y1[1], X[2], Y1[2], border = "white", col = adjustcolor("red",alpha=0.2))
polygon(c(pcyc5$TCI5_cycl, rev(pcyc5$TCI5_cycl)), 
        c(pcyc5$lower,rev(pcyc5$upper)), col="darkorange",
        border=NA)
lines(pcyc5$TCI5_cycl, pcyc5$response,  lwd=1)
quantile(Centdat_scores_ndist$TCI5_cycl,probs=c(0.50,0.95))
abline(v=1, col="black", lty=2)
abline(v=3, col="black", lty=2)


psal=listmods[[4]]

plot(psal$log.sal, psal$response, type="n", lwd=3, ylim=(c(0,1)),
     main="",ylab="Performance",
     xlab = "log(Salinity Index)",cex.axis=2,cex.lab=2)
X=c(-0.1,6.5)
rect(X[1], Y1[1], X[2], Y1[2], border = "white", col = adjustcolor("red",alpha=0.2))
polygon(c(psal$log.sal, rev(psal$log.sal)), 
        c(psal$lower,rev(psal$upper)), col="cadetblue1",
        border=NA)
lines(psal$log.sal, psal$response,  lwd=1)
quantile(Centdat_scores_ndist$log.sal,probs=c(0.50,0.95))
abline(v=0, col="black", lty=2)
abline(v=4, col="black", lty=2)
dev.off()


rm(list=setdiff(ls(), c("CPFtable","Centdat_scores_ndist","backtransform.est","devexp")))

##----------------------------------
###ACROPORA Branching-Best Mod
##----------------------------------
CPFtable$predictors[CPFtable$Group=="ACB"]
#[1] lagCycl+TCI5_COTS

mod=gam(ACBX_score~ s(lagCycl, k = 3, bs = "cr")+ s(TCI5_COTS, k = 3, bs = "cr")+
          s(time, bs = "cr",k = 5)+ s(FULLREEF_ID, bs = "re")+ s(REEF_SITE_NO,bs = "re") + 
          s(LATITUDE, k = 5), family = betar(link = "logit"), 
        data=Centdat_scores_ndist,method="REML")

### remove each predictor at a time
mp_cyc=update(mod, ~.-s(lagCycl, k = 3, bs = "cr"),sp=mod$sp[-1])
mp_cots=update(mod, ~.-s(TCI5_COTS, k = 3, bs = "cr"),sp=mod$sp[-2])


## Contribution to Explained deviance
listmods=setNames(lapply(ls(pattern="\\mp_"), get), paste(ls(pattern="\\mp_")))
listmods=listmods[c(1:5)]
names(listmods) <- c('mcyc', 'mcots')

pred_dev= as.data.frame(matrix(nrow = 2, ncol = 2))
colnames(pred_dev)=c("pred","devex")
for (i in 1:2){
  dev=devexp(listmods[[i]],mod)*100
  pred_dev[i, 1] <- names(listmods)[i]
  pred_dev[i, 2] <- dev
}
pred_dev[3,1]="total"
pred_dev[3,2]=pred_dev%>%summarize_if(is.numeric, sum, na.rm=TRUE)
pred_dev[]

rm(list=setdiff(ls(), c("CPFtable","Centdat_scores_ndist","mod","devexp")))



## Predict Effects
tdcots <- expand.grid(TCI5_COTS=seq(min(Centdat_scores_ndist$TCI5_COTS),max(Centdat_scores_ndist$TCI5_COTS),length.out = 50),
                      lagCycl=median(mod$model$lagCycl),
                      time=mean(mod$model$time),
                      REEF_SITE_NO=(mod$model$REEF_SITE_NO),
                      FULLREEF_ID=(mod$model$FULLREEF_ID),
                      LATITUDE=mean(mod$model$LATITUDE))%>%
  distinct()%>%
  glimpse()

fits <- predict.gam(mod, newdata=tdcots, type='response', se.fit=T, exclude= c("s(REEF_SITE_NO)","s(FULLREEF_ID)"))

pcots = tdcots %>%data.frame(fits)%>%
  group_by(TCI5_COTS)%>%
  summarise(response=mean(fit), lower=mean(fit-(1.96*se.fit)),
            upper=mean(fit+(1.96*se.fit)))%>%
  ungroup()



tdc <- expand.grid(lagCycl=seq(min(Centdat_scores_ndist$lagCycl),max(Centdat_scores_ndist$lagCycl),length.out = 50),
                   TCI5_COTS=median(mod$model$TCI5_COTS),
                   time=mean(mod$model$time),
                   REEF_SITE_NO=(mod$model$REEF_SITE_NO),
                   FULLREEF_ID=(mod$model$FULLREEF_ID),
                   LATITUDE=mean(mod$model$LATITUDE))%>%
  distinct()%>%
  glimpse()

fits <- predict.gam(mod, newdata=tdc, type='response', se.fit=T, exclude= c("s(REEF_SITE_NO)","s(FULLREEF_ID)"))

pcyc = tdc %>%data.frame(fits)%>%
  group_by(lagCycl)%>% 
  summarise(response=mean(fit), lower=mean(fit-(1.96*se.fit)),
            upper=mean(fit+(1.96*se.fit)))%>%
  ungroup()



listmods=setNames(lapply(ls(pattern="\\p"), get), paste(ls(pattern="\\p")))## 
names(listmods) <- c('pcots','pcyc')

for (i in names(listmods)){
  listmods[[i]]=mutate(listmods[[i]], response=(backtransform.est(response,nrow(Centdat_scores_ndist))),
                       lower=(backtransform.est(lower,nrow(Centdat_scores_ndist))),
                       upper=(backtransform.est(upper,nrow(Centdat_scores_ndist))))
}

# Plot effects in order of Importance
png(file="acbx_cpf_pred.png",
    width=900, height=250,res =100)
par(mfrow=c(1,4))

pcots=listmods[[1]]
plot(pcots$TCI5_COTS, pcots$response, type="n", lwd=3, ylim=(c(0,1)),
     main="",ylab="Performance",xlab = "5yr-COTS events",cex.axis=2,cex.lab=2)
X <- c(0,5)
Y1 <- c(0,0.41)
Y2 <- c(0.4,0.6)

rect(X[1], Y1[1], X[2], Y1[2], border = "white", col = adjustcolor("red",alpha=0.2))
polygon(c(pcots$TCI5_COTS, rev(pcots$TCI5_COTS)), 
        c(pcots$lower,rev(pcots$upper)), col="magenta",
        border=NA)
lines(pcots$TCI5_COTS, pcots$response,  lwd=1)
quantile(Centdat_scores_ndist$TCI5_COTS,probs=c(0.50,0.95))
abline(v=0, col="black", lty=2)
abline(v=2, col="black", lty=2)
dev.off()


##----------------------------------
##Other Branching-Best Mod
##----------------------------------
CPFtable$predictors[CPFtable$Group=="CBR"]
#[1] lagCycl+lagDHW+log.DINriv+sqrt.DINe+TCI5_COTS

mod=gam(CBRN_score~s(lagCycl, k = 3, bs = "cr")+s(lagDHW, k = 3, bs = "cr")+s(log.DINriv, k = 3, bs = "cr")+s(sqrt.DINe, k = 3, bs = "cr")+  s(TCI5_COTS, k = 3, bs = "cr")+
          s(time, bs = "cr", k = 5)+s(FULLREEF_ID, bs = "re")+s(REEF_SITE_NO, bs = "re")+ s(LATITUDE, k = 5), 
        family = betar(link = "logit"),data=Centdat_scores_ndist,method="REML")


### remove each predictor at a time
mp_cyc=update(mod, ~.-s(lagCycl, k = 3, bs = "cr"),sp=mod$sp[-1])
mp_dhw=update(mod, ~.-s(lagDHW, k = 3, bs = "cr"),sp=mod$sp[-2])## Alternative model with zone instead
mp_dinr=update(mod, ~.-s(log.DINriv, k = 3, bs = "cr"),sp=mod$sp[-3])
mp_dine=update(mod, ~.-s(sqrt.DINe, k = 3, bs = "cr"),sp=mod$sp[-4])
mp_cots=update(mod, ~.-s(TCI5_COTS, k = 3, bs = "cr"),sp=mod$sp[-5])


## Contribution to Explained deviance
listmods=setNames(lapply(ls(pattern="\\mp_"), get), paste(ls(pattern="\\mp_")))
listmods=listmods[c(1:5)]
names(listmods) <- c('mcyc', 'mcots','mp_dhw','mp_dine','mp_dinr')

pred_dev= as.data.frame(matrix(nrow = 5, ncol = 2))
colnames(pred_dev)=c("pred","devex")
for (i in 1:5){
  dev=devexp(listmods[[i]],mod)*100
  pred_dev[i, 1] <- names(listmods)[i]
  pred_dev[i, 2] <- dev
}
pred_dev[6,1]="total"
pred_dev[6,2]=pred_dev%>%summarize_if(is.numeric, sum, na.rm=TRUE)
pred_dev[]

rm(list=setdiff(ls(), c("CPFtable","Centdat_scores_ndist","mod","devexp")))


##Predict Effects
####################
tdcots <- expand.grid(TCI5_COTS=seq(min(Centdat_scores_ndist$TCI5_COTS),max(Centdat_scores_ndist$TCI5_COTS),length.out = 50),
                      log.DINriv=median(mod$model$log.DINriv),
                      lagCycl=median(mod$model$lagCycl),
                      sqrt.DINe=median(mod$model$sqrt.DINe),
                      lagDHW=median(mod$model$lagDHW),
                      #zone=levels(mod$model$zone),##Alternative model
                      time=mean(mod$model$time),
                      REEF_SITE_NO=(mod$model$REEF_SITE_NO),
                      FULLREEF_ID=(mod$model$FULLREEF_ID),
                      LATITUDE=mean(mod$model$LATITUDE))%>%
  distinct()%>%
  glimpse()

fits <- predict.gam(mod, newdata=tdcots, type='response', se.fit=T, exclude= c("s(REEF_SITE_NO)","s(FULLREEF_ID)"))

pcots = tdcots %>%data.frame(fits)%>%
  group_by(TCI5_COTS)%>%
  summarise(response=mean(fit), lower=mean(fit-(1.96*se.fit)),
            upper=mean(fit+(1.96*se.fit)))%>%
  ungroup()


####################
tddin <- expand.grid(log.DINriv=seq(min(Centdat_scores_ndist$log.DINriv),max(Centdat_scores_ndist$log.DINriv),length.out = 50),
                     TCI5_COTS=median(mod$model$TCI5_COTS),
                     lagCycl=median(mod$model$lagCycl),
                     sqrt.DINe=median(mod$model$sqrt.DINe),
                     lagDHW=median(mod$model$lagDHW),
                     #zone=levels(mod$model$zone),##Alternative model
                     time=mean(mod$model$time),
                     REEF_SITE_NO=(mod$model$REEF_SITE_NO),
                     FULLREEF_ID=(mod$model$FULLREEF_ID),
                     LATITUDE=mean(mod$model$LATITUDE))%>%
  distinct()%>%
  glimpse()

## Exclude reef sites for prediction to reduce CI 
fits <- predict.gam(mod, newdata=tddin, type='response', se.fit=T, exclude= c("s(REEF_SITE_NO)","s(FULLREEF_ID)"))

pdinr = tddin %>%data.frame(fits)%>%
  group_by(log.DINriv)%>% 
  summarise(response=mean(fit), lower=mean(fit-(1.96*se.fit)),
            upper=mean(fit+(1.96*se.fit)))%>%
  ungroup()



####################
tdcyc <- expand.grid(lagCycl=seq(min(Centdat_scores_ndist$lagCycl),max(Centdat_scores_ndist$lagCycl),length.out = 50),
                     TCI5_COTS=median(mod$model$TCI5_COTS),
                     log.DINriv=median(mod$model$log.DINriv),
                     sqrt.DINe=median(mod$model$sqrt.DINe),
                     lagDHW=median(mod$model$lagDHW),
                     #zone=levels(mod$model$zone),##Alternative model
                     time=mean(mod$model$time),
                     REEF_SITE_NO=(mod$model$REEF_SITE_NO),
                     FULLREEF_ID=(mod$model$FULLREEF_ID),
                     LATITUDE=mean(mod$model$LATITUDE),
                     LONGITUDE=mean(mod$model$LONGITUDE))%>%
  distinct()%>%
  glimpse()

fits <- predict.gam(mod, newdata=tdcyc, type='response', se.fit=T, exclude= c("s(REEF_SITE_NO)","s(FULLREEF_ID)"))

pcyc = tdcyc %>%data.frame(fits)%>%
  group_by(lagCycl)%>% 
  summarise(response=mean(fit), lower=mean(fit-(1.96*se.fit)),
            upper=mean(fit+(1.96*se.fit)))%>%
  ungroup()

####################
tddine <- expand.grid(sqrt.DINe=seq(min(Centdat_scores_ndist$sqrt.DINe),max(Centdat_scores_ndist$sqrt.DINe),length.out = 50),
                      TCI5_COTS=median(mod$model$TCI5_COTS),
                      log.DINriv=median(mod$model$log.DINriv),
                      lagCycl=median(mod$model$lagCycl),
                      lagDHW=median(mod$model$lagDHW),
                      #zone=levels(mod$model$zone),##Alternative model
                      time=mean(mod$model$time),
                      REEF_SITE_NO=(mod$model$REEF_SITE_NO),
                      FULLREEF_ID=(mod$model$FULLREEF_ID),
                      LATITUDE=mean(mod$model$LATITUDE),
                      LONGITUDE=mean(mod$model$LONGITUDE))%>%
  distinct()%>%
  glimpse()

fits <- predict.gam(mod, newdata=tddine, type='response', se.fit=T, exclude= c("s(REEF_SITE_NO)","s(FULLREEF_ID)"))

pdine = tddine %>%data.frame(fits)%>%
  group_by(sqrt.DINe)%>% 
  summarise(response=mean(fit), lower=mean(fit-(1.96*se.fit)),
            upper=mean(fit+(1.96*se.fit)))%>%
  ungroup()


listmods=setNames(lapply(ls(pattern="\\p"), get), paste(ls(pattern="\\p")))## 
names(listmods) <- c('pcots','pcyc','pdine','pdinr')

for (i in names(listmods)){
  listmods[[i]]=mutate(listmods[[i]], response=(backtransform.est(response,nrow(Centdat_scores_ndist))),
                       lower=(backtransform.est(lower,nrow(Centdat_scores_ndist))),
                       upper=(backtransform.est(upper,nrow(Centdat_scores_ndist))))
}


# Plot effects in order of Importance
png(file="cbrn_cpf_pred.png",
    width=900, height=250,res =100)
par(mfrow=c(1,4))

pcots=listmods[[1]]
plot(pcots$TCI5_COTS, pcots$response, type="n", lwd=3, ylim=(c(0,1)),
     main="",ylab="Performance",xlab = "5yr-COTS events",cex.axis=2,cex.lab=2)
X <- c(0,5)
Y1 <- c(0,0.41)
Y2 <- c(0.4,0.6)
rect(X[1], Y1[1], X[2], Y1[2], border = "white", col = adjustcolor("red",alpha=0.2))
polygon(c(pcots$TCI5_COTS, rev(pcots$TCI5_COTS)), 
        c(pcots$lower,rev(pcots$upper)), col="magenta",
        border=NA)
lines(pcots$TCI5_COTS, pcots$response,  lwd=1)
quantile(Centdat_scores_ndist$TCI5_COTS,probs=c(0.50,0.95))
abline(v=0, col="black", lty=2)
abline(v=2, col="black", lty=2)

pcyc=listmods[[2]]
plot(pcyc$lagCycl^2, pcyc$response, type="n", lwd=3, ylim=(c(0,1)),
     main="",ylab="Performance",
     xlab = "Lag-Cyclone exposure",cex.axis=2,cex.lab=2)
X=c(0,14)
rect(X[1], Y1[1], X[2], Y1[2], border = "white", col = adjustcolor("red",alpha=0.2))
polygon(c(pcyc$lagCycl^2, rev(pcyc$lagCycl^2)), 
        c(pcyc$lower,rev(pcyc$upper)), col="darkorange",
        border=NA)
lines(pcyc$lagCycl^2, pcyc$response,  lwd=1)
quantile(Centdat_scores_ndist$lagCycl^2,probs=c(0.50,0.95))
abline(v=0, col="black", lty=2)
abline(v=5, col="black", lty=2)

pdinr=listmods[[4]]
plot(exp(pdinr$log.DINriv)-1, pdinr$response, type="n", lwd=3, ylim=(c(0,1)),
     main="",ylab="Performance",
     xlab = "River DIN (mg m3)",cex.axis=2,cex.lab=2)
X=c(0,210)
rect(X[1], Y1[1], X[2], Y1[2], border = "white", col = adjustcolor("red",alpha=0.2))
polygon(c(exp(pdinr$log.DINriv)-1, rev(exp(pdinr$log.DINriv)-1)), 
        c(pdinr$lower,rev(pdinr$upper)), col="darkgoldenrod",
        border=NA)
lines(exp(pdinr$log.DINriv)-1, pdinr$response,  lwd=1)
abline(v=1.8, col="black", lty=2)
quantile(exp(Centdat_scores_ndist$log.DINriv)-1,probs=c(0.50,0.95))
abline(v=0.08, col="black", lty=2)
abline(v=55, col="black", lty=2)


pdine=listmods[[3]]
plot(pdine$sqrt.DINe^2, pdine$response, type="n", lwd=3, ylim=(c(0,1)),
     main="",ylab="Performance",
     xlab = "DIN (mg m3)",cex.axis=2,cex.lab=2)
X=c(0,2.49)
rect(X[1], Y1[1], X[2], Y1[2], border = "white", col = adjustcolor("red",alpha=0.2))
polygon(c(pdine$sqrt.DINe^2, rev(pdine$sqrt.DINe^2)), 
        c(pdine$lower,rev(pdine$upper)), col="darkgoldenrod",
        border=NA)
lines(pdine$sqrt.DINe^2, pdine$response,  lwd=1)
quantile(Centdat_scores_ndist$sqrt.DINe^2,probs=c(0.50,0.95))
abline(v=0.58, col="black", lty=2)
abline(v=1.27, col="black", lty=2)
dev.off()


rm(list=setdiff(ls(), c("CPFtable","Centdat_scores_ndist","backtransform.est","devexp")))


##----------------------------------
##MSE-Best Mod
##----------------------------------
CPFtable$predictors[CPFtable$Group=="MSE"]
#[1] lagCycl+log.sal+TCI5_cycl

mod=gam(MSE_score~s(lagCycl, k = 3, bs = "cr")+s(log.sal, k = 3, bs = "cr")+ s(TCI5_cycl, k = 3, bs = "cr")+
          s(time, bs = "cr", k = 5)+s(FULLREEF_ID, bs = "re")+s(REEF_SITE_NO, bs = "re")+ 
          s(LATITUDE, k = 5), family=betar(link="logit"), data=Centdat_scores_ndist,method="REML")

##Remove each at a time
mp_cyc=update(mod, ~.-s(lagCycl, k = 3, bs = "cr"), sp=mod$sp[-1])
mp_sal=update(mod, ~.-s(log.sal, k = 3, bs = "cr"), sp=mod$sp[-2])
mp_cyc5=update(mod, ~.-s(TCI5_cycl, k = 3, bs = "cr"), sp=mod$sp[-3])


## Contribution to Explained deviance
listmods=setNames(lapply(ls(pattern="\\mp_"), get), paste(ls(pattern="\\mp_")))
listmods=listmods[c(1:5)]
names(listmods) <- c('mcyc', 'mp_cyc5','mp_sal')

pred_dev= as.data.frame(matrix(nrow = 3, ncol = 2))
colnames(pred_dev)=c("pred","devex")
for (i in 1:3){
  dev=devexp(listmods[[i]],mod)*100
  pred_dev[i, 1] <- names(listmods)[i]
  pred_dev[i, 2] <- dev
}
pred_dev[4,1]="total"
pred_dev[4,2]=pred_dev%>%summarize_if(is.numeric, sum, na.rm=TRUE)
pred_dev[]



rm(list=setdiff(ls(), c("CPFtable","Centdat_scores_ndist","mod","devexp")))


# predict effects

tdcyc <- expand.grid(lagCycl=seq(min(Centdat_MSE$lagCycl),max(Centdat_MSE$lagCycl),length.out = 50),
                     TCI5_cycl=median(mod$model$TCI5_cycl),
                     log.sal=median(mod$model$log.sal),
                     time=median(mod$model$time),
                     REEF_SITE_NO=(mod$model$REEF_SITE_NO),
                     FULLREEF_ID=(mod$model$FULLREEF_ID),
                     LATITUDE=mean(mod$model$LATITUDE))%>%
  distinct()%>%
  glimpse()

fits <- predict.gam(mod, newdata=tdcyc, type='response', se.fit=T, exclude= c("s(REEF_SITE_NO)","s(FULLREEF_ID)"))

pcyc = tdcyc %>%data.frame(fits1)%>%
  group_by(lagCycl)%>% 
  summarise(response=mean(fit), lower=mean(fit-(1.96*se.fit)),
            upper=mean(fit+(1.96*se.fit)))%>%
  ungroup()


####################
tdcyc5 <- expand.grid(TCI5_cycl=seq(min(Centdat_MSE$TCI5_cycl),max(Centdat_MSE$TCI5_cycl),length.out = 50),
                      lagCycl=median(mod$model$lagCycl),
                      log.sal=median(mod$model$log.sal),
                      time=median(mod$model$time),
                      REEF_SITE_NO=(mod$model$REEF_SITE_NO),
                      FULLREEF_ID=(mod$model$FULLREEF_ID),
                      LATITUDE=mean(mod$model$LATITUDE))%>%
  distinct()%>%
  glimpse()

fits <- predict.gam(mod, newdata=tdcyc5, type='response', se.fit=T, exclude= c("s(REEF_SITE_NO)","s(FULLREEF_ID)"))

pcyc5 = tdcyc5 %>%data.frame(fits)%>%
  group_by(TCI5_cycl)%>% 
  summarise(response=mean(fit), lower=mean(fit-(1.96*se.fit)),
            upper=mean(fit+(1.96*se.fit)))%>%
  ungroup()

####################
tdsal <- expand.grid(log.sal=seq(min(Centdat_MSE$log.sal),max(Centdat_MSE$log.sal),length.out = 50),
                     lagCycl=mean(mod$model$lagCycl),
                     TCI5_cycl=mean(mod$model$TCI5_cycl),
                     time=mean(mod$model$time),
                     REEF_SITE_NO=(mod$model$REEF_SITE_NO),
                     FULLREEF_ID=(mod$model$FULLREEF_ID),
                     LATITUDE=mean(mod$model$LATITUDE),
                     LONGITUDE=mean(mod$model$LONGITUDE))%>%
  distinct()%>%
  glimpse()

fits <- predict.gam(mod, newdata=tdsal , type='response', se.fit=T, exclude= c("s(REEF_SITE_NO)","s(FULLREEF_ID)"))

psal = tdsal %>%data.frame(fits)%>%
  group_by(log.sal)%>%
  summarise(response=mean(fit), lower=mean(fit-(1.96*se.fit)),
            upper=mean(fit+(1.96*se.fit)))%>%
  ungroup()

# Plot effects in order of Importance

listmods=setNames(lapply(ls(pattern="\\p"), get), paste(ls(pattern="\\p")))## 
names(listmods) <- c('pcyc','pcyc5','psal')

for (i in names(listmods)){
  listmods[[i]]=mutate(listmods[[i]], response=(backtransform.est(response,nrow(Centdat_scores_ndist))),
                       lower=(backtransform.est(lower,nrow(Centdat_scores_ndist))),
                       upper=(backtransform.est(upper,nrow(Centdat_scores_ndist))))
}



png(file="mse_cpf_pred.png",
    width=900, height=250,res =100)
par(mfrow=c(1,4))

pcyc=listmods[[1]]
plot(pcyc$lagCycl^2, pcyc$response, type="n", lwd=3, ylim=(c(0,1)),
     main="",ylab="Performance",
     xlab = "Lag-Cyclone exposure",cex.axis=2,cex.lab=2)
X=c(0,14)
Y1 <- c(0,0.41)
rect(X[1], Y1[1], X[2], Y1[2], border = "white", col = adjustcolor("red",alpha=0.2))
polygon(c(pcyc$lagCycl^2, rev(pcyc$lagCycl^2)), 
        c(pcyc$lower,rev(pcyc$upper)), col="darkorange",
        border=NA)
lines(pcyc$lagCycl^2, pcyc$response,  lwd=1)
quantile(Centdat_MSE$lagCycl^2,probs=c(0.50,0.95))
abline(v=0, col="black", lty=2)
abline(v=5, col="black", lty=2)

pcyc5=listmods[[2]]
plot(pcyc5$TCI5_cycl, pcyc5$response, type="n", lwd=3, ylim=(c(0,1)),
     main="",ylab="Performance",
     xlab = "Lagged-storm events",cex.axis=2,cex.lab=2)
X=c(0,4.1)
Y1 <- c(0,0.41)
rect(X[1], Y1[1], X[2], Y1[2], border = "white", col = adjustcolor("red",alpha=0.2))
polygon(c(pcyc5$TCI5_cycl, rev(pcyc5$TCI5_cycl)), 
        c(pcyc5$lower,rev(pcyc5$upper)), col="darkorange",
        border=NA)
lines(pcyc5$TCI5_cycl, pcyc5$response,  lwd=1)
quantile(Centdat_MSE$TCI5_cycl,probs=c(0.50,0.95))
abline(v=1, col="black", lty=2)
abline(v=3, col="black", lty=2)
dev.off()

