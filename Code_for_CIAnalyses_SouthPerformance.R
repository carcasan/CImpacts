##------------------------------------------------------------------------##
## This script reproduces the analyses done to Performance variables      ##
## using Full-subset analysis for each coral group on Southern reefs      ##
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

Southdat_scores_ndist%<>% mutate(
  ACBX_score=transform01(ACBX_score),
  ACTO_score=transform01(ACTO_score),
  CBRN_score=transform01(CBRN_score),
  MSE_score=transform01(MSE_score),
  Total_score=transform01(Change.score.mean))%>%
  drop_na()

##-----------------------------------------------------------------
##Full Subset analyses run separately for each coral group-
##-----------------------------------------------------------------


South.preds=c("lagDHW","lagCOTS","lagCycl","TCI5","sqrt.sal","sqrt.MUD","sqrt.DINe","Chla","mPAR","log.DINriv")
lnvars=c("lagCIacute","TCI5_COTS","TCI5_cycl","TCI5_dhw")

factor.vars= c("zone","lagCOTSp")
Southdat_scores_ndist$lagCOTSp=as.character(Southdat_scores_ndist$lagCOTSp)


ModelHC=gam(HC_score~s(lagDHW,bs='cr',k=3)+ s(time,bs='cr',k=5)+s(LATITUDE, k=5)+ 
              s(FULLREEF_ID,bs='re')+s(REEF_SITE_NO,bs='re'),family=betar(link="logit"),data=Southdat_scores_ndist)


South_modelHC=generate.model.set(use.dat=Southdat_scores_ndist,
                                 max.predictors=5, 
                                 test.fit=ModelHC, k=3, 
                                 pred.vars.cont=South.preds,
                                 linear.vars=lnvars,
                                 pred.vars.fact = factor.vars,
                                 factor.smooth.interactions=F,
                                 smooth.smooth.interactions=F,
                                 cov.cutoff=0.28,                
                                 null.terms="s(time,bs='cr',k=5)+s(FULLREEF_ID,bs='re')+s(REEF_SITE_NO,bs='re')+s(LATITUDE, k=5)")


STotal.list=fit.model.set(STotal_model.set)## FITS mod
##extract Variable Importance Scores
TOTAL.var.imp=STotal.list$variable.importance$aic$variable.weights.raw
# Model selection table and best models
mod.table.STotal=STotal.list$mod.data.out
mod.table.STotal=mod.table.STotal[order(mod.table.STotal$AICc),]
STotal.less.2AICc=mod.table.STotal[which(mod.table.STotal$delta.AICc<2),]

##-------------------
## Acropora tabular
##-------------------

Model.ACTO=gam(ACTO_score~s(lagDHW,bs='cr',k=3)+s(time, bs = "cr", k = 5) +
                 s(LATITUDE, k = 5)+ s(FULLREEF_ID,bs = "re")+ 
                 s(REEF_SITE_NO, bs = "re"),family=betar(link="logit"),data=Southdat_scores_ndist)


SActo_model.set=generate.model.set(use.dat=Southdat_scores_ndist,
                                   max.predictors=5, 
                                   test.fit=Model.ACTO, k=3, 
                                   linear.vars=lnvars,
                                   pred.vars.cont=South.preds,
                                   pred.vars.fact =factor.vars,
                                   factor.smooth.interactions= F,
                                   cyclic.vars=NA,
                                   smooth.smooth.interactions=F,
                                   cov.cutoff=0.28,              
                                   null.terms="s(time,bs='cr',k=5)+s(LATITUDE,k=5)+s(FULLREEF_ID,bs ='re')+s(REEF_SITE_NO,bs='re')")

SActo.list=fit.model.set(SActo_model.set)## FITS mod
##extract Variable Importance Scores
ACTO.var.imp=SActo.list$variable.importance$aic$variable.weights.raw
# Model selection table and best models
mod.table.SActo=SActo.list$mod.data.out
mod.table.SActo=mod.table.SActo[order(mod.table.SActo$AICc),]
SActo.less.2AICc=mod.table.SActo[which(mod.table.SActo$delta.AICc<2),]

##-------------------
## Acropora branching
##-------------------


Model.ACBX=gam(ACBX_score~s(lagDHW,bs='cr',k=3)+s(time, bs = "cr", k = 5) +
                 s(LATITUDE, k = 5)+
                 s(FULLREEF_ID,bs = "re")+ 
                 s(REEF_SITE_NO, bs = "re"),family=betar(link="logit"),data=Southdat_scores_ndist)

SAcbr_model.set=generate.model.set(use.dat=Southdat_scores_ndist,
                                   max.predictors=5, 
                                   test.fit=Model.ACBX, k=3, 
                                   pred.vars.cont=South.preds,
                                   pred.vars.fact =factor.vars,
                                   factor.smooth.interactions= F,
                                   cyclic.vars=NA,
                                   smooth.smooth.interactions=F,
                                   cov.cutoff=0.28,              
                                   null.terms="s(time,bs='cr',k=5)+s(LATITUDE,k=5)+s(FULLREEF_ID,bs ='re')+s(REEF_SITE_NO,bs='re')")


SAcbr.list=fit.model.set(SAcbr_model.set)## FITS mod
##extract Variable Importance Scores
ACBR.var.imp=SAcbr.list$variable.importance$aic$variable.weights.raw
# Model selection table and best models
mod.table.SAcbr=SAcbr.list$mod.data.out
mod.table.SAcbr=mod.table.SAcbr[order(mod.table.SAcbr$AICc),]
SAcbr.less.2AICc=mod.table.SAcbr[which(mod.table.SAcbr$delta.AICc<2),]

##-------------------
## Other branching corals
##-------------------

Model.CBRN=gam(CBRN_score~s(lagDHW,bs='cr',k=3)+s(time, bs = "cr", k = 5) +
                 s(LATITUDE, k = 5)+
                 s(FULLREEF_ID,bs = "re")+ 
                 s(REEF_SITE_NO, bs = "re"),family=betar(link="logit"),data=Southdat_scores_ndist)


Scbrn_model.set=generate.model.set(use.dat=Southdat_scores_ndist,
                                         max.predictors=5, 
                                         test.fit=Model.CBRN, k=3, 
                                         pred.vars.cont=South.preds,
                                         pred.vars.fact =factor.vars,
                                         factor.smooth.interactions= F,
                                         cyclic.vars=NA,
                                         smooth.smooth.interactions=F,
                                         cov.cutoff=0.28,              
                                         null.terms="s(time,bs='cr',k=5)+s(LATITUDE,k=5)+s(FULLREEF_ID,bs ='re')+s(REEF_SITE_NO,bs='re')")

Scbrn.list=fit.model.set(Scbrn_model.set)## FITS mod
##extract Variable Importance Scores
CBRN.var.imp=Scbrn.list$variable.importance$aic$variable.weights.raw
# Model selection table and best models
mod.table.Scbrn=Scbrn.list$mod.data.out
mod.table.Scbrn=mod.table.Scbrn[order(mod.table.Scbrn$AICc),]
Scbrn.less.2AICc=mod.table.Scbrn[which(mod.table.Scbrn$delta.AICc<2),]

##----------------------------
## Massive-sub-massive-encrusting
##----------------------------

Model.MSE=gam(MSE_score~s(lagDHW,bs='cr',k=3)+s(time, bs = "cr", k = 5) +
                s(LATITUDE, k = 5)+
                s(FULLREEF_ID,bs = "re")+ 
                s(REEF_SITE_NO, bs = "re"),family=betar(link="logit"),data=Southdat_scores_ndist)


SMse_model.set=generate.model.set(use.dat=Southdat_scores_ndist,
                                  max.predictors=5, 
                                  test.fit=Model.MSE, k=3, 
                                  pred.vars.cont=South.preds,
                                  pred.vars.fact =factor.vars,
                                  factor.smooth.interactions= F,
                                  cyclic.vars=NA,
                                  smooth.smooth.interactions=F,
                                  cov.cutoff=0.28,              
                                  null.terms="s(time,bs='cr',k=5)+s(LATITUDE,k=5)+s(FULLREEF_ID,bs ='re')+s(REEF_SITE_NO,bs='re')")

SMse.list=fit.model.set(SMse_model.set)## FITS mod
##extract Variable Importance Scores
MSE.var.imp=SMse.list$variable.importance$aic$variable.weights.raw
# Model selection table and best models
mod.table.SMse=SMse.list$mod.data.out
mod.table.SMse=mod.table.SMse[order(mod.table.SMse$AICc),]
SMse.less.2AICc=mod.table.SMse[which(mod.table.SMse$delta.AICc<2),]


#-----------------------------------
##Gather data from Models for Variable Importance Scores
#-----------------------------------

ALLcorals.vis<- do.call("cbind",mget(ls(pattern = "var.imp*")))
SPF_vis=as.data.frame(ALLcorals.vis)
names(SPF_vis)<- substr(colnames(SPF_vis), 1, 4)
SPF_vis$predictor <- rownames(SPF_vis)

SPF_vis$predictor=dplyr::recode(SPF_vis$predictor,
                                TCI5_COTS="5yr-outbreaks",
                                sqrt.DINe="Environ DIN",
                                mPAR="PAR",
                                log.DINriv="River DIN",
                                sqrt.sal="Salinity Index",
                                sqrt.MUD="Sediments",
                                TCI5="5yr-acute events",
                                TCI5_cycl="5yr-storms",
                                TCI5_dhw="5yr-bleaching risk",
                                lagCOTSp="lagOutbreaks",
                                lagCIacute="lagAcute")


rownames(SPF_vis)=SPF_vis$predictor
##Order columns and remove predictor column
SPF_vis=SPF_vis[c(2,1,3,4,5)]
##Convert to matrix
SPF_vis <- data.matrix(SPF_vis)
##Transpose for heatmap
SPF_vis<-t(SPF_vis)

ordpred=c("lagCOTS","lagCycl","lagDHW","lagOutbreaks","lagAcute","5yr-outbreaks","5yr-storms",
           "5yr-bleaching risk","5yr-acute events","Chla","Environ DIN","PAR","River DIN",
           "Salinity Index","Sediments","zone")


##----------------------------------------------------
##Heatmap of Variable Importance Scores: Figure S7.A
##----------------------------------------------------

dev.new()
ggsave(file="VarImp_SPF.jpeg",width = 9, height = 6, units = c("in"),dpi = 300,
       heatmap.2(SPF_vis[,ordpred],notecex=0.3,  dendrogram ="none",
                 col=colorRampPalette(c("white","#D7FADE", "#A1FC9C", "#17D136", "#2B9C09"))(15),
                 main="South performance\n (1992-2017)",
                 trace="none",key.title = "",keysize=0.5,
                 notecol="black",key=T,
                 sepcolor = "black",margins=c(11,11), lhei=c(2,6),lwid=c(2,6),cexRow = 2,cexCol = 2,
                 Rowv=FALSE,Colv=FALSE))
dev.off()

#-----------------------------------
##Gather data for  Best Models
#-----------------------------------

SPF<- do.call("rbind",mget(ls(pattern = "2AICc*")))
SPFtable=SPF[c(1,3,9,11,5,7)]
SPFtable$Group <- substr(rownames(SPFtable), 1, 3)

##Table S3.Central
SPFtable%<>%as.data.frame(row.names = NULL)%>%
  rename(predictors=modname)%>%
  mutate(Model="SPF")%>%
  arrange(match(Group, c("Total","ACT", "ACB", "CBR","MSE", desc(delta.AICc))))%>%
  dplyr::select(8,7,everything())%>%
  glimpse()


#-----------------------------------------------------------------------
##Re-fit Best mod and Estimate Variance Contribution to select variables 
#-----------------------------------------------------------------------

SPFtable$predictors[SPFtable$Group=="Total"]

mod=gam(Total_score~s(lagCOTS, k =3, bs = "cr") + s(lagCycl, k = 3, bs = "cr") +
          s(lagDHW, k = 3, bs = "cr") + s(sqrt.sal, k = 3, bs = "cr") + 
          s(time, bs = "cr",k = 5) + s(FULLREEF_ID, bs = "re") + s(REEF_SITE_NO,bs = "re") + 
          s(LATITUDE, k = 5), family = betar(link = "logit"), 
        data=Southdat_scores_ndist,method="REML")

### remove each predictor at a time
ms_cots=update(mod, ~.-s(lagCOTS, k = 3, bs = "cr"),sp=mod$sp[-1])
ms_cyc=update(mod, ~.-s(lagCycl, k = 3, bs = "cr"),sp=mod$sp[-2])
ms_dhw=update(mod, ~.-s(lagDHW, k = 3, bs = "cr"),sp=mod$sp[-3])
ms_sal=update(mod, ~.-s(sqrt.sal, k = 3, bs = "cr"),sp=mod$sp[-4])


listmods=setNames(lapply(ls(pattern="\\ms_"), get), paste(ls(pattern="\\ms_")))
listmods=listmods[c(1:4)]
names(listmods) <- c('mcots', 'mcyc','mdhw','msal')


#function to estimate Explained variance per term 
devexp<-function(x,y){
  ((summary(y)$dev.expl)-(summary(x)$dev.expl))/(summary(y)$dev.expl)
}


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

rm(list=setdiff(ls(), c("SPFtable","Southdat_scores_ndist","mod","devexp")))


##---------------------------------------
##Get Model predictions and Plot effects 
##---------------------------------------

## pred Cycl
tdcycl <- expand.grid(lagCycl=seq(min(Southdat_scores_ndist$lagCycl),max(Southdat_scores_ndist$lagCycl),length.out = 50),
                      lagCOTS=median(mod$model$lagCOTS),
                      lagDHW=median(mod$model$lagDHW),
                      sqrt.MUD=median(mod$model$sqrt.MUD),
                      time=mean(mod$model$time),
                      REEF_SITE_NO=(mod$model$REEF_SITE_NO),
                      FULLREEF_ID=(mod$model$FULLREEF_ID),
                      LATITUDE=mean(mod$model$LATITUDE))%>%
  distinct()%>%
  glimpse()

## Exclude reef sites for prediction to reduce CI 
fits1 <- predict.gam(mod, newdata=tdcycl, type='response', se.fit=T, exclude= c('s(REEF_SITE_NO)','s(FULLREEF_ID)'))

pcyl = tdcycl%>%data.frame(fits1)%>%
  group_by(lagCycl)%>% #only change here
  summarise(response=mean(fit), lower=mean(fit-(1.96*se.fit)),
            upper=mean(fit+(1.96*se.fit)))%>%
  ungroup()

##PredCots
tdcots <- expand.grid(lagCOTS=seq(min(Southdat_scores_ndist$lagCOTS),max(Southdat_scores_ndist$lagCOTS),length.out = 50),
                      lagCycl=median(mod$model$lagCycl),
                      lagDHW=median(mod$model$lagDHW),
                      sqrt.MUD=median(mod$model$sqrt.MUD),
                      time=mean(mod$model$time),
                      REEF_SITE_NO=(mod$model$REEF_SITE_NO),
                      FULLREEF_ID=(mod$model$FULLREEF_ID),
                      LATITUDE=mean(mod$model$LATITUDE))%>%
  distinct()%>%
  glimpse()

## Exclude reef sites for prediction to reduce CI 
fits2 <- predict.gam(mod, newdata=tdcots, type='response', se.fit=T, exclude= c('s(REEF_SITE_NO)','s(FULLREEF_ID)'))

pcots = tdcots%>%data.frame(fits2)%>%
  group_by(lagCOTS)%>% #only change here
  summarise(response=mean(fit), lower=mean(fit-(1.96*se.fit)),
            upper=mean(fit+(1.96*se.fit)))%>%
  ungroup()

##PredDHW
tdhw <- expand.grid(lagDHW=seq(min(Southdat_scores_ndist$lagDHW),max(Southdat_scores_ndist$lagDHW),length.out = 50),
                    lagCycl=median(mod$model$lagCycl),
                    lagCOTS=median(mod$model$lagCOTS),
                    sqrt.MUD=median(mod$model$sqrt.MUD),
                    time=mean(mod$model$time),
                    REEF_SITE_NO=(mod$model$REEF_SITE_NO),
                    FULLREEF_ID=(mod$model$FULLREEF_ID),
                    LATITUDE=mean(mod$model$LATITUDE))%>%
  distinct()%>%
  glimpse()

## Exclude reef sites for prediction to reduce CI 
fits3 <- predict.gam(mod, newdata=tdhw, type='response', se.fit=T, exclude= c('s(REEF_SITE_NO)','s(FULLREEF_ID)'))

pdhw = tdhw%>%data.frame(fits3)%>%
  group_by(lagDHW)%>% #only change here
  summarise(response=mean(fit), lower=mean(fit-(1.96*se.fit)),
            upper=mean(fit+(1.96*se.fit)))%>%
  ungroup()


##PredSal
tdsal <- expand.grid(sqrt.sal=seq(min(Southdat_scores_ndist$sqrt.sal),max(Southdat_scores_ndist$sqrt.sal),length.out = 50),
                     lagCycl=median(mod$model$lagCycl),
                     lagCOTS=median(mod$model$lagCOTS),
                     lagDHW=median(mod$model$lagDHW),
                     time=mean(mod$model$time),
                     REEF_SITE_NO=(mod$model$REEF_SITE_NO),
                     FULLREEF_ID=(mod$model$FULLREEF_ID),
                     LATITUDE=mean(mod$model$LATITUDE))%>%
  distinct()%>%
  glimpse()

## Exclude reef sites for prediction to reduce CI 
fits5 <- predict.gam(mod, newdata=tdsal, type='response', se.fit=T, exclude= c('s(REEF_SITE_NO)','s(FULLREEF_ID)'))

psal = tdsal%>%data.frame(fits5)%>%
  group_by(sqrt.sal)%>% 
  summarise(response=mean(fit), lower=mean(fit-(1.96*se.fit)),
            upper=mean(fit+(1.96*se.fit)))%>%
  ungroup()


##Backtransform response for plotting
backtransform.est <- function(x, n) {
  y <- (x * n - 0.5) / (n - 1)
  return(y)
}

listmods=setNames(lapply(ls(pattern="\\p"), get), paste(ls(pattern="\\p")))## 
names(listmods) <- c('pcots','pcyl', 'pdhw','psal')#alphabetical


for (i in names(listmods)){
  listmods[[i]]=mutate(listmods[[i]], response=(backtransform.est(response,nrow(Southdat_scores_ndist))),
                       lower=(backtransform.est(lower,nrow(Southdat_scores_ndist))),
                       upper=(backtransform.est(upper,nrow(Southdat_scores_ndist))))
}

# Plot effects in order of Importance
png(file="total_spf_pred.png",
    width=900, height=250,res =100)
par(mfrow=c(1,4))


pcots=listmods[[1]]
plot(pcots$lagCOTS, pcots$response, type="n", lwd=3, ylim=(c(0,1)),
     main="",ylab="Performance",xlab = "lagCOTS density",cex.axis=2,cex.lab=2)
X <- c(0,4)
Y1 <- c(0,0.41)
rect(X[1], Y1[1], X[2], Y1[2], border = "white", col = adjustcolor("red",alpha=0.2))
polygon(c(pcots$lagCOTS, rev(pcots$lagCOTS)), 
        c(pcots$lower,rev(pcots$upper)), col="magenta",
        border=NA)
lines(pcots$lagCOTS, pcots$response,  lwd=1)
quantile(Southdat_scores_ndist$lagCOTS,probs=c(0.50,0.95))
abline(v=0, col="black", lty=2)
abline(v=0.5, col="black", lty=2)
text(x=3, y=0.85, "25%", cex=3)
mtext(side=1, line=5,"CoTS per tow\n1yr-lag", cex=2)


pcyc=listmods[[2]]
plot(pcyc$lagCycl^2, pcyc$response, type="n", lwd=3, ylim=(c(0,1)),
     main="",ylab="",xlab = "", cex.axis=2.1)
X=c(0,9)
Y1 <- c(0,0.41)

rect(X[1], Y1[1], X[2], Y1[2], border = "white", col = adjustcolor("red",alpha=0.2))
polygon(c(pcyc$lagCycl^2, rev(pcyc$lagCycl^2)), 
        c(pcyc$lower,rev(pcyc$upper)), col="darkorange",
        border=NA)
lines(pcyc$lagCycl^2, pcyc$response,  lwd=1)
quantile(Southdat_scores_ndist$lagCycl^2,probs=c(0.50,0.95))
abline(v=0, col="black", lty=2)
abline(v=0, col="black", lty=2)
text(x=6, y=0.85, "38%", cex=3)


pdhw=listmods[[3]]
plot(pdhw$lagDHW^2, pdhw$response, type="n", lwd=3, ylim=(c(0,1)),
     main="",ylab="",xlab = "", cex.axis=2.1)
X=c(0,11)
Y1 <- c(0,0.41)
rect(X[1], Y1[1], X[2], Y1[2], border = "white", col = adjustcolor("red",alpha=0.2))
polygon(c(pdhw$lagDHW^2, rev(pdhw$lagDHW^2)), 
        c(pdhw$lower,rev(pdhw$upper)), col="red",
        border=NA)
lines(pdhw$lagDHW^2, pdhw$response,  lwd=1)
quantile(Southdat_scores_ndist$lagDHW^2,probs=c(0.50,0.95))
abline(v=0.15, col="black", lty=2)
abline(v=3.7, col="black", lty=2)
text(x=6, y=0.85, "32%", cex=3)
dev.off()


rm(list=setdiff(ls(), c("SPFtable","Southdat_scores_ndist","backtransform.est","devexp")))




