#MLM HW4
# April 6th, 2016
library(LukeMLM)
library(lattice)
library(lme4)
# Just follow the hint
data(nlsy)
nwrk <- nlsy[nlsy$alcday >= 0 & nlsy$smkday >=0 & nlsy$mrjday >= 0,]
nwrk <- nwrk[1:2000,]
nwrk$pubid <- as.factor(nwrk$pubid)
nwrk$nonwhite
# the new data set is nwrk
summary(nwrk)
mrjday<-data.frame(nwrk$mrjday)

# From the hint given in the HW, it seems we can only catergorize race to white  and non-white
gmod1<-glmer(mrjday~+sex97+nonwhite+alcday+smkday+(1|pubid),data=nwrk, family=poisson)
gmod2<-glmer(mrjday~sex97*nonwhite*alcday + sex97*nonwhite*smkday + (1|pubid),data=nwrk, family=poisson)
summary(gmod1)
summary(gmod2)

###################
# Prediction graph 
#####################

## Note that: family=poisson, the link function is log, the inverse of log function is exp()
summary(nwrk$smkday)
smkday<-nwrk$smkday
mrjday<-nwrk$mrjday
summary(mrjday)
alcday<-nwrk$alcday
## (1) Manual type: 

x<-seq(0,30,0.5)
# intercept, sex, race, alchol 
y_fw<-exp(-2.723 + x*0.025)

y_mw<- exp(-2.72341 -0.391916+0+ (0.070433-0.060123)+(0.02511+0.0165)*x)

y_fnw<- exp(-2.723410+0-2.285982+(0.070433-0.087252 )+(0.025112+0.016454)*x)
y_mnw<- exp(-2.72341-0.391916-2.285982+(0.070433-0.087252+0.167951)-0.391916-0.060123
               +(0.02511+0.0165)*x)

plot(x,y_fw,type="l",col="blue",lty=1,xlim=c(0,30),ylim=c(0,0.3),main="Predicted  marijuana use changes (manual type)",
     xlab="number of cigarettes per day",ylab="marijuana use change")

points(x,y_mw, type="l",col="brown1",lty=1)
points(x,y_fnw, type="l",col="olivedrab",lty=1)
points(x,y_mnw, type="l",col="cyan",lty=1)
legend(c(0,9),c(0.25,0.3),c("White Female","White Male","Non-white Female","Non-white Male"),col=c("blue","brown1","olivedrab","cyan"),lty=c(1,1,1,1))
# legend means the description box size

## (2) Use predic function to draw prediction graph

summary(alcday)
newdata_wf <- data.frame(smkday = seq (0,30,.3),alcday= seq (0,30,.3), nonwhite = "White",sex97="Female")
summary(newdata_wf)
newdata_wm <- data.frame(smkday = seq (0,30,.3),alcday= seq (0,30,.3),nonwhite = "White",sex97="Male")
newdata_nwf <- data.frame (smkday = seq(0,30,.3),alcday= seq (0,30,.3),nonwhite = "Non-white",sex97="Female")
newdata_nwm <- data.frame (smkday = seq (0,30,.3),alcday= seq (0,30,.3),nonwhite = "Non-white",sex97="Male")

y_wf <-exp(predict(gmod2,newdata_wf,re.form=NA))
y_wm <-exp(predict(gmod2,newdata_wm,re.form=NA))
y_nwf<- exp(predict(gmod2,newdata_nwf,re.form=NA))
y_nwm<- exp(predict(gmod2,newdata_nwm,re.form=NA))

plot (newdata_wf$smkday,y_wf,type="l",col="blue",lty=1,
      xlim= c (0,30),ylim= c(0,0.5),
      main="Predicted marijuana use change (with prediction function)",xlab="cigarette use",
      ylab=" marijuana use changes")
points (newdata_wm$smkday,y_wm,type="l",lty=1,col="brown1")
points (newdata_nwf$smkday,y_nwf,type="l",lty=1,col="olivedrab")
points (newdata_nwm$smkday,y_nwm,type="l",lty=1,col="cyan")
legend(c(0,14),c(0.35,0.5),c("White Female","White Male","Non-white Female","Non-white Male"),col=c("blue","brown1","olivedrab","cyan"),lty=c(1,1,1,1))

#install.packages("stargazer")
library(stargazer)

stargazer(gmod1,gmod2, title="Results", align=TRUE)

################################################
#           Determine Overdispersion           #
#################################################
# Follow lecture notes:lab8b
# Since it is count data, we use family=poisson, and its link function is log.
# For the Poisson model, the variance is assumed to be equal to the mean. 


mean(nwrk$mrjday)   #[1] 1.7505
var(nwrk$mrjday)    #[1] 36.51151

# Reference link: http://glmm.wikidot.com/faq
# glmer() won't estimate overdispersion, 
#so we can use the fowlloing code to calculate the ratio


overdisp_fun <- function(model) {
  ## number of variance parameters in
  ## an n-by-n variance-covariance matrix
  vpars <- function(m) {
    nrow (m)*( nrow (m)+1)/2
  }
  model.df <- sum ( sapply ( VarCorr (model),vpars))+ length ( fixef (model))
  rdf <- nrow ( model.frame (model))-model.df
  rp <- residuals (model,type="pearson")
  Pearson.chisq <- sum (rp^2)
  prat <- Pearson.chisq/rdf
  pval <- pchisq (Pearson.chisq, df=rdf, lower.tail=FALSE)
  c (chisq=Pearson.chisq,ratio=prat,rdf=rdf,p=pval)
}

overdisp_fun(gmod1) # the main-effect model overdispersion index : 2.930726
overdisp_fun(gmod2) # ratio=2.974468

####################################
##   Zero-inflated (and zero-truncated) models to deal with excess zero counts
##   problems
#############################
install.packages ("R2admb")
install.package ("reshape")
install.packages ("glmmADMB",
                  repos= c ("http://glmmadmb.r-forge.r-project.org/repos",
                            getOption ("repos")),
                  type="source")
install.packages ("coefplot2",repos="http://www.math.mcmaster.ca/bolker/R", type="source")

library(glmmADMB) 
# See very clear description here: http://glmmadmb.r-forge.r-project.org/glmmADMB.html

nwrk$pubid<-as.factor(nwrk$pubid)

#(1) Basic Poisson model
fit_poiss<-glmmadmb(mrjday~alcday+smkday+sex97+nonwhite+(1|pubid),
                       data=nwrk,
                       zeroInflation=FALSE,
                       family="poisson")
summary (fit_poiss)

#(2) zero-inflated Poisson Model---"zeroInflation=TRUE": ignore overdispersion problem
fit_zipoiss<-glmmadmb(mrjday~alcday+smkday+sex97+nonwhite+(1|pubid),
                    data=nwrk,
                    zeroInflation=TRUE,
                    family="poisson")
summary (fit_zipoiss)

# (3) Negative-Binomial models: deal with overdispersion

fit_nb<-glmmadmb(mrjday~alcday+smkday+sex97+nonwhite+(1|pubid),
                      data=nwrk,
                      zeroInflation=FALSE,
                      family="nbinom")
summary (fit_nb)

# (4) Zero-Inflated Negative-Bionomial Model:i.e 
# deal with excess zero counts and overdispersion

fit_zinb<-glmmadmb(mrjday~alcday+smkday+sex97+nonwhite+(1|pubid),
                 data=nwrk,
                 zeroInflation=TRUE,
                 family="nbinom")
summary (fit_zinb)

install.packages("bbmle")
library(bbmle)

AICtab(fit_poiss,fit_zipoiss,fit_nb,fit_zinb)
install.packages("texreg")
library("texreg")
library(stargazer)
screenreg(list(fit_poiss,fit_nb,fit_zipoiss,fit_zinb))

library (coefplot2)
install.packages("coda")
library(coda)
## Loading required package: coda
# compare the parameters
vn <- c ("sex97","nonwhite","smkday","alcday")
coefplot2(list(Poisson=fit_poiss,
                 Zero_Inflated=fit_zipoiss,
                 Neg_bin=fit_nb,
                 ZINB=fit_zinb),
                 varnames=vn,
                legend=TRUE)
          
# Comparison
compare<-list("ML_Poisson"=fit_poiss,
              "NB_Poisson"=fit_nb,
              "ZI_Poisson"=fit_zipoiss,
              "ZINB"=fit_zinb)
sapply(compare, function(x) coef(x)[1:5])

cbind("ML_Poisson"=sqrt(diag(vcov(fit_poiss))),
      "NB_Poisson"=sqrt(diag(vcov(fit_nb))),
      "ZI_Poisson"=sqrt(diag(vcov(fit_zipoiss))),
      "ZINB"=sqrt(diag(vcov(fit_zinb))),
      sapply(compare[-1], function(x) sqrt(diag(vcov(x)))[1:5]))

      
rbind(logLik=sapply(compare, function(x) round(logLik(x),
                                               digits=0)),
      Df=sapply(compare, function(x) attr(logLik(x),"df")))     


