# Data and analyses are included in the paper "Life history evolution and phenotypic plasticity in parasitic eyebrights (Euphrasia, Orobanchaceae)", bioRxiv doi: https://doi.org/10.1101/362400. In Review.  

# This script codes the generalised linear mixed effect models used in four different data sets.
# 1. "Manyhosts.csv" contains morphological measurements from one Euphrasia population grown with many hosts
# 2. "Manyspecies.csv" contains morphological measurements from many species grown on a single host
# 3. "Earlylate.csv" contains repeated growth measures at different times of year
# 4. "Wildcommon.csv" Comparison between common garden grown plants and wild collected plants

# Composed by Max Brown, 05-08-2018
# © Max Brown, 2018
# Minor updates and upload to: https://github.com/Euphrasiologist/phenotypic_plasticity_euphrasia in October 2019


# libraries needed

library(plyr) 
library(dplyr)
library(tidyr)
library(lme4)
library(MCMCglmm)
library(lattice)
library(Hmisc); library(RcmdrMisc)

library(data.table)
library(broom)
library(multcomp)
devtools::install_github("wilkelab/cowplot")
library(cowplot)
library(emmeans)

# custom functions used

# 1. Repeatability (variance explained) for poisson MCMCglmm
# needs model and random effect of interest (or residual variance = "units") stored in the VCV object
# modified from Greg Albery (pers comms) and https://github.com/rforge/rptr/tree/master/R
# The estimates are on the link scale, not the original (Poisson log) scale

MCMCReppois<-function(mod, y = "variable"){
  var.a      <- mod$VCV[,y]
  var.e      <- mod$VCV[,"units"]
  beta0      <- sapply(1:dim(mod$Sol)[1],function(z) mean(as.matrix(mod$X)%*%as.matrix(mod$Sol[z,1:ncol(mod$X)])))
  postR.link <- var.a/(var.a + var.e+log(1/exp(beta0)+1))
  R.link     <- posterior.mode( postR.link )
  # CI's
  CI.link    <- coda::HPDinterval(postR.link)[1,]
  # compile the list of results
  res 	   <- list(R.link=R.link, CI.link=CI.link)
  res.2 <- unlist(res)*100
  names(res.2) <- c("Link Scale Point Estimate", "CI Link Scale Lower", "CI Link Scale Upper")
  res.2 <- data.frame(res.2)
  colnames(res.2)[1] <- "%"
  return(res.2)
}

# 2. Repeatability (variance explained) for gaussian MCMCglmm
# needs model and random effect name from the VCV
# modified from Greg Albery (pers comms) and https://github.com/rforge/rptr/tree/master/R

MCMCRepnorm <- function(mod, y = "variable"){
  var.a      <- mod$VCV[,y]
  var.e      <- rowSums(mod$VCV)
  postR.link <- var.a/(var.e)
  # point estimate
  R.link     <- posterior.mode( postR.link )
  # CI's
  CI.link    <- coda::HPDinterval(postR.link)[1,]
  # list of results
  res 	   <- list(R.link=R.link, CI.link=CI.link)
  res.2 <- unlist(res)*100
  names(res.2) <- c("Point Estimate", "CI Lower", "CI Upper")
  res.2 <- data.frame(res.2)
  colnames(res.2)[1] <- "%"
  return(res.2)
}

# 3. Easy way to look at traces of MCMCglmm model output. First is fixed effect means, 
# second is the random effect variances. Note careful when specifying pr = T in MCMCglmm
# as the random effects will be saved in solutions, will lead to very slow (and useless)
# traces output.

traces<-function(mod){
  par.settings = list(strip.background=list(col="lightgrey"))
  trace.1<-xyplot(mod$Sol,
                  mean.list = mod$Sol,
                  panel = function(x, y, mean.list) {
                    panel.lines(x, y, col = "black")
                    panel.abline(h = mean(y),
                                 col.line = "red")
                  },
                  par.settings=par.settings)
  trace.2<-xyplot(mod$VCV,
                  mean.list = mod$VCV,
                  panel = function(x, y, mean.list) {
                    panel.lines(x, y)
                    panel.abline(h = mean(y),
                                 col.line = "red")
                  },
                  par.settings=par.settings)
  Newlist <- list(trace.1, trace.2)
  return(Newlist)
}

# 4. Models checked for violation of assumptions using viusual tools, histograms
# and overdisperion using 'overdisp_fun' http://bbolker.github.io/mixedmodels-misc/glmmFAQ.html

overdisp_fun <- function(model) {
  # number of variance parameters in 
  #   an n-by-n variance-covariance matrix
  vpars <- function(m) {
    nrow(m)*(nrow(m)+1)/2
  }
  model.df <- sum(sapply(VarCorr(model),vpars))+length(fixef(model))
  rdf <- nrow(model.frame(model))-model.df
  rp <- residuals(model,type="pearson")
  Pearson.chisq <- sum(rp^2)
  prat <- Pearson.chisq/rdf
  pval <- pchisq(Pearson.chisq, df=rdf, lower.tail=FALSE)
  c(chisq=Pearson.chisq,ratio=prat,rdf=rdf,p=pval)
}

# 5. plot tukey comparisons, modified from MCMCfixplot

g_lmfixplot <- function(model, tukey = FALSE , factor = mcp(...)){
  if(!tukey){
    if(any(attributes(model)$class == "glm")){
      
      rsum <- tidy(summary(glht(model, factor)))
      colnames(rsum)[1] <- "comparison"
      colnames(rsum)[6] <- "adj.p.value"
      rsum$conf.low <- unlist(c(rsum[,3] - 2*rsum[,4]))
      rsum$conf.high <- unlist(c(rsum[,3] + 2*rsum[,4]))
      
    } else rsum <- broom::tidy(TukeyHSD(aov(model)))
    
    rsum$ptext <- rep(NA, dim(rsum)[1])
    
    for(i in 1:dim(rsum)[1]){
      if(rsum$adj.p.value[i] <= 0.05 & rsum$adj.p.value[i] > 0.01){
        rsum$ptext[i] <- "*"
      } else if(rsum$adj.p.value[i] <= 0.01 & rsum$adj.p.value[i] > 0.001){
        rsum$ptext[i] <- "**"
      } else if(rsum$adj.p.value[i] <= 0.001 & rsum$adj.p.value[i] >= 0){
        rsum$ptext[i] <- "***"
      } else rsum$ptext[i] <- NA
    }
    
    rsum$enhance <- !is.na(as.character(rsum$ptext))
    
    ggplot(rsum,aes(x = comparison , y = estimate, alpha = enhance))+
      geom_pointrange(aes(ymin=conf.low,
                          ymax=conf.high))+
      geom_text(aes(label = ptext), nudge_x = 0.2) +
      coord_flip()+
      theme_bw()+
      labs(y = "Point estimates and confidence intervals", x = "Variables")+
      theme(legend.position = "none")
  } else warning("Functionality not yet supported!")
    
  
}


# 6. plot tukey comparisons from lmer models
 # list(pairwise ~ factor)
g_lmerfixplot <- function(model, factor = list(...), reorder = FALSE){
    rsum <-  emmeans(model, factor, adjust = "tukey")[[2]]
    rsum <- as.data.frame(rsum)
    
    rsum$ptext <- rep(NA, dim(rsum)[1])
    
    for(i in 1:dim(rsum)[1]){
      if(rsum$p.value[i] <= 0.05 & rsum$p.value[i] > 0.01){
        rsum$ptext[i] <- "*"
      } else if(rsum$p.value[i] <= 0.01 & rsum$p.value[i] > 0.001){
        rsum$ptext[i] <- "**"
      } else if(rsum$p.value[i] <= 0.001 & rsum$p.value[i] >= 0){
        rsum$ptext[i] <- "***"
      } else rsum$ptext[i] <- NA
    }
    
    rsum$enhance <- !is.na(as.character(rsum$ptext))
    
    if(reorder){
      ggplot(rsum,aes(x = reorder(contrast, -estimate) , y = estimate, alpha = enhance))+
        geom_pointrange(aes(ymin=estimate-2*SE,
                            ymax=estimate+2*SE))+
        geom_text(aes(label = ptext), nudge_x = 0.2) +
        coord_flip()+
        theme_bw()+
        labs(y = "Point estimates and confidence intervals", x = "Variables")+
        theme(legend.position = "none")
    } else
      ggplot(rsum,aes(x = contrast , y = estimate, alpha = enhance))+
      geom_pointrange(aes(ymin=estimate-2*SE,
                          ymax=estimate+2*SE))+
      geom_text(aes(label = ptext), nudge_x = 0.2) +
      coord_flip()+
      theme_bw()+
      labs(y = "Point estimates and confidence intervals", x = "Variables")+
      theme(legend.position = "none")
    
  }

# 7. Calculating P-values from annoying gaussian lmer objects
calc_pvals <- function(lmermod){
  if(attributes(lmermod)$class == "lmerMod"){
    dat<-summary(lmermod)$coefficients
    
    dat <- data.table::as.data.table(as.data.frame(dat), keep.rownames = T)
    colnames(dat)[1] <- "Levels"
    
    dat$pval_upperdf <- 1-pt(q = dat$`t value`, df = stats::nobs(lmermod))
    dat$pval_lowerdf <- 1-pt(q = dat$`t value`, df = summary(lmermod)$ngrps)
    
    return(dat)
  } else warning("Model not of class lmerMod")

}

# 8. specify the number of decimal places (https://stackoverflow.com/questions/3443687/formatting-decimal-places-in-r)

specify_decimal <- function(x, k) trimws(format(round(x, k), nsmall=k))

##### Phenotypic plasticity experiment #####

# read in the data
NBerwickme<- read.csv("./Data/manyhosts.csv")
setDT(NBerwickme)

# calculate means and sd's
PPmeanSD<-NBerwickme[, .(`Mean Height` = mean(Height),
               `SE Height` = sd(Height)/sqrt(.N),
               `Mean Number of Branches` = mean(No_branches),
               `SE Number of Branches` = sd(No_branches)/sqrt(.N),
               `Mean Node to Flower` = mean(Nodes_to_flower_inc_coty),
               `SE Node to Flower` = sd(Nodes_to_flower_inc_coty)/sqrt(.N),
               `Mean Number of Leaf Teeth` = mean(No_teeth_lower_floral_leaf),
               `SE Number of Leaf Teeth` = sd(No_teeth_lower_floral_leaf)/sqrt(.N),
               `Mean Corolla Length` = mean(Standard_corolla_l, na.rm = TRUE),
               `SE Corolla Length` = sd(Standard_corolla_l, na.rm = TRUE)/sqrt(.N),
               `Mean Cauline:Internode Ratio` = mean(Cauline_ratio),
               `SE Cauline:Internode Ratio` = sd(Cauline_ratio)/sqrt(.N),
               `Mean Julian Days to Flower` = mean(Julian_days_to_flower, na.rm = TRUE),
               `SE Julian Days to Flower` = sd(Julian_days_to_flower, na.rm = TRUE)/sqrt(.N)
), by = .(Host)]
# write to file
write.csv(x = cbind(PPmeanSD[,1], specify_decimal(PPmeanSD[,-1], 2)), file = "./Output/Phenotypic_plasticity/Means_SE_PP.csv")

# correlations
# subset the seven variables of interest
cormanyhost<-NBerwickme[,c(8, 9, 11, 12, 16, 21, 24)]
# convert to matrix
cormanyhost <- as.matrix(cormanyhost)
# compute correlation matrix
cormanyhost1 <- rcorr.adjust(x = cormanyhost, use = "pairwise.complete.obs")

# calculation of r
cormanyhost1$R$r
write.csv(x = cormanyhost1$R$r, file = "./Output/Phenotypic_plasticity/CorrelationsR.csv")
# associated p-values (adjusted)
cormanyhost1$P
write.csv(x = cormanyhost1$P, file = "./Output/Phenotypic_plasticity/CorrelationsP.csv")
# Run with host as a fixed effect. #

NBerwickme<- read.csv("./Data/Manyhosts.csv")
# change the names
levels(NBerwickme$Host)[c(2,3,5,6)] <- c("Equisetum arvense", "Festuca rubra", "Marchantia polymorpha", "No host")

## 1.1 Height

### Distribution of height across all individuals to see if it's bimodal ###

lattice::densityplot(NBerwickme$Height) #or
ggplot(NBerwickme, aes(x = Height)) +geom_jitter(aes(y =0), height = 0.001)+
  theme_bw()+geom_density()

# relevel so no host is baseline
NBerwickme$Host <- relevel(x = NBerwickme$Host, ref = "No host")

# models are specified below, height is log(height)
height1 <- lm(log(Height) ~ Host, data=NBerwickme)
# write to csv the coefficients of the model and the anova
write.csv(x = tidy(height1), 
          file = "./Output/Phenotypic_plasticity/PP_height.csv")
write.csv(anova(height1), file =  "./Output/Phenotypic_plasticity/PP_height_anova.csv") #effect of host is highly significant
# tukey tests here
tukeyheight <- tidy(TukeyHSD(aov(height1)))
write.csv(tukeyheight, file =  "./Output/Phenotypic_plasticity/PP_height_tukey.csv")
# and a plot for the tukey tests
g_lmfixplot(height1)

## 1.2 Nodes_to_flower_inc_coty

nodes1 <- glm(Nodes_to_flower_inc_coty ~ Host, 
              data=NBerwickme, family = "poisson")
summary(nodes1)
write.csv(x = tidy(nodes1), 
          file = "./Output/Phenotypic_plasticity/PP_nodes.csv")
write.csv(anova(nodes1, test = "Chisq"), file =  "./Output/Phenotypic_plasticity/PP_nodes_anova.csv") #effect of host is highly significant

tukeynodes <- tidy(TukeyHSD(aov(nodes1)))

tukeynodes <- tidy(summary(glht(nodes1, mcp(Host="Tukey"))))
write.csv(tukeynodes, file =  "./Output/Phenotypic_plasticity/PP_nodes_tukey.csv")

g_lmfixplot(nodes1, factor = mcp(Host = "Tukey"))

## 1.3 Standard_corolla_l

corol1 <- lm(Standard_corolla_l ~ Host, 
             data=NBerwickme[c(-20, -37),], na.action = na.omit)

summary(corol1)
write.csv(x = tidy(corol1), 
          file = "./Output/Phenotypic_plasticity/PP_corolla.csv")
write.csv(anova(corol1), file =  "./Output/Phenotypic_plasticity/PP_corolla_anova.csv") #effect of host is highly significant

tukeycorol <- tidy(TukeyHSD(aov(corol1)))
write.csv(tukeycorol, file =  "./Output/Phenotypic_plasticity/PP_corolla_tukey.csv")

g_lmfixplot(corol1)

## 1.4 Cauline_ratio

ratio1 <- lm(Cauline_ratio ~ Host, 
             data=NBerwickme[-54,], na.action = na.omit)
summary(ratio1)

write.csv(x = tidy(ratio1), 
          file = "./Output/Phenotypic_plasticity/PP_cauline_ratio.csv")
write.csv(anova(ratio1), file =  "./Output/Phenotypic_plasticity/PP_cauline_ratio_anova.csv") # effect of host is highly significant

tukeyratio <- tidy(TukeyHSD(aov(ratio1)))
write.csv(tukeyratio, file =  "./Output/Phenotypic_plasticity/PP_cauline_ratio_tukey.csv")

g_lmfixplot(ratio1)

## 1.5 Julian_days_to_flower


julian1 <- glm(Julian_days_to_flower ~ Host, 
               data=NBerwickme, na.action = na.omit, family = "poisson")

summary(julian1)

write.csv(x = tidy(julian1), 
          file = "./Output/Phenotypic_plasticity/PP_julian.csv")
write.csv(anova(julian1, test = "Chisq"), file =  "./Output/Phenotypic_plasticity/PP_julian_anova.csv") #effect of host is highly significant

tukeyjulian <-  tidy(summary(glht(julian1, mcp(Host="Tukey"))))
write.csv(tukeyjulian, file =  "./Output/Phenotypic_plasticity/PP_julian_tukey.csv")

g_lmfixplot(julian1, tukey = FALSE, mcp(Host = "Tukey"))

## 1.6 No_branches. Not a good model
branches1 <- glm(No_branches ~ Host, 
                 data=NBerwickme, family = "quasipoisson", na.action = na.omit)

summary(branches1)

write.csv(x = tidy(branches1), 
          file = "./Output/Phenotypic_plasticity/PP_branches.csv")
write.csv(anova(branches1, test = "F"), file =  "/Users/mbrown/Dropbox/Euphrasia 2016 common garden hosts vs pops/Experiments 2016/Manuscript/AJB/V2/Analyses/Phenotypic_plasticity/PP_branches_anova.csv") #effect of host is highly significant

tukeybranches <- tidy(summary(glht(branches1, mcp(Host="Tukey"))))
write.csv(tukeybranches, file =  "./Output/Phenotypic_plasticity/PP_branches_tukey.csv")


g_lmfixplot(branches1, mcp(Host = "Tukey"), tukey = FALSE)

## 1.7 No_teeth_lower_floral_leaf

teeth1 <- glm(No_teeth_lower_floral_leaf ~ Host, 
              data=NBerwickme, family = "poisson", na.action = na.omit)

write.csv(x = tidy(teeth1), 
          file = "./Output/Phenotypic_plasticity/PP_teeth.csv")
write.csv(anova(teeth1, test = "Chisq"), file =  "./Output/Phenotypic_plasticity/PP_teeth_anova.csv") #effect of host is highly significant

tukeyteeth <- tidy(summary(glht(teeth1, mcp(Host="Tukey"))))
write.csv(tukeyteeth, file =  "./Output/Phenotypic_plasticity/PP_teeth_tukey.csv")

g_lmfixplot(teeth1,  mcp(Host = "Tukey"), tukey = FALSE)

# Principal component analysis for many hosts

# subset data but include hosts
PCANB <- NBerwickme[, c(5, 8, 9, 11, 12, 16, 21, 24)]
# remove NA values
PCANB2<-PCANB[complete.cases(PCANB), ]
# compute PCA
Nberwickpca<- prcomp(PCANB2[,-1], 
                     center = TRUE,
                     scale. = TRUE)
# summary gives proportion of variance
summary(Nberwickpca)
# below gives variance component of each variable
aload <- abs(Nberwickpca$rotation)
percentvarhosts<-sweep(aload, 2, colSums(aload), "/")
# tidy output
percentvarhosts2 <- summary(Nberwickpca)
percentvarhosts2 <- percentvarhosts2$importance
percentvarhosts3<-rbind(percentvarhosts, percentvarhosts2)

write.csv(x = specify_decimal(percentvarhosts3, 3), file = "./Output/PCA/Hosts_PCA.csv")

one.pca.2 <- as.data.frame(Nberwickpca$x)
one.pca.2$Host <- PCANB2$Host

PCA3<- ggplot(one.pca.2, aes(x = PC1, y= PC2,
                      group = Host,
                      colour=Host))+
  geom_point(aes(group = Host, colour=Host), size =2)+
  stat_ellipse(aes(group = Host, colour=Host), alpha = 1)+
  ylab(label = "PC2 (18.5%)")+
  xlab(label = "PC1 (51.8%)")+
  theme_bw()
 

##### Species differences experiment #####

# load the data
manysp<- read.csv("./Data/Manyspecies.csv")
setDT(manysp)
## WITHOUT HYBRIDS ##
manysp2 <- manysp[!grepl(pattern = " x ", fixed = TRUE, x = manysp$Expert_ID_living),]
manysp2$Expert_ID_living <- factor(manysp2$Expert_ID_living)

# table of means and standard errors for each species for each trait
SPMeanSD <- manysp[, .(`Mean Height` = mean(Height),
            `SE Height` = sd(Height)/sqrt(.N),
            `Mean Number of Branches` = mean(No_branches),
            `SE Number of Branches` = sd(No_branches)/sqrt(.N),
            `Mean Node to Flower` = mean(Nodes_to_flower_inc_coty),
            `SE Node to Flower` = sd(Nodes_to_flower_inc_coty)/sqrt(.N),
            `Mean Number of Leaf Teeth` = mean(No_teeth_lower_floral_leaf),
            `SE Number of Leaf Teeth` = sd(No_teeth_lower_floral_leaf)/sqrt(.N),
            `Mean Corolla Length` = mean(Standard_corolla_l, na.rm = TRUE),
            `SE Corolla Length` = sd(Standard_corolla_l, na.rm = TRUE)/sqrt(.N),
            `Mean Cauline:Internode Ratio` = mean(Cauline_ratio),
            `SE Cauline:Internode Ratio` = sd(Cauline_ratio)/sqrt(.N),
            `Mean Julian Days to Flower` = mean(Julian_days_to_flower, na.rm = TRUE),
            `SE Julian Days to Flower` = sd(Julian_days_to_flower, na.rm = TRUE)/sqrt(.N)
), by = .(Expert_ID_living)]
# the fold variation
SPMeanSD[!grepl(" x ", Expert_ID_living), lapply(.SD, function(x) max(round(x), na.rm = TRUE)/min(round(x), na.rm = TRUE)), .SDcols = names(SPMeanSD)[!grepl("SE|Expert", names(SPMeanSD))]]
# scaled variance
SPMeanSD[!grepl(" x ", Expert_ID_living), lapply(.SD, function(x) var(scale(x, center = FALSE), na.rm = TRUE)), .SDcols = names(SPMeanSD)[!grepl("SE|Expert", names(SPMeanSD))]]

write.csv(x = cbind(SPMeanSD[,1], specify_decimal(SPMeanSD[,-1], 2)), 
          file = "./Output/Species_differences/Means_SE_SD.csv")



# correlations
# subset the seven variables of interest
cormanysp<-manysp[,c(6,7,9,10,12,16,17)]
# convert to matrix
cormanysp <- as.matrix(cormanysp)
# compute correlation matrix
cormanysp1 <- rcorr.adjust(x = cormanysp, use = "pairwise.complete.obs")
# calculation of r
cormanysp1$R$r
write.csv(x = cormanysp1$R$r, file = "./Output/Species_differences/CorrelationsR.csv")
# associated p-values
cormanysp1$P 
write.csv(x = cormanysp1$P, file = "./Output/Species_differences/CorrelationsP.csv")

# 1. Height
# full model
manysp2$E4EandSp <- as.factor(paste(manysp2$E4E, manysp2$Expert_ID_living, sep = "-"))

manysp2$Expert_ID_living <- relevel(x = manysp2$Expert_ID_living, ref = "E. arctica")

modmanysp1<- lmer(log(Height) ~ Expert_ID_living + (1 | E4EandSp),
                  data=manysp2)

# omitting fixed effects
modmanysp1.1<- lmer(log(Height) ~ 1 + (1 | E4EandSp),
                    data=manysp2)
# omitting random effects
modmanysp1.2<- lm(log(Height) ~ Expert_ID_living, data=manysp2)

# LRT's
write.csv(x = anova(modmanysp1.1, modmanysp1), 
          file = "./Output/Species_differences/Height/LRT_SP.csv") # effect of species
write.csv(x = anova(modmanysp1, modmanysp1.2), 
          file = "./Output/Species_differences/Height/LRT_POP.csv") # effect of population


# tukey graphs
write.csv(x = emmeans(object = modmanysp1, specs = list(pairwise ~ Expert_ID_living))[2],
          file = "./Output/Species_differences/SD_height.csv")
# want all six traits in plot side by side so save this one for later...
g_lmerfixplot(model = modmanysp1, factor = list(pairwise ~ Expert_ID_living), reorder = T)

# MCMCglmm
# note parameter expanded prior
# variance structure is the population that each Euphrasia came from
prior.sp.1 <-list(R=list(V=diag(1), nu=0.002), 
                  G=list(G1=list(V=diag(1), nu=1, alpha.mu=rep(0,1), alpha.V=diag(1)*1000)))

modmanysp1.1<- MCMCglmm(log(Height) ~ Expert_ID_living, 
                      random = ~E4E:Expert_ID_living,
                      data=manysp2, prior = prior.sp.1, 
                      nitt = 13000*7, 
                      burnin = 3000*7, 
                      thin = 10*7)
# summary output
summary(modmanysp1.1)
# check traces
traces(modmanysp1.1)
# variance explained

aod::wald.test(cov(modmanysp1.1$Sol[,2:5, drop = FALSE]),
               colMeans(modmanysp1.1$Sol[,2:5, drop = FALSE]),
               Terms = 1:4)$result

write.csv(x = data.table(UNITS = MCMCRepnorm(modmanysp1.1, "units"),
                         POPULATION = MCMCRepnorm(modmanysp1.1, y = "E4E:Expert_ID_living")), 
          file = "./Output/Species_differences/Height/Variance_Explained_Pop_Height.csv")


# 2. Nodes to flower
# full model 
modmanysp2 <- glmer(Nodes_to_flower_inc_coty ~ Expert_ID_living + (1 | E4E:Expert_ID_living),
                    data = manysp2, family = "poisson")
# omitting fixed effects
modmanysp2.1 <- glmer(Nodes_to_flower_inc_coty ~ 1 + (1 | E4E:Expert_ID_living),
                      data = manysp2, family = "poisson")
# omitting random effects
modmanysp2.2 <- glm(Nodes_to_flower_inc_coty ~ Expert_ID_living,
                    data = manysp2, family = "poisson")
# LRT's
anova(modmanysp2.1, modmanysp2)
anova(modmanysp2, modmanysp2.2)

# LRT's
write.csv(x = anova(modmanysp2.1, modmanysp2), 
          file = "./Output/Species_differences/Nodes/LRT_SP.csv") # effect of species
write.csv(x = anova(modmanysp2, modmanysp2.2), 
          file = "./Output/Species_differences/Nodes/LRT_POP.csv") # effect of population


write.csv(x = emmeans(object = modmanysp2, specs = list(pairwise ~ Expert_ID_living))[2],
          file = "./Output/Species_differences/SD_nodes.csv")
# again save this plot
g_lmerfixplot(model = modmanysp2, factor = list(pairwise ~ Expert_ID_living), reorder = T)



# MCMCglmm

prior.sp2 <-list(R=list(V=diag(1), nu=0.002), 
                 G=list(G1=list(V=diag(1), nu=0.002)))

modmanysp2<- MCMCglmm(Nodes_to_flower_inc_coty ~ Expert_ID_living, 
                      random = ~E4E:Expert_ID_living,
                      data=manysp2, prior = prior.sp2, 
                      family = "poisson",
                      nitt = 13000*8, burnin = 3000*8, thin = 10*8)
# summary output
summary(modmanysp2)
# check traces
traces(modmanysp2)
# variance explained
MCMCReppois(modmanysp2, y = "E4E:Expert_ID_living")

write.csv(x = data.table(UNITS = MCMCReppois(modmanysp2, "units"),
                         POPULATION = MCMCReppois(modmanysp2, y = "E4E:Expert_ID_living")), 
          file = "./Output/Species_differences/Nodes/Variance_Explained_Pop_Nodes.csv")


# 3. Corolla length

# full model
modmanysp3 <- lmer(Standard_corolla_l ~ Expert_ID_living +
                     (1 | E4E:Expert_ID_living), data = manysp2, na.action = na.omit)
# omitting fixed effects
modmanysp3.1 <- lmer(Standard_corolla_l ~ 1 +
                       (1 | E4E:Expert_ID_living), data = manysp2, na.action = na.omit)
# omitting random effects
modmanysp3.2 <- lm(Standard_corolla_l ~ Expert_ID_living, data = manysp2, na.action = na.omit)

# LRT's
anova(modmanysp3.1, modmanysp3)
anova(modmanysp3, modmanysp3.2)

write.csv(x = anova(modmanysp3.1, modmanysp3), 
          file = "./Output/Species_differences/Corolla/LRT_SP.csv") # effect of species
write.csv(x = anova(modmanysp3, modmanysp3.2), 
          file = "./Output/Species_differences/Corolla/LRT_POP.csv") # effect of population


write.csv(x = emmeans(object = modmanysp3, specs = list(pairwise ~ Expert_ID_living))[2],
          file = "./Output/Species_differences/SD_corolla.csv")
g_lmerfixplot(model = modmanysp3, factor = list(pairwise ~ Expert_ID_living), reorder = T)

# MCMCglmm

prior.sp3 <-list(R=list(V=diag(1), nu=0.002), 
                 G=list(G1=list(V=diag(1), nu=0.002)))
modmanysp3<- MCMCglmm(Standard_corolla_l ~ Expert_ID_living, 
                      random = ~E4E:Expert_ID_living,
                      data=manysp2, prior = prior.sp3,
                      nitt = 13000*7, burnin = 3000*7, thin = 10*7, pr = TRUE)
# summary output
summary(modmanysp3)
# check traces
traces(modmanysp3)
# variance explained
MCMCRepnorm(modmanysp3, "E4E:Expert_ID_living")

write.csv(x = data.table(UNITS = MCMCRepnorm(modmanysp3, "units"),
                         POPULATION = MCMCRepnorm(modmanysp3, y = "E4E:Expert_ID_living")), 
          file = "./Output/Species_differences/Corolla/Variance_Explained_Pop_Corolla.csv")

MCMCranef(modmanysp3)
# we can see that Euphrasia arctica shows MASSIVE variation
apply(modmanysp3$Sol, 2, posterior.mode) %>% data.frame() %>% data.table(keep.rownames = T) %>% arrange(desc(.))

# 4. Cauline internode ratio

# full model 
modmanysp4 <- lmer(log(Cauline_ratio) ~ Expert_ID_living +
                     (1 | E4E:Expert_ID_living), data = manysp2)
# omitting fixed effects
modmanysp4.1 <- lmer(log(Cauline_ratio) ~ 1 +
                       (1 | E4E:Expert_ID_living), data = manysp2)
# omitting random effects
modmanysp4.2 <- lm(log(Cauline_ratio) ~ Expert_ID_living, data = manysp2)

# LRT's
anova(modmanysp4, modmanysp4.1)
anova(modmanysp4, modmanysp4.2)

write.csv(x = anova(modmanysp4.1, modmanysp4), 
          file = "./Output/Species_differences/Cauline_ratio/LRT_SP.csv") # effect of species
write.csv(x = anova(modmanysp4, modmanysp4.2), 
          file = "./Output/Species_differences/Cauline_ratio/LRT_POP.csv") # effect of population


write.csv(x = emmeans(object = modmanysp4, specs = list(pairwise ~ Expert_ID_living))[2],
          file = "./Output/Species_differences/Cauline_ratio/SD_cauline_ratio.csv")
g_lmerfixplot(model = modmanysp4, factor = list(pairwise ~ Expert_ID_living), reorder = T)

# MCMCglmm

prior.sp4 <-list(R=list(V=diag(1), nu=0.002), 
                 G=list(G1=list(V=diag(1), nu=0.002)))
modmanysp4 <- MCMCglmm(Cauline_ratio ~ Expert_ID_living, 
                       random = ~E4E:Expert_ID_living,
                       data=manysp2, prior = prior.sp4,
                       nitt = 13000*7, burnin = 3000*7, thin = 10*7, pr =T)
# check summary 
summary(modmanysp4)
# check traces
traces(modmanysp4)
# variance explained
MCMCRepnorm(mod = modmanysp4, y = "E4E:Expert_ID_living")

write.csv(x = data.table(UNITS = MCMCRepnorm(modmanysp4, "units"),
                         POPULATION = MCMCRepnorm(modmanysp4, y = "E4E:Expert_ID_living")), 
          file = "./Output/Species_differences/Cauline_ratio/Variance_Explained_Pop_Cauline_ratio.csv")

MCMCranef(mod = modmanysp4)
# 5. Days to flower

# full model 
modmanysp5 <- glmer(Julian_days_to_flower ~ Expert_ID_living + 
                      (1 | E4E:Expert_ID_living), data = manysp2[complete.cases(manysp2),], family = "poisson")
# omitting fixed effects
modmanysp5.1 <- glmer(Julian_days_to_flower ~ 1 + 
                        (1 | E4E:Expert_ID_living), data = manysp2[complete.cases(manysp2),], family = "poisson")
# omitting random effects
modmanysp5.2 <- glm(Julian_days_to_flower ~ Expert_ID_living, data = manysp2[complete.cases(manysp2),], family = "poisson")

# LRT's
anova(modmanysp5, modmanysp5.1)
anova(modmanysp5, modmanysp5.2)

write.csv(x = anova(modmanysp5.1, modmanysp5), 
          file = "./Output/Species_differences/Julian_days/LRT_SP.csv") # effect of species
write.csv(x = anova(modmanysp5, modmanysp5.2), 
          file = "./Output/Species_differences/Julian_days/LRT_POP.csv") # effect of population


write.csv(x = emmeans(object = modmanysp5, specs = list(pairwise ~ Expert_ID_living))[2],
          file = "./Output/Species_differences/Julian_days/SD_julian.csv")
g_lmerfixplot(model = modmanysp5, factor = list(pairwise ~ Expert_ID_living), reorder = T)


# MCMCglmm

prior.sp5 <-list(R=list(V=diag(1), nu=0.002), 
                 G=list(G1=list(V=diag(1), nu=0.02)))
modmanysp5<- MCMCglmm(Julian_days_to_flower ~ Expert_ID_living, 
                      random = ~E4E:Expert_ID_living,
                      data=manysp2, prior = prior.sp5,
                      family = "poisson",
                      nitt = 13000*7, burnin = 3000*7, thin = 10*7, pr = TRUE)
# summary output
summary(modmanysp5)
# check traces
traces(modmanysp5)
# variance explained
MCMCReppois(modmanysp5, "E4E:Expert_ID_living")

write.csv(x = data.table(UNITS = MCMCReppois(modmanysp5, "units"),
                         POPULATION = MCMCReppois(modmanysp5, y = "E4E:Expert_ID_living")), 
          file = "./Output/Species_differences/Julian_days/Variance_Explained_Pop_Cauline_ratio.csv")
MCMCranef(modmanysp5)

# 6. Number of branches

# MCMCglmm only

prior.sp6 <-list(R=list(V=diag(1), nu=0.002), 
                 G=list(G1=list(V=diag(1), nu=0.002)))
modmanysp6<- MCMCglmm(No_branches ~ Expert_ID_living, 
                      random = ~E4E:Expert_ID_living,
                      data=manysp2, prior = prior.sp6,
                      family = "poisson",
                      nitt = 13000*7, burnin = 3000*7, thin = 10*7)
# check summary
summary(modmanysp6)
# check traces
traces(modmanysp6)
# variance explained
MCMCReppois(modmanysp6, "E4E:Expert_ID_living")

# omitting fixed effects
modmanysp6.1<- MCMCglmm(No_branches ~ 1, 
                        random = ~E4E:Expert_ID_living,
                        data=manysp2, prior = prior.sp6,
                        family = "poisson",
                        nitt = 13000*7, burnin = 3000*7, thin = 10*7)
# no random effect variance prior
prior.sp6.1 <-list(R=list(V=diag(1), nu=0.002))
# omitting random effects
modmanysp6.2<- MCMCglmm(No_branches ~ Expert_ID_living, 
                        data=manysp2, prior = prior.sp6.1,
                        family = "poisson",
                        nitt = 13000*7, burnin = 3000*7, thin = 10*7)
# comparing deviance information criteria
modmanysp6$DIC
modmanysp6.1$DIC
modmanysp6.2$DIC


# 7. Number of teeth on the lower floral leaf
manysp2$Obs <- as.factor(1:nrow(manysp2))
# full model
modmanysp7 <- glmer(No_teeth_lower_floral_leaf ~ Expert_ID_living + 
                      (1 | E4E:Expert_ID_living), data = manysp2, family = "poisson")

g_lmerfixplot(modmanysp7, factor = list(pairwise ~ Expert_ID_living), reorder = T)
# omitting fixed effects
modmanysp7.1 <- glmer(No_teeth_lower_floral_leaf ~ 1 + 
                        (1 | E4E:Expert_ID_living) + (1 | Obs), data = manysp2, family = "poisson")
# omitting random effects
modmanysp7.2 <- glm(No_teeth_lower_floral_leaf ~ Expert_ID_living, data = manysp2, family = "poisson")

# LRT's
anova(modmanysp7, modmanysp7.1)
anova(modmanysp7, modmanysp7.2)

write.csv(x = anova(modmanysp7.1, modmanysp7), 
          file = "./Output/Species_differences/Teeth/LRT_SP.csv") # effect of species
write.csv(x = anova(modmanysp7, modmanysp7.2), 
          file = "./Output/Species_differences/Teeth/LRT_POP.csv") # effect of population


write.csv(x = emmeans(object = modmanysp7, specs = list(pairwise ~ Expert_ID_living))[2],
          file = "./Output/Species_differences/SD_teeth.csv")
g_lmerfixplot(model = modmanysp7, factor = list(pairwise ~ Expert_ID_living), reorder = T)

g_lmfixplot(modmanysp7.2, mcp(Expert_ID_living = "Tukey"))
# MCMCglmm

prior.sp7 <-list(R=list(V=diag(1), nu=0.002), 
                 G=list(G1=list(V=diag(1), nu=1, alpha.mu=rep(0,1), alpha.V=diag(1)*1000)))
modmanysp7 <- MCMCglmm(No_teeth_lower_floral_leaf ~ Expert_ID_living, 
                       random = ~E4E:Expert_ID_living,
                       data=manysp2, prior = prior.sp7,
                       family = "poisson",
                       nitt = 13000*7, burnin = 3000*7, thin = 10*7)
# summary output
summary(modmanysp7)
# check traces
traces(modmanysp7)
# variance explained
MCMCReppois(mod = modmanysp7, y = "E4E:Expert_ID_living")

write.csv(x = data.table(UNITS = MCMCReppois(modmanysp7, "units"),
                         POPULATION = MCMCReppois(modmanysp7, y = "E4E:Expert_ID_living")), 
          file = "./Output/Species_differences/Teeth/Variance_Explained_Pop_Teeth.csv")


# cowplot all of these tukey graphs together...

b1 <- g_lmerfixplot(model = modmanysp1, factor = list(pairwise ~ Expert_ID_living)) + ggtitle("Height differences") + theme(axis.title.x = element_blank(),
                                                                                                                            plot.title = element_text(size = 10))
b2 <- g_lmerfixplot(model = modmanysp2, factor = list(pairwise ~ Expert_ID_living)) + ggtitle("Node to flower differences") + theme(axis.text.y = element_blank(),
                                                                                                                   axis.ticks.y = element_blank(),
                                                                                                                   axis.title.y = element_blank(),
                                                                                                                   axis.title.x = element_blank(),
                                                                                                                   plot.title = element_text(size = 10))
b3 <- g_lmerfixplot(model = modmanysp3, factor = list(pairwise ~ Expert_ID_living)) + ggtitle("Corolla length differences")+ theme(axis.text.y = element_blank(),
                                                                                                                             axis.ticks.y = element_blank(),
                                                                                                                             axis.title.y = element_blank(),
                                                                                                                             axis.title.x = element_blank(),
                                                                                                                             plot.title = element_text(size = 10))
b4 <- g_lmerfixplot(model = modmanysp4, factor = list(pairwise ~ Expert_ID_living)) + ggtitle("Cauline:internode ratio differences")+ theme(axis.text.y = element_blank(),
                                                                                                                                      axis.ticks.y = element_blank(),
                                                                                                                                      axis.title.y = element_blank(),
                                                                                                                                      axis.title.x = element_blank(),
                                                                                                                                      plot.title = element_text(size = 10))

rsum <-  emmeans(modmanysp4, list(pairwise ~ Expert_ID_living), adjust = "tukey")[[2]]
rsum <- as.data.table(rsum)

rsum2 <- emmeans(modmanysp5, list(pairwise ~ Expert_ID_living), adjust = "tukey")[[2]]
rsum2 <- as.data.table(rsum2)

rsum <- rsum2[rsum, on = "contrast"]


rsum$ptext <- rep(NA, dim(rsum)[1])
rsum$enhance <- !is.na(as.character(rsum$ptext))

b5 <- ggplot(rsum,aes(x = contrast , y = estimate, alpha = enhance))+
  geom_pointrange(aes(ymin=estimate-2*SE,
                      ymax=estimate+2*SE))+
  geom_text(aes(label = ptext), nudge_x = 0.2) +
  coord_flip()+
  theme_bw()+
  labs(y = "Point estimates and confidence intervals")+
  theme(legend.position = "none")+ ggtitle("Julian days to flower differences")+ theme(axis.text.y = element_blank(),
                                                                                         axis.ticks.y = element_blank(),
                                                                                       axis.title.y = element_blank(),
                                                                                       axis.title.x = element_blank(),
                                                                                       plot.title = element_text(size = 10))

b6 <- g_lmerfixplot(model = modmanysp7, factor = list(pairwise ~ Expert_ID_living)) + ggtitle("Number of teeth on lower floral leaf differences") + theme(axis.text.y = element_blank(),
                                                                                                                                                   axis.ticks.y = element_blank(),
                                                                                                                                                   axis.title.y = element_blank(),
                                                                                                                                                   axis.title.x = element_blank(),
                                                                                                                                                   plot.title = element_text(size = 10))
cowplot::plot_grid(b1, b2, b3, b4, b5, b6, rel_widths = c(1, rep(0.5, 5)), nrow = 1, label_x = "Point estimates")



# Principal Component Analysis many species

# subset data
PCAone<-manysp[,-c(1,2,5, 8, 11, 13, 14, 15, 17, 18)]
# remove NA's
PCAone2 <- PCAone[complete.cases(PCAone),]

# compute PCA
one.pca<- prcomp(PCAone2[,-c(1,2)], 
                 center = TRUE,
                 scale. = TRUE)
# percentage variance of each variable
manyspaload <- abs(one.pca$rotation)
PCAmanysp<-sweep(manyspaload, 2, colSums(aload), "/")

# summary and cumulative variance of components
summarypcaone  <- summary(one.pca)
summarypcaone <- summarypcaone$importance
data.frame(summarypcaone)
# tidy data frame
PCAmanysp2<-rbind(PCAmanysp, summarypcaone)

write.csv(x = specify_decimal(PCAmanysp2, 3), file = "./Output/PCA/Many_sp_incl_hybrids.csv")

one.pca.1 <- as.data.frame(one.pca$x)
one.pca.1$Pop <- PCAone2$E4E
one.pca.1$Species <- PCAone2$Expert_ID_living

PCA1<- ggplot(one.pca.1, aes(x = PC1, y= PC2,
                      group = Species,
                      colour = Species))+
  geom_point(aes(group = Species, colour =Species))+
  stat_ellipse(aes(group = Species, colour = Species))+
  ylab(label = "PC2 (20.1%)")+
  xlab(label = "PC1 (50.3%)")+
  theme_bw()

# minus hybrids

# compute PCA
PCAone3 <- PCAone2[!grepl("x", Expert_ID_living),]
one.pca2<- prcomp(PCAone3[,-c(1,2)], 
                 center = TRUE,
                 scale. = TRUE)
# percentage variance of each variable
manyspaload2 <- abs(one.pca2$rotation)
PCAmanysphyb<-sweep(manyspaload2, 2, colSums(aload), "/")

# summary and cumulative variance of components
summarypcatwo  <- summary(one.pca2)
summarypcatwo <- summarypcatwo$importance
data.frame(summarypcatwo)
# tidy data frame
PCAmanysphyb2<-rbind(PCAmanysphyb, summarypcatwo)

write.csv(specify_decimal(PCAmanysphyb2, 3),file = "./Output/PCA/Many_sp_no_hybrids_PCA.csv")

one.pca.2 <- as.data.frame(one.pca2$x)
one.pca.2$Pop <- PCAone3$E4E
one.pca.2$Species <- PCAone3$Expert_ID_living

PCA2 <- ggplot(one.pca.2, aes(x = PC1, y= PC2,
                      group = Species,
                      colour = Species))+
  geom_point(aes(group = Species, colour =Species))+
  stat_ellipse(aes(group = Species, colour = Species))+
  ylab(label = "PC2 (20.6%)")+
  xlab(label = "PC1 (52.8%)")+
  theme_bw()



##### Early growth vs late growth #####
# MCMCglmm analysis for 
# load the data
earlylate <- read.csv("./Data/Earlylate.csv")
setDT(earlylate)

write.csv(earlylate[, .(Early_height = specify_decimal(mean(Early.season.growth, na.rm = TRUE), 1),
              SE = specify_decimal(sd(Early.season.growth, na.rm = TRUE)/sqrt(.N), 1),
              End_height = specify_decimal(mean(Height_end_season, na.rm = TRUE), 1),
              SE = specify_decimal(sd(Height_end_season, na.rm = TRUE)/sqrt(.N), 1)), by = "Host"], 
          file = "./Output/Early_late_means_sd/Early_late_means_sd.csv")

# 1. Height end of season as a function of height at first flowering
prior.height <- list(R=list(V=diag(1), nu=0.002), 
                     G=list(G1=list(V=diag(1), nu=0.002)))
mcmcheight <- MCMCglmm(Height_end_season ~ Height, random = ~Host, data=earlylate,
                       prior = prior.height,
                       nitt = 13000*7,
                       burnin=3000*7,
                       thin=10*7)
# summary output
summary(mcmcheight)
# check traces
traces(mcmcheight)

# correlation between height at end of season and height at first flowering
rcorr(earlylate$Height, earlylate$Height_end_season)$r[2]

# 2. Height at end of season as a function of early season growth

prior.height2 <- list(R=list(V=diag(1), nu=0.002), 
                      G=list(G1=list(V=diag(1), nu=0.002)))
mcmcheight2 <- MCMCglmm(Height_end_season ~ Early.season.growth, random = ~Host, 
                        data=earlylate[complete.cases(earlylate$Early.season.growth),],
                        prior = prior.height2,
                        nitt = 13000*7,
                        burnin=3000*7,
                        thin=10*7)
# summary output
summary(mcmcheight2)
# check traces
traces(mcmcheight2)

# correlation between height at end of season and height at first flowering
rcorr(earlylate$Early.season.growth, earlylate$Height_end_season)$r[2]


# 3. Days to flower as a function of height at end of season

prior.height3 <- list(R=list(V=diag(1), nu=0.002), 
                      G=list(G1=list(V=diag(1), nu=0.002)))
mcmcheight3 <- MCMCglmm(Julian.days.to.flower ~ Height_end_season, random = ~Host, 
                        data=earlylate[complete.cases(earlylate$Height_end_season),],
                        prior = prior.height3,
                        nitt = 13000*7,
                        burnin=3000*7,
                        thin=10*7)
# summary output
summary(mcmcheight3)
# check traces
traces(mcmcheight3)

# correlation between height at end of season and height at first flowering
rcorr(earlylate$Julian.days.to.flower, earlylate$Height_end_season)$r[2]


# 4. Number of branches as a function of height at end of season
prior.height4 <- list(R=list(V=diag(1), nu=0.002), 
                      G=list(G1=list(V=diag(1), nu=0.002)))
mcmcheight4 <- MCMCglmm(No..branches ~ Height_end_season, random = ~Host, 
                        data=earlylate[complete.cases(earlylate$Height_end_season),],
                        prior = prior.height4,
                        family = "poisson",
                        nitt = 13000*7,
                        burnin=3000*7,
                        thin=10*7)
# summary output
summary(mcmcheight4)
# check traces
traces(mcmcheight4)

# correlation between height at end of season and height at first flowering
rcorr(earlylate$No..branches, earlylate$Height_end_season)$r[2]


##### Common garden vs Wild #####

# read in data
pope4ecom <- read.csv("./Data/Wildcommon.csv")
tidy.pop <- gather(pope4ecom, "Variable", "Value", 4:11) ## merge into a long format data set
tidy.pop.n <- setDT(tidy.pop)[Variable == "Nodes.to.flower.inc.coty"]

# 1. Nodes to flower

tidy.pop.n <- filter(tidy.pop, tidy.pop$Variable == "Nodes.to.flower.inc.coty")
ggplot(tidy.pop.n, aes(x = Value))+
  geom_histogram(aes(fill = Experiment), alpha=0.3, position = position_dodge())

prior.mcmcwild1<-list(R=list(V=diag(1), nu=0.002), 
                      G=list(G1=list(V=diag(1), nu=1, alpha.mu=rep(0,1), alpha.V=diag(1)*1000),
                             #G2=list(V=diag(1), nu=1, alpha.mu=rep(0,1), alpha.V=diag(1)*1000),
                             G3=list(V=diag(1), nu=1, alpha.mu=rep(0,1), alpha.V=diag(1)*1000)))

tidy.pop.n$Exp.E4E<-paste(tidy.pop.n$Experiment, tidy.pop.n$E4E, sep = "-") 
# no. of pops
length(unique(tidy.pop.n$Exp.E4E))

mcmc.cor.n<- MCMCglmm(Value ~ Experiment, random=~Species + E4E, 
                      family = "poisson",
                      data=tidy.pop.n,
                      nitt = 13000*10,
                      burnin = 3000*10,
                      thin = 10*10,
                      prior = prior.mcmcwild1)
# check traces
traces(mcmc.cor.n)
# summary output
write.csv(x = specify_decimal(summary(mcmc.cor.n)$solutions, 3), 
          file = "./Output/Common_wild/Nodes/CW_nodes.csv")
# variance explained
write.csv(x = data.table(UNITS = MCMCReppois(mcmc.cor.n, "units"),
                         SPECIES = MCMCReppois(mcmc.cor.n, y = "Species"),
                         #EXP.E4E = MCMCReppois(mcmc.cor.n, y = "Exp.E4E"),
                         E4E = MCMCReppois(mcmc.cor.n, y = "E4E")), 
          file = "./Output/Common_wild/Nodes/CW_nodes_variances.csv")

# 2. Standard corolla length

tidy.pop.c <- filter(tidy.pop, tidy.pop$Variable == "Standard.corolla.l")
# no. of pops
length(unique(tidy.pop.c$E4E))

tidy.pop.c$Exp.E4E<-paste(tidy.pop.c$Experiment, tidy.pop.c$E4E, sep = "-") 


prior.mcmcwild2<-list(R=list(V=diag(1), nu=0.002), 
                      G=list(G1=list(V=diag(1), nu=1, alpha.mu=rep(0,1), alpha.V=diag(1)*1000),
                             G1=list(V=diag(1), nu=1, alpha.mu=rep(0,1), alpha.V=diag(1)*1000)))
mcmc.cor.c<- MCMCglmm(Value ~ Experiment, random=~E4E + Species, 
                      data=tidy.pop.c,
                      nitt = 13000*10,
                      burnin = 3000*10,
                      thin = 10*10,
                      prior = prior.mcmcwild2)
# check traces
traces(mcmc.cor.c)
# summary output
summary(mcmc.cor.c)
# variance explained
MCMCRepnorm(mcmc.cor.c, "Species")

write.csv(x = specify_decimal(summary(mcmc.cor.c)$solutions, 3), 
          file = "./Output/Common_wild/Corolla/CW_corolla.csv")
# variance explained
write.csv(x = data.table(UNITS = MCMCRepnorm(mcmc.cor.c, "units"),
                         SPECIES = MCMCRepnorm(mcmc.cor.c, y = "Species"),
                         E4E = MCMCRepnorm(mcmc.cor.c, y = "E4E")), 
          file = "./Output/Common_wild/Corolla/CW_corolla_variances.csv")

# 3. Cauline internode ratio

tidy.pop.lr <- filter(tidy.pop, tidy.pop$Variable == "Cauline.internode.leaf.ratio")
# no. of pops
length(unique(tidy.pop.lr$E4E))

tidy.pop.lr$Exp.E4E<-paste(tidy.pop.lr$Experiment, tidy.pop.lr$E4E, sep = "-")

prior.mcmcwild3<-list(R=list(V=diag(1), nu=0.002), 
                      G=list(G1=list(V=diag(1), nu=1, alpha.mu=rep(0,1), alpha.V=diag(1)*1000),
                             G1=list(V=diag(1), nu=1, alpha.mu=rep(0,1), alpha.V=diag(1)*1000)))
mcmc.cor.lr<- MCMCglmm(Value ~ Experiment, random=~E4E+Species, 
                       data=tidy.pop.lr,
                       nitt = 13000*8,
                       burnin = 3000*8,
                       thin = 10*8,
                       prior = prior.mcmcwild3)
plot(mcmc.cor.lr)
# check traces
traces(mcmc.cor.lr)
# cummary output
summary(mcmc.cor.lr)
# variance explained
MCMCRepnorm(mcmc.cor.lr, "Species")

write.csv(x = specify_decimal(summary(mcmc.cor.lr)$solutions, 3), 
          file = "./Output/Common_wild/Cauline_ratio/CW_cauline_ratio.csv")
# variance explained
write.csv(x = data.table(UNITS = MCMCRepnorm(mcmc.cor.lr, "units"),
                         SPECIES = MCMCRepnorm(mcmc.cor.lr, y = "Species"),
                         E4E = MCMCRepnorm(mcmc.cor.lr, y = "E4E")), 
          file = "./Output/Common_wild/Cauline_ratio/CW_cauline_ratio_variances.csv")

# 4. Number of teeth on the lower floral leaf

tidy.pop.t <- filter(tidy.pop, tidy.pop$Variable == "No..teeth.lower.floral.leaf")

tidy.pop.t$Exp.E4E<-paste(tidy.pop.t$Experiment, tidy.pop.t$E4E, sep = "-")


prior.mcmcwild4<-list(R=list(V=diag(1), nu=0.002), 
                      G=list(G1=list(V=diag(1), nu=1, alpha.mu=rep(0,1), alpha.V=diag(1)*1000),
                             G1=list(V=diag(1), nu=1, alpha.mu=rep(0,1), alpha.V=diag(1)*1000)))
mcmc.cor.t<- MCMCglmm(Value ~ Experiment, random=~E4E+Species, 
                      data=tidy.pop.t,
                      family = "poisson",
                      nitt = 13000*8,
                      burnin = 3000*8,
                      thin = 10*8,
                      prior = prior.mcmcwild4)

# check traces
traces(mcmc.cor.t)
# summary output
summary(mcmc.cor.t)
# variance explained
MCMCReppois(mcmc.cor.t, "Species")

write.csv(x = specify_decimal(summary(mcmc.cor.t)$solutions, 3), 
          file = "./Output/Common_wild/Teeth/CW_teeth.csv")
# variance explained
write.csv(x = data.table(UNITS = MCMCReppois(mcmc.cor.t, "units"),
                         SPECIES = MCMCReppois(mcmc.cor.t, y = "Species"),
                         E4E = MCMCReppois(mcmc.cor.t, y = "E4E")), 
          file = "./Output/Common_wild/Teeth/CW_teeth_variances.csv")



# 5. Number of branches

tidy.pop.b <- filter(tidy.pop, tidy.pop$Variable == "No.of.branches")

tidy.pop.b$Exp.E4E<-paste(tidy.pop.b$Experiment, tidy.pop.b$E4E, sep = "-")

prior.mcmcwild5<-list(R=list(V=diag(1), nu=0.002), 
                      G=list(G1=list(V=diag(1), nu=1, alpha.mu=rep(0,1), alpha.V=diag(1)*1000),
                             G1=list(V=diag(1), nu=1, alpha.mu=rep(0,1), alpha.V=diag(1)*1000)))

mcmc.cor.b<- MCMCglmm(Value ~ Experiment, random=~E4E + Species, 
                      data=tidy.pop.b,
                      family = "poisson",
                      nitt = 13000*10,
                      burnin = 3000*10,
                      thin = 10*10,
                      prior = prior.mcmcwild5)
# check traces
traces(mcmc.cor.b)
# summary output
summary(mcmc.cor.b)
# variance explained
MCMCReppois(mcmc.cor.b, "Species")

write.csv(x = specify_decimal(summary(mcmc.cor.b)$solutions, 3), 
          file = "./Output/Common_wild/Branches/CW_branches.csv")
# variance explained
write.csv(x = data.table(UNITS = MCMCReppois(mcmc.cor.b, "units"),
                         SPECIES = MCMCReppois(mcmc.cor.b, y = "Species"),
                         E4E = MCMCReppois(mcmc.cor.b, y = "E4E")), 
          file = "./Output/Common_wild/Branches/CW_branches_variances.csv")

setDT(tidy.pop.b)
setDT(tidy.pop.t)
setDT(tidy.pop.lr)
setDT(tidy.pop.c)
setDT(tidy.pop.n)

a1<-ggplot(tidy.pop.n[, .(mean = mean(Value),
                      sem = sd(Value)/sqrt(.N)), by = c("E4E", "Species", "Experiment")], aes(x = Experiment, y = mean, ymin = mean-sem, ymax=mean+sem, group = E4E))+
  geom_errorbar(position = position_dodge(0.6))+
  geom_line(aes(group=E4E), position = position_dodge(0.6), lty =3, alpha = 0.7)+
  geom_point(aes(colour = Species), position = position_dodge(0.6), size=4)+
  theme_bw()+
  theme(legend.position = "none")+
  ylab(label = "Mean Node to Flower")

a2<-ggplot(tidy.pop.c[, .(mean = mean(Value),
                      sem = sd(Value)/sqrt(.N)), by = c("E4E", "Species", "Experiment")], aes(x = Experiment, y = mean, ymin = mean-sem, ymax=mean+sem, group = E4E))+
  geom_errorbar(position = position_dodge(0.6))+
  geom_line(aes(group=E4E), position = position_dodge(0.6), lty =3, alpha = 0.7)+
  geom_point(aes(colour = Species), position = position_dodge(0.6), size=4)+
  theme_bw()+
  theme(legend.position = "none")+
  ylab(label = "Mean Corolla Length (mm)")

a3<-ggplot(tidy.pop.lr[, .(mean = mean(Value),
                       sem = sd(Value)/sqrt(.N)), by = c("E4E", "Species", "Experiment")][!is.na(mean),], aes(x = Experiment, y = mean, ymin = mean-sem, ymax=mean+sem, group = E4E))+
  geom_errorbar(position = position_dodge(0.6))+
  geom_line(aes(group=E4E), position = position_dodge(0.6), lty =3, alpha = 0.7)+
  geom_point(aes(colour = Species), position = position_dodge(0.6), size=4)+
  theme_bw()+
  theme(legend.position = "none")+
  ylab(label = "Mean Cauline:Internode Ratio")

a4<-ggplot(tidy.pop.t[, .(mean = mean(Value),
                      sem = sd(Value)/sqrt(.N)), by = c("E4E", "Species", "Experiment")], aes(x = Experiment, y = mean, ymin = mean-sem, ymax=mean+sem, group = E4E))+
  geom_errorbar(position = position_dodge(0.6))+
  geom_line(aes(group=E4E), position = position_dodge(0.6), lty =3, alpha = 0.7)+
  geom_point(aes(colour = Species), position = position_dodge(0.6), size=4)+
  theme_bw()+
  theme(legend.position = "none")+
  ylab(label = "Mean Number of Leaf Teeth")

a5<-ggplot(tidy.pop.t[, .(mean = mean(Value),
                      sem = sd(Value)/sqrt(.N)), by = c("E4E", "Species", "Experiment")], aes(x = Experiment, y = mean, ymin = mean-sem, ymax=mean+sem, group = E4E))+
  geom_errorbar(position = position_dodge(0.6))+
  geom_line(aes(group=E4E), position = position_dodge(0.6), lty =3, alpha = 0.7)+
  geom_point(aes(colour = Species), position = position_dodge(0.6), size=4)+
  theme_bw()+
  ylab(label = "Mean Number of teeth")

a5.1 <- get_legend(a5)

a6<-ggplot(tidy.pop.b[, .(mean = mean(Value),
                          sem = sd(Value)/sqrt(.N)), by = c("E4E", "Species", "Experiment")], aes(x = Experiment, y = mean, ymin = mean-sem, ymax=mean+sem, group = E4E))+
  geom_errorbar(position = position_dodge(0.6))+
  geom_line(aes(group=E4E), position = position_dodge(0.6), lty =3, alpha = 0.7)+
  geom_point(aes(colour = Species), position = position_dodge(0.6), size=4)+
  theme_bw()+
  theme(legend.position = "none")+
  ylab(label = "Mean Number of Branches")

# ditch multiple populations and branches!!
plot_grid(... = a1,a2,a3,a4,a5.1, nrow = 1, labels = "AUTO", rel_widths = c(3,3,3,3,1.5))

##### Additional analysis; removing Pinus and Marchantia #####

`%!in%` <- function(x,y)!('%in%'(x,y))

NBerwickme_add <- setDT(NBerwickme)[Host %!in% c("Marchantia polymorpha", "Pinus sylvestris")]
# relevel so no host is baseline
NBerwickme_add$Host <- relevel(x = NBerwickme_add$Host, ref = "No host")



# models are specified below, height is log(height)
height1.1 <- lm(log(Height) ~ Host, data=NBerwickme_add)
# write to csv the coefficients of the model and the anova
write.csv(x = tidy(height1.1), 
          file = "./Output/Phenotypic_plasticity/PP_height_minusPM.csv")
write.csv(anova(height1.1), file =  "./Output/Phenotypic_plasticity/PP_height_anova_minusPM.csv") #effect of host is highly significant
# tukey tests here
tukeyheight1.1 <- tidy(TukeyHSD(aov(height1.1)))
write.csv(tukeyheight1.1, file =  "./Output/Phenotypic_plasticity/PP_height_tukey_minusPM.csv")
# and a plot for the tukey tests
g_lmfixplot(height1.1)

## 1.2 Nodes_to_flower_inc_coty

nodes1.1 <- glm(Nodes_to_flower_inc_coty ~ Host, 
              data=NBerwickme_add, family = "poisson")
write.csv(x = tidy(nodes1.1), 
          file = "./Output/Phenotypic_plasticity/PP_nodes_mins_PM.csv")
write.csv(anova(nodes1.1, test = "Chisq"), file =  "./Output/Phenotypic_plasticity/PP_nodes_anova_minusPM.csv") 

tukeynodes1.1 <- tidy(summary(glht(nodes1.1, mcp(Host="Tukey"))))
write.csv(tukeynodes1.1, file =  "./Output/Phenotypic_plasticity/PP_nodes_tukey_minus_PM.csv")

g_lmfixplot(nodes1.1, factor = mcp(Host = "Tukey"))

## 1.3 Standard_corolla_l

corol1.1 <- lm(Standard_corolla_l ~ Host, 
             data=NBerwickme_add[c(-20, -37),], na.action = na.omit)

summary(corol1.1)
write.csv(x = tidy(corol1.1), 
          file = "./Output/Phenotypic_plasticity/PP_corolla_minus_PM.csv")
write.csv(anova(corol1.1), file =  "./Output/Phenotypic_plasticity/PP_corolla_anova_minusPM.csv") #effect of host is highly significant

tukeycorol1.1 <- tidy(TukeyHSD(aov(corol1.1)))
write.csv(tukeycorol1.1, file =  "./Output/Phenotypic_plasticity/PP_corolla_tukey_minus_PM.csv")

g_lmfixplot(corol1.1)

## 1.4 Cauline_ratio

ratio1.1 <- lm(Cauline_ratio ~ Host, 
             data=NBerwickme_add[-54,], na.action = na.omit)
summary(ratio1.1)

write.csv(x = tidy(ratio1.1), 
          file = "./Output/Phenotypic_plasticity/PP_cauline_ratio_minus_PM.csv")
write.csv(anova(ratio1.1), file =  "./Output/Phenotypic_plasticity/PP_cauline_ratio_anova_minus_PM.csv") # effect of host is highly significant

tukeyratio1.1 <- tidy(TukeyHSD(aov(ratio1.1)))
write.csv(tukeyratio1.1, file =  "./Output/Phenotypic_plasticity/PP_cauline_ratio_tukey_minus_PM.csv")

g_lmfixplot(ratio1.1)

## 1.5 Julian_days_to_flower


julian1.1 <- glm(Julian_days_to_flower ~ Host, 
               data=NBerwickme_add, na.action = na.omit, family = "poisson")

summary(julian1.1)

write.csv(x = tidy(julian1.1), 
          file = "./Output/Phenotypic_plasticity/PP_julian_minusPM.csv")
write.csv(anova(julian1.1, test = "Chisq"), file =  "./Output/Phenotypic_plasticity/PP_julian_anova_minus_PM.csv") #effect of host is highly significant

tukeyjulian1.1 <-  tidy(summary(glht(julian1.1, mcp(Host="Tukey"))))
write.csv(tukeyjulian1.1, file =  "./Output/Phenotypic_plasticity/PP_julian_tukey_minus_PM.csv")

g_lmfixplot(julian1.1, tukey = FALSE, mcp(Host = "Tukey"))

## 1.6 No_branches. Not a good model (omit...)
## 1.7 No_teeth_lower_floral_leaf

teeth1.1 <- glm(No_teeth_lower_floral_leaf ~ Host, 
              data=NBerwickme_add, family = "poisson", na.action = na.omit)

write.csv(x = tidy(teeth1.1), 
          file = "./Output/Phenotypic_plasticity/PP_teeth_minusPM.csv")
write.csv(anova(teeth1.1, test = "Chisq"), file =  "./Output/Phenotypic_plasticity/PP_teeth_anova_minus_PM.csv") #effect of host is highly significant

tukeyteeth1.1 <- tidy(summary(glht(teeth1.1, mcp(Host="Tukey"))))
write.csv(tukeyteeth1.1, file =  "./Output/Phenotypic_plasticity/PP_teeth_tukey_minus_PM.csv")

g_lmfixplot(teeth1.1,  mcp(Host = "Tukey"), tukey = FALSE)

# cowplot the PCA's

cowplot::plot_grid(PCA1, PCA2, PCA3)