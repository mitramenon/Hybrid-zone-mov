############Morphoogical data for clines###############################
##Author: menonm2 (menonm2@mymail.vcu.edu)
##Date modified: 2-Oct-2018
##NOTE: First several steps are data prep. If your data is in the appropriate format (mean, var, sample size per population), skip to last step
####################################################################

library(MASS)
library(rcompanion)
library(geosphere)
library(hzar)

#load data 
mexican<-read.table("Data Mexico-USA P. strobiformis.txt",header=T,sep="\t")
PIFL<-read.table("PIFL_coneData.txt",header = T,sep="\t")
Chapter1<-read.table("Pops_55_chapter1.txt",header=T,sep = "\t")
Chapter1<-Chapter1[!(duplicated(Chapter1$RenamedPop)), ]

#1. check normality prior to cline fit
cone<-c(mexican[ ,5],PIFL[ ,8])
seed<-c(mexican[ ,7],PIFL[ ,9])

qqnorm(cone)
qqline(cone,col=2)
qqnorm(seed)
qqline(seed,col=2)

plotNormalHistogram(cone)
plotNormalHistogram(seed)
##mostly normal fit

######STEP 2. DATA PREP :
This step was needed to prepare the dataset for input into HZAR, here we estimate mean, var and N per population###
#################################################

####SWWP
colnames(mexican)[1]<-"RenamedPop"
commonUS<-merge(Chapter1,mexican,by="RenamedPop") #determine shared pops first to preserve these
commonUS_Locs<-commonUS[!(duplicated(commonUS$RenamedPop)), ]

#split by pops and estimate mean per pop
commonUS<-split(commonUS,f = commonUS$RenamedPop,drop=TRUE)
commonUS_mean<-lapply(commonUS, function(x) return(apply(x[ ,c(11,13)],2,mean,na.rm=TRUE)))
commonUS_var<-lapply(commonUS, function(x) return(apply(x[ ,c(11,13)],2,var,na.rm=TRUE)))

#create data for use, needs to include number of inds per pop
commonUS_Ph<-as.data.frame(matrix(nrow=length(commonUS),ncol = 6))
for (i in 1:length(commonUS)){
  commonUS_Ph[i,1]<-names(commonUS)[i]
  commonUS_Ph[i,2]<-commonUS_mean[[i]][1]
  commonUS_Ph[i,3]<-commonUS_mean[[i]][2]
  commonUS_Ph[i,4]<-commonUS_var[[i]][2]
  commonUS_Ph[i,5]<-commonUS_var[[i]][2]
  commonUS_Ph[i,6]<-nrow(commonUS[[i]])
}
colnames(commonUS_Ph)<-c("Pop","Cone_lt_mean","seedWt_10_mean","Cone_lt_var",
                         "seedWt_10_var","N")
commonUS_Ph<-merge(commonUS_Ph,commonUS_Locs[ ,2:5],by="Pop")


mexican<-mexican[!(mexican$RenamedPop%in%commonUS_Locs$RenamedPop), ]
mexicanLocs<-mexican[!(duplicated(mexican$RenamedPop)),]
SWWP<-split(mexican,f = mexican$RenamedPop,drop=TRUE)

#estimate mean and var and N per matching pop in the non-hybrid range
SWWP_mean<-lapply(SWWP, function(x) return(apply(x[ ,c(5,7)],2,mean,na.rm=TRUE)))
SWWP_var<-lapply(SWWP, function(x) return(apply(x[ ,c(5,7)],2,var,na.rm=TRUE)))
SWWP_Ph<-as.data.frame(matrix(nrow=length(SWWP),ncol = 6))
for (i in 1:length(SWWP_mean)){
  SWWP_Ph[i,1]<-names(SWWP_mean)[i]
  SWWP_Ph[i,2]<-SWWP_mean[[i]][1]
  SWWP_Ph[i,3]<-SWWP_mean[[i]][2]
  SWWP_Ph[i,4]<-SWWP_var[[i]][1]
  SWWP_Ph[i,5]<-SWWP_var[[i]][2]
  SWWP_Ph[i,6]<-nrow(SWWP[[i]])
}

colnames(SWWP_Ph)<-c("RenamedPop","Cone_lt_mean","seedWt_10_mean","Cone_lt_var",
                     "seedWt_10_var","N")
SWWP_Ph<-merge(SWWP_Ph,mexicanLocs[ ,1:4],by="RenamedPop")
SWWP_Ph<-SWWP_Ph[complete.cases(SWWP_Ph), ]

## data prep LP
PIFLlocs<-PIFL[!(duplicated(PIFL$Pop)), ]
LP<-split(PIFL,stratum = PIFL$Pop,drop=TRUE)
LP_mean<-lapply(LP, function(x) return(apply(x[ ,c(8,9)],2,mean)))
LP_var<-lapply(LP, function(x) return(apply(x[ ,c(8,9)],2,var)))
LP_Ph<-as.data.frame(matrix(nrow=length(LP),ncol = 6))
for (i in 1:length(LP)){
  LP_Ph[i,1]<-names(LP)[i]
  LP_Ph[i,2]<-LP_mean[[i]][1]
  LP_Ph[i,3]<-LP_mean[[i]][2]
  LP_Ph[i,4]<-LP_var[[i]][2]
  LP_Ph[i,5]<-LP_var[[i]][2]
  LP_Ph[i,6]<-nrow(LP[[i]])
}
colnames(LP_Ph)<-c("Pop","Cone_lt_mean","seedWt_10_mean","Cone_lt_var",
                   "seedWt_10_var","N")
LP_Ph<-merge(LP_Ph,PIFLlocs[ ,c(1,5:6)],by="Pop")

#########STEP 3#################################
#loop to find nearest neighbours for SWWP, this step was needed in our case to increase sample size because several populations did not have
#exact match between morphological dataset and genetic dataset
#We tested the robustness of the choice of nearest neighbour distance by conducting downstream analyses for varying nearest neighbour distance.
#We also assessed decay in relatdness and genetic distance as a function of increasing neighbour distance to obtain appropriate value
###########################################

#determine pops from Hybrid zone that does not have exact matches between phenotypic and genotypic dataset
Chpt1_diff<-Chapter1[!(Chapter1$RenamedPop%in%commonUS_Locs$RenamedPop), ]
neighbours<-as.data.frame(matrix(nrow = nrow(Chpt1_diff),ncol = nrow(SWWP_Ph))) #matrix of all distances

minDis<-5 #minimum dist (km) from target Gpop to keep a PhPop

Phenotype<-as.data.frame(matrix(nrow = nrow(Chpt1_diff),ncol=8))
colnames(Phenotype)<-c("RenamedPop","Cone_lt_mean","seedWt_10_mean","Cone_lt_var",
                       "seedWt_10_var","N","Phen_N","PhName")
Phenotype[ ,1]<-Chpt1_diff$RenamedPop

for (i in 1:nrow(neighbours)){
  num<-0
  for (c in 1:ncol(neighbours)){
    neighbours[i,c]<-(distVincentyEllipsoid(p1 = Chpt1_diff[i,c(4,5)], p2 =SWWP_Ph[c,c(9,8)]))/1000
    if (neighbours[i,c] <= minDis){
      #print(SWWP_Ph[c,1])
      num<-num+1
      Phenotype[i,7]<-num
      Phenotype[i,2]<-mean(c(Phenotype[i,2],SWWP_Ph[c,2]),na.rm=TRUE)
      Phenotype[i,3]<-mean(c(Phenotype[i,3],SWWP_Ph[c,3]),na.rm=TRUE)
      Phenotype[i,4]<-mean(c(Phenotype[i,4],SWWP_Ph[c,4]),na.rm=TRUE)
      Phenotype[i,5]<-mean(c(Phenotype[i,5],SWWP_Ph[c,5]),na.rm=TRUE)
      Phenotype[i,6]<-sum(c(Phenotype[i,6],SWWP_Ph[i,6]),na.rm=TRUE)
      Phenotype[i,8]<-SWWP_Ph[c,1]
    }
  }
}



Phenotype<-merge(Phenotype,Chpt1_diff[ ,c(1,4,5)],by="RenamedPop")
Phenotype<-Phenotype[complete.cases(Phenotype), ]
colnames(LP_Ph)[1]<-"Pop"
colnames(LP_Ph)[7]<-"lat"
colnames(LP_Ph)[8]<-"long"
colnames(Phenotype)[1]<-"Pop"

##FINAL DATAFRAME PREP
traits<-rbind(Phenotype[ ,c(1:6,10,9)],commonUS_Ph[ ,c(1:6,9,8)],LP_Ph)
write.table(SppPhenotype,file="Traits_5km.txt",sep="\t",row.names = F,quote=F)


##############################################################################################
######DATA MANIPULATION DONE AND SAVED AS Traits_5km#############

#######Performing cline fits in HZAR################################
################################################################################################
SppPhenotype<-read.table("Traits_5km.txt",header=T,sep="\t")
SppPhenotype<-SppPhenotype[order(SppPhenotype$lat), ]
distN <- NULL
for(i in 1:nrow(SppPhenotype)) {
  distN[i] <- distVincentyEllipsoid(p1 = SppPhenotype[1,c(8,7)], p2 =SppPhenotype[i,c(8,7)]) #longitude then latitude
}


## Make each model run off a separate seed ,one for each cline model
mainSeed<-list(A=c(978,544,99,596,528,124),B=c(544,99,596,528,124,978),C=c(99,596,528,124,978,544),
               D=c(596,528,124,978,544,99), E=c(528,124,978,544,99,596))
## Blank out space in memory to hold morphological analysis
if(length(apropos("^SWWP$",ignore.case=FALSE)) == 0 ||!is.list(SWWP) ) 
  SWWP <- list()

##For cone length data
SWWP$Cone <- list();
## Space to hold the observed data
SWWP$Cone$obs <- list();
## Space to hold the models to fit
SWWP$Cone$models <- list();
## Space to hold the compiled fit requests
SWWP$Cone$fitRs <- list();
## Space to hold the output data chains
SWWP$Cone$runs <- list();
## Space to hold the analysed data
SWWP$Cone$analysis <- list();

##For seed wt data
SWWP$Seed <- list();
SWWP$Seed$obs <- list();
SWWP$Seed$models <- list();
SWWP$Seed$fitRs <- list();
SWWP$Seed$runs <- list();
SWWP$Seed$analysis <- list();


SWWP$Cone$obs<-hzar.doNormalData1DPops(distance = distN,muObs = SppPhenotype$Cone_lt_mean,
                                       varObs = SppPhenotype$Cone_lt_var,nEff = SppPhenotype$N)
SWWP$Seed$obs<-hzar.doNormalData1DPops(distance = distN,muObs = SppPhenotype$seedWt_10_mean,
                                       varObs = SppPhenotype$seedWt_10_var,nEff = SppPhenotype$N)

###NOTE: Analysis from this point on is only for Cone, for seed data or any other morpho data simply replace the name of the list element
SWWP$Cone$models$M1<-hzar.makeCline1DNormal(data = SWWP$Cone$obs,tails = "none")
SWWP$Cone$models$M2<-hzar.makeCline1DNormal(data = SWWP$Cone$obs,tails = "right")
SWWP$Cone$models$M3<-hzar.makeCline1DNormal(data = SWWP$Cone$obs,tails = "left")
SWWP$Cone$models$M4<-hzar.makeCline1DNormal(data = SWWP$Cone$obs,tails = "mirror")
SWWP$Cone$models$M5<-hzar.makeCline1DNormal(data = SWWP$Cone$obs,tails = "both")

#change tuning parameters of complex models, needed based on preassesment of model fits
hzar.meta.tune(SWWP$Cone$models$M5)<-1.2

maxD<- 2100000 #maximum distance of the 1D transect for the hybrid zone
SWWP$Cone$models <- sapply(SWWP$Cone$models,
                            hzar.model.addBoxReq,
                            0 , maxD,
                            simplify=FALSE)
                            
##for gaussian (cline1DNormal) use fitRequest.gC
SWWP$Cone$fitRs$init <- sapply(SWWP$Cone$models,hzar.first.fitRequest.gC,
                               obsData=SWWP$Cone$obs,verbose=FALSE,simplify=FALSE)

SWWP$Cone$fitRs$init <- sapply(SWWP$Cone$fitRs$init,
                               function(mdl) {
                                 mdl$mcmcParam$chainLength <-
                                   1e5; 
                                 mdl$mcmcParam$burnin <-
                                   1e4; 
                                 mdl },
                               simplify=FALSE)

SWWP$Cone$fitRs$init$M1$mcmcParam$seed[[1]] <-mainSeed$A
SWWP$Cone$fitRs$init$M2$mcmcParam$seed[[1]] <-mainSeed$B
SWWP$Cone$fitRs$init$M3$mcmcParam$seed[[1]] <-mainSeed$C
SWWP$Cone$fitRs$init$M4$mcmcParam$seed[[1]] <-mainSeed$D
SWWP$Cone$fitRs$init$M5$mcmcParam$seed[[1]] <-mainSeed$E

#RUN the chains for each model, 6 per model equal to the number of seeds
SWWP$Cone$runs$init<-lapply(SWWP$Cone$fitRs$init,hzar.doFit)

#plot trace
plot(hzar.mcmc.bindLL(SWWP$Cone$runs$init$M3))

#intiates parallel runs of chain (from each seed), for all models
SWWP$Cone$fitRs$chains <-lapply(SWWP$Cone$runs$init,hzar.next.fitRequest)

#Repeating fit 3 replicate times using 6 different seeds per model
SWWP$Cone$fitRs$chains <-hzar.multiFitRequest(SWWP$Cone$fitRs$chains,each=3,baseSeed=NULL)
#now running this fit
SWWP$Cone$runs$chains <-  hzar.doChain.multi(SWWP$Cone$fitRs$chains,
                                              doPar=TRUE,
                                              inOrder=FALSE,
                                              count=3)

########test for convergence of each model across their 3 chains, if not converged change tuning parameter of the model/run longer
#hzar.mcmc.bindLL(x[[1]]) to x[[3]] to check values across each rep

summary(do.call(mcmc.list,lapply(SWWP$Cone$runs$chains[13:15],
                                 function(x) hzar.mcmc.bindLL(x[[3]]) )) )
#visualize
plot(do.call(mcmc.list,lapply(SWWP$Cone$runs$chains[4:6],
                              function(x) hzar.mcmc.bindLL(x[[3]]))) )           
                              
#compile analysis
SWWP$Cone$analysis$initDGs <- list( )

SWWP$Cone$analysis$initDGs<-lapply(SWWP$Cone$runs$init,hzar.dataGroup.add)
names(SWWP$Cone$analysis$initDGs)<-names(SWWP$Cone$runs$init)
SWWP$Cone$analysis$oDG <-hzar.make.obsDataGroup(SWWP$Cone$analysis$initDGs)
SWWP$Cone$analysis$oDG <-hzar.copyModelLabels(SWWP$Cone$analysis$initDGs,
                                            SWWP$Cone$analysis$oDG)

## Convert all runs to hzar.dataGroup objects, 
# adding them to the hzar.obsDataGroup object.
SWWP$Cone$analysis$oDG <-hzar.make.obsDataGroup(lapply(SWWP$Cone$runs$chains,
                          hzar.dataGroup.add),SWWP$Cone$analysis$oDG)
print(summary(SWWP$Cone$analysis$oDG$data.groups))

#can only compare to null model graphically
hzar.plot.cline(SWWP$Cone$analysis$oDG)


##Final best fit model based on AICc
print(SWWP$Cone$analysis$AICcTable <-
        hzar.AICc.hzar.obsDataGroup(SWWP$Cone$analysis$oDG))
print(SWWP$Cone$analysis$model.name <-rownames(SWWP$Cone$analysis$AICcTable
                    )[[ which.min(SWWP$Cone$analysis$AICcTable$AICc )]])

SWWP$Cone$analysis$model.selected <-
  SWWP$Cone$analysis$oDG$data.groups [[SWWP$Cone$analysis$model.name]]

##get variation around the best 2LL estimate space (similar to C.I)
print(hzar.getLLCutParam(SWWP$Cone$analysis$model.selected,
                          names(SWWP$Cone$analysis$model.selected$data.param)))
print(hzar.get.ML.cline(SWWP$Cone$analysis$model.selected))

##Plot the cline and the 2LL around it
hzar.plot.fzCline(SWWP$Cone$analysis$model.selected)                             
                            
