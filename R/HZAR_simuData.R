######################################################################
######Geographic cline analysis using Qscore estimated via faststructure##########
######Script estimates cline centre and width for each generation output per MC run per Model
######Author: menonm2 (menonm2@mymail.vcu.edu)
######Last modified: 9-Nov-2018
##########################################################################

    library(hzar)
    library(data.table)
    library(geosphere)

PhPops<-read.table("~/eckertlab/Mitra/HZ_mov/Traits_5km.txt",header=T,sep="\t")

###Generate a vector of the generations###
mid<-c(500,520)
In<-seq(570,1020,by=50)
In<-c(mid,In)

###Generate vectors for all the structure output file and structure input file which contains the coordinates
meanQ<-paste0("strFiles/ind",In,".2.meanQ")
coords<-paste0("strFiles/ind",In,".str")

##Generating null list equal to the number of generations it needs to analyse
files<-vector("list",length(meanQ))

for (f in 1:length(files)){

files[[f]][1]<-meanQ[f]
files[[f]][2]<-coords[f]

}

##################################################################################
########Long function to conduct cline anlysis#############################
#assumes that structure output is from K=2
#Function takes a single input, Infile is a list with two elements:
##the first being structure output and second being coordinate information
##Function output the ML estimate of centre, slope and the name of the best fit model per generation
####################################################################################
##***********START FUNCTION********************##

estCline<-function(Infile){

    K2<-read.table(Infile[[1]],sep=" ")
    K2<-K2[ ,-2]
    colnames(K2)<-c("cluster1","cluster2")
    
    ID<-fread(Infile[[2]],sep=" ",header=F,data.table=F)
    ID<-ID[!(duplicated(ID$V4)),c(1:4,6)]
    colnames(ID)<-c("patchID","X","Y","Ind","Pop")
    dim(ID)

  #estimate meanQ by population ID (here indicated by patchID)
    K2<-cbind(ID,K2)
    K2pops<-split(K2,f=K2$patchID)
    K2mean<-lapply(K2pops, function(df) return(mean(df[ ,7])))
    K2mean<-do.call(rbind,K2mean)
    K2mean<-cbind(ID[!(duplicated(ID$patchID)), ],K2mean)

    sampleN<-as.data.frame(table(K2$patchID))
    colnames(sampleN)[1]<-"patchID"

    K2mean<-merge(K2mean,sampleN,by="patchID")
    #order by increaseing Y values on the X-Y landscape
    K2mean<-K2mean[order(K2mean$Y), ]
    head(K2mean)

  #estimate geographical distance using X and Y
     distN <- NULL
      for(d in 1:nrow(K2mean)) {
        distN[d] <- dist(rbind(K2mean[1,c(2,3)],K2mean[d,c(2,3)]))
      }

  #This step is important to conduct cline analyses, it determines the max extent of the landscape. Could be increased by a few more units
  #Tried increasing this value but did not have much effect on the overall pattern
  
      maxD<-round(max(distN))
  
  #Setup seeds to generate replicate models
  
      rotateModelSeeds <- function(fitL)
        hzar.multiFitRequest(fitL,
                             rotateSeed=TRUE, skip=50,
                             baseChannel=NULL, each=1,
                             baseSeed=c(596,528,124,978,544,99))

      if(length(apropos("^SWWP$",ignore.case=FALSE)) == 0 ||
         !is.list(SWWP) ) SWWP <- list()

  #setup null list object to hold elements

       SWWP$cnppd <- list();
      SWWP$cnppd$obs <- list();
      SWWP$cnppd$models <- list();
      SWWP$cnppd$fitRs <- list();
      SWWP$cnppd$runs <- list();
      SWWP$cnppd$analysis <- list();

  #Fit cline, remember to provide the number of samples per pop here
      SWWP$cnppd$obs <-hzar.doMolecularData1DPops(distN,K2mean$K2mean,K2mean$Freq)

     SWWP.loadcnppdmodel <- function(scaling,tails,
                                      id=paste(scaling,tails,sep="."))
      SWWP$cnppd$models[[id]] <<- hzar.makeCline1DFreq(SWWP$cnppd$obs, scaling, tails)

      SWWP.loadcnppdmodel("none" ,"none"  ,"typN");
      SWWP.loadcnppdmodel("none" ,"left"  ,"typL");
      SWWP.loadcnppdmodel("none" ,"right" ,"typR");
      SWWP.loadcnppdmodel("none" ,"mirror","typM");
      SWWP.loadcnppdmodel("none" ,"both"  ,"typB");
      SWWP.loadcnppdmodel("fixed","none"  ,"fixN");
      SWWP.loadcnppdmodel("fixed","left"  ,"fixL");
      SWWP.loadcnppdmodel("fixed","right" ,"fixR");
      SWWP.loadcnppdmodel("fixed","mirror","fixM");
      SWWP.loadcnppdmodel("fixed","both"  ,"fixB");
      SWWP.loadcnppdmodel("free" ,"none"  ,"optN");
      SWWP.loadcnppdmodel("free" ,"left"  ,"optL");
      SWWP.loadcnppdmodel("free" ,"right" ,"optR");
      SWWP.loadcnppdmodel("free" ,"mirror","optM");
      SWWP.loadcnppdmodel("free" ,"both"  ,"optB");

      SWWP$cnppd$models <- sapply(SWWP$cnppd$models,
                                  hzar.model.addBoxReq,
                                  0 , maxD,
                                  simplify=FALSE)

      SWWP$cnppd$fitRs$init <- sapply(SWWP$cnppd$models,
                                      hzar.first.fitRequest.old.ML,
                                      obsData=SWWP$cnppd$obs,
                                      verbose=FALSE,
                                      simplify=FALSE)

      SWWP$cnppd$fitRs$init <- sapply(SWWP$cnppd$fitRs$init,
                                      function(mdl) {
                                        mdl$mcmcParam$chainLength <-
                                          1e5; #1e5 is ideal, but shorter might work too
                                        mdl$mcmcParam$burnin <-
                                          1e4; #1e4 is ideal, but shorter might work too
                                        mdl },
                                      simplify=FALSE)

      SWWP$cnppd$fitRs$init <-
        rotateModelSeeds(SWWP$cnppd$fitRs$init)

      SWWP$cnppd$runs$init <-
        hzar.doFit.multi(SWWP$cnppd$fitRs$init) #this will automate sequential fitting
      names(SWWP$cnppd$runs$init) <- names(SWWP$cnppd$fitRs$init)

      SWWP$cnppd$fitRs$chains <-
        lapply(SWWP$cnppd$runs$init,
               hzar.next.fitRequest)

      SWWP$cnppd$fitRs$chains <-
        hzar.multiFitRequest(SWWP$cnppd$fitRs$chains,
                             each=3,
                             baseSeed=NULL)

      SWWP$cnppd$runs$chains <-  hzar.doChain.multi(SWWP$cnppd$fitRs$chains,
                                                    doPar=TRUE,
                                                    inOrder=FALSE,
                                                    count=3)

      SWWP$cnppd$analysis$initDGs <- list(
        nullModel =  hzar.dataGroup.null(SWWP$cnppd$obs))

      SWWP$cnppd$analysis$initDGs <-
        c( SWWP$cnppd$analysis$initDGs,
           sapply(SWWP$cnppd$runs$init,
                  hzar.dataGroup.add,
                  simplify=FALSE))
                  
      print(names(SWWP$cnppd$analysis$initDGs))

      SWWP$cnppd$analysis$oDG <-
        hzar.make.obsDataGroup(SWWP$cnppd$analysis$initDGs)
      SWWP$cnppd$analysis$oDG <-
        hzar.copyModelLabels(SWWP$cnppd$analysis$initDGs,
                             SWWP$cnppd$analysis$oDG)

      SWWP$cnppd$analysis$oDG <-
        hzar.make.obsDataGroup(lapply(SWWP$cnppd$runs$chains,
                                      hzar.dataGroup.add),
                                   SWWP$cnppd$analysis$oDG)
      
      #assesing best fit model in relation to the null and to others based on AICc
      oDGkey <- which(!(names(SWWP$cnppd$analysis$oDG$data.groups) %in%
                          "nullModel"))
      print(hzar.getLLCutParam(SWWP$cnppd$analysis$oDG$data.groups[oDGkey ],
                               c("center","width")));
      print(SWWP$cnppd$analysis$AICcTable <-
              hzar.AICc.hzar.obsDataGroup(SWWP$cnppd$analysis$oDG))
      print(SWWP$cnppd$analysis$model.name <-
              rownames(SWWP$cnppd$analysis$AICcTable
              )[[ which.min(SWWP$cnppd$analysis$AICcTable$AICc )]])

      SWWP$cnppd$analysis$model.selected <-
        SWWP$cnppd$analysis$oDG$data.groups[[SWWP$cnppd$analysis$model.name]]

  #Get ML estimates from the best fit model
     Best<-hzar.get.ML.cline(SWWP$cnppd$analysis$model.selected)

     center<-Best[[1]][1]
     width<-Best[[1]][2]
     type<-SWWP$cnppd$analysis$model.name
     params<-c(center,width,type)

     #cat("The parameters are", params)


     return(params) 
    }


   output<-lapply(files,function(x) return(estCline(x)))
   output<-do.call(rbind,output)
   output<-cbind(coords,output)


write.table(output,file="CDPop/data_July2018/ModelB_ItoIII_2_iii_Genes_1537378126/batchrun0mcrun0/Cline_Biii2_batch0.txt",row.names=F,quote=F,sep="\t")

