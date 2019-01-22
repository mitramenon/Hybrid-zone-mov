##modifying CDPop genetic output for hierf and estimating Hierarchial Fst########
###author: menonm2 (menonm2@mymail.vcu.edu)
##Date:20-July-2018####
########################################################

library(data.table)
library(hierfstat)
library(dplyr)

#loading files for each run, 1 model (several generations) per script
model<-list.files("SimuRuns/data_July2018/ModelB_ItoIII_RUN2_i_Genes_1533254936/batchrun0mcrun0/",pattern = ".csv")
phase1<-seq(from=50,to=500,by = 50)
mid<-c(504,510,514,520)
phase1b<-seq(570,1020,by=50)
phase1<-c(phase1,mid,phase1b)

phase1_files<-vector("list",length(phase1))
phase1_files<-paste0("SimuRuns/data_July2018/ModelB_ItoIII_RUN2_i_Genes_1533254936/batchrun0mcrun0/ind",phase1,".csv")

##function to read files, modify and run hierfstat###############
##START FUNCTION###########
#Since we are implementing Hierfstat, we assume that CDPop was run in a hierachial manner with SubPops and Populations
#In the outputs used here, SubPops are called PatchID and Populations are called SubPatchID. PatchID is nested within SubPatchID

myRead<-function(N,s,e){
  ************************************************************************
  ###Function takes three arguments:
  ##N: a dataframe that is outputed from CDMetaPop, one per generation per Model
  ##s: column at which the loci information begin
  ##e: column where the loci information ends
  **************************************************************************
  print(N)
  file<-fread(N,header=T,sep=",",data.table=F)
  loci<-file[ ,c(s:e)] #first several columns of CDPop genetic output lists mating and movement information, not neeeded for hierfstat
  
  ##combine every conseqtive column into one, starting from the second column
  combined<-data.frame( loci[1], mapply( paste0, loci[-1][c(T,F)], loci[-1][c(F,T)] ) )
  combined<-apply(combined[ ,-1],MARGIN = 2,as.character)
  
  #convert to hierf format
  formated<- apply(combined, 2, function(df) gsub('11','het',df))
  formated <- apply(formated, 2, function(df) gsub('02','maj',df))
  formated <- apply(formated, 2, function(df) gsub('20','min',df))
  
  formated<- apply(formated, 2, function(df) gsub('het','12',df))
  formated <- apply(formated, 2, function(df) gsub('maj','22',df))
  formated <- apply(formated, 2, function(df) gsub('min','11',df))
  
  formated<-apply(formated,2,as.numeric)
  formated<-as.data.frame(formated)
  
  print("estimating fst")

  Groups<-data.frame(file$SubPatchID,file$PatchID)
  diff<-varcomp.glob(Groups,formated)

  FCT<-diff$F[1,1]
  DF<-cbind(as.character(file$PatchID),formated)
  fst<-wc(DF,diploid=TRUE)

  out<-c(fst$FST,FCT)
  return(out)
  
}

##
FILES<-lapply(phase1_files,function(x) return(myRead(x)))
FILES<-do.call(rbind,FILES)
FILES<-cbind(phase1,FILES)

write.table(FILES,file="SimuRuns/outputs/ModelB_i_run2_batch0.txt",row.names = F,quote = F,sep="\t")




