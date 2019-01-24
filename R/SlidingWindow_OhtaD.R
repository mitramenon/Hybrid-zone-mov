###This script is used to analyse the empirical dataset used in this study but the function estOhta was used for analysing simulation output too
###For the simulation output, loci in the top10 percentile of allele frequency difference were first selected as the TopCandidates (See associated paper)

library(data.table)
library(ohtadstats)

#load in genetic and phenotypic data
Pilotdata<-fread("good_snps.recode.vcf.gz.012",sep="\t",data.table=F)
Inds<-read.table("good_snps.recode.vcf.gz.012.indv")
 PhPops<-read.table("Traits_5km.txt",header=T,sep="\t")

 Pilotpops<-read.table("Pops.txt",header=T,sep="\t")
 head(Pilotpops)

 pos<-read.table("good_snps.recode.vcf.gz.012.pos")

 pos$location<-paste(pos[ ,1],pos[ ,2],sep=":")
 colnames(Pilotdata)<-pos$location
 Pilotdata<-as.matrix(Pilotdata)
 Pilotdata[Pilotdata==-1]<-NA
 Pilotdata<-cbind(Inds,Pilotpops,Pilotdata)

#Utilize the SNPs in top10 percentile (nearly diagnostic) for conducting analysis

 FstTop<-read.table("HierachialFst_forCline.txt",header=T,sep="\t")
 TopSNPs<-Pilotdata[ ,colnames(Pilotdata)%in%FstTop$location]
 TopSNPs<-cbind(Pilotdata[ ,1:6],TopSNPs) #add back the pop and coordinate information

#Since we need to generate windows to estimate change in DIS/DST, we will order the pops by latitude. 
#This will depend on the gradient of hybridization, if the gradient is east-west, you'll want to order by Longitude

 TopSNPs<-TopSNPs[order(TopSNPs$lat), ]
 PopsUnq<-TopSNPs[!(duplicated(TopSNPs$RenamedPop)), ]
 PopsUnq<-PopsUnq$RenamedPop
 
 #####generate a list to hold all the pop names with overlap of one##########
names<-NULL

e<-4 #determines the window size
s<-1
l<-1

while (s < 34){
names[[l]]<-PopsUnq[s:e]
s<-e+1
e<-e+3 #overlap every 4th pop betwe consequtive elements of the list
l<-l+1
}

length(names)
names

#####Function to get Ohta's D#############################
##Function takes single input and outputs mean Dis and mean Dst##
############################################################

estOhta<-function(df){

#Input is a vector (df) with list of pops to perform variance partioning for

  set<-TopSNPs[TopSNPs$RenamedPop%in%df, ] #subset to retain only the 4 pops
  
  nam<-set$RenamedPop #keep IDs
  set<-set[ ,-c(1:7)]
  set<-as.matrix(set)
  
  rownames(set)<-nam  #this is important to add to the matrix eventually, otherwise Ohta's estimator will not work
  
  OHTD<-dwrapper(set, tot_maf = 0.01, pop_maf = 0.01)
  dst<-mean(OHTD$d2st_mat,na.rm=TRUE)
  dis<-mean(OHTD$d2is_mat,na.rm=TRUE)
  #dstP<-mean(OHTD$dp2st_mat,na.rm=TRUE)
  #disP<-mean(OHTD$dp2is_mat,na.rm=TRUE)
  
  out<-c(dst,dis)

  return(out)

  }

#############END OF FUNCTION##############################################################

D_stat<-lapply(names,function(x) return(estOhta(x)))
D_stat<-do.call(rbind,D_stat)
colnames(D_stat)<-c("Dst","Dis")

write.table(D_stat,file="Dstat_overlap1_set4.txt",row.names=F,quote=F,sep="\t")
