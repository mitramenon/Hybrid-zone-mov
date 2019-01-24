#Script assumes that you are in a directory with all the required generation files, if not be sure to change paths while generating vectors


library(data.table)

##Generate vector to hold all generation names for which output is available
phase1<-seq(from=50,to=500,by = 50)
mid<-c(504,510,514,520)
phase1b<-seq(570,1020,by=50)
phase1<-c(phase1,mid,phase1b)

#Vector to hold names of all generation files
phase1_files<-vector("list",length(phase1))
phase1_files<-paste0("ind",phase1,".csv")

##################Function to convert CDPop outfile to strFormat and write it out#############################
##############################################################################################################

strConvert<-function(N,s,e){

#****Function takes 3 inputs######
#N has to be a vector contianing the names of all files (.csv format) to be processed
#s is the column number at which loci info begins
#e is the last column with loci info
#***********************************

print(N)

gen<-fread(N,header=T,sep=",",data.table=F)
Info<-gen[ ,c(1:5,17)]
loci<-gen[ ,s:e]

#Combine every subsequent allele column per loci to get coded as 0,1,2
combined<-data.frame( loci[1], mapply( paste0, loci[-1][c(T,F)], loci[-1][c(F,T)] ) )
combined<-apply(combined[ ,-1],MARGIN = 2,as.character)

#remove monomorphic snps
mono <- numeric(ncol(combined))

for (i in 1:ncol(combined)) {
  mono[i] <- length(table(combined[,i])) #techinically hets can be included but very rarely will we have all hets at a loci
}
names(mono)<-colnames(combined)
remove<-mono[mono==1]
snp_reformat<-combined[,!(colnames(combined)%in%remove)]

#convert to structure format, two rows per individual, corresponding to the diploid genotype format for STRUCTURE
formated<- apply(combined, 2, function(df) gsub('11','1/2',df))
formated <- apply(formated, 2, function(df) gsub('02','1/1',df))
formated <- apply(formated, 2, function(df) gsub('20','2/2',df))
formated[is.na(formated)]<-'-9'

Info.expand<-Info[rep(row.names(Info), each=2), ]

## split each row at a time into two rows by the / and then rbind all
strFormat<-apply(formated,1, function(x) return(strsplit(x,"/")))
strFormat2<-lapply(strFormat,function(x) return(t(do.call(rbind,x))))
strFinal<-do.call(rbind,strFormat2)

strFinal<-cbind(Info.expand,strFinal)

#Generate name for the Outfile based on name of the Infile
outfile<-paste0("strFiles/",strsplit(N,".csv")[[1]],".str")

write.table(strFinal,file=outfile,sep=" ",col.names =F,quote=F,row.names = F)

}

#########END OF FUNCTION#####################################################################3

#Perform conversion
lapply(phase1_files,function(x) return(strConvert(x)))
