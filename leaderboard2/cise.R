library(RWeka)
library(cluster)
args<-commandArgs(TRUE)
vcfdat = read.table(args[1],sep='\t',comment.char='#')
batdat = read.table(args[2],sep='\t',comment.char='#')
datacol=as.integer(args[3]) + 9
namecol=9
tumour_stat = data.frame(do.call(rbind, strsplit(as.vector(vcfdat[,datacol]), split = ":", fixed = TRUE)))
colnames(tumour_stat) = strsplit(as.vector(unique(vcfdat[,namecol])),':')[[1]]
TDP = as.integer(as.vector(tumour_stat[,'DP']))
TAD1 <-  as.integer(unlist(lapply(strsplit(as.vector(tumour_stat[,'AD']),','),'[[',1)))
TAD2 <-  as.integer(unlist(lapply(strsplit(as.vector(tumour_stat[,'AD']),','),'[[',2)))
BQ <- as.integer(as.vector(tumour_stat[,'BQ']))
if(datacol==11){datacol=10}else{datacol=11}
tumour_stat = data.frame(do.call(rbind, strsplit(as.vector(vcfdat[,datacol]), split = ":", fixed = TRUE)))
colnames(tumour_stat) = strsplit(as.vector(unique(vcfdat[,namecol])),':')[[1]]
NDP <- as.integer(as.vector(tumour_stat[,'DP']))
NAD1 <-  as.integer(unlist(lapply(strsplit(as.vector(tumour_stat[,'AD']),','),'[[',1)))
NAD2 <-  as.integer(unlist(lapply(strsplit(as.vector(tumour_stat[,'AD']),','),'[[',2)))
chrvcf <-as.vector(vcfdat[,1])
posvcf<-as.integer(as.vector(vcfdat[,2]))
vcfref<- as.vector(vcfdat[,4])
vcfalt<- as.vector(vcfdat[,5])
chrbat <- as.vector(batdat[2:length(batdat[,1]),1])
strbat <- as.integer(as.vector(batdat[2:length(batdat[,1]),2]))
endbat <- as.integer(as.vector(batdat[2:length(batdat[,1]),3]))
BAFbat <- as.double(as.vector(batdat[2:length(batdat[,1]),4]))
LRbat <- as.double(as.vector(batdat[2:length(batdat[,1]),6]))
BAF=NULL
LR=NULL
A=NULL
T=NULL
C=NULL
G=NULL
tA=NULL
tT=NULL
tC=NULL
tG=NULL
BAF[length(TDP)]=NA
LR[length(TDP)]=NA
l = 1
print(length(TDP))
for(i in 1 : length(TDP)){
  for(j in 1:(length(batdat[,1])-1)){
    if(vcfref[i]=="A"){
      A[i]=1
    }else{
      A[i]=0
    }
    if(vcfref[i]=="T"){
      T[i]=1
    }else{
      T[i]=0
    }
    if(vcfref[i]=="C"){
      C[i]=1
    }else{
      C[i]=0
    }
    if(vcfref[i]=="G"){
      G[i]=1
    }else{
      G[i]=0
    }
    if(vcfalt[i]=="A"){
      tA[i]=1
    }else{
      tA[i]=0
    }
    if(vcfalt[i]=="T"){
      tT[i]=1
    }else{
      tT[i]=0
    }
    if(vcfalt[i]=="C"){
      tC[i]=1
    }else{
      tC[i]=0
    }
    if(vcfalt[i]=="G"){
      tG[i]=1
    }else{
      tG[i]=0
    }
    if(chrvcf[i] == chrbat[j] && posvcf[i]>=strbat[j] && posvcf[i]<=endbat[j]){
      BAF[i]=BAFbat[j]
      LR[i]=LRbat[j]
    }
  }
}
nNAD2=NAD2
BBAF=BAF
TFA= TAD2/TDP
NFA= NAD2/NDP
CNV=TDP/NDP
MR=(TAD2+NAD2)/(TAD1+TAD2+NAD1+NAD2)
ttt= data.frame(NAD2,TDP,NDP,MR,CNV,BQ,A,T,G,C,tA,tT,tC,tG,TFA)
ttt=as.data.frame(scale(ttt))
LR=scale(LR)
BAF=scale(BAF)
LR=as.factor(LR)
BAF = as.factor(BAF)
ttt=data.frame(BAF,LR,ttt)
write.arff(ttt,file="tempfile.arff")
system("export WEKA_HOME=/d/wekafiles")
out = system("java -cp /usr/share/java/weka.jar -mx1024m weka.attributeSelection.CfsSubsetEval -P 1 -E 1 -i tempfile.arff -s \"weka.attributeSelection.BestFirst -D 1 -N 5
\"",TRUE)
a=grep("Selected attributes:",out)
john=strsplit(out[a]," ",fixed=T);
john1=strsplit(john[[1]]," ");
if(grepl("4",john1[3])){
}else{
    john1[3]=paste(john1[3],",4",sep="")
}
if(grepl("12",john1[3])){
}else{
    john1[3]=paste(john1[3],",12",sep="")
}
call=paste(john1[3],",17\"",sep="")
call=paste("\"",call,sep="")
call=paste("java -cp /usr/share/java/weka.jar -mx1024m weka.filters.unsupervised.attribute.Remove -i tempfile.arff -o temp.arff -V -R",call,sep=" ")
system(call)
out=system("java -cp /usr/share/java/weka.jar -mx1024m weka.clusterers.EM -t temp.arff -I 100 -N -1 -X 10 -max -1 -ll-cv 1.0E-6 -ll-iter 1.0E-6 -M 1.0E-6 -K 100 -num-slots 23 -S 100",TRUE)
clusterNum = 0
 
for(c in 1 : length(out)) {
   if(out[c] == "Clustered Instances") {
     c = c + 1
     for(p in c : length(out)) {
       if(out[p] != "") {
         clusterNum = clusterNum + 1;
       }else {
         break
       }
    }
     break
   }
}
print (clusterNum)
pam(ttt,clusterNum)
runPam= data.frame(TFA,BQ,NAD2,TDP,MR)
AAAA=pam(runPam,clusterNum)
cellularity=NULL
clust=AAAA$clustering
print(clust)
print(length(clust))

fuckkkkk=as.data.frame(clust)
for(j in 1 : clusterNum){
  m0 = NULL
  B0 = NULL
  mk0=1
  kk0=1
  for(i in 1 : length(TDP)){
    if(clust[i]==j){
      a1=TAD2[i]/(TAD2[i]+TAD1[i])
      a2=TAD1[i]/(TAD2[i]+TAD1[i])
      b1=nNAD2[i]/(nNAD2[i]+NAD1[i])
      b2=NAD1[i]/(nNAD2[i]+NAD1[i])
      a=c(a1,a2)
      b=c(b1,b2)
      qqq=rbind(a,b)
      if(is.na(BBAF[i])){
        m0[kk0]=daisy(qqq,matrix("manhattan"))
      }else{
        m0[kk0]=daisy(qqq,matrix("manhattan"))*(1.5-BBAF[i])
      }
      kk0=kk0+1
    }
  }
  cellularity[j]=mean(m0)
}
mcell=min(cellularity)
for(i in 1 : clusterNum ){
   if(cellularity[i] == mcell){
     cellularity[i]=0
     tempi = i
   }
}
countNum=NULL
for(j in 1: clusterNum){
    countNum[j]=0
    for(i in 1:length(TDP)){
        if(clust[i]==j)
            countNum[j]=countNum[j]+1
    }
}

if(tempi!=clusterNum){
    cellularity[tempi]=cellularity[clusterNum]
    cellularity[clusterNum]=0
    
    tempnum=countNum[tempi]
    countNum[tempi]=countNum[clusterNum]
    countNum[clusterNum]=tempnum
    for(i in 1:length(clust)){
        if(clust[i]==tempi){
            clust[i]=clusterNum
        }else if (clust[i]==clusterNum){
            clust[i]=tempi
        }
    }
 }
 write.table(max(cellularity),"subchallenge1A.txt",row.names=F,col.names=F,quote=F,sep="\t")
 write.table(cbind(1:clusterNum,countNum,cellularity),"subchallenge1C.txt",row.names=F,col.names=F,quote=F,sep="\t")
 write.table(clusterNum-1,"subchallenge1B.txt",row.names=F,col.names=F,quote=F,sep="\t")
 write.table(clust,"subchallenge2A.txt",row.names=F,col.names=F,quote=F,sep="\t")