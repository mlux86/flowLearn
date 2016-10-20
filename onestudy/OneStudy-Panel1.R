library(flowCore)
library(flowDensity)
#library(pryr)
source('helperfunc.R', echo=TRUE)
#mem.use <-c()
#mem.use <- c(mem.use,mem_used())
#ptm <- proc.time()
panel <- "ONE_1"
#type is "" for all datas except the CNTRP that has whole blood samples and PBMCs
type <-c("","WBC","PBMCs")[2]
centre <-c("DuraClone-T1D","German","CNTRP","CNTRP-PhaseII")[3]
#Files are in Bioinformatics folder
path <- paste("/mnt/f/FCS data/OneStudy",centre,panel,sep="/")
if (centre=="CNTRP")
  path <- paste(path,type,sep="/")

dir.create(paste("~/Downloads/One-Study",centre,"Results",type,panel,sep="/"),recursive=T)
res <- paste("~/Downloads/One-Study",centre,"Results",type,panel,sep="/")
#path to FCS files
fcs.path <- list.files(path,full.names=T,pattern=".LMD",recursive = T)

#Read the second dataset in the LMD file as it's not transformed or compensated
#In here we only keep the markers names from the first data set in order to use it for the second dataset within each FCS file
#Reading files
#------------------------------------------------------------------------------------

f <- read.FCS(fcs.path[1],dataset = 1)  
markers <- unlist(lapply(as.vector(f@parameters@data[,2]), function(x) unlist(strsplit(x," "))[1]))
markers <- markers[-c(1:3)]
fs <- read.flowSet(fcs.path,dataset = 2)
fs<- fsApply(fs, function(x) {
x@parameters@data[,2] <- c("FS-A","FS-W","SS-A",markers);
return(x)
})

sampleNames(fs)<-unlist(lapply(fcs.path, function(x) return(gsub("/","-",unlist(strsplit(unlist(strsplit(x,"_1/"))[2]," "))[1]))))
#Reading compensation files
print(sampleNames(fs))
#mem.use <- c(mem.use,mem_used())
#-----------------------------------------------------------------------------------
if(centre=="German")
{
  comp.path <- list.files(path,full.names=T,pattern=".rds")
  names(comp.path)<-sampleNames(fs)
}
if(centre=="CNTRP")
{
#   if ( type=="WBC")
#   {
     comp.path <- list.files(paste("/mnt/f/FCS data/OneStudy/CNTRP/Comp/",type,sep = ""),full.names=T,pattern=".csv")
     names(comp.path) <- unlist(lapply(comp.path, function(x) unlist(strsplit(x,"-"))[2]))
     sub.centre <- unlist(lapply(sampleNames(fs), function(x) return(unlist(strsplit(x,"-"))[2])))
     names(sub.centre) <-sampleNames(fs)
#   }else{
#     comp.path <- list.files("/mnt/f/FCS data/OneStudy/CNTRP/Comp/PBMCs",full.names=T,pattern=".csv")
#   }
}

  
#Marker names
#---------------------------------------------------------------------------------------
#****************************************************************************************
#                    Finding markers
#****************************************************************************************
channels.ind <- c(1:ncol(fs[[1]]))
names(channels.ind) <- c("FS-A","FS-W","SS-A",markers)

#Doublet removal
#----------------------------------------------------------------------------------------
#*****************************************************************************************
#                      Removing doublets
#*****************************************************************************************


fsw.gate <-fsApply(fs,function(x) return(tail(getPeaks(frame = x,chans = channels.ind["FS-W"],tinypeak.removal=.4)$Peaks,1)+4*sd(exprs(x)[,channels.ind["FS-W"]])))
names(fsw.gate) <- sampleNames(fs)
dev.set(2)
par(mfrow=c(3,4)) 
Singlet <- fsApply(fs, function(x) {
  print(identifier(x))
  x <- removeMargins(x, chans = c(1:3))
  singlets <- flowDensity(x,channels=channels.ind[c(2,1)],position=c(F,NA),gates=c(fsw.gate[identifier(x)],NA))
  #plot(density(exprs(x)[,2]))
# plotDens(x,channels=channels.ind[c(2,1)],main=identifier(x))#,nbin=3000
 # lines(singlets@filter, col=1, lty=3, lwd=3)
  return(singlets)
})
fs.sngl <-lapply(Singlet, function(x) getflowFrame(x))
fs.sngl <- as(object=fs.sngl, Class="flowSet")

#Compensating and transforming files
#----------------------------------------------------------------------------------------
#*****************************************************************************************
#                        /Compensation/
#*****************************************************************************************
#Using the compensation matrix provided

if (centre=="German")
{
  fs.sngl <- fsApply(fs.sngl, function(x) {
    comp <-  readRDS(comp.path[identifier(x)])[[1]][1]$KalComp
    colnames(comp) <- colnames(x)[channels.ind[4:length(channels.ind)]]
    rownames(comp)<-NULL
    return(compensate(x, comp))})
}
if(centre=="DuraClone-T1D")
{
  comp<- read.csv("/mnt/f/FCS data/OneStudy/CNTRP/Comp/WBC/Compensation-YVR-.csv",check.names = F)
  comp <- comp[,-1]
  fs.sngl <- fsApply(fs.sngl, function(x) return(compensate(x, comp)))
}
if(centre=="CNTRP")
{
  fs.sngl <- fsApply(fs.sngl, function(x) {
    comp <-  read.csv(comp.path[sub.centre[identifier(x)]],check.names = F)
    comp <- comp[,-1]
    return(compensate(x, comp))})
}
if(centre=="CNTRP-PhaseII")
{
  comp<- read.csv("/mnt/f/FCS data/OneStudy/CNTRP-PhaseII/Compensation-PhaseII.csv",check.names = F)
  comp <- comp[,-1]
  fs.sngl <- fsApply(fs.sngl, function(x) return(compensate(x, comp)))
}
#*****************************************************************************************
#                      /Transformation/
#*****************************************************************************************

if (centre=="CNTRP")
{
  temp <- lapply(unique(sub.centre), function(x){
    return(estimate.logicle(fs.sngl[which(sub.centre==x)],med=F,trans.chans=grep(colnames(fs[[1]]),pattern="FL"),estimate=T))
  })
  fs.temp <-temp[[1]]
  if ( length(temp) >1){
  for ( i in 2:length(temp))
    fs.temp <- rbind2(fs.temp,temp[[i]])
  }
  fs.sngl <- fs.temp
}else{
  fs.sngl <-estimate.logicle(fs=fs.sngl,med=F,trans.chans=grep(colnames(fs[[1]]),pattern="FL"),estimate=T)
}
#Gating CD45, and Granulocytes
#----------------------------------------------------------------------------------------
#*****************************************************************************************
#                      Gating CD45
#*****************************************************************************************

cd45.gate <- fsApply(fs.sngl,deGate,channels.ind["CD45"],percentile=NA, tinypeak.removal=1/40,upper=F,alpha=.01)
#This was for WBC and T1D and worked perfectly but it didn't work on PBMCs,
#The above threshold need to be checked for these datasets as well
##cd45.gate <-averageGates(fsApply(fs.sngl,deGate,channels.ind["CD45"]),sd.coeff = 2)
names(cd45.gate) <- sampleNames(fs.sngl)
dev.set(2)
par(mfrow=c(5,5))
CD45 <- fsApply(fs.sngl, function(x) {
  print(identifier(x))
  cd45 <- flowDensity(x, channels.ind[c("CD45","SS-A")],position=c(T,NA),gates=c(cd45.gate[identifier(x)],NA))
  #plotDens(x,channels=channels.ind[c("CD45","SS-A")],main=identifier(x))#,nbin=3000
  #lines(cd45@filter, col=1, lty=3, lwd=3)
  return(cd45)
})

fs.45 <-lapply(CD45, function(x) getflowFrame(x))
fs.45 <- as(object=fs.45, Class="flowSet")

ss.gate <- fsApply(fs.45,deGate,channels.ind["SS-A"],percentile=.999,tinypeak.removal=1/5)
gran.gate <- averageGates(fsApply(fs.45,deGate,channels.ind["CD14"]))
names(ss.gate) <- names(gran.gate) <- sampleNames(fs.45)

dev.set(2)
par(mfrow=c(3,3))
PBMC <- fsApply(fs.45, function(x){
 temp <- flowDensity(x,channels.ind[c("CD14","SS-A")],position = c(F,F),gates=c(gran.gate[identifier(x)], ss.gate[identifier(x)]))
 ss.lo<- deGate(getflowFrame(temp),channels.ind["SS-A"],percentile = NA,upper=T,alpha=.01)
 ss.lo <-ifelse(ss.lo < deGate(getflowFrame(temp),channels.ind["SS-A"],percentile = .5,use.percentile=T),no = ss.lo,
        yes = deGate(getflowFrame(temp),channels.ind["SS-A"],upper=T,use.upper =T,tinypeak.removal=.9,alpha=.05) )
 s.gate <- min(ss.lo, ss.gate[identifier(x)]*.95)
 gran <- flowDensity(x,channels.ind[c("CD14","SS-A")],position = c(F,T),gates=c(gran.gate[identifier(x)], s.gate))
  not.gran <- notSubFrame(x,channels.ind[c("CD14","SS-A")],position = c(F,T),gates =gran@gates)
  #plotDens(x,channels.ind[c("CD14","SS-A")],main="CD45+")
 # lines(gran@filter,lty=3,lwd=3)
  not.gran@filter <- gran@filter
  return(not.gran)
})

fs.pbmc<-lapply(PBMC, function(x) getflowFrame(x))
fs.pbmc <- as(object=fs.pbmc, Class="flowSet")

#The old version of gating granulocytes
# dev.set(2)
# par(mfrow=c(3,4))
# PBMC<- fsApply(fs.45, function(x){
#   x2<- rotate.data(x,chans = c(1,3),theta =pi/5 )$data
#   fs.gate <- min(deGate(x2,1,upper=F,tinypeak.removal = .9),0)
#   temp <- flowDensity(x2,c(1,3),position = c(T,NA),gates=c(fs.gate, NA))
#   ss.gate <- deGate(temp@flow.frame,3)
#   not.gran <- flowDensity(x2,c(1,3),position = c(T,F),gates=c(fs.gate, ss.gate))
#   plotDens(x2,c(1,3),main="CD45+")
#   lines(not.gran@filter,lty=3,lwd=3)
#   not.gran@filter <- rotate.data(not.gran@filter,chans = c(1,3),theta =-pi/5 )$data
#   not.gran@flow.frame <- rotate.data(getflowFrame(not.gran),chans = c(1,3),theta =-pi/5 )$data
#   plotDens(x,c(1,3),main="CD45+")
#   lines(not.gran@filter,lty=3,lwd=3)
#   return(not.gran)
# })
# fs.pbmc<-lapply(PBMC, function(x) getflowFrame(x))
# fs.pbmc <- as(object=fs.pbmc, Class="flowSet")

#Gating lymphocytes
#----------------------------------------------------------------------------------------
#*****************************************************************************************
#                      Gating Lymph
#*****************************************************************************************
dev.set(2)
par(mfrow=c(3,4))
fs.hi <- fsApply(fs.pbmc,deGate,channels.ind["FS-A"],upper=T,percentile=NA,tinypeak.removal=.9)
Lymph <- fsApply(fs.pbmc, function(x) {
  print(identifier(x))
  x2 <- rotate.data(x,chans = channels.ind[c("FS-A","SS-A")],theta = -pi/6)$data
  peaks<- getPeaks(x2,channels.ind["FS-A"],tinypeak.removal = .05)$Peaks
  fs.gate <-ifelse(test = length(peaks)>2,yes = deGate(x2,channels.ind["FS-A"],tinypeak.removal = .05,all.cuts =T)[2],no=peaks[1]+sd(exprs(x2)[,channels.ind["FS-A"]]))
  fs.hi <- deGate(x2,channels.ind["FS-A"],upper=T,percentile=NA,tinypeak.removal=.9)*1.05
  #fs.p <- deGate(x2,1,use.percentile=T,tinypeak.removal=.5,percentile=.84)
  fs.gate <- min(fs.gate, fs.hi)
  ss.gate <- deGate(x,channels.ind["SS-A"], tinypeak.removal=.9,percentile=NA,upper=T,alpha=.01)
  temp <- flowDensity(x2,channels.ind[c("FS-A","SS-A")],position=c(F,NA), gates=c(fs.gate,NA))
  no.mono <- rotate.data(getflowFrame(temp),chans = channels.ind[c("FS-A","SS-A")],theta = pi/6)$data                    
  lymph <- flowDensity(no.mono, channels.ind[c("FS-A","SS-A")],position=c(NA,F), gates=c(NA,ss.gate*1.2))
  #plotDens(x,channels=channels.ind[c("FS-A","SS-A")],main=identifier(x))#,nbin=3000
  #lines(lymph@filter, col=1, lty=3, lwd=3)
  lymph@proportion <- lymph@cell.count/nrow(x)*100
  return(lymph)
  
})
fs.lymph <-lapply(Lymph, function(x) getflowFrame(x))
fs.lymph <- as(object=fs.lymph, Class="flowSet")
#*****************************************************************************************
#                      Gating Bcells
#*****************************************************************************************
dev.set(2)
par(mfrow=c(3,4))
Bcells <- fsApply(fs.lymph, function(x) {
  print(identifier(x))
#  cd19.gate <- deGate(x, channels.ind["CD19"], tinypeak.removal=.1,use.upper=T,upper=T,percentile=NA,alpha = .001)
  # cd19.gate <- tail(deGate(x, channels.ind["CD19"], tinypeak.removal=1/60,all.cuts=T),1)
  temp <-  flowDensity(x, channels.ind[c("CD19","CD3")],position=c(NA,F))
  cd19.gate <- deGate(temp@flow.frame,channels.ind["CD19"],percentile=NA, upper=T)*1.1
  quad.1 <- flowDensity(x, channels.ind[c("CD19","CD3")],position=c(F,F),gates= c(cd19.gate,NA))
  quad.2 <- flowDensity(x, channels.ind[c("CD19","CD3")],position=c(T,F),gates= c(cd19.gate,NA))#,use.percentile=c(T,T), percentile=c(.999,.999))
  quad.3 <- flowDensity(x, channels.ind[c("CD19","CD3")],position=c(F,T),gates= c(cd19.gate,NA))
  quad.4 <- flowDensity(x, channels.ind[c("CD19","CD3")],position=c(T,T),gates= c(cd19.gate,NA))
 # plotDens(x,channels=channels.ind[c("CD19","CD3")],main=identifier(x))#,nbin=3000
 # lines(quad.3@filter, col=1, lty=3, lwd=3)
  return(list(quad1=quad.1,quad2=quad.2,quad3=quad.3, quad4=quad.4))
})
fs.tcell <- lapply(Bcells, function(x) getflowFrame(x$quad3))
fs.tcell <- as(object=fs.tcell, Class="flowSet")
#*****************************************************************************************
#                      Gating CD4
#*****************************************************************************************
dev.set(2)
par(mfrow=c(3,4))
CD4.8 <- fsApply(fs.tcell, function(x) {
  print(identifier(x))
  cd4.gate <- deGate(x,channels.ind["CD4"])
  temp <- flowDensity(x, channels.ind[c("CD4","CD8")],position=c(T,NA),gates=c(cd4.gate,NA))
  cd8.gate <- deGate(temp@flow.frame,channels.ind["CD8"],use.upper=T,upper=T,percentile=NA,tinypeak.removal = .9)*1.05
  temp.2 <- flowDensity(x, channels.ind[c("CD4","CD8")],position=c(F,NA),gates=c(cd4.gate,NA))
  cd8.gate.lo <- min(deGate(temp.2@flow.frame,channels.ind["CD8"],percentile=.08),deGate(temp.2@flow.frame,9,percentile=NA,tinypeak.removal =.5, upper=F))
  cd4<-flowDensity(x, channels.ind[c("CD4","CD8")],position=c(T,F),gate=c(cd4.gate,cd8.gate))
  dn<-flowDensity(x, channels.ind[c("CD4","CD8")],position=c(F,F),gate=c(cd4.gate*.95,cd8.gate.lo*.98))
  cd8 <- flowDensity(x, channels.ind[c("CD4","CD8")],position=c(F,T),gates=c(cd4.gate,cd8.gate.lo))
 # plotDens(x,channels=channels.ind[c("CD4","CD8")],main=identifier(x),ylim=c(0,5))#,nbin=3000
 # lines(cd4@filter, col=1, lty=3, lwd=3)
 # lines(cd8@filter, col=1, lty=3, lwd=3)
 # lines(dn@filter, col=1, lty=3, lwd=3)
  return(list(cd4=cd4,cd8=cd8, dn=dn))
})
#*****************************************************************************************
#                      Gating CD14
#*****************************************************************************************
cd14.gate <- averageGates(fsApply(fs.pbmc,deGate, channels.ind["CD14"],percentile=NA,tinypeak.removal=.9, upper=T),sd.coeff = 1.5)
names(cd14.gate) <- sampleNames(fs.lymph)
dev.set(2)
par(mfrow=c(3,4))
#It was pi/11 for WBC
CD14 <- fsApply(fs.pbmc, function(x) {
  print(identifier(x))
  x2 <- rotate.data(x,chans = channels.ind[c("CD16","CD14")],theta = pi/14)$data
  cd14.hi <- deGate(x2,channels.ind["CD14"],percentile=NA,tinypeak.removal = .99,upper=T,alpha=.05)
  cd14.gate <- deGate(x2,channels.ind["CD14"],percentile=NA)
  cd14.gate <- min(cd14.gate,cd14.hi)
  cd14 <- flowDensity(x2, channels.ind[c("CD16","CD14")],position=c(NA,T),gates=c(NA,cd14.gate))
  cd14@flow.frame <- rotate.data(getflowFrame(cd14),chans = channels.ind[c("CD16","CD14")],theta = -pi/14)$data
  cd14@filter <- rotate.data(cd14@filter,theta=-pi/14)$data
  #plotDens(x,channels=channels.ind[c("CD16","CD14")],main=identifier(x))#,nbin=3000
  #lines(cd14@filter, col=1, lty=3, lwd=3)
  return(cd14)
})
fs.14 <-lapply(CD14, function(x) getflowFrame(x))
fs.14 <- as(object=fs.14, Class="flowSet")
#*****************************************************************************************
#                      Gating CD14
#*****************************************************************************************
cd16.hi <- fsApply(fs.14,deGate,channels.ind["CD16"],upper=T,percentile=NA,tinypeak.removal=.8,alpha=.08) #.2 for WBC alpha=.01
names(cd16.hi) <- sampleNames(fs.14)
dev.set(2)
par(mfrow=c(5,5))
CD14.hi <- fsApply(fs.14, function(x) {
  print(identifier(x))
  temp <- flowDensity(x, channels.ind[c("CD16","CD14")],position=c(F,NA),gates=c(cd16.hi[identifier(x)],NA))
  cd14.gate <- deGate(temp@flow.frame, channels.ind["CD14"],upper=F,percentile=NA,tinypeak.removal = .99)# max of deGate(temp@flow.frame,  channels.ind["CD14"],upper=F,percentile=NA,tinypeak.removal = .2))
  quad.1 <- flowDensity(x, channels.ind[c("CD16","CD14")],position=c(F,F),gates=c(cd16.hi[identifier(x)],cd14.gate))
  quad.2 <- flowDensity(x, channels.ind[c("CD16","CD14")],position=c(T,F),gates=quad.1@gates)
  quad.3 <- flowDensity(x, channels.ind[c("CD16","CD14")],position=c(F,T),gates=quad.1@gates)
  quad.4 <- flowDensity(x, channels.ind[c("CD16","CD14")],position=c(T,T),gates=quad.1@gates)
 # plotDens(x,channels=channels.ind[c("CD16","CD14")],main=identifier(x))#,nbin=3000
 # abline(v=quad.1@gates[1],h=quad.1@gates[2])
  return(list(quad1=quad.1,quad2=quad.2,quad3=quad.3, quad4=quad.4))
})

#*****************************************************************************************
#                      Gating CD56
#*****************************************************************************************
cd3.gate <- fsApply(fs.lymph,deGate,channels.ind["CD3"])
names(cd3.gate) <- sampleNames(fs.lymph)
dev.set(2)
par(mfrow=c(3,4))
CD56 <- fsApply(fs.lymph, function(x) {
  print(identifier(x))
  temp <- flowDensity(x,channels=channels.ind[c('CD56','CD3')],position=c(NA,F),gates=c(NA, cd3.gate[identifier(x)]))
  cd56.gate <- deGate(temp@flow.frame,channels.ind['CD56'],percentile=NA)
  temp <- flowDensity(x,channels=channels.ind[c('CD56','CD3')],position=c(F,T), gates=c(cd56.gate, cd3.gate[identifier(x)]))
  cd56.gate.lo <- deGate(temp@flow.frame,channels.ind['CD56'],upper=T,percentile=NA,tinypeak.removal = .99,alpha=.01)
  quad.2 <- flowDensity(x,channels=channels.ind[c('CD56','CD3')],position=c(T,F), gates=c(cd56.gate, cd3.gate[identifier(x)]))
  quad.4 <- flowDensity(x,channels=channels.ind[c('CD56','CD3')],position=c(T,T), gates=c(cd56.gate.lo, cd3.gate[identifier(x)]))
  cd56.gate.hi <- deGate(quad.4@flow.frame,channels.ind['CD56'],percentile =.99995,use.percentile=T)
  #plotDens(x,channels=channels.ind[c('CD56','CD3')],main="Lymphocytes")#,nbin=3000
 # abline(v=cd56.gate.hi,h=quad.2@gates[2])
  #lines(quad.2@filter,type="l")
  #lines(quad.4@filter,type="l")
  print(cd56.gate.hi)
  return(list(nk=quad.2, nkt=quad.4,gate=cd56.gate.hi))
})
fs.nk <-lapply(CD56, function(x) getflowFrame(x$nk))
fs.nk <- as(object=fs.nk, Class="flowSet")
#*****************************************************************************************
#                      Gating CD16/56
#*****************************************************************************************
#cd16.gate <- averageGates(fsApply(fs.nk,deGate,channels.ind["CD16"],upper=F,percentile=NA,tinypeak.removal=.2),sd.coeff=1)*.8
#names(cd16.gate)   <-sampleNames(fs.nk)
dev.set(2)
par(mfrow=c(3,4))
CD56.nk <- fsApply(fs.nk, function(x) {
  print(identifier(x))
  #x2 <- rotate.data(x,channels.ind[c("CD56","CD16")],theta = pi/18)$data
  cd56.upper<- deGate(x,channels.ind["CD56"],percentile=NA,upper=T,use.upper=T,tinypeak.removal = .5,alpha=.2)
  #It was tail
  cd56.gate <- getPeaks(x, channels.ind["CD56"],tinypeak.removal = .1)$Peaks[1]+1.5*sd(exprs(x)[,channels.ind["CD56"]])
  cd56.gate <- min(cd56.gate,cd56.upper)
  cd56.hi <-  flowDensity(x,channels=channels.ind[c("CD56","CD16")],position=c(T,NA),gates=c(cd56.gate,NA))
  #cd56.hi@flow.frame <- rotate.data(getflowFrame(cd56.hi),chans = channels.ind[c("CD56","CD16")],theta = -pi/18)$data
  #cd56.hi@filter <- rotate.data(cd56.hi@filter,theta = -pi/18)$data
  #cd56.hi@proportion <- cd56.hi@cell.count/nrow(x)*100
  #plotDens(x,channels=channels.ind[c("CD56","CD16")],main=identifier(x))#,nbin=3000
 # plot(density(exprs(x)[,5]))
 #abline(v=c(cd56.upper,cd56.gate),col=c(1,3))
 # lines(cd56.hi@filter,type="l",lwd=2)
  return(cd56.hi)
})
fs.56 <- lapply(CD56.nk, function(x) getflowFrame(x))
fs.56 <- as(object=fs.56, Class="flowSet")
#*****************************************************************************************
#                      Gating CD16/64
#*****************************************************************************************
cd64.gate <- fsApply(fs.14,deGate,channels.ind["CD64"],tinypeak.removal=.5,percentile=NA,upper=F)
names(cd64.gate)  <- sampleNames(fs.14)
dev.set(2)
par(mfrow=c(3,4))
CD64.14 <- fsApply(fs.14, function(x) {
  print(identifier(x))
  cd64 <- flowDensity(x,channels=channels.ind[c("CD64","CD16")],position=c(T,T),gates=c(cd64.gate[identifier(x)],cd16.hi[identifier(x)]))
  #plotDens(x,channels=channels.ind[c("CD64","CD16")],main=identifier(x))#,nbin=3000
  #abline(v=cd64@gates[1],h=cd64@gates[2])
  return(cd64)
})
#****************************************************************************************
#****************************************************************************************
dev.set(2)
par(mfrow=c(3,4))
CD56.16 <- fsApply(fs.56, function(x) {
  cd16.neg <- flowDensity(x,channels=channels.ind[c("CD56","CD16")],position=c(NA,F),gates=c(NA,cd16.hi[identifier(x)]))
  #plotDens(x,channels=channels.ind[c("CD56","CD16")],main=identifier(x),cex=3,col=4)#,nbin=3000
  #abline(h=cd16.neg@gates[2])

  return(cd16.neg)
})

#*****************************************************************************************
#                     Plotting
#*****************************************************************************************

cell.prop <- c()
for ( i in 1:length(fs))
{
  png(paste(res,"/",identifier(fs[[i]]),".png",sep=""),width=1700,height=1300,pointsize=23)
  par(mfrow=c(3,4), mar=c(4,4,2,2),oma=c(1,1,4,1))
  plotDens(fs[[i]],channels=channels.ind[c("FS-W","FS-A")],main="All")#,nbin=3000
  lines(Singlet[[i]]@filter, col=1, lty=3, lwd=3)
  legend("bottomright","Singlets",cex=.8)
  plotDens(fs.sngl[[i]],channels=channels.ind[c("CD45","SS-A")],main="Singlets")#,nbin=3000
  lines(CD45[[i]]@filter, col=1, lty=3, lwd=3)
  legend("bottomright","CD45+",cex=.8)
  plotDens(fs.45[[i]],channels=channels.ind[c("CD14","SS-A")],main="CD45+")#,nbin=3000
  lines(PBMC[[i]]@filter, col=1, lty=3, lwd=3)
  legend("bottomright","PBMCs",cex=.8)
  plotDens(fs.pbmc[[i]],channels=channels.ind[c("FS-A","SS-A")],main="PBMCs")#,nbin=3000
  lines(Lymph[[i]]@filter, col=1, lty=3, lwd=3)
  legend("topleft","Lymphocytes",cex=.8)
  plotDens(fs.lymph[[i]],channels=channels.ind[c("CD19","CD3")],main="Lymphocytes")#,nbin=3000
  abline(v=Bcells[[i]]$quad1@gates[1],h=Bcells[[i]]$quad1@gates[2])
  legend("topleft","Tcells",cex=.8)
  legend("bottomright","Bcells",cex=.8)
  plotDens(fs.tcell[[i]],channels=channels.ind[c("CD4","CD8")],main="Tcells")#,nbin=3000
  lines(CD4.8[[i]]$cd4@filter, col=1, lty=3, lwd=3)
  lines(CD4.8[[i]]$cd8@filter, col=1, lty=3, lwd=3)
  lines(CD4.8[[i]]$dn@filter, col=1, lty=3, lwd=3)
  legend(x = CD4.8[[i]]$cd8@gates[1],y=0,"CD4",cex=.8)
  legend(x = CD4.8[[i]]$cd8@gates[1],y=3,"CD8",cex=.8)
  legend("bottomlef","DN Tcells",cex=.8)
  plotDens(fs.pbmc[[i]],channels=channels.ind[c("CD16","CD14")],main="PBMCs")#,nbin=3000
  lines(CD14[[i]]@filter, col=1, lty=3, lwd=3)
  legend("topleft","CD14+",cex=.8)
  plotDens(fs.14[[i]],channels=channels.ind[c("CD16","CD14")],main="CD 14+")#,nbin=3000
  abline(v=CD14.hi[[i]]$quad1@gates[1],h=CD14.hi[[i]]$quad1@gates[2])
  legend("bottomleft","CD14+CD16-",cex=.8)
  legend("bottomright","CD14+CD16+",cex=.8)
  legend("topleft","CD14++CD16-",cex=.8)
  legend("topright","CD14++CD16+",cex=.8)
  plotDens(fs.lymph[[i]],channels=channels.ind[c("CD56","CD3")],main="Lymphocytes")#,nbin=3000
  lines(CD56[[i]]$nk@filter, col=1, lty=3, lwd=3)
  lines(CD56[[i]]$nkt@filter, col=1, lty=3, lwd=3)
  legend("bottomright","NK",cex=.8)
  legend("topright","NKT",cex=.8)
  plotDens(fs.nk[[i]],channels=channels.ind[c("CD56","CD16")],main="NK")#,nbin=3000
  lines(CD56.nk[[i]]@filter,col=1,lty=3,lwd=3)
  legend("topright","CD56Bright",cex=.8)
  plotDens(fs.14[[i]],channels=channels.ind[c("CD16","CD64")],main="CD 14+")#,nbin=3000
  abline(v=CD64.14[[i]]@gates[1],h=CD64.14[[i]]@gates[2])
  legend("topright","CD64++CD16+",cex=.8)
  plotDens(fs.56[[i]],channels=channels.ind[c("CD56","CD16")],main="CD56 Bright",col=4,cex=3)#,nbin=3000
  abline(h=CD56.16[[i]]@gates[2])
  legend("bottomright","CD56Bright CD16-",cex=.8)
  mtext(outer = T, side=3,text=sampleNames(fs)[i],cex=.8,font=4)
  dev.off()
  cell.prop <- rbind(cell.prop, c(Lymph[[i]]@proportion,
                                  Bcells[[i]]$quad2@proportion,Bcells[[i]]$quad3@proportion,
                                  CD4.8[[i]]$dn@proportion,CD4.8[[i]]$cd4@proportion,CD4.8[[i]]$cd8@proportion,
                                  CD14[[i]]@proportion, CD14.hi[[i]]$quad2@proportion,CD14.hi[[i]]$quad3@proportion,
                                  CD14.hi[[i]]$quad4@proportion,CD56[[i]]$nk@proportion,CD56[[i]]$nkt@proportion,
                                  CD56.nk[[i]]@proportion, CD64.14[[i]]@proportion,CD56.16[[i]]@proportion))
}
colnames(cell.prop) <- c("Lymphocytes","Bcells","Tcells","DN T cells","CD4+ T cells","CD8+ T cells","CD14+",
                         "CD14+CD16+","CD14++CD16-","CD14++CD16+","NK cells","NKT cells","CD56bright NK cells","CD64++CD16+","CD56bright CD16-")
rownames(cell.prop) <- sampleNames(fs)
write.csv(cell.prop,paste(res,"Proportions.csv",sep="/"))
#proc.time()-ptm