# these are `rohproc` R pacakge functions (see original package on GitHub: https://github.com/CeballosGene/rohproc)

# FUNCTIONS ROHi
roh_island<-function(pop,chr,p1,p2){
  names(pop)<-tolower(names(pop))
  a<-pop[pop$chr==chr,]
  island<-subset(a,pos1<=p1 & pos2>=p2)
  n<-length(unique(island$iid))/length(unique(pop$iid))
  return(n)
}
poisson.roh_island<-function(pop,chr,p1,p2){
  names(pop)<-tolower(names(pop))
  a<-pop[pop$chr==chr,]
  island<-subset(a,pos1<=p1 & pos2>=p2)
  n<-length(unique(island$iid))
  return(n)
}

# Islands of ROH
get_RHOi<-function(POP,ChroNumber,population){
  SizeWindow=10000
  if(ChroNumber==1){lenChro=250000000}
  if(ChroNumber==2){lenChro=250000000}
  if(ChroNumber==3){lenChro=200000000}
  if(ChroNumber==4){lenChro=191000000}
  if(ChroNumber==5){lenChro=182000000}
  if(ChroNumber==6){lenChro=171000000}
  if(ChroNumber==7){lenChro=160000000}
  if(ChroNumber==8){lenChro=146000000}
  if(ChroNumber==9){lenChro=139000000}
  if(ChroNumber==10){lenChro=133900000}
  if(ChroNumber==11){lenChro=136000000}
  if(ChroNumber==12){lenChro=134000000}
  if(ChroNumber==13){lenChro=115000000}
  if(ChroNumber==14){lenChro=108000000}
  if(ChroNumber==15){lenChro=102000000}
  if(ChroNumber==16){lenChro=91000000}
  if(ChroNumber==17){lenChro=84000000}
  if(ChroNumber==18){lenChro=81000000}
  if(ChroNumber==19){lenChro=59000000}
  if(ChroNumber==20){lenChro=64000000}
  if(ChroNumber==21){lenChro=49000000}
  if(ChroNumber==22){lenChro=52000000}
  data.n = apply(data.frame(seq(0,lenChro-SizeWindow,SizeWindow), seq(SizeWindow,lenChro,SizeWindow)), MARGIN=1,
                 function(x,y,z,a) poisson.roh_island(POP, ChroNumber, x[1], x[2]))
  data.p = apply(data.frame(seq(0,lenChro-SizeWindow,SizeWindow), seq(SizeWindow,lenChro,SizeWindow)), MARGIN=1,
                 function(x,y,z,a) roh_island(POP, ChroNumber, x[1], x[2]))
  av<-mean(data.n)
  pval<-ppois(data.n,lambda=av,lower=FALSE)
  x<-seq(1:round(lenChro/SizeWindow))
  pos<-(x*SizeWindow)
  prop<-data.p*100
  data<-data.frame(cbind(x,data.n,prop,pos,pval))
  data.p<-subset(data,data$pval<=0.05/(lenChro/SizeWindow))
  data.re<-data.p |> dplyr::group_by(new=cumsum(c(1,diff(x)!=1))) |>
    dplyr::summarise(pos1=min(pos),pos2=max(pos),n.ind=mean(data.n),per.ind=mean(prop))
  Chr<-rep(ChroNumber,length(data.re$pos1))
  ROHi<-data.frame(cbind(Chr,data.re))
  ROHi<-dplyr::mutate(ROHi,len=(pos2-pos1)/1000000)
  ROHi<-dplyr::mutate(ROHi,pop=rep(population,length(ROHi$Chr)))
  ROHi<-dplyr::select(ROHi,Chr,pos1,pos2,len,n.ind,per.ind,pop)
  colnames(ROHi)<-c("Chr","Start","End","Length","N_Individuals","P_Individuals","Population")
  return(ROHi)
}
# Summarize by population
rohi_summ_pop<-function(mypath){
  files_list <- list.files(path=mypath, full.names=TRUE)
  dat <- data.frame()
  for (i in 1:length(files_list)) {
    dat <- rbind(dat, read.csv((files_list[i]),header=TRUE))
  }
  dat<-dat|>
    dplyr::group_by(Population)|>
    dplyr::summarise(Number=length(Chr),
                     mean_length=mean(Length),
                     sd_length=sd(Length),
                     median_Length=median(Length),
                     iqr_Length=IQR(Length),
                     max_Length=max(Length),
                     mean_N_Individuals=mean(N_Individuals),
                     sd_N_Individuals=sd(N_Individuals),
                     median_N_Individuals=median(N_Individuals),
                     iqr_N_Individuals=IQR(N_Individuals),
                     max_N_Individuals=max(N_Individuals),
                     mean_P_Individuals=mean(P_Individuals),
                     sd_P_Individuals=sd(P_Individuals),
                     median_P_Individuals=median(P_Individuals),
                     iqr_P_Individuals=IQR(P_Individuals),
                     max_P_Individuals=max(P_Individuals))
  out<-as.data.frame(dat)
  return(out)
}