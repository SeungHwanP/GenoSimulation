system("mkdir -p /storage0/lab/psh/simulation_n10K_p10K_continuous/5.plink")
setwd("/storage0/lab/psh/simulation_n10K_p10K_continuous/5.plink")

rm(list=ls())
library(data.table)

# make function to make ped-map file
pid=NULL;mid=NULL;sex=NULL;phe=NULL
make_ped <- function(name,pid=NULL,mid=NULL,sex=NULL,phe=NULL){
  setwd(mldose_mlinfo)
  print("Read file")
  dose <- fread(paste0("ukb_imp_",name,".mldose"));n <- nrow(dose)
  info <- fread(paste0("ukb_imp_",name,".mlinfo"))
  fid <- dose[,1];iid <- dose[,2]
  dose <- as.matrix(dose)
  
  print("make file")
  allele1 <- round(dose[,-c(1:2)])
  allele2 <- 2-allele1
  tmp <- lapply(1:nrow(allele1),
                function(x) rep(as.vector(t(as.matrix(info[,c("Al1","Al2")]))),as.vector(t(matrix(c(allele1[x,],allele2[x,]),ncol=2)))))
  rm(allele1,allele2,info)
  tmp <- unlist(tmp)
  tmp <- matrix(tmp,nrow=n,byrow=T)
  
  if(is.null(pid)){ pid <- 0 }
  if(is.null(mid)){ mid <- 0 }
  if(is.null(sex)){ sex <- 1 }
  if(is.null(phe)){ phe <- 1 }
  tmp <- as.data.frame(cbind(fid,iid,pid,mid,sex,phe,tmp))
  tmp[1:10,1:10]
  
  print("save file")
  setwd(saving_point)
  fwrite(tmp,paste0("ukb_imp_",name,".ped"),quote=F,row.names=F,col.names = F)
}

make_map <- function(name){
  setwd(dose_info)
  info = as.matrix(fread(paste0("ukb_imp_",name,".info")))
  setwd(mldose_mlinfo)
  info2 <- fread(paste0("ukb_imp_",name,".mlinfo"))
  chr_bp <- matrix(as.numeric(info[,c(1,4)]),ncol=2)
  
  tmp <- floor(chr_bp[,2]/1000000)
  chr_bp <- data.frame(cbind(chr_bp[,1],info[,2],tmp,chr_bp[,2]))
  colnames(chr_bp) <- NULL
  
  setwd(saving_point)
  fwrite(chr_bp,paste0("ukb_imp_",name,".map"),quote=F,row.names=F,col.names = F)
}

# run function to make ped-map file
dose_info <- "/storage0/lab/psh/simulation_n10K_p10K_continuous/3.dose_info/"
mldose_mlinfo <- "/storage0/lab/psh/simulation_n10K_p10K_continuous/4.mldose_mlinfo"
saving_point <- "/storage0/lab/psh/simulation_n10K_p10K_continuous/5.plink"

for(case in c("training","testing")){
  make_map(case)
  make_ped(case)
}
