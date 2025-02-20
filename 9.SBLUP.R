system("mkdir -p /storage0/lab/psh/simulation_n10K_p10K_continuous/9.SBLUP")
#step1 : make ma
rm(list=ls())
library(data.table)
seed <- 1
nsnp <- 1
phe <- 1
for(seed in 1:100){
  system(paste0("mkdir -p /storage0/lab/psh/simulation_n10K_p10K_continuous/9.SBLUP/",seed,"/ma"))
  for(nsnp in c(1,seq(10,100,10))){
    for(phe in 1:2){
      bolt <- fread(paste0("zcat /storage0/lab/psh/simulation_n10K_p10K_continuous/7.bolt/",seed,"/Bolt_p",nsnp,"00_",phe,".bgen.stats.gz"))
      
      res <- data.table(SNP=bolt$SNP,A1=bolt$ALLELE1,A0=bolt$ALLELE0,freq=bolt$A1FREQ,
                        b=bolt$BETA,se=bolt$SE,p=bolt$P_BOLT_LMM,N=10000)
      
      fwrite(res,paste0("/storage0/lab/psh/simulation_n10K_p10K_continuous/9.SBLUP/",seed,"/ma/p",nsnp,"00_",phe,".ma"),row.names=F,quote=F,sep="\t",col.names=T)
    }
  }
}


#step2 : run SBLUP
rm(list=ls())
library(data.table)
library(glue)
seed <- 1
nsnp <- 1
phe <- 1
prev <- 0.2
for(seed in 1:50){
  for(nsnp in c(1,seq(10,100,10))){
    for(phe in 2){
      h2 <- system(glue(paste0("cat /storage0/lab/psh/simulation_n10K_p10K_continuous/8.GRM,GC/{seed}/GC_p{nsnp}00_{phe}.hsq |grep 'V(G)/Vp'")),intern=T)
      h2 <- as.numeric(strsplit(h2,"\t")[[1]][2])
      h2 <- round(((1/h2)-1)*10000)
      sblup <- glue("/storage0/program/gcta_1.93.2beta/gcta64 --bfile /storage0/lab/psh/simulation_n10K_p10K_continuous/5.plink/train --cojo-file /storage0/lab/psh/simulation_n10K_p10K_continuous/9.SBLUP/{seed}/ma/p{nsnp}00_{phe}.ma --cojo-sblup {h2} --cojo-wind 1000 --out /storage0/lab/psh/simulation_n10K_p10K_continuous/9.SBLUP/{seed}/SBLUP_p{nsnp}00_{phe}")
      system(sblup)
    }
  }
}
