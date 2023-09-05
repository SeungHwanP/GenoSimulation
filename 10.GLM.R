rm(list=ls())
library(glue)
library(MASS)
library(glmnet)
library(data.table)

dose <- fread("/storage/lab/psh/simulation_n10K_p10K_continous/4.plink/train_raw.raw")[,-c(2:6)]

nsnps <- c(0.01,seq(0.1,1,0.1))*10000
ps <- c(0.01,0.001,0.0001)
for(seed in 1:100){
 print(paste0("Thresholding pvalue < ",ps[p]))
 res <- list()
 for(nsnp in 1:11){
   print(Sys.time())
   print(nsnps[nsnp])
   pheno <- fread(paste0("/storage/lab/psh/simulation_n10K_p10K_continous/5.phenotype/",seed,"/1.pheno_p",nsnps[nsnp],"_train.phe"))[,-2]
   pheno <- as.data.frame(merge(dose[,"FID"],pheno,by="FID",sort=F))
   res_tmp <- list()
   for(phe in 1:2){
     print(phe)
     gwas <- fread(glue("zcat /storage/lab/psh/simulation_n10K_p10K_continous/7.bolt/{seed}/Bolt_p{nsnps[nsnp]}_{phe}.bgen.stats.gz"))
     print("lasso")
     cv1 <- cv.glmnet(x=as.matrix(data[,-1],),y=pheno[,(phe+1)],alpha=1,type.measure = "mse")
     print("EN0.5")
     cv2 <- cv.glmnet(x=as.matrix(data[,-1],),y=pheno[,(phe+1)],alpha=0.5,type.measure = "mse")
     print("EN0.1")
     cv3 <- cv.glmnet(x=as.matrix(data[,-1],),y=pheno[,(phe+1)],alpha=0.1,type.measure = "mse")
     print("ridge")
     cv4 <- cv.glmnet(x=as.matrix(data[,-1],),y=pheno[,(phe+1)],alpha=0,type.measure = "mse")
     res_tmp[[phe]] <- list(lasso=cv1,'EN0.5'=cv2,'EN0.1'=cv3,ridge=cv4);rm(cv1,cv2,cv3,cv4)
   }
   res[[nsnp]] <- list(phe1=res_tmp[[1]],phe2=res_tmp[[2]])
   save(res,file=paste0("/storage/lab/psh/simulation_n10K_p10K_continous/10.GLM/",seed,"_p",ps[p],"_tmp.Rdata"))
 }
 system(paste0("mv /storage/lab/psh/simulation_n10K_p10K_continous/10.GLM/",seed,"_p",ps[p],"_tmp.Rdata /storage/lab/psh/simulation_n10K_p10K_continous/10.GLM/",seed,"_p",ps[p],".Rdata"))
}