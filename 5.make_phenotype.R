rm(list=ls())
library(glue)
library(MASS)
library(data.table)
nsnp <- 100
h2_1 <- 0.6;h2_2 <- 0.3
crr <- 0.4
SNP <- fread("/storage/lab/psh/simulation_n10K_p10K_continous/4.plink/train.bim")

for(seed in 1:100){
  system(paste0("mkdir -p /storage/lab/psh/simulation_n10K_p10K_continous/5.phenotype/",seed))
  for(nsnp in c(0.01,seq(0.1,1,0.1))*10000){
    h2= c(h2_1,h2_2); dd = h2 / nsnp
    rho12 = crr * sqrt(dd[1]) * sqrt(dd[2])
    beta <- data.frame(SNP=SNP$V2,A1=SNP$V5,beta1=0,beat2=0)
    set.seed(seed)
    num <- sort(sample(nrow(SNP),nsnp))
    covar = diag(dd);   covar[2,1] = covar[1,2] = rho12
    beta[num,3:4] = mvrnorm(n=nsnp,rep(0,2),covar)
    fwrite(beta,paste0("/storage/lab/psh/simulation_n10K_p10K_continous/5.phenotype/",seed,"/1.beta_",nsnp,".txt"),row.names=F,quote=F,sep="\t",col.names=T)
    fwrite(beta[beta$beta1 !=0,]$SNP,paste0("/storage/lab/psh/simulation_n10K_p10K_continous/5.phenotype/",seed,"/0.SNPlist_",nsnp,".txt"),row.names=F,quote=F,sep="\t",col.names=F)
    WT <- paste0("/storage/lab/psh/simulation_n10K_p10K_continous/5.phenotype/",seed,"/1.beta_",nsnp,".txt")
    N0 <- paste0("/storage/lab/psh/simulation_n10K_p10K_continous/5.phenotype/",seed,"/0.SNPlist_",nsnp,".txt")
    #xbeta-train
    out <- paste0("/storage/lab/psh/simulation_n10K_p10K_continous/5.phenotype/",seed,"/0.phe_wo_error_p",nsnp,"_train")
    code <- glue("/storage/program/plink_linux_x86_64_20200616/plink1.9 --bfile /storage/lab/psh/simulation_n10K_p10K_continous/4.plink/train --extract {N0} --recode A --a1-allele {WT} 2 1 --out {out}")
    system(code)
    
    dose <- as.data.frame(fread(glue("{out}.raw"))[,-c(3:6)])
    xbeta <- data.frame(dose[,1:2])
    dose_col <- colnames(dose)[-c(1:2)]
    dose_col <- unlist(lapply(1:length(dose_col),function(x) paste(strsplit(dose_col,"_",fixed=T)[[x]][-4],collapse = "_")))
    dose_col <- data.frame(SNP=dose_col)
    sel_col <- merge(dose_col,beta,by="SNP",sort=F,all.x=T)
    dose <- as.matrix(scale(dose[,paste0(sel_col$SNP,"_",sel_col$A1)]))
    xbeta <- cbind(xbeta,dose %*% as.matrix(sel_col[,-c(1:2)]))
    
    #xbeta-test
    out <- paste0("/storage/lab/psh/simulation_n10K_p10K_continous/5.phenotype/",seed,"/0.phe_wo_error_p",nsnp,"_test")
    code <- glue("/storage/program/plink_linux_x86_64_20200616/plink1.9 --bfile /storage/lab/psh/simulation_n10K_p10K_continous/4.plink/test --extract {N0} --recode A --a1-allele {WT} 2 1 --out {out}")
    system(code)
    
    dose <- as.data.frame(fread(glue("{out}.raw"))[,-c(3:6)])
    xbeta_test <- data.frame(dose[,1:2])
    dose_col <- colnames(dose)[-c(1:2)]
    dose_col <- unlist(lapply(1:length(dose_col),function(x) paste(strsplit(dose_col,"_",fixed=T)[[x]][-4],collapse = "_")))
    dose_col <- data.frame(SNP=dose_col)
    sel_col <- merge(dose_col,beta,by="SNP",sort=F,all.x=T)
    dose <- as.matrix(scale(dose[,paste0(sel_col$SNP,"_",sel_col$A1)]))
    xbeta_test <- cbind(xbeta_test,dose %*% as.matrix(sel_col[,-c(1:2)]))
    colnames(xbeta) <- colnames(xbeta_test) <- c("FID","IID","phe1","phe2")
    y <- xbeta;y_test <- xbeta_test
    
    err.vec = NULL
    set.seed(seed)
    k <- 1
    for (k in 1:2) err.vec = c(err.vec, rnorm(nrow(xbeta),sd=sqrt((1-h2[k]))))
    err = matrix(err.vec, ncol=2)
    y[,3:4] <- y[,3:4]+err
    y_test[,3:4] <- y_test[,3:4]+err
    
    cat("crr =",crr,",",cor(beta[,3],beta[,4]),"\n") # genetic correlation
    cat("h1 =",h2[1],",",var(xbeta[,3])/var(y[,3]),"\n") # heritability for trait 1
    cat("h2 =",h2[2],",",var(xbeta[,4])/var(y[,4]),"\n") # heritability for trait 2
    
    #test
    cat("h1 =",h2[1],",",var(xbeta_test[,3])/var(y_test[,3]),"\n") # heritability for trait 1
    cat("h2 =",h2[2],",",var(xbeta_test[,4])/var(y_test[,4]),"\n") # heritability for trait 2
    
    fwrite(NULL,paste0(paste0("/storage/lab/psh/simulation_n10K_p10K_continous/5.phenotype/",seed,"/1..p",nsnp,"-h1 =",h2[1],",",round(var(xbeta[,3])/var(y[,3]),4)),paste0(", h2 =",h2[2],",",round(var(xbeta[,4])/var(y[,4]),4)),"-1.train"))
    fwrite(NULL,paste0(paste0("/storage/lab/psh/simulation_n10K_p10K_continous/5.phenotype/",seed,"/1..p",nsnp,"-h1 =",h2[1],",",round(var(xbeta_test[,3])/var(y_test[,3]),4)),paste0(", h2 =",h2[2],",",round(var(xbeta_test[,4])/var(y_test[,4]),4)),"-2.test"))
    
    fwrite(y, paste0("/storage/lab/psh/simulation_n10K_p10K_continous/5.phenotype/",seed,"/1.pheno_p",nsnp,"_train.phe"), col.names=T, row.names=F, quote=F, sep="\t")
    fwrite(y_test, paste0("/storage/lab/psh/simulation_n10K_p10K_continous/5.phenotype/",seed,"/1.pheno_p",nsnp,"_test.phe"), col.names=T, row.names=F, quote=F, sep="\t")
  }
}