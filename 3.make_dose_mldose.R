rm(list=ls())
library(data.table)

for(case in c("train","test")){
  info <- dose <- mlinfo <- NULL
  for(chr in 1:22){
    print(paste0(case,", CHR=",chr))
    print("Load")
    gen <- fread(paste0("/storage/lab/psh/simulation_n10K_p933K/3.gen/",case,"_Chr",chr,".gen"))
    res <- gen[,c(1:6)]
    dos1 <- gen[,-c(1:6)];rm(gen)
    
    print("Make1")
    n <- ncol(dos1)/3
    dos1 <- as.matrix(dos1)
    dos1 = Reduce(rbind,lapply(1:n,function(x) 3-apply(dos1[,c((x*3-2):(x*3)) ],1,which.max)))
    
    info <- rbind(info,res)
    dose <- cbind(dose,dos1)
    
    print("Make2")
    info_new = res[,c(2,5,6,4,4)]
    colnames(info_new) = c("SNP","Al1","Al2","Freq1","MAF")
    info_new[,4] = apply(dos1,2,mean)/2
    info_new[,5] = 1-as.numeric(info_new$Freq1)
    info_new$Quality = 0.8; info_new$Rsq = 0.8
    
    mlinfo <- rbind(mlinfo,info_new)
  }
  id <- fread(paste0("/storage/lab/psh/simulation_n10K_p933K/3.gen/",case,"_Chr1.sample"))[-1,1:2];length(unique(id$ID_1))
  mldose <- cbind(id,dose)

  fwrite(info,paste0("/storage/lab/psh/simulation_n10K_p10K_continuous/3.dose/",case,".info"),row.names=F, col.names=F,quote=F,sep="\t")
  fwrite(dose,paste0("/storage/lab/psh/simulation_n10K_p10K_continuous/3.dose/",case,".dose"),row.names=F, col.names=F,quote=F,sep="\t")
  
  fwrite(mlinfo, paste0("/storage/lab/psh/simulation_n10K_p10K_continuous/3.dose/",case,".mlinfo"), col.names=T, row.names=F, quote=F, sep="\t")
  fwrite(mldose, paste0("/storage/lab/psh/simulation_n10K_p10K_continuous/3.dose/",case,".mldose"), col.names=F, row.names=F, quote=F, sep="\t")  
}
