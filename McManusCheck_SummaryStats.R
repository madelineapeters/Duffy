sum.stat.raw = function(){
  
  library(vcfR)
  library(pegas)
  library(dplyr)
  library(tidyr)
  
  vcf = read.vcfR(paste("~/SLiM/testSampleVCF",x,".txt",sep=""), verbose = FALSE )
  poly.list = as.numeric(vcf@fix[,2])
  
  DNA.in = vcfR2DNAbin(vcf,extract.indels=FALSE,verbose=FALSE)
  
  haplos = haplotype(DNA.in)
  haploNum = dim(haplotype(DNA.in))[1]
  FSF = site.spectrum(DNA.in)
  
  SegSites = sum(FSF[FSF>0])
  
  pi.fun = function(){
    
    poly.list = as.numeric(vcf@fix[,2])
    tidy.DNA = extract_gt_tidy(vcf,verbose=FALSE)
    tidy.DNA$Key = poly.list
    
    tidy.DNA$POS1 = sapply(1:nrow(tidy.DNA),FUN=function(x){
      strsplit(tidy.DNA$gt_GT[x],"|")[[1]][1]
    })
    tidy.DNA$POS2 = sapply(1:nrow(tidy.DNA),FUN=function(x){
      strsplit(tidy.DNA$gt_GT[x],"|")[[1]][3]
    })
    
    Haplo1.mat = select(tidy.DNA,-gt_GT,-gt_GT_alleles,-POS2) 
    Haplo1.mat = matrix(Haplo1.mat$POS1,ncol=length(poly.list),nrow=101,byrow=TRUE) 
    rownames(Haplo1.mat) = paste(1:nrow(Haplo1.mat),"a",sep="")
    
    Haplo2.mat = select(tidy.DNA,-gt_GT,-gt_GT_alleles,-POS1) 
    Haplo2.mat = matrix(Haplo2.mat$POS2,ncol=length(poly.list),nrow=101,byrow=TRUE) 
    rownames(Haplo2.mat) = paste(1:nrow(Haplo2.mat),"b",sep="")
    
    Haplo.mat = rbind(Haplo1.mat,Haplo2.mat)
    
    Haplo.df = as.data.frame(Haplo.mat)
    Haplo.unite = unite(as.data.frame(Haplo.mat),"Haplo",1:ncol(Haplo.mat),sep="")
    Haplo.tab = (table(Haplo.unite$Haplo))
    Haplo.seq = names(Haplo.tab)
    
    n = nrow(Haplo.df)
    
    pi.list = c()
    
    for (i in 1:(length(Haplo.seq)-1)){
      for (j in (i+1):length(Haplo.seq)){
        
        seq1 = c(strsplit(as.character(Haplo.seq[i]),split=""))
        seq2 = c(strsplit(as.character(Haplo.seq[j]),split=""))
        
        numDiff = sum(abs(c(sapply(seq1,FUN=as.numeric)) - c(sapply(seq2,FUN=as.numeric))))
        
        freq1 = as.numeric(Haplo.tab[i])/n
        freq2 = as.numeric(Haplo.tab[j])/n
        
        pi.part = freq1*freq2*numDiff    
        
        pi.list = append(pi.list,pi.part)
        
      }
      
    }
    
    pi = sum(pi.list)

    return(pi)
    
  }
  pi = pi.fun()
  
  TajD = tajima.test(DNA.in)$D
  
  Singletons = site.spectrum(DNA.in)[1]
  Doubletons = site.spectrum(DNA.in)[2]
  Fixed = 5001 - sum(site.spectrum(DNA.in))
  
  #Assume mutation lost
  Frequency = 0
  
  fixed.txt = read.delim(paste("~/SLiM/fixedMutations",x,".txt",sep=""))
  n.mut = nrow(fixed.txt)
  POS.list = vector()
  for (r in 2:n.mut){
    txt.row = as.character(fixed.txt[r,1])
    POS.list = append(POS.list,as.integer(strsplit(txt.row," ")[[1]][4]))
  }
  
  if (3501 %in% POS.list){
    Frequency = 1
  } 
  if (Frequency != 1){
    if ("3501" %in% poly.list){
      POS = which(poly.list == "3501")
      temp = strsplit(vcf@fix[POS,8],";")$INFO[7]
      temp = strsplit(temp,split=NULL)[[1]][4]
      temp = as.integer(temp)
      Frequency = 2/2*(ncol(vcf@gt)-1)
    }  
  }
  
  #Homozygosity statistics (EHH, H statistics)
  Homozygosity.fun = function(){
    poly.list = as.numeric(vcf@fix[,2])
    tidy.DNA = extract_gt_tidy(vcf,verbose=FALSE)
    tidy.DNA$Key = poly.list
    
    tidy.DNA$POS1 = sapply(1:nrow(tidy.DNA),FUN=function(x){
      strsplit(tidy.DNA$gt_GT[x],"|")[[1]][1]
    })
    tidy.DNA$POS2 = sapply(1:nrow(tidy.DNA),FUN=function(x){
      strsplit(tidy.DNA$gt_GT[x],"|")[[1]][3]
    })
    
    Haplo1.mat = select(tidy.DNA,-gt_GT,-gt_GT_alleles,-POS2) 
    Haplo1.mat = matrix(Haplo1.mat$POS1,ncol=length(poly.list),nrow=101,byrow=TRUE) 
    rownames(Haplo1.mat) = paste(1:nrow(Haplo1.mat),"a",sep="")
    
    Haplo2.mat = select(tidy.DNA,-gt_GT,-gt_GT_alleles,-POS1) 
    Haplo2.mat = matrix(Haplo2.mat$POS2,ncol=length(poly.list),nrow=101,byrow=TRUE) 
    rownames(Haplo2.mat) = paste(1:nrow(Haplo2.mat),"b",sep="")
    
    Haplo.mat = rbind(Haplo1.mat,Haplo2.mat)
    
    Haplo.df = as.data.frame(Haplo.mat)
    Haplo.unite = unite(as.data.frame(Haplo.mat),"Haplo",1:ncol(Haplo.mat),sep="")
    Haplo.tab = (table(Haplo.unite$Haplo))
    Haplo.seq = names(Haplo.tab)
    
    #epsilon.i is the number of differences between haplotype and ancestral sequence
    #i is number of times haplotype occurs in sample
    n = nrow(Haplo.df)
    
    theta.h = (2/(n*(n-1)))*sum(sapply(1:length(Haplo.seq),FUN=function(x){
      
      alleles = c(strsplit(Haplo.seq[x],split="")[[1]])
      epsilon.i = length(alleles[alleles=="1"])
      i.squared = as.numeric(Haplo.tab[x])^2
      return(epsilon.i*i.squared)
      
    }))
    
    theta.pi = (2/(n*(n-1)))*sum(sapply(1:length(Haplo.seq),FUN=function(x){
      
      alleles = c(strsplit(Haplo.seq[x],split="")[[1]])
      epsilon.i = length(alleles[alleles=="1"])
      i = as.numeric(Haplo.tab[x])
      return(epsilon.i*i*(n-i))
      
    }))
    
    H = theta.pi - theta.h
    
    EHH.val = 0
    for (i in 1:haploNum){
      
      nh = Haplo.tab[i]
      num = nh*(nh-1)/2
      denom = n*(n-1)/2
      EHH.val = EHH.val + as.numeric(num/denom)
      
    }
    
    Haplo.list = sort(Haplo.tab,decreasing=TRUE)
    H1 = sum((Haplo.list/n)^2)
    H2 = as.numeric(H1 - (Haplo.list[1]/n)^2)
    H2.H1 = H2/H1
    H12 = as.numeric(H1 + 2*(Haplo.list[1]/n)*(Haplo.list[2]/n))
    
    return(c(EHH.val,H1,H2,H12,H2.H1,H))
    
  }
  Homozygosity.list = Homozygosity.fun()

  H1 = Homozygosity.list[2]
  H2 = Homozygosity.list[3]
  H12 = Homozygosity.list[4]
  H2.H1 = Homozygosity.list[5]
  
  #Fay and Wu's H
  H = Homozygosity.list[6]
  
  #EHH and IHH
  EHH.fun = function(Frequency){
    
    library(rehh)
    library(caTools)
    
    hap = data2haplohh(hap_file=paste("~/SLiM/testSampleVCF",x,".txt",sep=""),polarize_vcf = FALSE, verbose = FALSE)
    
    if (3501 %in% hap@positions){
      
      if ((1000 %in% hap@positions)&(6001 %in% hap@positions)){
        
        hap@positions = hap@positions
        
      } else if ((1000 %in% hap@positions)) {
        
        hap@positions = as.numeric(c(hap@positions,6001))
        
        mat.6001 = matrix(0,nrow=nrow(hap@haplo),ncol=1)
        mat.temp2 = cbind(hap@haplo,mat.6001)
        mode(mat.temp2) = "integer"
        hap@haplo = mat.temp2
        
      } else if ((6001 %in% hap@positions)) {
        
        hap@positions = as.numeric(c(1000,hap@positions))
        
        mat.1000 = matrix(0,nrow=nrow(hap@haplo),ncol=1)
        mat.temp2 = cbind(mat.1000,hap@haplo)
        mode(mat.temp2) = "integer"
        hap@haplo = mat.temp2
        
      } else {
        
        hap@positions = as.numeric(c(1000,hap@positions,6001))
        
        mat.1000 = matrix(0,nrow=nrow(hap@haplo),ncol=1)
        mat.6001 = matrix(0,nrow=nrow(hap@haplo),ncol=1)
        mat.temp1 = cbind(mat.1000,hap@haplo)
        mat.temp2 = cbind(mat.temp1,mat.6001)
        mode(mat.temp2) = "integer"
        hap@haplo = mat.temp2
        
      }
      
      slice = which(hap@positions == 3501)
      
      ehh.obj = calc_ehh(hap,mrk=slice,limehh=0.05)
      
    } else {
      
      if ((1000 %in% hap@positions)&(6001 %in% hap@positions)){
        
        hap@positions = as.numeric(c(hap@positions[hap@positions<3501],3501,hap@positions[hap@positions>3501]))
        
        mat.3501 = matrix(Frequency,nrow=nrow(hap@haplo),ncol=1)
        
        mat.temp1 = cbind(hap@haplo[,1:length(hap@positions[hap@positions<3501])],mat.3501)
        
        mat.temp2 = cbind(mat.temp1,hap@haplo[,(length(hap@positions[hap@positions<3501])+1):ncol(hap@haplo)])
        
        mode(mat.temp2) = "integer"
        hap@haplo = mat.temp2
        
        hap@positions = as.numeric(c(hap@positions[hap@positions<3501],3501,hap@positions[hap@positions>3501]))
        
      } else if ((1000 %in% hap@positions)) {
        
        hap@positions = as.numeric(c(hap@positions[hap@positions<3501],3501,hap@positions[hap@positions>3501],6001))
        
        mat.3501 = matrix(Frequency,nrow=nrow(hap@haplo),ncol=1)
        mat.6001 = matrix(0,nrow=nrow(hap@haplo),ncol=1)
        
        mat.temp1 = cbind(hap@haplo[,1:length(hap@positions[hap@positions<3501])],mat.3501)
        
        mat.temp2 = cbind(hap@haplo[,(length(hap@positions[hap@positions<3501])+1):ncol(hap@haplo)],mat.6001)
        mat.temp2 = cbind(mat.temp1,mat.temp2)
        
        mode(mat.temp2) = "integer"
        hap@haplo = mat.temp2
        
        hap@positions = as.numeric(c(hap@positions[hap@positions<3501],3501,hap@positions[hap@positions>3501],6001))
        
      } else if ((6001 %in% hap@positions)) {
        
        mat.1000 = matrix(0,nrow=nrow(hap@haplo),ncol=1)
        mat.3501 = matrix(Frequency,nrow=nrow(hap@haplo),ncol=1)
        
        mat.temp1 = cbind(mat.1000,hap@haplo[,1:length(hap@positions[hap@positions<3501])])
        
        mat.temp2 = cbind(mat.3501,hap@haplo[,(length(hap@positions[hap@positions<3501])+1):ncol(hap@haplo)])
        mat.temp2 = cbind(mat.temp1,mat.temp2)
        
        mode(mat.temp2) = "integer"
        hap@haplo = mat.temp2
        
        hap@positions = as.numeric(c(1000,hap@positions[hap@positions<3501],3501,hap@positions[hap@positions>3501]))
        
      } else {
        
        mat.1000 = matrix(0,nrow=nrow(hap@haplo),ncol=1)
        mat.3501 = matrix(Frequency,nrow=nrow(hap@haplo),ncol=1)
        mat.6001 = matrix(0,nrow=nrow(hap@haplo),ncol=1)
        
        mat.temp1 = cbind(mat.1000,hap@haplo[,1:length(hap@positions[hap@positions<3501])])
        mat.temp1 = cbind(mat.temp1,mat.3501)
        
        mat.temp2 = cbind(hap@haplo[,(length(hap@positions[hap@positions<3501])+1):ncol(hap@haplo)],mat.6001)
        mat.temp2 = cbind(mat.temp1,mat.temp2)
        
        mode(mat.temp2) = "integer"
        hap@haplo = mat.temp2
        
        hap@positions = as.numeric(c(1000,hap@positions[hap@positions<3501],3501,hap@positions[hap@positions>3501],6001))
        
      }
      
      slice = which(hap@positions == 3501)
      
      ehh.obj = calc_ehh(hap,mrk=slice,limehh=0.05)
      
    }
    
    EHH_A_1000 = ehh.obj$ehh$EHH_A[1]
    EHH_A_6001 = ehh.obj$ehh$EHH_A[nrow(ehh.obj$ehh)]
    
    EHH_A = mean(c(EHH_A_1000,EHH_A_6001))
    
    EHH_D_1000 = ehh.obj$ehh$EHH_D[1]
    EHH_D_6001 = ehh.obj$ehh$EHH_D[nrow(ehh.obj$ehh)]
    
    EHH_D = mean(c(EHH_D_1000,EHH_D_6001))
    
    IHH_A = trapz(ehh.obj$ehh$POSITION,ehh.obj$ehh$EHH_A)
    IHH_D = trapz(ehh.obj$ehh$POSITION,ehh.obj$ehh$EHH_D)
    
    return(c(EHH_A,EHH_D,IHH_A,IHH_D))
  }
  EHH.list = EHH.fun(Frequency)
  EHH_A = EHH.list[1]
  EHH_D = EHH.list[2]
  IHH_A = EHH.list[3]
  IHH_D = EHH.list[4]
  
  file.remove(paste("~/SLiM/testSampleVCF",x,".txt",sep=""))
  file.remove(paste("~/SLiM/fixedMutations",x,".txt",sep=""))
  
  return(c(pi,SegSites,TajD,H,Fixed,Singletons,Doubletons,H1,H2,H12,H2.H1,haploNum,Frequency,Singletons/Fixed,EHH_A,EHH_D,IHH_A,IHH_D))
  
}

