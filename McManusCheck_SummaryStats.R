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

  EHH = NA #Homozygosity.list[1]
  H1 = Homozygosity.list[2]
  H2 = Homozygosity.list[3]
  H12 = Homozygosity.list[4]
  H2.H1 = Homozygosity.list[5]
  
  #Fay and Wu's H
  H = Homozygosity.list[6]
  
  #iHH
  iHH = NA
  
  file.remove(paste("~/SLiM/testSampleVCF",x,".txt",sep=""))
  file.remove(paste("~/SLiM/fixedMutations",x,".txt",sep=""))
  
  return(c(pi,SegSites,TajD,H,iHH,Fixed,Singletons,Doubletons,H1,H2,H12,H2.H1,haploNum,Frequency,EHH,Singletons/Fixed))
  
}

