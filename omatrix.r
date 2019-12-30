
# the script to calculate the emissional probability matrix of offspring 
# Notes:
# p: initial probability of parents for each marker
# snp: SNP marker genotype with 'hmp' format, the rows (number of markers) is equal to 'nummarker'
# oid: offspring line index used to estimate the IBD structure
# nummarker: the number of test SNP markers, the value is equal to the rows of 'snp' file
# G: generations that the offsprings decented from the parents
# numparent: the number of ancestor parents, default as 24
# rou: correlations between any pairs of flanking markers, that estimated with the offspring-LD level after corrected by parent-LD level

omatrix=function(p,snp,oid=1,nummarker,G,numparent=24,rou){
	
    # Z<-array(0,c(dim(snp)[1],numparent,2))
    P<-matrix(0,nrow=dim(snp)[1],ncol=numparent)        
    Q<-matrix(0,nrow=dim(snp)[1],ncol=numparent) 
	
    if(snp[1,(11+numparent+oid)]==paste(p[2,2],p[2,2],sep="")){
      P[1,]<-as.numeric(p[2,3:(numparent+2)])
    } else if(snp[1,(11+numparent+oid)]==paste(p[3,2],p[3,2],sep="")){
      P[1,]<-as.numeric(p[3,3:(numparent+2)])} else
     {P[1,]<-as.numeric(p[4,3:(numparent+2)])} 
 
      
    if(snp[nummarker,(11+numparent+oid)]==paste(p[2+4*(nummarker-1),2],p[2+4*(nummarker-1),2],sep="")){
      Q[nummarker,]<-as.numeric(p[2+4*(nummarker-1),3:(numparent+2)])
    } else if(snp[nummarker,(11+numparent+oid)]==paste(p[3+4*(nummarker-1),2],p[3+4*(nummarker-1),2],sep="")){
      Q[nummarker,]<-as.numeric(p[3+4*(nummarker-1),3:(numparent+2)])} else
     {Q[nummarker,]<-as.numeric(p[4+4*(nummarker-1),3:(numparent+2)])} 
  
      
    for (i in 1:(nummarker-1)){
	pai<-matrix(0,nrow=1,ncol=numparent)
	Qpai<-matrix(0,nrow=1,ncol=numparent)
                
# determined by generations
    if (rou[i]>1E-6)     
    {lamda<-G*rou[i]/50} else  
    {lamda<-G*(2E-8)}
   #  lamda<-G*rou[i]/1E4   

    r<-matrix((1-exp(-lamda))/numparent,nrow=numparent,ncol=numparent)
   diag(r)<-exp(-lamda)+(1-exp(-lamda))/numparent
  
   if (rou[nummarker-i]>1E-6)
    {Qlamda<-G*(rou[nummarker-i]/50)} else
    {Qlamda<-G*(2E-8)}
   #  Qlamda<-G*rou[nummarker-i]/1E4   

    Qr<-matrix((1-exp(-Qlamda))/numparent,nrow=numparent,ncol=numparent)
   diag(Qr)<-exp(-Qlamda)+(1-exp(-Qlamda))/numparent


    if(snp[i+1,(11+numparent+oid)]==paste(p[2+4*i,2],p[2+4*i,2],sep=""))
      {pai<-as.numeric(p[2+4*i,3:(numparent+2)])} else if
      (snp[i+1,(11+numparent+oid)]==paste(p[3+4*i,2],p[3+4*i,2],sep=""))
     {pai<-as.numeric(p[3+4*i,3:(numparent+2)])} else
      {pai<-as.numeric(p[4+4*i,3:(numparent+2)])}

     if(snp[nummarker-i,(11+numparent+oid)]==paste(p[2+4*(nummarker-1-i),2],p[2+4*(nummarker-1-i),2],sep=""))
      {Qpai<-as.numeric(p[2+4*(nummarker-1-i),3:(numparent+2)])} else if
      (snp[nummarker-i,(11+numparent+oid)]==paste(p[3+4*(nummarker-1-i),2],p[3+4*(nummarker-1-i),2],sep=""))
     {Qpai<-as.numeric(p[3+4*(nummarker-1-i),3:(numparent+2)])} else
      {Qpai<-as.numeric(p[4+4*(nummarker-1-i),3:(numparent+2)])}
 
    
    f=apply(r,2,function(m){
    m*pai/as.vector(m%*%pai)})

    P[i+1,]<-P[i,]%*%t(f)
    
    Qf=apply(Qr,2,function(m){
    m*Qpai/as.vector(m%*%Qpai)})
    Q[nummarker-i,]<-Q[nummarker-i+1,]%*%t(Qf)
     
}
# Z[,,1]<-P
# Z[,,2]<-Q
Z=matrix((P+Q)/2,nrow=nrow(P))
return(Z)
}
 
