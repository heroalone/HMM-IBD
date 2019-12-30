# the script to calculate the initial probability matrix of parents
# Notes:
# k: the number of tested SNP markers
# snp_parent: the parents' genotype with the first four cols as marker infomation
# numparent: the number of parents, the default as 24


pmatrix=function(k,snp_parent,numparent=24,p1=0.97){
	pres=lapply(1:k,function(i){
   # p1=0.97               
	marker<-vector(length=0)
	
	marker[1]<-"marker"
	marker[2]<-snp_parent[i,1]
	marker[3]<-3
    marker[4]<-snp_parent[i,3]
    marker[5:(2+numparent)]<-""
    
    allele1<-vector(length=0)
	allele1[1]<-"allele"
	allele1[2]<-"NN"
	allele1[3:(2+numparent)]<-round(1/numparent, 3)	
	
	allele2<-vector(length=0)
	allele3<-vector(length=0)
	
	a2<-unlist(strsplit(snp_parent[i,2],"/"))[1]
	a3<-unlist(strsplit(snp_parent[i,2],"/"))[2]
	aa2<-paste(a2,a2,sep = "")
	aa3<-paste(a3,a3,sep = "")
	allele2[1]<-"allele"
	allele2[2]<-a2
	
	allele3[1]<-"allele"
	allele3[2]<-a3
	
	n1<-(snp_parent[i,5:(4+numparent)]==aa2)
	n2<-(snp_parent[i,5:(4+numparent)]==aa3)
	

if (sum(n1+n2)==numparent) 
    {allele2[2+which(snp_parent[i,5:(4+numparent)]==aa2)]<-p1/sum(n1)
     allele2[2+which(snp_parent[i,5:(4+numparent)]!=aa2)]<-(1-p1)/sum(n2)	
     allele3[2+which(snp_parent[i,5:(4+numparent)]==aa3)]<-p1/sum(n2)
     allele3[2+which(snp_parent[i,5:(4+numparent)]!=aa3)]<-(1-p1)/sum(n1)}
    if (sum(n1+n2)!=numparent)
      { p<-which(n1+n2==0)
      	n1[1,p]<-1/2
        n2[1,p]<-1/2
        allele2_m<-n1*p1/sum(n1)
        allele3_m<-n2*p1/sum(n2)
        allele2_m[(allele2_m==0)]<-(1-p1)/sum(allele2_m==0)
        allele3_m[(allele3_m==0)]<-(1-p1)/sum(allele3_m==0)
        allele2[3:(2+numparent)]<-allele2_m
        allele3[3:(2+numparent)]<-allele3_m
      }

	pres<-data.frame(marker,allele2,allele3,allele1)
	
	return(t(pres))
})
return(do.call(rbind,pres))
}
