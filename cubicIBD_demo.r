# cubicIBD_demo.r

rm(list=ls())

# load the script to calculate the initial probability matrix of parents
source('pmatrix.r')
# load the script to calculate the emissional probability matrix of offspring 
source('omatrix.r')


# calculate the initial probs. for parents 
 # read hmp genotype 

snp = read.delim('./data/demo_geno_chr10.hmp.txt',stringsAs=F,check.names=F)
pp<-pmatrix(nrow(snp),snp[,c(1:4,12:35)])  

# read pseudo map 
mapchr = read.delim('./data/demo_map.txt')
rou = diff(mapchr$pos_g)
rou = ifelse(rou<0,0,rou)

# calculate the offspring probs 
Z<-omatrix(pp,snp,oid=1,nrow(snp),G=9,24,rou)

write.table(Z,'./res/posterior_probs_chr10.txt',sep='\t',quote=FALSE,row.names=FALSE,col.names=FALSE)
cat('===== IBD probability calculation done! =====\n')

