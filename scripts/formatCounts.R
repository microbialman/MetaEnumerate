library(stringr)

#load the counts table
counts=read.table(commandArgs(trailingOnly=T)[2],sep="\t",header=T)

names=colnames(counts)[7:ncol(counts)]
m=str_match(names,"^sample_mappings\\.dir\\.(\\S+)\\.(\\S+)\\.mapped\\.bam$")
colnames(counts)[7:ncol(counts)]=m[,3]

#remove extra columns and write out the raw counts
raw=counts[,c(1,7:ncol(counts))]
colnames(raw)[1]=commandArgs(trailingOnly=T)[1]
write.table(raw,commandArgs(trailingOnly=T)[4],sep="\t",row.names = F)

#generate TPM normalised table
#divide count by feature length
lendiv=apply(counts[,7:ncol(counts)],2,function(x) x/(counts[,6]/1000))
#get per million scale per sample
if(nrow(counts)>1){
permil=colSums(lendiv)/1000000}else{
permil=lendiv/1000000
}
#get tpm from rpk
if(nrow(counts)>1){
tpm=t(apply(lendiv,1,function(x) x/permil))}else{
  tpm=t(lendiv/permil)
}
#add back feature names
tpmdat=cbind(raw[,1],data.frame(tpm))
colnames(tpmdat)[1]=commandArgs(trailingOnly=T)[1]
write.table(tpmdat,commandArgs(trailingOnly=T)[3],sep="\t",row.names = F)
