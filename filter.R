args=commandArgs(T)
infile=args[1]
outfile=args[2]
thr=as.numeric(args[3])
mis=as.numeric(args[4])
pop=args[5]

dat=read.table(infile,row.names = 1,col.names = c('ID','a','b','h','miss'))
dat$missRate=dat$miss/rowSums(dat)
muchMis=dat$missRate>mis

if(pop == "RIL2")pop ="F2" 
if( grepl ("RIL",pop)){
pvalue=sapply(rownames(dat[!muchMis,]),function(x) chisq.test(c(dat[x,1],dat[x,2]),p = c(0.5,0.5))$p.value )
}else if(pop=="F2"){
pvalue=sapply(rownames(dat[!muchMis,]),function(x) chisq.test(c(dat[x,1],dat[x,2],dat[x,3]),p = c(0.25,0.25,0.5))$p.value )
}

valid=names(pvalue[pvalue>=thr])

out=dat
out$pvalue=NA
out[names(pvalue),'pvalue']=pvalue

write.table(valid,file=paste0(outfile,'.valid'),quote=F,row.names = F,col.names = F,sep="\t")
write.csv(out,file=paste0(outfile,'.all.pvalues.csv'),quote=F,row.names=T)
