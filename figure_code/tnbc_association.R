setwd('./TNBC/')

tnbcTable= read.delim('tnbc_association_matrix.txt',header=TRUE,sep='\t')


pValueMatrix = matrix(nrow = length(3:ncol(tnbcTable)),ncol=2)
colnames(pValueMatrix)= c('p-value','bonferroni corrected p-value')
rownames(pValueMatrix)= colnames(tnbcTable[3:ncol(tnbcTable)])
for(i in 3:ncol(tnbcTable)){
	
	print(colnames(tnbcTable)[i])
	m = matrix(nrow=2,ncol=2,dimnames = list(c("Sensitive", "Insensitive"),c("Feature", "No_Feature")))
	m[1,1] = length(intersect(which(tnbcTable[,2]==1),which(tnbcTable[,i]==1)))
	m[1,2] = length(intersect(which(tnbcTable[,2]==1),which(tnbcTable[,i]==0)))
	m[2,1] = length(intersect(which(tnbcTable[,2]==0),which(tnbcTable[,i]==1)))
	m[2,2] = length(intersect(which(tnbcTable[,2]==0),which(tnbcTable[,i]==0)))
	pValue = fisher.test(m,alternative ='greater')$p.value
	correctedPValue = pValue * length(3:ncol(tnbcTable))
	pValueMatrix[(i-2),] = c(pValue,correctedPValue)
         
            
       
	
}

write.table(pValueMatrix,'tnbc_association_fishers.txt',row.names=TRUE,col.names=TRUE,sep='\t')