#150511_tnbc_brd4_figures.R
setwd('./TNBC/')


#==========================================================
#=====================DATA INPUT FILES=====================
#==========================================================

#load up the signal tables
hcc1395_table = read.table('tables/HG19_ENRICHED_HCC1395_BRD4_DMSO_-0_+0_signal.txt',header=TRUE,sep='\t')

mda436_table = read.table('tables/HG19_ENRICHED_MDA436_BRD4_DMSO_-0_+0_signal.txt',header=TRUE,sep='\t')

sum1315_table = read.table('tables/HG19_ENRICHED_SUM1315_BRD4_DMSO_-0_+0_signal.txt',header=TRUE,sep='\t')

sum149_table = read.table('tables/HG19_ENRICHED_SUM149_BRD4_DMSO_-0_+0_signal.txt',header=TRUE,sep='\t')

sum149R_table = read.table('tables/HG19_ENRICHED_SUM149R_BRD4_DMSO_-0_+0_signal.txt',header=TRUE,sep='\t')

sum159_table = read.table('tables/HG19_ENRICHED_SUM159_BRD4_DMSO_R1_-0_+0_signal.txt',header=TRUE,sep='\t')

sum159R_table = read.table('tables/HG19_ENRICHED_SUM159R_BRD4_DMSO_-0_+0_signal.txt',header=TRUE,sep='\t')

sum185_table = read.table('tables/HG19_ENRICHED_SUM185_BRD4_DMSO_-0_+0_signal.txt',header=TRUE,sep='\t')

t47d_table = read.table('tables/HG19_ENRICHED_T47D_BRD4_DMSO_-0_+0_signal.txt',header=TRUE,sep='\t')

#==========================================================
#=============ALL REGIONS FOLD CHANGE BOX PLOT=============
#==========================================================



hcc1395_fold = log2(hcc1395_table[,4]/hcc1395_table[,3])
mda436_fold = log2(mda436_table[,4]/mda436_table[,3])
sum1315_fold = log2(sum1315_table[,5]/sum1315_table[,3])
sum149_fold = log2(sum149_table[,4]/sum149_table[,3])
sum149R_fold = log2(sum149R_table[,4]/sum149R_table[,3])
sum159_fold = log2(sum159_table[,4]/sum159_table[,3])
sum159R_fold = log2(sum159R_table[,5]/sum159R_table[,3])
sum185_fold = log2(sum185_table[,4]/sum185_table[,3])
t47d_fold = log2(t47d_table[,4]/t47d_table[,3])

boxplot(hcc1395_fold,mda436_fold,sum1315_fold,sum149_fold,sum149R_fold,sum159_fold,sum159R_fold,sum185_fold,t47d_fold,cex=0)

pdf(file='figures/150511_tnbc_brd4_fold_box.pdf',width = 8,height = 5)

namesList = c('SUM149R','SUM159R','HCC1395','SUM185','T47D','MDA436','SUM149','SUM1315','SUM159')
boxplot(sum149R_fold,sum159R_fold,hcc1395_fold,sum185_fold,t47d_fold,mda436_fold,sum149_fold,sum1315_fold,sum159_fold,cex=0,las=2,ylim =c(-14,2),names=namesList,ylab='Log2 fold change in BRD4 +/- JQ1',main='BRD4 bound regions in TNBC')
abline(h=0)
dev.off()
#==========================================================
#===========TOP 1000 REGIONS FOLD CHANGE BOX PLOT==========
#==========================================================


foldMatrix = matrix(nrow=1000,ncol=9)



topFold <- function(sigTable,dmsoCol,jq1Col,top=3000){
	
	dmsoOrder = order(sigTable[,dmsoCol],decreasing=TRUE)
	foldVector = log2(sigTable[dmsoOrder[1:top],jq1Col]/sigTable[dmsoOrder[1:top],dmsoCol])
	
	return(foldVector)
	
	
}

hcc1395_fold = topFold(hcc1395_table,3,4)
mda436_fold = topFold(mda436_table,3,4)
sum1315_fold = topFold(sum1315_table,3,5)
sum149_fold = topFold(sum149_table,3,4)
sum149R_fold = topFold(sum149R_table,3,4)
sum159_fold = topFold(sum159_table,3,4)
sum159R_fold = topFold(sum159R_table,3,5)
sum185_fold = topFold(sum185_table,3,4)
t47d_fold = topFold(t47d_table,3,4)

boxplot(hcc1395_fold,mda436_fold,sum1315_fold,sum149_fold,sum149R_fold,sum159_fold,sum159R_fold,sum185_fold,t47d_fold,cex=0,las=2,ylim =c(-10,2))


boxplot(sum149R_fold,sum159R_fold,sum185_fold,t47d_fold,sum149_fold,hcc1395_fold,mda436_fold,sum1315_fold,sum159_fold,cex=0,las=2,ylim =c(-10,2))
abline(h=0)



