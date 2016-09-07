setwd('./TNBC/')

#=====================================================================
#===============================DATA INPUT============================
#=====================================================================



#chemseq
chemseqTable = read.delim('./TNBC/dynamic/TNBC_CHEMSEQ/output/HG19_SUM159_JQ1_CHEM_SEQ_R1_SUM159R_JQ1_CHEM_SEQ_R1_merged_MERGED_SUPERS_RANK_TABLE.txt',sep='\t')


#sum159_brd4
sum159_brd4 = read.delim('mappedFolder/n1/HG19_SUM159_JQ1_CHEM_SEQ_R1_SUM159R_JQ1_CHEM_SEQ_R1_merged_MERGED_SUPERS_RANK/HG19_SUM159_JQ1_CHEM_SEQ_R1_SUM159R_JQ1_CHEM_SEQ_R1_merged_MERGED_SUPERS_RANK_SUM159_BRD4_DMSO_R1.gff')

sum159R_brd4 = read.delim('mappedFolder/n1/HG19_SUM159_JQ1_CHEM_SEQ_R1_SUM159R_JQ1_CHEM_SEQ_R1_merged_MERGED_SUPERS_RANK/HG19_SUM159_JQ1_CHEM_SEQ_R1_SUM159R_JQ1_CHEM_SEQ_R1_merged_MERGED_SUPERS_RANK_SUM159R_BRD4_DMSO_R2.gff')

m = cbind(chemseqTable[,7],sum159_brd4[,3],sum159R_brd4[,3])
rownames(m) = chemseqTable[,1]
colnames(m) = c('LOG2_CHANGE','SUM159_BRD4','SUM159R_BRD4')

#median normalizing brd4

medianFoldVector = log2(m[,3]/median(m[,3])/m[,2]/median(m[,2]))



#===================================================================
#===========================GAIN LOSS REGIONS=========================
#===================================================================
gainedRows = which(m[,1] >= 1)
conservedRows = intersect(which(m[,1]<= 0.25),which(m[,1]>= -0.25))
lostRows = which(m[,1]<= -1)



par(mfrow=c(1,3))
boxplot(m[gainedRows,c(2,3)],cex=0)
boxplot(m[conservedRows,c(2,3)],cex=0)

boxplot(m[lostRows,c(2,3)],cex=0)

pdf(file='figures/141016_chemseq_change_brd4.pdf',width = 6,height =4)

boxplot(medianFoldVector[gainedRows],medianFoldVector[conservedRows],medianFoldVector[lostRows],cex=0,names=c('gained','conserved','lost'),ylab='Log2 fold change in median normalized BRD4 signal',main='BRD4 change at regions of gained/conserved/lost JQ1 binding')
abline(h=0)
dev.off()
