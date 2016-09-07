setwd('./TNBC/')

#=====================================================================
#===============================DATA INPUT============================
#=====================================================================


expTable = read.table('cufflinks/cuffnorm/output/TNBC_all_fpkm_exprs.txt',sep='\t',header=TRUE)

#chemseq
chemseqTable = read.delim('./TNBC/dynamic/TNBC_CHEMSEQ/output/HG19_SUM159_JQ1_CHEM_SEQ_R1_SUM159R_JQ1_CHEM_SEQ_R1_merged_MERGED_SUPERS_GENE_TABLE.txt',sep='\t')

#k27ac

#k27
k27Table=read.delim('./TNBC/dynamic/TNBC_H3K27AC_BASELINE/output/HG19_SUM159_H3K27AC_DMSO_SUM159R_H3K27AC_DMSO_merged_MERGED_SUPERS_GENE_TABLE.txt',sep='\t')

#k27_r1
k27R1Table=read.delim('./TNBC/dynamic/TNBC_H3K27AC_R1_BASELINE/output/HG19_SUM159_H3K27AC_DMSO_R1_SUM159R_H3K27AC_DMSO_R1_merged_MERGED_SUPERS_GENE_TABLE.txt',sep='\t')


brd4Table = read.delim('./TNBC/dynamic/TNBC_BASELINE/output/HG19_SUM159_BRD4_DMSO_SUM159R_BRD4_DMSO_merged_MERGED_SUPERS_GENE_TABLE.txt')
#now we need to make a mean matrix?
sum159R_3H = log2(apply(expTable[,grep('SUM159R_DMSO_3H',colnames(expTable))],1,mean)/apply(expTable[,grep('SUM159_DMSO_3H',colnames(expTable))],1,mean))


sum159R_24H = log2(apply(expTable[,grep('SUM159R_DMSO_24H',colnames(expTable))],1,mean)/apply(expTable[,grep('SUM159_DMSO_24H',colnames(expTable))],1,mean))



m = as.matrix(cbind(sum159R_3H,sum159R_24H))


fullM = cbind(expTable[,grep('SUM159_DMSO_24H',colnames(expTable))],expTable[,grep('SUM159R_DMSO_24H',colnames(expTable))])

#fullM = cbind(expTable[,grep('SUM159_DMSO_3H',colnames(expTable))],expTable[,grep('SUM159R_DMSO_3H',colnames(expTable))])


fullM = cbind(apply(expTable[,grep('SUM159_DMSO_3H',colnames(expTable))],1,mean),apply(expTable[,grep('SUM159_DMSO_24H',colnames(expTable))],1,mean),apply(expTable[,grep('SUM159R_DMSO_3H',colnames(expTable))],1,mean),apply(expTable[,grep('SUM159R_DMSO_24H',colnames(expTable))],1,mean))
colnames(fullM) = c('SUM159_DMSO_3H','SUM159_DMSO_24H','SUM159R_DMSO_3H','SUM159R_DMSO_24H')

gseaTable = cbind(apply(expTable[,grep('SUM159_DMSO_24H',colnames(expTable))],1,mean),apply(expTable[,grep('SUM159R_DMSO_24H',colnames(expTable))],1,mean))
colnames(gseaTable) = c('SUM159_DMSO_24H','SUM159R_DMSO_24H')
write.table(gseaTable,'./tables/GSEA_SUM159_SUM159R.txt',quote=FALSE,sep='\t')

#===================================================================
#===========================GAIN LOSS GENES=========================
#===================================================================


seTable = chemseqTable
#getting expressed genes
expressedGenes = rownames(fullM)[which(apply(fullM,1,max) >= 1)]


#gain genes 


gainedGenes = intersect(expressedGenes,as.character(seTable[which(seTable[,6]>=1),1]))

lostGenes = intersect(expressedGenes,as.character(seTable[which(seTable[,6]<= -1),1]))

conservedGenes = intersect(as.character(seTable[which(seTable[,6]<= 0.25),1]),as.character(seTable[which(seTable[,6] >= -0.25),1]))
conservedGenes = intersect(expressedGenes,conservedGenes)

gainedRows = c()
for(gene in gainedGenes){
	
	gainedRows = c(gainedRows,which(rownames(fullM)==gene))
	
}

lostRows = c()
for(gene in lostGenes){
	
	lostRows = c(lostRows,which(rownames(fullM)==gene))
	
}

conservedRows = c()
for(gene in conservedGenes){
	
	conservedRows = c(conservedRows,which(rownames(fullM)==gene))
	
}

#===================================================================
#========================EXPRESSION BOX PLOTS=======================
#===================================================================


#par(mfrow=c(1,2))
#boxplot(log2(fullM[gainedRows,]),cex=0,las=2)

#boxplot(log2(fullM[lostRows,]),cex=0,las=2)

pdf(file='figures/141016_k27_change_expression.pdf',width = 6,height =4)
boxplot(m[gainedRows,2],m[conservedRows,2],m[lostRows,2],cex=0,col=c('red','grey','blue'),cex=0,ylim =c(-2.5,2.5),names=c('gained','conserved','lost'),ylab='Log2 fold change in exprssion',main='Chromatin changes based off of H3K27ac Chip-Seq')

abline(h=0)
dev.off()






#===================================================================
#========================ORDERING BY FOLD CHANGE====================
#===================================================================

#at 24h


changeRows = which(abs(m[,2])>1)



#===================================================================
#============================CLUSTERING=============================
#===================================================================

d_row = dist(fullM)

hc_row= hclust(d_row)

d_col= dist(t(fullM))
hc_col = hclust(d_col)


#===================================================================
#=========================MAKING HEATMAPS===========================
#===================================================================
clusterMatrix = fullM[changeRows[order(m[changeRows,2])],]

#clusterMatrix = fullM[hc_row$order,hc_col$order]
#clusterMatrix = fullM[order(m[,2]),]

for(i in 1:nrow(clusterMatrix)){
	
	clusterMatrix[i,] = log2(clusterMatrix[i,]/median(as.numeric(clusterMatrix[i,])))
	
}

#Set the color spectrum
colorSpectrum <- colorRampPalette(c("blue","white","red"))(100)

#setting a color data range
minValue <- quantile(clusterMatrix,na.rm=TRUE,prob=0.025,names=FALSE)
maxValue <- quantile(clusterMatrix,na.rm=TRUE,prob=0.975,names=FALSE)
#color_cuts <- seq(-1,1,length=100)
minValue = -2
maxValue = 2
color_cuts <- seq(minValue,maxValue,length=100)

color_cuts <- c(min(clusterMatrix), color_cuts,max(clusterMatrix))

#add one extra min color to even out sampling
colorSpectrum <- c(colorSpectrum[1],colorSpectrum)

#Making png

#clusterPNGFile = paste(outputFolder,genome,'_',analysisName,'_2d_cluster.png',sep='')
png(filename = 'figures/141015_parentalResistantExpressionHeatmap.png',width = 800,height =800)
layout(matrix(data=c(1,1,1,1,1,2,2),ncol= 7))


image(1:ncol(clusterMatrix),1:nrow(clusterMatrix),t(clusterMatrix),breaks=color_cuts,col=colorSpectrum,xaxt="n",yaxt="n",ylab='',xlab='')


image(1:2,color_cuts[2:101],t(matrix(data=color_cuts[2:101],ncol=2,nrow=100)),breaks=color_cuts,col=colorSpectrum,xaxt="n",xlab="",ylab="Fold vs. median")
dev.off()


#===================================================================
#=========================MAKING HEATMAPS===========================
#===================================================================
