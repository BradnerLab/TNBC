setwd('./TNBC/')

#=====================================================================
#===============================DATA INPUT============================
#=====================================================================


expTable = read.table('cufflinks/cuffnorm/output/TNBC_all_fpkm_exprs.txt',sep='\t',header=TRUE)

seTable = read.delim('./TNBC/rose_stitch/SUM159_JQ1_CHEM_SEQ_R1_ROSE/SUM159_JQ1_CHEM_SEQ_R1_peaks_SuperEnhancers_ENHANCER_TO_TOP_GENE.txt',sep='\t')
#now we need to make a mean matrix?

sum159_3H = log2(apply(expTable[,grep('SUM159_JQ1_3H',colnames(expTable))],1,mean)/apply(expTable[,grep('SUM159_DMSO_3H',colnames(expTable))],1,mean))
sum159_12H = log2(apply(expTable[,grep('SUM159_JQ1_12H',colnames(expTable))],1,mean)/apply(expTable[,grep('SUM159_DMSO_12H',colnames(expTable))],1,mean))
sum159_24H = log2(apply(expTable[,grep('SUM159_JQ1_24H',colnames(expTable))],1,mean)/apply(expTable[,grep('SUM159_DMSO_24H',colnames(expTable))],1,mean))


sum159R_24H = log2(apply(expTable[,grep('SUM159R_JQ1_24H',colnames(expTable))],1,mean)/apply(expTable[,grep('SUM159R_DMSO_24H',colnames(expTable))],1,mean))



m = as.matrix(cbind(sum159_24H,sum159R_24H))

meanMatrix = as.matrix(cbind(apply(expTable[,grep('SUM159_DMSO_24H',colnames(expTable))],1,mean),apply(expTable[,grep('SUM159_JQ1_24H',colnames(expTable))],1,mean),apply(expTable[,grep('SUM159R_DMSO_24H',colnames(expTable))],1,mean),apply(expTable[,grep('SUM159R_JQ1_24H',colnames(expTable))],1,mean)))
colnames(meanMatrix) = c('SUM159_DMSO_24H','SUM159_JQ1_24H','SUM159R_DMSO_24H','SUM159R_JQ1_24H')



fullM = cbind(expTable[,grep('SUM159_DMSO_3H',colnames(expTable))],expTable[,grep('SUM159_DMSO_24H',colnames(expTable))],expTable[,grep('SUM159_JQ1_3H',colnames(expTable))],expTable[,grep('SUM159_JQ1_24H',colnames(expTable))],expTable[,grep('SUM159R_DMSO_3H',colnames(expTable))],expTable[,grep('SUM159R_DMSO_24H',colnames(expTable))],expTable[,grep('SUM159R_JQ1_3H',colnames(expTable))],expTable[,grep('SUM159R_JQ1_24H',colnames(expTable))])



fullM = cbind(expTable[,grep('SUM159_DMSO_24H',colnames(expTable))],expTable[,grep('SUM159_JQ1_24H',colnames(expTable))],expTable[,grep('SUM159R_DMSO_24H',colnames(expTable))],expTable[,grep('SUM159R_JQ1_24H',colnames(expTable))])



#===================================================================
#========================ORDERING BY FOLD CHANGE====================
#===================================================================

#at 24h


changeRows = which(abs(m[,1])>1)

upRows = which(m[,1]>1)
downRows = which(m[,1] < -1)


par(mfrow=c(2,1))

boxplot(log2(meanMatrix[upRows,]),cex=0,ylim =c(-8,10))
boxplot(log2(meanMatrix[downRows,]),cex=0,ylim =c(-8,10))

pdf(file='figures/141019_exp_resistance_jq1_box.pdf',width  =4,height = 8)
par(mfrow=c(2,1))
boxplot(meanMatrix[upRows,],cex=0,ylim =c(0,50))
boxplot(meanMatrix[downRows,],cex=0,ylim =c(0,50))
dev.off()

pdf(file='figures/141019_exp_fold_change_resistance_jq1_box.pdf',width  =4,height = 8)
par(mfrow=c(2,1))
boxplot(m[upRows,],cex=0,ylim =c(-3,3))
boxplot(m[downRows,],cex=0,ylim = c(-3,3))
dev.off()
#===================================================================
#=========================MAKING HEATMAPS===========================
#===================================================================
clusterMatrix = fullM[changeRows[order(m[changeRows,1])],]
for(i in 1:nrow(clusterMatrix)){
	
	clusterMatrix[i,] = log2(clusterMatrix[i,]/median(as.numeric(clusterMatrix[i,])))
	
}

#Set the color spectrum
colorSpectrum <- colorRampPalette(c("blue","white","red"))(100)

#setting a color data range
minValue <- quantile(clusterMatrix,na.rm=TRUE,prob=0.025,names=FALSE)
maxValue <- quantile(clusterMatrix,na.rm=TRUE,prob=0.975,names=FALSE)
color_cuts <- seq(minValue,maxValue,length=100)
#color_cuts <- seq(-1,1,length=100)
color_cuts <- c(min(clusterMatrix), color_cuts,max(clusterMatrix))

#add one extra min color to even out sampling
colorSpectrum <- c(colorSpectrum[1],colorSpectrum)

#Making png

#clusterPNGFile = paste(outputFolder,genome,'_',analysisName,'_2d_cluster.png',sep='')
png(filename = 'figures/141019_expressionHeatmap.png',width = 800,height =800)
layout(matrix(data=c(1,1,1,1,1,2,2),ncol= 7))


image(1:ncol(clusterMatrix),1:nrow(clusterMatrix),t(clusterMatrix),breaks=color_cuts,col=colorSpectrum,xaxt="n",yaxt="n",ylab='',xlab='')


image(1:2,color_cuts[2:101],t(matrix(data=color_cuts[2:101],ncol=2,nrow=100)),breaks=color_cuts,col=colorSpectrum,xaxt="n",xlab="",ylab="Fold vs. median")
dev.off()


#===================================================================
#=========================MAKING HEATMAPS===========================
#===================================================================
