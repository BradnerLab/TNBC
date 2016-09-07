setwd('./TNBC/')

#=====================================================================
#===============================DATA INPUT============================
#=====================================================================


expTable = read.table('cufflinks/cuffnorm/output/TNBC_all_fpkm_exprs.txt',sep='\t',header=TRUE)

seTable = read.delim('rose_stitch/SUM159_JQ1_CHEM_SEQ_R1_ROSE/SUM159_JQ1_CHEM_SEQ_R1_peaks_SuperEnhancers_ENHANCER_TO_TOP_GENE.txt',sep='\t')
#now we need to make a mean matrix?

sum159R_3H = log2(apply(expTable[,grep('SUM159R_DMSO_3H',colnames(expTable))],1,mean)/apply(expTable[,grep('SUM159_DMSO_3H',colnames(expTable))],1,mean))


sum159R_24H = log2(apply(expTable[,grep('SUM159R_DMSO_24H',colnames(expTable))],1,mean)/apply(expTable[,grep('SUM159_DMSO_24H',colnames(expTable))],1,mean))



m = as.matrix(cbind(sum159R_3H,sum159R_24H))


fullM = cbind(expTable[,grep('SUM159_DMSO_24H',colnames(expTable))],expTable[,grep('SUM159R_DMSO_24H',colnames(expTable))])



gseaTable = cbind(apply(expTable[,grep('SUM159_DMSO_24H',colnames(expTable))],1,mean),apply(expTable[,grep('SUM159R_DMSO_24H',colnames(expTable))],1,mean))
colnames(gseaTable) = c('SUM159_DMSO_24H','SUM159R_DMSO_24H')
#write.table(gseaTable,'./tables/GSEA_SUM159_SUM159R.txt',quote=FALSE,sep='\t')
#===================================================================
#========================ORDERING BY FOLD CHANGE====================
#===================================================================

#at 24h

changeRows = which(abs(m[,2])>1)


#===================================================================
#======================MAKING TABLE FOR METACORE====================
#===================================================================

pVector = c()

for(i in 1:nrow(fullM)){
	if(length(unique(as.numeric(fullM[i,c(1,2,3,4)]))) == 1){
		print('balls')
		pVector = c(pVector,1)
	}else{
	pVector = c(pVector,t.test(fullM[i,c(1,2)],fullM[i,c(3,4)])$p.value)
	}
	
}

sigRows = which(pVector < 0.01)


diffRows = intersect(sigRows,changeRows)

diffM = cbind(expTable[,grep('SUM159_DMSO_24H',colnames(expTable))],expTable[,grep('SUM159R_DMSO_24H',colnames(expTable))],m[,3],pVector)

colnames(diffM) = c('SUM159_DMSO_24H_1','SUM159_DMSO_24H_2','SUM159R_DMSO_24H_1','SUM159R_DMSO_24H_2','LOG2_FOLD_SUM159R_SUM159_24H','PVALUE')

diffM  = diffM[diffRows,]
diffMOrdered = diffM[order(diffM[,5]),]

write.table(diffMOrdered,file='tables/HG19_SUM159R_VS_SUM159_DMSO_24H.txt',quote=FALSE,sep='\t')

#===================================================================
#============================CLUSTERING=============================
#===================================================================



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
