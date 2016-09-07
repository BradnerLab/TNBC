

#150901_ic50_heatmap.R
setwd('./TNBC/')
#==========================================================================
#=========================DATA INPUT=======================================
#==========================================================================

ic50 = read.delim('./tnbc_ic50_revision.txt',header=TRUE,sep='\t')

ic50_matrix = as.matrix(ic50[2:nrow(ic50),2:ncol(ic50)])

ic50_matrix = apply(ic50_matrix,c(1,2),as.numeric)
rownames(ic50_matrix) = ic50[2:nrow(ic50),1]

classList = as.vector(ic50[1,2:ncol(ic50)])


#==========================================================================
#=============================CLUSTERING===================================
#==========================================================================

d_row = dist(ic50_matrix)

hc_row= hclust(d_row)

d_col= dist(t(ic50_matrix))
hc_col = hclust(d_col)
#==========================================================================
#==============================HEATMAP=====================================
#==========================================================================



clusterMatrix = log10(ic50_matrix[,order(ic50_matrix[1,])])

colOrder = order(clusterMatrix[2,])
rowOrder = c(7,1,4,3,6,5,2)
#Set the color spectrum
colorSpectrum <- colorRampPalette(c("red","red","pink","black"))(100)

#setting a color data range
#minValue <- quantile(clusterMatrix,na.rm=TRUE,prob=0.025,names=FALSE)
#maxValue <- quantile(clusterMatrix,na.rm=TRUE,prob=0.975,names=FALSE)
minValue = -2
maxValue = 1.3
color_cuts <- seq(minValue,maxValue,length=100)
color_cuts <- c(-2, color_cuts,max(clusterMatrix))

#add one extra min color to even out sampling
colorSpectrum <- c(colorSpectrum[1],colorSpectrum)
pdf(file='./TNBC_IC50_HEATMAP_REVISION.pdf',width = 10,height = 10)

layout(matrix(data=c(1,1,1,1,1,1,2),ncol= 7))

par(mar=c(1,6,8,1))
image(1:ncol(clusterMatrix),1:nrow(clusterMatrix),t(clusterMatrix[rowOrder,colOrder]),breaks=color_cuts,col=colorSpectrum,xaxt="n",yaxt="n",ylab='',xlab='')
axis(3,1:ncol(clusterMatrix),colnames(clusterMatrix)[colOrder],las=2)
axis(2,1:nrow(clusterMatrix),rownames(clusterMatrix)[rowOrder],las=2)

par(mar=c(1,4,8,1))

image(1:2,color_cuts[2:101],t(matrix(data=color_cuts[2:101],ncol=2,nrow=100)),breaks=color_cuts,col=colorSpectrum,xaxt="n",xlab="",ylab="Log10 uM")
dev.off()
