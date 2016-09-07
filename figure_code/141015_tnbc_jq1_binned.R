setwd('./TNBC/')

#===================================================================
#========================HELPER FUNCTIONS===========================
#===================================================================

error.bar <- function(x, y, upper, lower=upper, length=0.1,...){
	if(length(x) != length(y) | length(y) !=length(lower) | length(lower) != length(upper))
	stop("vectors must be same length")
	arrows(x,y+upper, x, y-lower, angle=90, code=3, length=length, ...)
	}

#=========================================================================
#===========================DATA INPUT====================================
#=========================================================================


chemseqMapped = read.table('mappedFolder/n1/HG19_ENRICHED_SUM159_JQ1_CHEM_SEQ_R1_-500_+500/HG19_ENRICHED_SUM159_JQ1_CHEM_SEQ_R1_-500_+500_SUM159_JQ1_CHEM_SEQ_R1.gff',sep='\t',header=TRUE)

brd4Mapped = read.table('mappedFolder/n1/HG19_ENRICHED_SUM159_JQ1_CHEM_SEQ_R1_-500_+500/HG19_ENRICHED_SUM159_JQ1_CHEM_SEQ_R1_-500_+500_SUM159_BRD4_DMSO_R1.gff',sep='\t',header=TRUE)

brd4JQ1Mapped = read.table('mappedFolder/n1/HG19_ENRICHED_SUM159_JQ1_CHEM_SEQ_R1_-500_+500/HG19_ENRICHED_SUM159_JQ1_CHEM_SEQ_R1_-500_+500_SUM159_BRD4_JQ1_R1.gff',sep='\t',header=TRUE)

wceMapped = read.table('mappedFolder/n1/HG19_ENRICHED_SUM159_JQ1_CHEM_SEQ_R1_-500_+500/HG19_ENRICHED_SUM159_JQ1_CHEM_SEQ_R1_-500_+500_SUM159_BRD4_DMSO_WCE_R1.gff',sep='\t',header=TRUE)

h3k27acMapped = read.table('mappedFolder/n1/HG19_ENRICHED_SUM159_JQ1_CHEM_SEQ_R1_-500_+500/HG19_ENRICHED_SUM159_JQ1_CHEM_SEQ_R1_-500_+500_SUM159_H3K27AC_DMSO_R1.gff',sep='\t',header=TRUE)

h3k27acJQ1Mapped = read.table('mappedFolder/n1/HG19_ENRICHED_SUM159_JQ1_CHEM_SEQ_R1_-500_+500/HG19_ENRICHED_SUM159_JQ1_CHEM_SEQ_R1_-500_+500_SUM159_H3K27AC_JQ1_R1.gff',sep='\t',header=TRUE)
#m = cbind(chemseqMapped[,3]-wceMapped[,3],brd4Mapped[,3]-wceMapped[,3],brd4JQ1Mapped[,3]-wceMapped[,3],h3k27acMapped[,3]-wceMapped[,3],h3k27acJQ1Mapped[,3]-wceMapped[,3])

m = cbind(chemseqMapped[,3],brd4Mapped[,3],brd4JQ1Mapped[,3],h3k27acMapped[,3],h3k27acJQ1Mapped[,3])


m[which(m < 0)] <- 0

rownames(m) = chemseqMapped[,2]
colnames(m) = c('SUM159_JQ1_CHEM_SEQ_R1','SUM159_BRD4_DSMO_R1','SUM159_BRD4_JQ1_R1','SUM159_H3K27AC_DMSO_R1','SUM159_H3K27AC_JQ1_R1')
 


#=========================================================================
#=====================SIMPLE CORRELATION SCATTER==========================
#=========================================================================



pdf(file='./figures/140924_tnbc_jq1_scatter.pdf',width = 8,height =4)
par(mfrow=c(1,2))
#for BRD4
plot(m[,1],m[,2],xlim =c(0,3),ylim = c(0,3),col=rgb(1,0,0,0.1),pch=19,cex=0.5,xlab='BIO-JQ1',ylab='BRD4')
abline(lm(m[,2]~m[,1]))
#lines(loess.smooth(m[,1],m[,2]),col='black',lwd=2)
cor(m[,1],m[,2])



#for h3k27ac
plot(m[,1],m[,4],xlim =c(0,8),ylim = c(0,8),col=rgb(0,0,1,0.1),pch=19,cex=0.5,xlab='BIO-JQ1',ylab='H3K27ac')
abline(lm(m[,4]~m[,1]))

cor(m[,1],m[,4])
dev.off()


#=========================================================================
#===================DELTA BRD4 as a function of chemseq===================
#=========================================================================
pdf(file='./figures/140924_tnbc_jq1_delta_scatter.pdf',width = 8,height =4)

par(mfrow=c(1,2))
plot(m[,1],log2(m[,3]/m[,2]),ylim = c(-3,1),xlim = c(0,2),col=rgb(1,0,0,0.1),pch=19,cex=0.5,xlab='BIO-JQ1',ylab='Log2 fold BRD4 +/- JQ1')
#abline(lm(log2(m[,3]/m[,2])~m[,1]))
lines(loess.smooth(m[,1],log2(m[,3]/m[,2])),col='black',lwd=2)
abline(h=0,lty=2)
cor(m[,1],log2(m[,3]/m[,2]))



plot(m[,1],log2(m[,5]/m[,4]),ylim = c(-3,3),xlim = c(0,2),col=rgb(0,0,1,0.1),pch=19,cex=0.5,xlab='BIO-JQ1',ylab='Log2 fold H3K27ac +/- JQ1')
#abline(lm(log2(m[,3]/m[,2])~m[,1]))
lines(loess.smooth(m[,1],log2(m[,5]/m[,4])),col='black',lwd=2)
abline(h=0,lty=2)

cor(m[,1],log2(m[,5]/m[,4]))
dev.off()


#=========================================================================
#========================DELTA BRD4 BINNED APPROACH=======================
#=========================================================================


bioJQ1Order = order(m[,1])

nBin = 10
binSize = nrow(m)%/%nBin
bins = seq(1,nrow(m),binSize)


boxMatrixBRD4 = matrix(nrow=binSize,ncol=nBin)
boxMatrixH3K27AC = matrix(nrow=binSize,ncol=nBin)
boxMatrixJQ1 = matrix(nrow=binSize,ncol=nBin)

for(i in 1:(nBin)){
	
	binRows = seq(bins[i],bins[(i+1)]-1,1)
	boxMatrixBRD4[,i] = log2(m[bioJQ1Order[binRows],3]/m[bioJQ1Order[binRows],2])
	boxMatrixH3K27AC[,i] = log2(m[bioJQ1Order[binRows],5]/m[bioJQ1Order[binRows],4])
	boxMatrixJQ1[,i] = m[bioJQ1Order[binRows],1]
	
	
}

pdf(file='./figures//140924_tnbc_jq1_delta_bin_box.pdf',width = 8,height =4)

par(mfrow=c(1,2))
boxplot(boxMatrixBRD4,cex=0,col='red',ylim =c(-3.5,0.5),ylab='Log2 fold change in BRD4',xlab='Regions ranked by increasing JQ1 binding')
boxplot(boxMatrixH3K27AC,cex=0,col='blue',ylim =c(-3,2),ylab= 'Log2 fold change in H3K27ac',xlab='Regions ranked by increasing JQ1 binding')

dev.off()
#=========================================================================
#========================DELTA BRD4 RESAMPLE CI APPROACH==================
#=========================================================================
brd4MeanVector = c()
brd4CIUpper = c()
brd4CILower = c()
nSample = 1000

for(i in 1:nBin){

	sampleMatrix = matrix(nrow=binSize,ncol=nSample)
	for(j in 1:nSample){
		
		sampleMatrix[,j] = sample(boxMatrixBRD4[,i],binSize,TRUE)
		
	}
	
	meanVector = apply(sampleMatrix,2,mean)
	brd4MeanVector = c(brd4MeanVector,mean(meanVector))
	brd4CILower = c(brd4CILower,mean(meanVector)-quantile(meanVector,c(0.025)))
	brd4CIUpper = c(brd4CIUpper,quantile(meanVector,c(0.975))-mean(meanVector))
	
		
}


h3k27MeanVector = c()
h3k27CIUpper = c()
h3k27CILower = c()

for(i in 1:nBin){

	sampleMatrix = matrix(nrow=binSize,ncol=nSample)
	for(j in 1:nSample){
		
		sampleMatrix[,j] = sample(boxMatrixH3K27AC[,i],binSize,TRUE)
		
	}
	
	meanVector = apply(sampleMatrix,2,mean)
	h3k27MeanVector = c(h3k27MeanVector,mean(meanVector))
	h3k27CILower = c(h3k27CILower,mean(meanVector)-quantile(meanVector,c(0.025)))
	h3k27CIUpper = c(h3k27CIUpper,quantile(meanVector,c(0.975))-mean(meanVector))
	
		
}

par(mfrow=c(1,2))
xVector = apply(boxMatrixJQ1,2,mean)
plot(xVector,brd4MeanVector,type='b',col='red',xlim = c(0.4,1.3),ylim =c(-2.2,-1),xlab='JQ1 binding (rpm/bp)',ylab='Log2 change in BRD4',cex=0)
error.bar(xVector,brd4MeanVector,brd4CIUpper,brd4CILower,col='red')

xVector = apply(boxMatrixJQ1,2,mean)
plot(xVector,h3k27MeanVector,type='b',col='blue',xlim = c(0.4,1.3),ylim =c(-1,0.2),xlab='JQ1 binding (rpm/bp)',ylab='Log2 change in BRD4',cex=0)
error.bar(xVector,h3k27MeanVector,h3k27CIUpper,h3k27CILower,col='blue')


