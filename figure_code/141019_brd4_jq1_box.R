#141016_brd4_jq1_box.R


setwd('./TNBC/')



#========================================================
#==================DATA INPUT============================
#========================================================


sum159_brd4_dmso = read.delim('mappedFolder/n1/HG19_ENRICHED_SUM159_JQ1_CHEM_SEQ_R1_-500_+500/HG19_ENRICHED_SUM159_JQ1_CHEM_SEQ_R1_-500_+500_SUM159_BRD4_DMSO_R1.gff')
sum159_brd4_jq1 = read.delim('mappedFolder/n1/HG19_ENRICHED_SUM159_JQ1_CHEM_SEQ_R1_-500_+500/HG19_ENRICHED_SUM159_JQ1_CHEM_SEQ_R1_-500_+500_SUM159_BRD4_JQ1_R1.gff')

sum159R_brd4_dmso = read.delim('mappedFolder/n1/HG19_ENRICHED_SUM159R_JQ1_CHEM_SEQ_R1_-500_+500/HG19_ENRICHED_SUM159R_JQ1_CHEM_SEQ_R1_-500_+500_SUM159R_BRD4_DMSO.gff')
sum159R_brd4_jq1 = read.delim('mappedFolder/n1/HG19_ENRICHED_SUM159R_JQ1_CHEM_SEQ_R1_-500_+500/HG19_ENRICHED_SUM159R_JQ1_CHEM_SEQ_R1_-500_+500_SUM159R_BRD4_JQ1_R1.gff')

#=======
sum159_h3k27ac_dmso = read.delim('mappedFolder/n1/HG19_ENRICHED_SUM159_JQ1_CHEM_SEQ_R1_-500_+500/HG19_ENRICHED_SUM159_JQ1_CHEM_SEQ_R1_-500_+500_SUM159_H3K27AC_DMSO_R1.gff')
sum159_h3k27ac_jq1 = read.delim('mappedFolder/n1/HG19_ENRICHED_SUM159_JQ1_CHEM_SEQ_R1_-500_+500/HG19_ENRICHED_SUM159_JQ1_CHEM_SEQ_R1_-500_+500_SUM159_H3K27AC_JQ1_R1.gff')


sum159R_h3k27ac_dmso = read.delim('mappedFolder/n1/HG19_ENRICHED_SUM159R_JQ1_CHEM_SEQ_R1_-500_+500/HG19_ENRICHED_SUM159R_JQ1_CHEM_SEQ_R1_-500_+500_SUM159R_H3K27AC_DMSO_R1.gff')
sum159R_h3k27ac_jq1 = read.delim('mappedFolder/n1/HG19_ENRICHED_SUM159R_JQ1_CHEM_SEQ_R1_-500_+500/HG19_ENRICHED_SUM159R_JQ1_CHEM_SEQ_R1_-500_+500_SUM159R_H3K27AC_JQ1_R1.gff')


#========================================================
#==================FOLD BOX PLOTS========================
#========================================================
#fold boxplots
sum159_brd4_fold = log2(sum159_brd4_jq1[,3]/sum159_brd4_dmso[,3])
sum159R_brd4_fold = log2(sum159R_brd4_jq1[,3]/sum159R_brd4_dmso[,3])


sum159_h3k27ac_fold = log2(sum159_h3k27ac_jq1[,3]/sum159_h3k27ac_dmso[,3])
sum159R_h3k27ac_fold = log2(sum159R_h3k27ac_jq1[,3]/sum159R_h3k27ac_dmso[,3])

pdf(file='figures/141019_brd4_h3k27ac_jq1_box.pdf',width = 4,height = 8)
par(mfrow=c(2,1))
boxplot(sum159_brd4_fold,sum159R_brd4_fold,cex=0,ylim =c(-3,2),names = c('SUM159','SUM159R'),ylab='Log2 fold change in BRD4 +/- JQ1')
abline(h=0)
boxplot(sum159_h3k27ac_fold,sum159R_h3k27ac_fold,cex=0,ylim =c(-3,2),names= c('SUM159','SUM159R'),ylab='Log2 fold change in H3K27AC +/- JQ1')
abline(h=0)
dev.off()



sum159Matrix = cbind(sum159_brd4_dmso[,3],sum159_brd4_jq1[,3])


sum159RMatrix = cbind(sum159R_brd4_dmso[,3],sum159R_brd4_jq1[,3])


par(mfrow=c(1,2))

boxplot(sum159Matrix,cex=0,ylim =c(0,2))
boxplot(sum159RMatrix,cex=0,ylim =c(0,2))

