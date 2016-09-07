#!/usr/bin/python
#pipeline_template.py

'''
The MIT License (MIT)

Copyright (c) 2013 Charles Lin

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
'''

#generic pipeline template for human data


#==========================================================================
#=============================DEPENDENCIES=================================
#==========================================================================

#need to clone the bradner lab pipeline scripts into the same folder
pipeline_dir = './pipeline/' #set location of pipeline_dfci.py and utils.py

import sys
sys.path.append(pipeline_dir)

import pipeline_dfci
import datetime
import random
import string
import os
import utils
#==========================================================================
#============================PARAMETERS====================================
#==========================================================================

#project folders
projectFolder = './' #PATH TO YOUR REPO

dataFile = './CHIP_DATA_TABLE.txt'
brd4_dataFile ='./BRD4_DATA_TABLE.txt'

sum159_rna_dataFile = './SUM159_RNA_TABLE.txt'
sum149_rna_dataFile = './SUM149_RNA_TABLE.txt'


genome ='hg19'
annotFile = './pipeline/annotation/hg19_refseq.ucsc'


#standard folder names
gffFolder ='%sgff/' % (projectFolder)
macsFolder = '%smacsFolder/' % (projectFolder)
macsEnrichedFolder = '%smacsEnriched/' % (projectFolder)
mappedEnrichedFolder = '%smappedEnriched/' % (projectFolder)
mappedFolder = '%smappedFolder/' % (projectFolder)
wiggleFolder = '%swiggles/' % (projectFolder)
metaFolder = '%smeta/' % (projectFolder)

#making folders
folderList = [gffFolder,macsFolder,macsEnrichedFolder,mappedEnrichedFolder,mappedFolder,wiggleFolder,metaFolder]

for folder in folderList:
    pipeline_dfci.formatFolder(folder,True)


#a list of good BRD4 datasets 
goodList = ['SUM159_BRD4_DMSO_3h','SUM159_BRD4_DMSO','SUM159_BRD4_DMSO_R1','SUM159_JQ1_CHEM_SEQ_R1','SUM159R_BRD4_DMSO_3h','SUM159R_BRD4_DMSO','SUM159R_JQ1_CHEM_SEQ_R1', 'SUM159R_BRD4_JQ1_3h','SUM159R_BRD4_JQ1_R1']



#==========================================================================
#=======================LOADING DATA ANNOTATION============================
#==========================================================================

##THIS SECTION LOADS A DATA TABLE.  MUST BE UNCOMMENTED FOR REST OF CODE TO WORK


#LOADING THE DATA TABLE
dataDict = pipeline_dfci.loadDataTable(dataFile)
print(dataDict.keys())

pipeline_dfci.summary(dataFile)

#==========================================================================
#=========================READ ME ON HOW TO USE============================
#==========================================================================


'''
tnbc_pipeline_scripts.py contains all python scripts for data analysis in Shu et al., 2016

Many of the functions and scripts here refer to functions and code from the bradner lab pipeline
https://github.com/BradnerLab/pipeline


Each section represents a block of code to be run sequentially in order to perform
analysis 
Recommend running one block at a time by uncommenting and running the code in python

DO NOT COMMENT ANYTHING ABOVE THIS BLOCK

PREREQUISITES:

1. Clone this repo and cd into the TNBC folder
-all paths are set relative to this folder
-recommended that you cd into the TNBC folder to run all analysis

2. please install bamliquidator https://github.com/BradnerLab/pipeline/wiki/bamliquidator
-set bamliquidator path globally as bamliquidator 
-set bamliquidator_batch path globally as bamliquidator_batch
-install macs1.4.2 and add to global path as macs14 (http://liulab.dfci.harvard.edu/MACS/Download.html)
-install bowtie2, samtools, cufflinks, cuffquant, R

3. Clone the bradnerlab pipeline https://github.com/BradnerLab/pipeline
-edit section of paths in pipeline_dfci.py to correctly point to the pipeline directory, samtools, bamliquidator, and bamliquidator_batch

4. Data tables
-A number of data tables are included to help organize files for analysis
-Each row represents a specific dataset. The UNIQUE_ID and NAME  corresponds to the sample name in GEO for Series GSE63584
-The BACKGROUND column represents the control for a given dataset
-ENRICHED_REGION and ENRICHED_MACS refer to locations of peak files generated by various peak calling algorithms
-For more information consult the loadDataTable function in the pipeline_dfci.py module

5. Obtaining raw data
-raw fastqs may be obtained here: http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE63584
-recommended that fastqs be deposited into the ./fastq folder

6. Generate sorted and indexed bam files for each sample.
-align raw fastqs to the hg19 genome. bowtie2 index can be found here: ftp://igenome:G3nom3s4u@ussd-ftp.illumina.com/Homo_sapiens/UCSC/hg19/Homo_sapiens_UCSC_hg19.tar.gz
-for alignment we used bowtie2 v2.2.1 with default parameters except for -k 1
-place bams in the ./bam folder or edit the "FILE_PATH" column in CHIP_DATA_TABLE.txt to point to the folder containing bams
-bams should be named to include the UNIQUE_ID as in the CHIP_DATA_TABLE.txt. The UNIQUE_ID is derived from the GEO sample name.
-bams shouled be sorted and indexed and include a corresponding .bai file with the same name
-NOTE: code will not work if multiple bams in the directory share the same UNIQUE_ID

7. Running python analysis
-cd into the repo directory
-Run each block of code by uncommenting and then running python ./tnbc_pipeline_scripts.py
-All paths are set relativeto the TNBC folder 

8. Running R scripts
-all R scripts are found in ./figure_code
-additional R scripts are used for expression analysis and to generate plots of ChIP-Seq data.
-output figures can be found in ./figures
'''






#==========================================================================
#========================CALL MACS FOR SUM159==============================
#==========================================================================

##THIS SECTION CALLS THE MACS ERROR MODEL


#namesList = dataDict.keys()

#print(namesList)
#pipeline_dfci.callMacs(dataFile,macsFolder,namesList,overwrite=False,pvalue='1e-9',useBackground=True)


#==========================================================================
#==================FORMAT MACS OUTPUT FOR SUM159===========================
#==========================================================================

##THIS SECTION FORMATS THE OUTPUT FROM MACS, CREATES THE MACSENRICHED FOLDER AND MOVES WIGGLES TO THE DESTINATION

#pipeline_dfci.formatMacsOutput(dataFile,macsFolder,macsEnrichedFolder,wiggleFolder,wigLink='/ark/wiggles/',useBackground=False)


#==========================================================================
#==============FIGURE 1F RUNNING ROSE WITH AUTOSTITCH======================
#==========================================================================

##just on the jq1 chemseq

# namesList = ['SUM159_JQ1_CHEM_SEQ_R1']
# print(namesList)

# parentFolder = '%srose_stitch' % (projectFolder)
# parentFolder = utils.formatFolder(parentFolder,True)

# maskFile ='./tables/hg19_encode_blacklist.bed' 
# bashFileName = '%srose_stitch/tnbc_brd4_rose.sh' %(projectFolder)

# pipeline_dfci.callRose2(dataFile,macsEnrichedFolder,parentFolder,namesList,[],'',2500,'',bashFileName,maskFile)




#==========================================================================
#===================FIGURE 1D MAKING GFFS FOR HEATMAP======================
#==========================================================================

# #+/- 10kb from center of sum159 jq1 chemseq enhancers
# enhancerGFFFile = '%srose_stitch/SUM159_JQ1_CHEM_SEQ_R1_ROSE/gff/SUM159_JQ1_CHEM_SEQ_R1_peaks_8KB_STITCHED_TSS_DISTAL.gff' % (projectFolder)

# enhancerGFF = utils.parseTable(enhancerGFFFile,'\t')
# newGFF= []
# for line in enhancerGFF:
#     start = min([int(line[3]),int(line[4])])
#     stop  = max([int(line[3]),int(line[4])])
#     center = (start+stop)/2
#     newLine = list(line)
#     newLine[3] = center -10000
#     newLine[4] = center + 10000
#     newGFF.append(newLine)

# utils.unParseTable(newGFF,'%sHG19_SUM159_JQ1_CHEM_SEQ_R1_peaks_8KB_STITCHED_TSS_DISTAL_-10KB_+10KB_CENTER.gff' % (gffFolder),'\t')


# #+/- 10kb from all TSS

# tssGFF = utils.parseTable('%sHG19_TSS_ALL_-5000_+5000.gff' % (gffFolder),'\t')

# newGFF = []
# for line in tssGFF:
#     start = min([int(line[3]),int(line[4])])
#     stop  = max([int(line[3]),int(line[4])])
#     newLine = list(line)
#     newLine[3] = start -5000
#     newLine[4] = stop + 5000
#     newGFF.append(newLine)

# utils.unParseTable(newGFF,'%sHG19_TSS_ALL_-10KB_+10KB.gff' % (gffFolder),'\t')



# # +/- 500bp from the chemseq regions

# namesList = ['SUM159_JQ1_CHEM_SEQ_R1','SUM159R_JQ1_CHEM_SEQ_R1']

# pipeline_dfci.makeEnrichedGFFs(dataFile,namesList,gffFolder,macsEnrichedFolder,macs=True,window=500)



#==========================================================================
#===========================FIGURE 1D HEATMAPS=============================
#==========================================================================

# # tssGFFFile = '%sHG19_TSS_ALL_-10KB_+10KB.gff' % (gffFolder)
# # enhancerGFFFile = '%sHG19_SUM159_JQ1_CHEM_SEQ_R1_peaks_8KB_STITCHED_TSS_DISTAL_-10KB_+10KB_CENTER.gff' % (gffFolder)
# cellTypeList = ['SUM159','SUM159R']
# gffList = [tssGFFFile,enhancerGFFFile]


# #n=200
# mappedFolder = utils.formatFolder('%smappedFolder/n200/' % (projectFolder),True)
# pipeline_dfci.mapBams(dataFile,cellTypeList,gffList,mappedFolder,nBin = 200,overWrite =False,nameList = goodDataList)


# #for TSS

# tssGFFFile = '%sHG19_TSS_ALL_-10KB_+10KB.gff' % (gffFolder)
# orderByName = 'SUM159_JQ1_CHEM_SEQ_R1'
# outputFolder = utils.formatFolder('%sheatmap' % (projectFolder),True)
# mappedFolder = utils.formatFolder('%smappedFolder/n200/' % (projectFolder),True)
# pipeline_dfci.callHeatPlotOrdered(dataFile,tssGFFFile,goodDataList,orderByName,'',outputFolder,mappedFolder,relative=True)



# #for ENHANCER
# enhancerGFFFile = '%sHG19_SUM159_JQ1_CHEM_SEQ_R1_peaks_8KB_STITCHED_TSS_DISTAL_-10KB_+10KB_CENTER.gff' % (gffFolder)
# orderByName = 'SUM159_JQ1_CHEM_SEQ_R1'
# outputFolder = utils.formatFolder('%sheatmap' % (projectFolder),True)
# mappedFolder = utils.formatFolder('%smappedFolder/n200/' % (projectFolder),True)
# pipeline_dfci.callHeatPlotOrdered(dataFile,enhancerGFFFile,goodDataList,orderByName,'',outputFolder,mappedFolder,relative=True)



#==========================================================================
#========================FIGURE 1E HIF1A PLOT==============================
#==========================================================================

# plotName = 'HIF1A_SE'
# geneID = 'NM_001530'
# regionString = 'chr14:+:61922788-62213660'
# outputFolder = utils.formatFolder('%sgenePlot/' % (projectFolder),True)
# plotList = ['SUM159_JQ1_CHEM_SEQ_R1','SUM159_BRD4_DMSO_R1','SUM159_BRD4_JQ1_R1']


# pipeline_dfci.callGenePlot(dataFile,geneID,plotName,annotFile,plotList,outputFolder,region=regionString,yScale = 'UNIFORM')




#==========================================================================
#=========EXTENDED DATA FIGURE 3 MAPPING BAMS AT CHEMSEQ REGIONS===========
#==========================================================================


# chemseqGFFFile = '%sHG19_ENRICHED_SUM159_JQ1_CHEM_SEQ_R1_-500_+500.gff' % (gffFolder)
# chemseqRGFFFile = '%sHG19_ENRICHED_SUM159R_JQ1_CHEM_SEQ_R1_-500_+500.gff' % (gffFolder)
# gffList = [chemseqGFFFile,chemseqRGFFFile]
# cellTypeList = ['SUM159','SUM159R']

# mappedFolder = utils.formatFolder('%smappedFolder/n1/' % (projectFolder),True)

# #n=1
# pipeline_dfci.mapBams(dataFile,cellTypeList,gffList,mappedFolder,nBin = 1,overWrite =False,nameList = [])






#==========================================================================
#========================FIGURE2C DYNAMIC ENHANCER=========================
#==========================================================================


# #analysis name
# analysisName = 'TNBC_CHEMSEQ'

# #dataset names
# name1 = 'SUM159_JQ1_CHEM_SEQ_R1'
# name2=  'SUM159R_JQ1_CHEM_SEQ_R1'

# namesString = '%s,%s' % (name1,name2)

# roseParentFolder = '%srose_stitch/' % (projectFolder)

# roseString = '%s%s_ROSE/,%s%s_ROSE/' % (roseParentFolder,name1,roseParentFolder,name2)


# #first set the output folder
# outputFolder = '%sdynamic/%s/' % (projectFolder,analysisName)
# outputFolder = utils.formatFolder(outputFolder,True)


# bashFileName = '%s%s_dynamic.sh' % (outputFolder,analysisName)
# bashFile = open(bashFileName,'w')

# bashFile.write('#!/usr/bin/bash\n')

# #for now change into pipelinedir just to be safe
# bashFile.write('cd %s\n' % pipeline_dir)

# dynamicCmd = 'python %sdynamicEnhancer.py -g hg19 -d %s -n %s -r %s -o %s -e super -m -p' % (pipeline_dir,dataFile,namesString,roseString,outputFolder)

# print dynamicCmd

# bashFile.write(dynamicCmd +'\n')

# bashFile.close()

# print("wrote commands to %s" % (bashFileName))


#==========================================================================
#====================FIGURE 2D MAPPING BAMS TO CHEMSEQ REGIONS=============
#==========================================================================

# #first make a gff of the ordered chemseq regions

# rankTable = utils.parseTable('%sdynamic/TNBC_CHEMSEQ/output/HG19_SUM159_JQ1_CHEM_SEQ_R1_SUM159R_JQ1_CHEM_SEQ_R1_merged_MERGED_SUPERS_RANK_TABLE.txt' % (projectFolder),'\t')

# rankGFF = []
# for line in rankTable[1:]:

#     gffLine = [line[1],line[0],line[0],line[2],line[3],'','.',line[6],line[0]]
#     rankGFF.append(gffLine)

# utils.unParseTable(rankGFF,'%sHG19_SUM159_JQ1_CHEM_SEQ_R1_SUM159R_JQ1_CHEM_SEQ_R1_merged_MERGED_SUPERS_RANK.gff' % (gffFolder),'\t')

# rankGFFFile = '%sHG19_SUM159_JQ1_CHEM_SEQ_R1_SUM159R_JQ1_CHEM_SEQ_R1_merged_MERGED_SUPERS_RANK.gff' % (gffFolder)


# gffList = [rankGFFFile]

# cellTypeList = ['SUM159','SUM159R']

# mappedFolder = utils.formatFolder('%smappedFolder/n1/' % (projectFolder),True)

# #n=1
# pipeline_dfci.mapBams(dataFile,cellTypeList,gffList,mappedFolder,nBin = 1,overWrite =False,nameList = goodDataList)




#==========================================================================
#=======================FIGURE 2F BCL-xL PLOTTING==========================
#==========================================================================


# plotName = 'BCL2L1_CHEMSEQ'
# geneID = 'NM_001191'
# regionString = 'chr20:-:30242478-30322514'
# outputFolder = utils.formatFolder('%sgenePlot/' % (projectFolder),True)
# plotList = ['SUM159_JQ1_CHEM_SEQ_R1','SUM159R_JQ1_CHEM_SEQ_R1']


# pipeline_dfci.callGenePlot(dataFile,geneID,plotName,annotFile,plotList,outputFolder,region=regionString,yScale = 'UNIFORM')


# plotName = 'BCL2L1_BRD4'
# geneID = 'NM_001191'
# regionString = 'chr20:-:30242478-30322514'
# outputFolder = utils.formatFolder('%sgenePlot/' % (projectFolder),True)
# plotList = ['SUM159_BRD4_DMSO_R1','SUM159R_BRD4_DMSO']


# pipeline_dfci.callGenePlot(dataFile,geneID,plotName,annotFile,plotList,outputFolder,region=regionString,yScale = 'UNIFORM')



# plotName = 'BCL2L1_H3K27AC'
# geneID = 'NM_001191'
# regionString = 'chr20:-:30242478-30322514'
# outputFolder = utils.formatFolder('%sgenePlot/' % (projectFolder),True)
# plotList = ['SUM159_H3K27AC_DMSO_R1','SUM159R_H3K27AC_DMSO_R1']


# pipeline_dfci.callGenePlot(dataFile,geneID,plotName,annotFile,plotList,outputFolder,region=regionString,yScale = 'UNIFORM')


#==========================================================================
#========================FIGURE 2G SOD2 PLOTTING===========================
#==========================================================================



# #NM_000636
# #SOD2

# plotName = 'SOD2_SUM159_BRD4'
# geneID = 'NM_000636'
# regionString = 'chr6:-:160085100-160141427'
# outputFolder = utils.formatFolder('%sgenePlot/' % (projectFolder),True)
# plotList = ['SUM159_BRD4_DMSO_R1','SUM159_BRD4_JQ1_R1']


# pipeline_dfci.callGenePlot(dataFile,geneID,plotName,annotFile,plotList,outputFolder,region=regionString,yScale = 'UNIFORM')


# plotName = 'SOD2_SUM159R_BRD4'
# geneID = 'NM_000636'
# regionString = 'chr6:-:160085100-160141427'
# outputFolder = utils.formatFolder('%sgenePlot/' % (projectFolder),True)
# plotList = ['SUM159R_BRD4_DMSO','SUM159R_BRD4_JQ1_R1']


# pipeline_dfci.callGenePlot(dataFile,geneID,plotName,annotFile,plotList,outputFolder,region=regionString,yScale = 'UNIFORM')



# plotName = 'SOD2_CHEMSEQ'
# geneID = 'NM_000636'
# regionString = 'chr6:-:160085100-160141427'
# outputFolder = utils.formatFolder('%sgenePlot/' % (projectFolder),True)
# plotList = ['SUM159_JQ1_CHEM_SEQ_R1','SUM159R_JQ1_CHEM_SEQ_R1']


# pipeline_dfci.callGenePlot(dataFile,geneID,plotName,annotFile,plotList,outputFolder,region=regionString,yScale = 'UNIFORM')




# plotName = 'SOD2_SUM159_H3K27AC'
# geneID = 'NM_000636'
# regionString = 'chr6:-:160085100-160141427'
# outputFolder = utils.formatFolder('%sgenePlot/' % (projectFolder),True)
# plotList = ['SUM159_H3K27AC_DMSO_R1','SUM159_H3K27AC_JQ1_R1']


# pipeline_dfci.callGenePlot(dataFile,geneID,plotName,annotFile,plotList,outputFolder,region=regionString,yScale = 'UNIFORM')


# plotName = 'SOD2_SUM159R_H3K27AC'
# geneID = 'NM_000636'
# regionString = 'chr6:-:160085100-160141427'
# outputFolder = utils.formatFolder('%sgenePlot/' % (projectFolder),True)
# plotList = ['SUM159R_H3K27AC_DMSO_R1','SUM159R_H3K27AC_JQ1_R1']


# pipeline_dfci.callGenePlot(dataFile,geneID,plotName,annotFile,plotList,outputFolder,region=regionString,yScale = 'UNIFORM')





#==========================================================================
#========================CALL MACS ACROSS TNBC=============================
#==========================================================================

##THIS SECTION CALLS THE MACS ERROR MODEL

#dataDict = pipeline_dfci.loadDataTable(brd4_dataFile)
#namesList = dataDict.keys()

#print(namesList)
#pipeline_dfci.callMacs(brd4_dataFile,macsFolder,namesList,overwrite=False,pvalue='1e-9',useBackground=True)


#==========================================================================
#==================FORMAT MACS OUTPUT FOR SUM159===========================
#==========================================================================

##THIS SECTION FORMATS THE OUTPUT FROM MACS, CREATES THE MACSENRICHED FOLDER AND MOVES WIGGLES TO THE DESTINATION

#pipeline_dfci.formatMacsOutput(brd4_dataFile,macsFolder,macsEnrichedFolder,wiggleFolder,wigLink='/ark/wiggles/',useBackground=False)


#==========================================================================
#==================EXTENDED DATA FIGURE 6C MAPPING=========================
#==========================================================================


#dataDict = pipeline_dfci.loadDataTable(brd4_dataFile)

#first make a gff for each peak

# namesList = [name for name in dataDict.keys() if name.count('BRD4_DMSO') == 1]


#pipeline_dfci.makeEnrichedGFFs(dataFile,namesList,gffFolder,macsEnrichedFolder,macs=True,window=0)

# #for each line, map signal to the set of peaks present in DMSO treated cells

# cellTypeList = [name.split('_')[0] for name in dataDict.keys()]
# cellTypeList = utils.uniquify(cellTypeList)

# print(cellTypeList)


# for cellType in cellTypeList:

#     namesList = [name for name in dataDict.keys() if name.count(cellType) == 1]

#     dmsoName = [name for name in namesList if name.count('BRD4_DMSO') == 1 and name.count('WCE') == 0][0]

#     enrichedGFFFile = '%sHG19_ENRICHED_%s_-0_+0.gff' % (gffFolder,dmsoName)
#     gffList = [enrichedGFFFile]
    
#     print(gffList)
#     print(namesList)

#     pipeline_dfci.mapBamsBatch(dataFile,gffList,mappedFolder,False,namesList,200)




#making signal table




# cellTypeList = [name.split('_')[0] + '_' for name in dataDict.keys()]
# cellTypeList = utils.uniquify(cellTypeList)
# print(cellTypeList)
# cellTypeList = ['SUM159R']
# for cellType in cellTypeList:

#     namesList = [name for name in dataDict.keys() if name.count(cellType) == 1]
#     namesList.sort()
#     dmsoName = [name for name in namesList if name.count('BRD4_DMSO') == 1 and name.count('WCE') == 0][0]

#     enrichedGFFFile = '%sHG19_ENRICHED_%s_-0_+0.gff' % (gffFolder,dmsoName)

    
#     print(enrichedGFFFile)
#     print(namesList)

#     output = '%stables/HG19_ENRICHED_%s_-0_+0_signal.txt' % (projectFolder,dmsoName)
#     print(output)

#     pipeline_dfci.makeSignalTable(dataFile,enrichedGFFFile,mappedFolder,namesList,False,output)


#==========================================================================
#=============EXTENDED DATA FIGURE 3D RUNNING ROSE ON SUM149===============
#==========================================================================

# #all samples at baseline

# dataDict = pipeline_dfci.loadDataTable(brd4_dataFile)
# namesList = ['SUM149_BRD4_DMSO']

# parentFolder = '%srose_stitch' % (projectFolder)
# parentFolder = utils.formatFolder(parentFolder,True)

# maskFile ='./tables/hg19_encode_blacklist.bed' 
# bashFileName = '%srose/sum149_brd4_rose.sh' %(projectFolder)

# pipeline_dfci.callRose2(dataFile,macsEnrichedFolder,parentFolder,namesList,[],'',2500,'',bashFileName,maskFile)





#==========================================================================
#==================EXTENDED DATA FIGURE 6A,B,D PLOTTING====================
#==========================================================================


# #NM_000636
# #SOD2


# figureGFF = [['chr6','SOD2','SOD2',160085100,160141427,'','-','','SOD2'],
#              ['chr20','BCL2L1','BCL2L1',30242478,30322514,'','-','','BCL2L1'],
#              ['chr14','HIF1A','HIF1A',61922788,62213660,'','+','','HIF1A'],
#              ]


# figureGFFFile = '%sHG19_TNBC_FIGURE_GENES.gff' % (gffFolder)
# utils.unParseTable(figureGFF,figureGFFFile,'\t')


# outputFolder = utils.formatFolder('%sgenePlot/' % (projectFolder),True)                        


# cellTypeList = [name.split('_')[0] + '_' for name in dataDict.keys()]
# cellTypeList = utils.uniquify(cellTypeList)
# print(cellTypeList)


# for cellType in cellTypeList:

#     namesList = [name for name in dataDict.keys() if name.count(cellType) == 1 and name.count('WCE') == 0 and name.count('BRD4') == 1]
#     namesList.sort()



#     plotName = 'TNBC_FIGURE_%sBRD4' % (cellType)
#     outputFolder = '%sgenePlot' % (projectFolder)

#     print(namesList)
#     print(plotName)
#     pipeline_dfci.callBatchPlot(dataFile,figureGFFFile,plotName,outputFolder,namesList,uniform=True,bed ='',plotType= 'MULTIPLE',extension=200,multiPage = False,debug=False,nameString = '')


#================================


# figureGFFFile = '%sHG19_TNBC_FIGURE_GENES.gff' % (gffFolder)

# outputFolder = utils.formatFolder('%sgenePlot/' % (projectFolder),True)                        

# namesList = [name for name in dataDict.keys() if name.count('H3K27AC') == 1]

# plotName = 'TNBC_FIGURE_H3K27AC'
# outputFolder = '%sgenePlot' % (projectFolder)

# print(namesList)
# print(plotName)
# pipeline_dfci.callBatchPlot(dataFile,figureGFFFile,plotName,outputFolder,namesList,uniform=True,bed ='',plotType= 'MULTIPLE',extension=200,multiPage = False,debug=False,nameString = '')





#================================


# figureGFFFile = '%sHG19_TNBC_FIGURE_GENES.gff' % (gffFolder)

# outputFolder = utils.formatFolder('%sgenePlot/' % (projectFolder),True)                        

# namesList = [name for name in dataDict.keys() if name.count('H3K27AC') == 1 and name.count('SUM149') == 1]

# plotName = 'TNBC_FIGURE_SUM149_H3K27AC'
# outputFolder = '%sgenePlot' % (projectFolder)

# print(namesList)
# print(plotName)
# pipeline_dfci.callBatchPlot(dataFile,figureGFFFile,plotName,outputFolder,namesList,uniform=True,bed ='',plotType= 'MULTIPLE',extension=200,multiPage = False,debug=False,nameString = '')






#==========================================================================
#=============MAKING RNA TABLES FOR SUM149 AND SUM159======================
#==========================================================================


# #FOR SUM159
# cufflinksFolder = utils.formatFolder('%scufflinks/' % (projectFolder),True)
# bashFileName = '%sSUM159_CUFFLINKS.sh' % (cufflinksFolder)
# gtfFile = './tables/genes.gtf' #this should point to the hg19 gene gtf from igenomes ftp://igenome:G3nom3s4u@ussd-ftp.illumina.com/Homo_sapiens/UCSC/hg19/Homo_sapiens_UCSC_hg19.tar.gz
# analysisName = 'SUM149'

# dataTable = utils.parseTable(sum159_rna_dataFile,'\t')

# groupList = [['SUM159R_DMSO_24H_1','SUM159R_DMSO_24H_2'],
#              ['SUM159R_DMSO_3H_1','SUM159R_DMSO_3H_2'],
#              ['SUM159R_JQ1_24H_1','SUM159R_JQ1_24H_2'],
#              ['SUM159R_JQ1_3H_1','SUM159R_JQ1_3H_2'],
#              ['SUM159_DMSO_12H_1','SUM159_DMSO_12H_2'],
#              ['SUM159_DMSO_24H_1','SUM159_DMSO_24H_2'],
#              ['SUM159_DMSO_3H_1','SUM159_DMSO_3H_2'],
#              ['SUM159_JQ1_12H_1','SUM159_JQ1_12H_2'],
#              ['SUM159_JQ1_24H_1','SUM159_JQ1_24H_2'],
#              ['SUM159_JQ1_3H_1','SUM159_JQ1_3H_2'],
#              ]

# pipeline_dfci.makeCuffTable(sum159_rna_dataFile,analysisName,gtfFile,cufflinksFolder,groupList,bashFileName)




# #FOR SUM149
# cufflinksFolder = utils.formatFolder('%scufflinks/' % (projectFolder),True)
# bashFileName = '%sSUM149_CUFFLINKS.sh' % (cufflinksFolder)
# gtfFile = './tables/genes.gtf' #this should point to the hg19 gene gtf from igenomes ftp://igenome:G3nom3s4u@ussd-ftp.illumina.com/Homo_sapiens/UCSC/hg19/Homo_sapiens_UCSC_hg19.tar.gz
# analysisName = 'SUM149'

# dataTable = utils.parseTable(sum149_rna_dataFile,'\t')


# groupList = [['SUM149_DMSO_12H_R1','SUM149_DMSO_12H_R2'],
#              ['SUM149_JQ1_12H_R1','SUM149_JQ1_12H_R2'],
#              ['SUM149R_DMSO_12H_R1','SUM149R_DMSO_12H_R2'],
#              ['SUM149R_JQ1_12H_R1','SUM149R_JQ1_12H_R2'],
#              ]

# pipeline_dfci.makeCuffTable(sum149_rna_dataFile,analysisName,gtfFile,cufflinksFolder,groupList,bashFileName)


