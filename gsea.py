# IMPLEMENTATION OF GSEA (GENE SET ENRICHMENT ANALYSIS)
# A method for interpreting gene expression data by focusing on gene sets.

import numpy as np
from collections import OrderedDict

def gsea(expressionProfiles, geneSets):
    #The main function to run the GSEA algorithm.
    #INPUT: expressionProfiles: a .txt file with the gene expression data
    #       geneSets: a .txt file with gene sets (pathways)
    #OUTPUT: a .txt file called 'results.txt' with calculated results of the GSEA analysis, formatted as:
    #       pathway name | NES | p-value | FDR value | pathway comments
    
    #step 1: rank order the genes to produce gene list L, according to the correlation of their expression profiles with phenotype C.
    #        There are several methods possible; as in the paper, we use signal-to-noise ratio.
    
    #read the file with the data with expression profiles.
    readData=profiles_reader(expressionProfiles)
    genes=readData[0]
    expressionData=readData[1]
    classes=readData[2]
    
    #separate data by phenotype.
    
    classSet=set(classes) #here we assume the data is separated in two classes.This set should have 2 elements.
    class1=classSet.pop()
    class2=classSet.pop()
    class1index=[]
    class2index=[]
    for i in range (0, len(classes)):
        if classes[i]==class1:
            class1index.append(i)
        elif classes[i]==class2:
            class2index.append(i)

    #for each gene, we calculate signal to noise ratio. We store genes and their stn ratios in an ordered dictionary.
    
    ratios={}
    for i in range(0,len(genes)):
        class1observations=expressionData[i, class1index]
        class2observations=expressionData[i, class2index]
        r=signal_to_noise(class1observations, class2observations)
        ratios[genes[i]] = r
    orderedRatios=OrderedDict(sorted(ratios.items(), key=lambda t: t[1], reverse=True))
    L=list(orderedRatios.keys())
    R=list(orderedRatios.values())
    
    #step 2: for each set S we calculate ES.
    
    #read the file with the pathways data.
    
    read_pathways=pathways_reader(geneSets)
    pathwayNames=read_pathways[0]
    pathwayComments=read_pathways[1]
    geneLists=read_pathways[2]
    
    #for each pathway, we calculate the enrichment score. 
    
    ES=[]
    for i in range (0, len(geneLists)):
        S=geneLists[i]
        enrichmentScore=enrichment_score(S, L, R)
        ES.append(enrichmentScore)
    
    #step 3: estimating significance: creating a vector ESnull with ES of 1000 fixed permutations of phenotype
    ESnull=[]
    for i in range (1000):      
        print ('permutation '+ str(i))
        tempClasses=classes
        np.random.shuffle(tempClasses)
        class1index=[]
        class2index=[]
        for j in range (0, len(tempClasses)):
            if tempClasses[j]==class1:
                class1index.append(j)
            elif tempClasses[j]==class2:
                class2index.append(j)
        
        
        ratios={}
        for j in range(0,len(genes)):
            class1observations=expressionData[j, class1index]
            class2observations=expressionData[j, class2index]
            r=signal_to_noise(class1observations, class2observations)
            ratios[genes[j]] = r
        
        orderedRatios=OrderedDict(sorted(ratios.items(), key=lambda t: t[1], reverse=True))
        
        L=list(orderedRatios.keys())
        R=list(orderedRatios.values())       
        
        ESnull1=[]
        for S in geneLists:
            enrichmentScore=enrichment_score(S, L, R)
            ESnull1.append(enrichmentScore)
        ESnull.append(ESnull1)
        
        
    ESnull=np.asarray(ESnull)
    
    #calculate nES, pValues and FDR q values
    
    significance=estimate_significance(ES, ESnull)
    nES=significance[0]
    pValues=significance[1]
    FDRvalues=significance[2]
    
    #step4: create a results.txt file, which stores the calculated data. Each line has the form:
    #       pathway name | NES | p-value | FDR value | pathway comments
    #       The results are sorted according to descending NES score.
    
    indices=sorted(range(len(nES)), key=lambda k: nES[k], reverse=True)
    
    results = open('results.txt','w') 
    
    results.write('{: <73} | {: >10} | {: >10} | {:>10} | {} \n'.format('PATHWAY NAME', 'NES', 'p-VALUE', 'FDR q-VALUE', 'PATHWAY COMMENTS \n'))
    for i in range (0, len(indices)):
        index=indices[i]
        results.write('{: <73} | {: >10.6f} | {: >10.6f} | {:>10.6f} | {} \n'.format(pathwayNames[index], nES[index], pValues[index], FDRvalues[index], pathwayComments[index]))
    results.close() 
        
    

def profiles_reader(expression_profiles):
    # Reads the file with gene expression profiles.
    # INPUT: expression_profiles: a .txt file, data is separated by '\t'.
    #                             First row: header, list of patients. We denote the number of patients by k.
    #                             Each consecutive row has a gene name in the first column and gene expression data in all following columns. We denote the number of genes by N.
    # OUTPUT: a list with three items:
    #         genes:              list (of length N) of all the genes.
    #         expression_data:    matrix (of dimension N x k) of all the gene expression data.
    #         classes:            list (of length k) of phenotype classes
    
    file=expression_profiles

    classes=open(file).readlines()[0].strip('\n').split('\t')[1::]
    NumCols = len(classes)+1

    NumRows = 0
    data = []
    genes = []

    with open(file) as f:
        for line in f.readlines()[1::]:
            NumRows = NumRows + 1
            line = line.split('\t')
            datarow = []
            for i in range(0, NumCols):
                if i == 0:
                    genes.append(line[i])
                else:
                    datarow.append(float(line[i]))
            data.append(datarow)
                
        expression_data = np.asarray(data)
        
    return (genes, expression_data, classes)
    

def pathways_reader(gene_sets):
    # Reads the file with sets of genes (pathways).
    # INPUT: expression_profiles: a .txt file, data is separated by '\t'.
    #                             Each row represets a pathway. It is of the form:
    #                             'pathway name \t pathway comment \t gene1 \t gene 2...'
    # OUTPUT: a list with three items
    #         pathway_names:      list of all the pathway names.
    #         pathway_comments:   list of pathway comments. 
    #         gene_lists:         list of lists: for each pathway, it contains a list of included genes.
    
    file=gene_sets
    
    
    pathway_names=[]
    pathway_comments=[]
    gene_lists=[]
    
    with open(file) as f:
        for line in f.readlines():
            line=line.split('\t')
            pathway_names.append(line[0])
            pathway_comments.append(line[1])
            gene_lists.append(line[2::])
    
    return [pathway_names, pathway_comments, gene_lists]
            
def signal_to_noise (class1, class2):
    #calculates correlation (signal-to-noise ratio).
    #INPUT : class1: list of values for the first class (first phenotype)
    #        class2: list of values for the second class (second phenotype)
    #OUTPUT: r:      singnal-to-noise ratio.
    r=(np.mean(class1)-np.mean(class2))/(np.std(class1)+np.std(class2))
    
    return r

def enrichment_score (S, genelist, ratios):
    #calculates the enrichment score for the set S
    #INPUT: S: list of genes (a pathway)
    #       genelist: list of genes we are interested in, ordered according to correlation with phenotype 
    #       ratios: list of correlation values
    #OUTPUT: ES: enrichment score
    
    p=1 #exponent to control the weigth; used 1 as in the paper.
    
    #check which elements of genelist are in S. We create a vector indicator as following:
    #indicator[i]=1 if genelist[i] in S, indicator[i]=0 if genelist[i] not in S
    
    indicator=[]
    for gene in genelist:
        if gene in S:
            indicator.append(1)
        else: indicator.append(0)
    
    NR=((np.abs(np.asarray(ratios)))**p).dot(np.asarray(indicator))
    
    
    #calculate the max and the min of the running sum

    runningSum=0
    minSum=0
    maxSum=0
    for i in range (0, len(indicator)):
        if indicator[i]==1:
            runningSum=runningSum + ((abs(ratios[i]))**p)/NR
            if runningSum > maxSum: maxSum=runningSum
        else:
            runningSum=runningSum-1/(len(genelist)-len(S))
            if runningSum < minSum: minSum=runningSum

            

    if max(abs(minSum), maxSum)==maxSum:
        ES=maxSum
    else: ES=minSum
   
    return ES

def normalize(ES, ESnull):
    #normalizes the ES and ESnull 
    #INPUT: ES: a vector of enrichment scores
    #       ESnull: an array of enrichment scores, calculated for all the sets with 1000 permutations of class labels
    #OUTPUT: a list with two values:
    #        NES: normalized vector of ES enrichment scores
    #        nESnull: normalized array ESnull of enrichment scores
    
    NES=[]
    nESnull=[]
    for i in range (0, len(ES)):
        negatives=[a for a in ESnull[::,i] if a <= 0]
        positives=[a for a in ESnull[::,i] if a>=0]
        negativesMean=np.mean(negatives)
        positivesMean=np.mean(positives)
        if ES[i] <=0:
            NESi=-ES[i]/negativesMean
        elif ES[i] >=0:
            NESi=ES[i]/positivesMean
        NES.append(NESi)
        
        normalizedNulls=[]
        for j in range (0, np.shape(ESnull)[0]):
            if ESnull[j,i] <= 0:
                a = - ESnull[j,i]/negativesMean
                normalizedNulls.append(a)
            elif ESnull[j,i] >=0:
                a= ESnull[j,i]/positivesMean
                normalizedNulls.append(a)
        nESnull.append(normalizedNulls)
    
    nESnull=np.asarray(nESnull).T
    return [NES, nESnull]

def estimate_significance (ES, ESnull):
    #Calculates the statistics for estimating significance.
    #INPUT: ES: a vector of enrichment scores, calculated for all sets(pathways)
    #       ESnull: an array of enrichment scores, calculated for all the sets(pathways) with 1000 permutations of class labels
    #OUTPUT:a list with three items:
    #       nES: a vector of normalized enrichment scores
    #       pValues: a vector of nominal p-values
    #       FDRValues: a vector of FDR q-values
    
    #p-value for each pathway
    pValues=[]
    for i in range (0, len(ES)):
        if ES[i] <=0:
            negatives=[x for x in ESnull[::, i] if x <= 0]
            pVal=(len([x for x in negatives if x <= ES[i]]))/(len(negatives))
            pValues.append(pVal)
        elif ES[i] >=0:
            positives=[x for x in ESnull[::, i] if x >= 0]
            pVal=(len([x for x in positives if x >= ES[i]]))/(len(positives))
            pValues.append(pVal)
            
    # FDR q-value
    normalized=normalize(ES, ESnull)
    nES=normalized[0]
    nESnull=normalized[1]
    
    FDRValues=[]
    flattenedNESnull=nESnull.flatten()
    for i in range (0, len(nES)):
        if nES[i] >=0:
            allPositive=len([x for x in flattenedNESnull if x >=0])
            greaterPositive=len([x for x in flattenedNESnull if x>=nES[i]])
            positiveNES=len([x for x in nES if x >=0])
            greaterNES=len([x for x in nES if x >= nES[i]])
            FDR=(greaterPositive/allPositive)/(greaterNES/positiveNES)
        else:
            allNegative=len([x for x in flattenedNESnull if x <=0])
            smallerNegative=len([x for x in flattenedNESnull if x<=nES[i]])
            negativeNES=len([x for x in nES if x <=0])
            smallerNES=len([x for x in nES if x <= nES[i]])
            FDR=(smallerNegative/allNegative)/(smallerNES/negativeNES)
        FDRValues.append(FDR)
        
    return [nES, pValues, FDRValues]