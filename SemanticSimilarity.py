#!/usr/bin/env python
# coding: utf-8


import numpy

#Let's define the Genes with their classes.
geneA = ['GO_0016020', 'GO_0003677'];
geneB = ['GO_0016021'];
geneC = ['GO_0003677'];

#Builds a dictionary of owlclasses which contain arrays of superclasses
def buildClassDictionary ():
    owlClasses = [];
    string = "";
    stringArr = [];
    superClasses = [];
    f = open("./superclasses.tsv", "r");
    owlIndex = 0;
    for line in f:
        owlClasses.append(line[0 : line.index('\t')]);
        for i in range(line.index('\t') + 1, len(line)):
            if i == line.index('\t') + 1:
                stringArr.append(owlClasses[owlIndex]);
            if (line[i] != ',') and (line[i] != ' ') and (line[i] != '\n'):
                string += (line[i]);
            if line[i] == ',' or line[i] == '\n':
                stringArr.append(string);
                string = "";
        owlIndex += 1;
        superClasses.append(stringArr);
        stringArr = [];
    dictionary = {};
    for index in range(0, len(owlClasses)):
        dictionary[owlClasses[index]] = superClasses[index];
    return dictionary;
    
owlDict = buildClassDictionary();

#Do Jaccard Similarity measure on GO Term 1 vs GO Term 2
def jaccardSimilarity (term1, term2, geneDict):
    #Start union by getting lengths.
    union = len(geneDict[term1]) + len(geneDict[term2]);
    
	#Get Intersection and Union of term1 and term2
    intersection = 0;
    for element1 in geneDict[term1]:
        for element2 in geneDict[term2]:
            if element1 == element2:
				#add element for intersection if equal; remove duplicate of union set if equal.
                intersection += 1;
                union -= 1;
    
	#Return the cardinality of intersection / cardinality of union.
    return intersection / union;

#Create a list of Jaccard scores for more easy computation later
def jaccardList (gene1, gene2, geneDict):
    jaccardScores = [];
    for go_A in gene1:
        for go_B in gene2:
            jaccardScores.append(jaccardSimilarity(go_A, go_B, geneDict));
    return jaccardScores;

#Create a dictionary of Jaccard dictinary to make best pairs analysis easier
def jaccardDict (gene1, gene2, geneDict):
    jaccardScores = {};
    for go_A in gene1:
        for go_B in gene2:
            jaccardScores[go_A] = {go_B : jaccardSimilarity(go_A, go_B, geneDict)};
    return jaccardScores;

#Get the information content of each inferred set from the GO Terms for use later in Resnik Similarity
def informationContent (inferredSets):
    ICDict = {};
    #set of all inferred sets combined (union).
    goSet = set();
    
    for gene in inferredSets:
        for goTerm in inferredSets[gene]:
            goSet.add(goTerm);
            
    goList = list(goSet);
    
    numerator = 0;
    denominator = len(inferredSets);
    for goTerm in goList:
        for gene in inferredSets:
            for elem in inferredSets[gene]:
                if goTerm == elem:
                    numerator += 1;
                    break;
        currentScore = numpy.math.log10((numerator / denominator));
        ICDict[goTerm] = currentScore * -1;
        numerator = 0;
    return ICDict;

#Find the lowest common parent node between two terms and compute the highest Information Content
def lowestCommonSubsumer (term1, term2, ICScores, geneDict):
    highestScore = 0;
    commonElements = [];
    lowestCommonDict = {};
    
    for elem1 in geneDict[term1]:
        for elem2 in geneDict[term2]:
            if elem1 == elem2:
                commonElements.append(elem1);
    for goTerm in commonElements:
        if ICScores[goTerm] > highestScore:
            highestScore = ICScores[goTerm];
    return highestScore;

#Compute Resnik Similarit measure
def resnikSimilarity (term1, term2, listOfGenes, geneNames, geneDict):
    ICScores = {};
    inferredSets = {};
    currentSet = set();
    counter = 0;
	
	#Get the inferred sets (implied set of superclasses) for each GO Term in the list of genes
    for gene in listOfGenes:
        for goTerm in gene:
            for superClass in geneDict[goTerm]:
                currentSet.add(superClass);
        inferredSets[geneNames[counter]] = list(currentSet);
        currentSet.clear();
        counter += 1;
    
	#Get information content of each GO value.
    ICScores = informationContent(inferredSets);
    lowestCommon = lowestCommonSubsumer(term1, term2, ICScores, geneDict);
    return lowestCommon;
   
#Create a dictionary of Resnik scores that can be used to more easily compute Best Pairs scores
def resnikDict (gene1, gene2, listOfGenes, geneNames, geneDict):
    resnikScores = {};
    for go_A in gene1:
        for go_B in gene2:
            resnikScores[go_A] = {go_B : resnikSimilarity(go_A, go_B, listOfGenes, geneNames, geneDict)};
    return resnikScores;
            
#Simply calculate average from a list of numbers
def average (numList):
    aggregateSum = 0;
    for num in numList:
        aggregateSum += num;
    return aggregateSum / len(numList);

#All pairs takes all semantic similarity scores (either Jaccard or Resnik) and gets the average
def allPairs (semanticScores):
    return average(semanticScores);

#Best pairs finds the best scores amongst each pairing of GO Terms and then finds the average
def bestPairs (gene1, gene2, scoreDict):
    maxList = [];
    currentMax = 0;
	
	#Iterate through the scores of semantic similarity measures as long as each gene list has GO Terms.
	#Then append the currentMax to a list of max scores for later computation.
    if len(gene1) > 1 and len(gene2) > 1:
        for go_A in scoreDict:
            for go_B in scoreDict[go_A]:
                if scoreDict[go_A][go_B] > currentMax:
                    currentMax = scoreDict[go_A][go_B];
            maxList.append(currentMax);
            currentMax = 0;
    else:
        for go_C in scoreDict:
            for go_D in scoreDict[go_C]: 
                maxList.append(scoreDict[go_C][go_D]);
	#Find average of all max scores found from each list.
    return average(maxList);


print("Jaccard All Pairs (GeneA and GeneC): " + str( allPairs(jaccardList(geneA, geneC, owlDict))) );
print("Jaccard Best Pairs (GeneA and GeneC): " + str( bestPairs(geneA, geneC, jaccardDict(geneA, geneC, owlDict))) );
geneList = [geneA, geneB, geneC];
geneNames = ["geneA", "geneB", "geneC"];
print("Resnik Best Pairs (GeneA and GeneB): " + str( bestPairs(geneA, geneB, resnikDict(geneA, geneB, geneList, geneNames, owlDict))) );





