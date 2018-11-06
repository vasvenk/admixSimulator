import numpy as np
import random
from math import exp

def readData(indFile, genoFile, snpFile):
    """
    Function for creating the relevant data types
    """
    snpData, snpDataCols = open(snpFile).read().split(), []
    for index in range(len(snpData) // 6):
        snpDataCols.append(snpData[6 * index: 6 * index + 6])
    snpDataCols.sort(key=lambda x: x[3])
    genoData = open(genoFile, encoding="ISO-8859-1")
    genoData = genoData.read().split()
    indData, indDataCols = open(indFile).read().split(), []
    for index in range(len(indData) // 3):
        indDataCols.append(indData[3 * index: 3 * index + 3])
    return {"ind": indDataCols, "geno": genoData, "snp": snpDataCols}

def mixedPopulation(genoA, genoB, alpha):
    #Returns a mixed segment
    allB, allA, newGeno = [x for x in range(len(genoB))], [x for x in range(len(genoA))], ""
    for ind in range(len(genoA)):
        if random.uniform(0, 1) < alpha:
            currPop, currGeno = allA, genoA
        else:
            currPop, currGeno = allB, genoB
        chosenInd = random.choice(currPop)
        currPop.remove(chosenInd)
        numRef = currGeno[chosenInd]
        newGeno += numRef
    return newGeno

def singlePopulation(geno):
    #Return a segment for a single population
    numInd = len(geno)
    newGeno, allInds = "", [x for x in range(numInd)]
    for ind in range(numInd):
        chosenInd = random.choice(allInds)
        allInds.remove(chosenInd)
        numRef = geno[chosenInd]
        newGeno += numRef
    return newGeno

def singleAdmixture(popA, popB, alpha, numGenerations):
    """
    Model for a single admixture event.
    Takes in dictionaries as created by readData for each population A and B.
    Alpha is the portion of population A.
    NumGenerations is the number of generations since admixture.
    """
    e = 2.7182
    admixedGeno, admixedSnp = [], []
    indA, genoA, snpA = popA["ind"], popA["geno"], popA["snp"]
    indB, genoB, snpB = popB["ind"], popB["geno"], popB["snp"]
    for currPos in range(1, len(snpA)):
        distance = float(snpA[currPos][2]) - float(snpA[currPos - 1][2])
        if np.random.uniform() < exp(- numGenerations * distance):
            #Choosing from pop A or B
            if np.random.uniform() < alpha:
                #Choose pop A
                newGeno = singlePopulation(genoA[currPos])
                admixedSnp.append(snpA[currPos])
                admixedGeno.append(newGeno)
            else:
                #Choose pop B
                newGeno = singlePopulation(genoB[currPos])
                admixedSnp.append(snpB[currPos])
                admixedGeno.append(newGeno)
        else:
            #Mixed population
            newGeno = mixedPopulation(genoA[currPos], genoB[currPos], alpha)
            admixedSnp.append(snpA[currPos])
            admixedGeno.append(newGeno)
    return admixedGeno, admixedSnp

"""Create CEU- spefic to my files"""
ceu = readData("../geneticData/CEU.phind", "../geneticData/CEU.phgeno", "../geneticData/CEU.phsnp")
ceuInd, ceuGeno, ceuSnp =  ceu["ind"], ceu["geno"], ceu["snp"]

yri = readData("../geneticData/YRI.phind", "../geneticData/YRI.phgeno",  "../geneticData/YRI.phsnp")
yriInd, yriGeno, yriSnp = yri["ind"], yri["geno"], yri["snp"]

admixedGeno, admixedSnp = singleAdmixture(ceu, yri, 0.5, 10)
