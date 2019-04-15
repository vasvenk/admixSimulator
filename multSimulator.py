import random
import sys
from math import exp


def readData(genoFile, snpFile, indFile):
    """
    Function for creating the relevant data types
    """
    snpData, snpDataCols = open(snpFile).read().split(), []
    for index in range(len(snpData) // 6):
        snpDataCols.append(snpData[6 * index: 6 * index + 6])
    genoData = open(genoFile)
    genoData = genoData.read().split()
    indData, indDataCols = open(indFile).read().split(), []
    for index in range(len(indData) // 3):
        indDataCols.append(indData[3 * index: 3 * index + 3])
    return {"ind": indDataCols, "geno": genoData, "snp": snpDataCols}


def initializeState(numInd, alpha, geno, numA, numB):
    """
    Initializes populations, assigning each admixed individual to an individual
    from a sample population based on the admixture proportion.
    """
    usedInd, freeInd, newGeno = {}, {'A': [x for x in range(numA)], 'B': [x for x in range(numB)]}, ""
    for ind in range(numInd):
        chosenPop = 'A' if random.uniform(0, 1) < alpha else 'B'
        random.shuffle(freeInd[chosenPop])
        sampleInd = freeInd[chosenPop].pop()
        usedInd[ind] = {'pop': chosenPop, 'indNum': sampleInd}
        newGeno += geno[chosenPop][0][sampleInd]
    return newGeno, usedInd, freeInd


def swapSample(usedInd, freeInd, ind, geno, currPos, alpha):
    chosenPop = 'A' if random.uniform(0, 1) < alpha else 'B'
    prevAncestor = usedInd[ind]
    random.shuffle(freeInd[chosenPop])
    newInd = freeInd[chosenPop].pop()
    usedInd[ind] = {'pop': chosenPop, 'indNum': ind}
    freeInd[prevAncestor['pop']].append(prevAncestor['indNum'])
    return geno[chosenPop][currPos][newInd], usedInd, freeInd


def haploidToDiploid(haploidGeno):
    diploidGeno = ''
    for index in range(0, len(haploidGeno), 2):
        currInd, nextInd = int(haploidGeno[index]), int(haploidGeno[index + 1])
        diploidInd = str(currInd + nextInd)
        diploidGeno += diploidInd
    return diploidGeno


def writeSnp(snp, popName):
    with open(popName + '.snp', 'w') as f:
        for entry in snp:
            f.write(str(entry[0]) + " " + str(entry[1]) + " " + str(entry[2]) + " " + str(entry[3]) + " " +
                    str(entry[4]) + " " + str(entry[5]) + "\n")


def writeGeno(geno, popName):
    with open(popName + '.geno', 'w') as f:
        for entry in geno:
            f.write(entry + "\n")


def writeInd(ind, popName):
    with open(popName + '.ind', 'w') as f:
        f.write(ind)

def writeAncestry(ancestry, popName):
    with open(popName + '.anc', 'w') as f:
        for entry in ancestry:
            f.write(entry + "\n")



def singleAdmixture(popA, popB, numGenerations, alpha, numInd, haploid, track, popName):
    """
    Model for a single admixture event.
    Takes in dictionaries as created by readData for each population A and B.
    Alpha is the portion of population A.
    NumGenerations is the number of generations since admixture.
    """
    admixedGeno, admixedSnp, ancestry = [], [], []
    ind = {'A': popA["ind"], 'B': popB["ind"]}
    geno = {'A': popA["geno"], 'B': popB["geno"]}
    snp = {'A': popA["snp"], 'B': popB["snp"]}
    numA, numB = len(ind['A']), len(ind['B'])

    if not haploid:
        numInd *= 2
        newGeno, currUsed, free = initializeState(numInd, alpha, geno, numA, numB)
        newGeno = haploidToDiploid(newGeno)
    else:
        newGeno, currUsed, free = initializeState(numInd, alpha, geno, numA, numB)
    admixedGeno.append(newGeno)
    admixedSnp.append(snp['A'][0])
    totalLength = len(snp['A'])
    currPercent, lastPercent = 0, -1

    for currPos in range(1, totalLength):
        distance = float(snp['A'][currPos][2]) - float(snp['A'][currPos - 1][2])
        currGeno = ""
        if snp['A'][currPos][1] == snp['A'][currPos - 1][1]:
            #Same Chromosome
            currAncestry = ""
            for ind in range(numInd):
                if random.uniform(0, 1) < 1 - exp(- numGenerations * distance):
                    #Recombination Event-we switch ancestors
                    newPos, currUsed, free = swapSample(currUsed, free, ind, geno, currPos, alpha)
                    currGeno += newPos
                else:
                    #No recombination- continue as usual
                    chosenPop, chosenInd = currUsed[ind]['pop'], currUsed[ind]['indNum']
                    currGeno += geno[chosenPop][currPos][chosenInd]

                if track:
                    currAncestry += currUsed[ind]['pop']
        else:
            #Next Chromosome
            currGeno, currUsed, free = initializeState(numInd, alpha, geno, numA, numB)

            if track:
                currAncestry = ""
                for ind in range(numInd):
                    currAncestry += currUsed[ind]['pop']

        if not haploid:
            currGeno = haploidToDiploid(currGeno)

        admixedSnp.append(snp['A'][currPos])
        admixedGeno.append(currGeno)

        if track:
            ancestry.append(currAncestry)

        currPercent = int((currPos / totalLength) * 100)
        if currPercent > lastPercent:
            print("Completed: " + str(currPercent))
            lastPercent = currPercent

    admixedInd, numIters = '', numInd if haploid else numInd // 2
    for indNum in range(numIters):
        admixedInd += popName + str(indNum) + " " + "U" + " " + popName + "\n"

    return admixedGeno, admixedSnp, admixedInd, ancestry

allGenerations = [5*x for x in range(20)]
admixtureProportion = 0.6
numInd = 10
baseFile = "/global/scratch/vasvenk/data/"
outFile = "/global/scratch/vasvenk/rfMixFiles/admixtures/CEU6/"
popA = readData(baseFile + "CEU.phgeno", baseFile + "CEU.phsnp", baseFile + "CEU.phind")
popB = readData(baseFile + "YRI.phgeno", baseFile + "YRI.phsnp", baseFile + "YRI.phind")
for numGenerations in allGenerations:
    popName = "CEUYRI" + str(numGenerations)
    admixedGeno, admixedSnp, admixedInd, ancestry = singleAdmixture(popA, popB, numGenerations, admixtureProportion, numInd, True, True, popName)
    writeGeno(admixedGeno, outFile + popName)
    writeSnp(admixedSnp, outFile + popName)
    writeInd(admixedInd, outFile + popName)
    writeAncestry(ancestry, outFile + popName)
