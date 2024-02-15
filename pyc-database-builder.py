import pandas as pd
import math
from copy import deepcopy

class Cation:
    #class for holding matched mass and charges for metal cations
    def __init__(self, name, isotope, metal, charge, mass):
        self.name = name
        self.isotope = str(isotope)
        self.metal = metal
        self.charge = charge
        self.mass = mass

#O: 15.994915
#easy to update library of metals and polyatomics for consideration
#[isotope label, neutral monoisotopic mass, list(oxidation states to consider)]
class metalLib():
    Mg = [24, "Mg", 23.985045, [2]]
    Ca = [40, "Ca", 39.962591, [2]]
    V = [51, "V", 50.943963, [2, 3]]
    VO = [51, "V", 66.938878, [2, 3]]
    VO2 = [51, "V", 82.933793, [1]]
    Cr = [52, "Cr", 51.940510, [2, 3]]
    Mn = [55, "Mn", 54.938046, [2, 3]]
    Fe = [56, "Fe", 55.934939, [2, 3]]
    Ni = [58, "Ni", 57.935347, [1, 2]]
    Co = [59, "Co", 58.933198, [2, 3]]
    Cu = [63, "Cu", 62.929599, [1, 2]]
    Zn = [64, "Zn", 63.929145, [2]]
    Se = [80, "Se", 79.916521, [2]]
    SeO = [80, "Se", 95.911436, [2]]
    Mo = [98, "Mo", 97.905405, [4, 6]]
    Ag = [107, "Ag", 106.905095, [1, 2]]
    Cd = [114, "Cd", 113.903361, [2]]
    Te = [130, "Te", 129.906229, [2]]
    Hg = [202, "Hg", 201.970632, [1, 2]]
    Tl = [205, "Tl", 204.974410, [1, 3]]
    Pb = [208, "Pg", 207.976641, [2]]

#reference mass changes used in calculations. mass subtractive modifications written as negatives
#hydrations/dehydrations included in masses listed
class massChanges():
    pyc2 = 482.114120
    yEC = 232.051779
    desE = -129.042593
    G = 57.021464
    S = 87.032028
    Q = 128.058578
    A = 71.037114
    E = 129.042593
    disulfide = -2.015650
    thiolate = -1.007276

class adductCalcs():
    m3H = ["M+3H", 1.007276 * 3, 3]
    m2HNa = ["M+2H+Na", 22.989218 + 1.007276 * 2, 3]
    mH2Na = ["M+H+2Na", 22.989218 * 2 + 1.007276, 3]
    m3Na = ["M+3Na", 22.989218 * 3, 3]
    m2H = ["M+2H", 1.007276 * 2, 2]
    mHNH4 = ["M+H+NH4", 1.007276 + 18.033823, 2]
    mHNa = ["M+H+Na", 1.007276 + 22.989218, 2]
    mHK = ["M+H+K", 1.007276 + 38.963158, 2]
    mACN2H = ["M+ACN+2H", 1.007276 * 2 + 41.026547, 2]
    m2Na = ["M+2Na", 22.989218 * 2, 2]
    m2ACN2H = ["M+2ACN+2H", 41.026547 * 2 + 1.007276 * 2, 2]
    m3ACN2H = ["M+3ACN+2H", 41.026547 * 3 + 1.007276 * 2, 2]
    mH = ["M+H", 1.007276, 1]
    mNH4 = ["M+NH4", 18.033823, 1]
    mNa = ["M+Na", 22.989218, 1]
    mK = ["M+K", 38.963158, 1]
    mACNH = ["M+ACN+H", 41.026547 + 1.007276, 1]
    m2NaMinH = ["M+2Na-H", 22.989218 * 2 - 1.007276, 1]
    mACNNa = ["M+ACN+Na", 41.026547 + 22.989218, 1]
    m2KMinH = ["M+2K-H", 38.963158 * 2 - 1.007276, 1]
    m2ACNH = ["M+2ACN+H", 41.026547 * 2 + 1.007276, 1]

#builds list of cation objects based off of members of metalLib
def makeCations(metalLib):
    metals = []
    cations = []
    #builds list of metals in library and excludes inherited methods
    for i in vars(metalLib):
        if "_" not in i:
            metals += [i]
    for metal in metals:
        element = metalLib.__dict__[metal]
        for oxidation in element[3]:
            #name follows "[isotope label][Me.charge]" syntax. neutral mass - (oxidation state * mass of e) used to calculate cation mass
            cations += [Cation(str(element[0])+metal+".+"+str(oxidation), element[0], element[1], oxidation, element[2] - (oxidation * 0.0005))]
    return cations

#fixes pythons broken rounding
def properRound(n, decimals=6):
    expoN = n * 10 ** decimals
    if abs(expoN) - abs(math.floor(expoN)) < 0.5:
        result = math.floor(expoN) / 10 ** decimals
    else:
        result = math.ceil(expoN) / 10 ** decimals
    if decimals <= 0:
        return int(result)
    return result

def progressBar(current, total, bar_length=20):
    fraction = current / total
    arrow = int(fraction * bar_length - 1) * '-' + '>'
    padding = int(bar_length - len(arrow)) * ' '
    ending = '\n' if current == total else '\r'
    print(f'Progress: [{arrow}{padding}] {int(fraction*100)}%', end=ending)

#makes starter table based on range of subunits (smallest:largest) to use
def getBasePyCs(smallest, largest, massDict):
    basePyClist = pd.DataFrame(columns = ["base.pyc", "cys.num", "monoisotopic.mass"])
    for i in range(smallest, largest+1):
        newPyC = (i - 2) * massDict.yEC + massDict.pyc2
        basePyClist.loc[len(basePyClist)] = ["PyC" + str(i), i, newPyC]
    return basePyClist

#calculates de-gamma-glutamylation metabolites, assuming presence of GGT
def calcGGTmetabs(PyCTable, massDict):
    GGTmetabs = pd.DataFrame(columns = ["seq.name", "seq.mass", "base.pyc", "cys.num", "ggt.metab"])
    for i in PyCTable.index:
        desE = PyCTable.loc[i, "monoisotopic.mass"] + massDict.desE
        GGTmetabs.loc[len(GGTmetabs)] = [PyCTable.loc[i,"base.pyc"], PyCTable.loc[i,"monoisotopic.mass"],
                                         PyCTable.loc[i,"base.pyc"], PyCTable.loc[i,"cys.num"], 0]
        GGTmetabs.loc[len(GGTmetabs)] = ["des.yE." + PyCTable.loc[i,"base.pyc"], desE,
                                         PyCTable.loc[i,"base.pyc"], PyCTable.loc[i,"cys.num"], 1]
    return GGTmetabs

#calculates various forms by carboxy-terminal AA. Defaults to all possibilities
def calcCarboxyTerm(PyCTable, massDict, options= ["G", "S", "Q", "A", "E"]):
    carboxyTerms = pd.DataFrame(columns = ["seq.name", "seq.mass", "base.pyc", "cys.num", "ggt.metab", "c.terminus"])
    basePyCcounter = 0
    for i in PyCTable.index:
        if len(carboxyTerms) != 0 and PyCTable.loc[i, "base.pyc"] != PyCTable.loc[i-1, "base.pyc"]:
            basePyCcounter = len(carboxyTerms)
        carboxyTerms.loc[len(carboxyTerms)] = [PyCTable.loc[i, "seq.name"], PyCTable.loc[i, "seq.mass"],
                                               PyCTable.loc[i, "base.pyc"], PyCTable.loc[i, "cys.num"],
                                               PyCTable.loc[i, "ggt.metab"], "X"]
        for terminus in options:
            if terminus == "E" and "des.yE." in PyCTable.loc[i, "seq.name"]:
                carboxyTerms.loc[basePyCcounter, "seq.name"] = carboxyTerms.loc[basePyCcounter, "seq.name"] + "; " + PyCTable.loc[i, "seq.name"] + "." + terminus
            else:
                carboxyTerms.loc[len(carboxyTerms)] = [PyCTable.loc[i, "seq.name"] + "." + terminus, PyCTable.loc[i, "seq.mass"] + massDict.__dict__[terminus],
                                                        PyCTable.loc[i, "base.pyc"], PyCTable.loc[i, "cys.num"],
                                                        PyCTable.loc[i, "ggt.metab"], terminus]
    return carboxyTerms

#determines maximum number of disulfides and creates all possibilities from 0 to maximum
def calcDisulfides(PyCTable, massDict):
    disulfides = pd.DataFrame(columns = ["seq.name", "seq.mass", "base.pyc", "cys.num", "ggt.metab", "c.terminus", "disulfides"])
    for i in PyCTable.index:
        maxSS = PyCTable.loc[i, "cys.num"] // 2
        for s in range(0, maxSS+1):
            disulfides.loc[len(disulfides)] = [PyCTable.loc[i, "seq.name"], PyCTable.loc[i, "seq.mass"] + massDict.disulfide * s, PyCTable.loc[i, "base.pyc"],
                                               PyCTable.loc[i, "cys.num"], PyCTable.loc[i, "ggt.metab"], PyCTable.loc[i, "c.terminus"], s]
    return disulfides

#determines all possible thiolates among available reduced cysteines. creates all possibilities and records formal charge
def calcThiolates(PyCTable, massDict):
    thiolates = pd.DataFrame(columns = ["seq.name", "seq.mass", "base.pyc", "cys.num", "ggt.metab", "c.terminus", "disulfides", "thiolates",  "charge"])
    for i in PyCTable.index:
        maxThiolates = PyCTable.loc[i, "cys.num"] - 2 * PyCTable.loc[i, "disulfides"]
        for t in range(0, maxThiolates+1):
            thiolates.loc[len(thiolates)] = [PyCTable.loc[i, "seq.name"], PyCTable.loc[i, "seq.mass"] + massDict.thiolate * t, PyCTable.loc[i, "base.pyc"],
                                             PyCTable.loc[i, "cys.num"], PyCTable.loc[i, "ggt.metab"], PyCTable.loc[i, "c.terminus"],
                                             PyCTable.loc[i, "disulfides"], t, t * -1]
    return thiolates


####Building out Dictionary of Cation Combinations#####

class cationGroup:
    def __init__(self):
        self.charge = 0
        self.cations = []
        self.metals = []
        self.mass = 0
    
    def info(self):
        print(f'Charge: {self.charge}, Mass: {self.mass}, Cations: {self.cations}, Metals: {self.metals}')
    
    def update(self, cation):
        self.charge += cation.charge
        self.cations += [cation.name]
        self.metals += [cation.metal]
        self.mass += cation.mass

def getCationCombos(minInputZ, maxInputZ, cations, maxOutputZ):
    comboDict = {}
    maxPos = maxInputZ - minInputZ + maxOutputZ
    for z in range(0, maxPos + 1):
        comboDict[str(z)] = []
    comboList = [cationGroup()]
    comboList += getCationCombosChild(comboList, maxPos, cations)
    print("Organizing cation combination dictionary by charge.")
    for i in range(0, len(comboList)):
        if comboList[i].charge <= maxPos:
            key = str(comboList[i].charge)
            comboDict[key] += [comboList[i]]
        progressBar(i, len(comboList) - 1)
    return comboDict

def getCationCombosChild(comboList, maxPos, cations, printProgress=False):
    updatedList = []
    for i in range(0, len(comboList)):
        for cation in cations:
            newCombo = deepcopy(comboList[i])
            newCombo.update(cation)
            updatedList += [newCombo]
        if printProgress:
            progressBar(i, len(comboList) - 1)
    rerunList = []
    minCharge = maxPos
    for combo in updatedList:
        if combo.charge < maxPos:
            rerunList += [combo]
        if combo.charge < minCharge:
            minCharge = combo.charge
    if len(rerunList) > 0:
        print(f'Finished with charge +{minCharge}. Max charge sought: +{maxPos}.')
        updatedList += getCationCombosChild(rerunList, maxPos, cations, True)
    return updatedList

def reduceCombos(comboList):
    uniqueMass = {}
    for combo in comboList:
        if str(combo.mass) not in uniqueMass:
            uniqueMass[str(combo.mass)] = [combo]
        else:
            uniqueMass[str(combo.mass)] += [combo]
    finalList = []
    for group in list(uniqueMass.keys()):
        trueUnique = {}
        for combo in uniqueMass[group]:
            if str(sorted(combo.cations)) not in trueUnique:
                trueUnique[str(sorted(combo.cations))] = [combo]
        for key in list(trueUnique.keys()):
            finalList += trueUnique[key]
    return finalList

def reduceComboDict(comboDict):
    for charge in list(comboDict.keys()):
        comboDict[charge] = reduceCombos(comboDict[charge])
    return comboDict

#took roughly 7-8 hours with all options applied to PyC2-5 and min=-1 and max=-2
def alignPyCCombos(PyCTable, comboDict, minCharge, maxCharge):
    conjugateTable = pd.DataFrame(columns = ["seq.name", "seq.mass", "base.pyc", "cys.num", "ggt.metab", "c.terminus", "disulfides", "thiolates",  "charge", "cations", "metals"])
    pycCharges = PyCTable['charge'].unique()
    matchDict = {}
    for initialCharge in pycCharges:
        comboChargeMatches = []
        for charge in list(comboDict.keys()):
            if initialCharge + int(charge) >= minCharge and initialCharge + int(charge) <= maxCharge:
                comboChargeMatches += [charge]
        matchDict[str(initialCharge)] = comboChargeMatches
    for key in list(matchDict.keys()):
        print(f'\nSearching for cation combinations for {key} charged phytochelatins.')
        for i in range(0, len(matchDict[key])):
            cationKey = matchDict[key][i]
            print(f'Matching cation combinations with {cationKey} charge, charge group {i+1} of {len(matchDict[key])}.')
            for j in range(0, len(comboDict[cationKey])):
                combo = comboDict[cationKey][j]
                pycSubset = PyCTable[PyCTable["charge"] == int(key)].copy()
                pycSubset["seq.mass"] += combo.mass
                pycSubset["charge"] += combo.charge
                pycSubset["cations"] = [combo.cations] * len(pycSubset)
                pycSubset["metals"] = [combo.metals] * len(pycSubset)
                conjugateTable = pd.concat([conjugateTable, pycSubset], ignore_index=True)
                if len(comboDict[cationKey]) > 1:
                    progressBar(j, len(comboDict[cationKey])-1)
    return conjugateTable

def createPseudocharges(conjugateTable,  massChanges, minCharge, maxCharge):
    pseudoCharged = pd.DataFrame(columns = conjugateTable.columns)
    for charge in range(minCharge, maxCharge+1):
        conjugateSubset = conjugateTable[conjugateTable["charge"] == charge].copy()
        if charge == 0:
            pseudoCharged = pd.concat([pseudoCharged, conjugateSubset])
        else:
            conjugateSubset["monoisotopic.mass"] += (charge * massChanges.thiolate)
            pseudoCharged = pd.concat([pseudoCharged, conjugateSubset])
    pseudoCharged["pseudocharge.modification"] = -pseudoCharged["charge"]
    return pseudoCharged

PyCTable = getBasePyCs(2, 5, massChanges)
PyCTable = calcGGTmetabs(PyCTable, massChanges)
PyCTable = calcCarboxyTerm(PyCTable, massChanges)
PyCTable = calcDisulfides(PyCTable, massChanges)
PyCTable = calcThiolates(PyCTable, massChanges)
cations = makeCations(metalLib)

minCharge = -1
maxCharge = 2

comboDict = getCationCombos(min(PyCTable['charge']), max(PyCTable['charge']), cations, maxCharge)
comboDict = reduceComboDict(comboDict)

conjugateTable = alignPyCCombos(PyCTable, comboDict, minCharge, maxCharge)
conjugateTable = conjugateTable.rename(columns={"seq.mass": "monoisotopic.mass"})

