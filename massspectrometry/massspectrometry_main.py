from massspectrometry.objs import Peptide as p
from massspectrometry.objs import SequenzDB as sdb
from massspectrometry.objs import AminoAcid as AA


def getCandidatesFor(peptideList, db, margin):
    candidatesList = list()
    for peptide in peptideList["Peptides"]:
        jobDict = dict()
        jobDict["Peptide"] = peptide
        candidates = list()
        unchargedMass = peptide["PepmassUncharged"]
        for entry in db:
            if float(entry["PepmassUncharged"]) > unchargedMass - margin and float(
                    entry["PepmassUncharged"]) < unchargedMass + margin:
                candidates.append(entry)
        jobDict["Candidates"] = candidates.copy()
        candidatesList.append(jobDict)
    return candidatesList


def getScoresForJobs(jobs, aaList):
    for j in jobs:
        score = 0
        for candidate in j["Candidates"]:
            candidate["IonMass"] = getBetaGammaIons(candidate, aaList)
        j["BestSequence"] = compareIonMassWithSpectrumMass(j, 0.04)


def compareIonMassWithSpectrumMass(candidateJob, threshold):
    peptide = candidateJob["Peptide"]
    candidates = candidateJob["Candidates"]

    bestScoreSeq = ""
    bestScore = float()

    for candidate in candidates:
        matches = 0
        for beta in candidate["IonMass"]["BetaIons"]:
            matches += countPositveIonHitsWithSpectrumAndThreshold(beta[1], peptide["m/z"], threshold)
        for gamma in candidate["IonMass"]["GammaIons"]:
            matches += countPositveIonHitsWithSpectrumAndThreshold(gamma[1], peptide["m/z"], threshold)

        #calculate score
        totalBetaGammaIons = len(candidate["IonMass"]["BetaIons"])+len(candidate["IonMass"]["GammaIons"])
        score = matches/totalBetaGammaIons

        if score > bestScore:
            bestScore = score
            bestScoreSeq = candidate["Sequence"]


    result = dict()
    result["Best Sequence"] = bestScoreSeq
    result["Score"] = bestScore
    return result


def countPositveIonHitsWithSpectrumAndThreshold(ionMass, mzList, threshold):
    positiveHits = 0
    for mz in mzList:
        if mz[0] + threshold > ionMass > mz[0] - threshold:
            positiveHits += 1
    return positiveHits


def getBetaGammaIons(candidate, aaList):
    result = dict()
    betaIons = list()
    gammaIons = list()
    betaIons = getBetaIons(candidate["Sequence"], aaList)
    gammaIons = getGammaIons(candidate["Sequence"], aaList)

    result["BetaIons"] = betaIons.copy()
    result["GammaIons"] = gammaIons.copy()
    return result


def getBetaIons(sequence, aaList):
    # from first until last-1 position
    result = list()
    for i in range(1, len(sequence)):
        partseq = sequence[0:i]
        mass = getBetaIonMass(partseq, aaList)
        result.append((partseq, mass))
    return result


def getGammaIons(sequence, aaList):
    result = list()
    for i in range(len(sequence) - 1, 0, -1):
        partseq = sequence[i:]
        mass = getGammaIonMass(partseq, aaList)
        result.append((partseq, mass))

    return result


def getBetaIonMass(sequence, aaList):
    mass = 0
    for s in sequence:
        mass += getAAMassValue(s, "Single Letter Code", "Mass", aaList)
    return mass + 1.007  # adding proton mass


def getGammaIonMass(sequence, aaList):
    mass = 0
    for s in sequence:
        mass += getAAMassValue(s, "Single Letter Code", "Mass", aaList)
    return mass + 1.007 + 18  # adding proton and h2o mass


def getAAMassValue(match, matchFeatureName, searchForValue, inAAList):
    for aa in inAAList:
        if match == aa[matchFeatureName]:
            res = aa[searchForValue]
            return float(res)
        # for key, value in aa.items():
        #    if key == searchForValue:
        #        return value


def printResults(jobs):
    for job in jobs:
        print(job["Peptide"]["Title"], " => Best database candidate: '",
              job["BestSequence"]["Best Sequence"],"' Score: ",
              job["BestSequence"]["Score"])

peptides = p.Peptide("spectraMS2.mgf").getPeptides()
db = sdb.SequenzDB("database.csv").getDataBase()
jobs = getCandidatesFor(peptides, db, 0.02)
aa = AA.AminoAcid("aminoacids.csv").getAminoAcids()
getScoresForJobs(jobs, aa)
printResults(jobs)

