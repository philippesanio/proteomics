from massspectrometry.objs import Peptide as p
from massspectrometry.objs import  SequenzDB as sdb

def getCandidatesFor(peptideList, db, margin):
    candidatesList = list()
    for peptide in peptideList["Peptides"]:
        jobDict = dict()
        jobDict["Peptide"]= peptide
        candidates = list()
        unchargedMass = peptide["PepmassUncharged"]
        for entry in db:
            if float(entry["PepmassUncharged"]) > unchargedMass - margin and float(entry["PepmassUncharged"]) < unchargedMass + margin:
                candidates.append(entry)
        jobDict["Candidates"] = candidates.copy()
        candidatesList.append(jobDict)
    return candidatesList

peptides = p.Peptide("spectraMS2.mgf").getPeptides()
db = sdb.SequenzDB("database.csv").getDataBase()
jobs = getCandidatesFor(peptides,db,0.02)

#calc for every possible canidate all b and y ion masses
#compare ion with masses aus spectrum (margine 0,04)
#count positive hits
#calc scoore = Sum(posHits)/Sum(b+y Ions)
#report best candidat with scoore + print sequence

print()