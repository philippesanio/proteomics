"""
@topic: Parsing a simple mass spectral file (.mgf)
@author: Philippe Sanio
@version: 1.0.0
"""

class Peptide:
    __peptides = dict()

    def __init__(self, fileName):
        self.__peptides["Peptides"] = list()
        self.__getPeptideFromFile(fileName)

    def __getPeptideFromFile(self, filename):
        peptide = dict()
        with open(filename) as reader:
            for untreatedLine in reader:
                line = untreatedLine.replace("\n", "")
                line.strip()
                if line.startswith("BEGIN IONS"):
                    peptide = dict()
                    peptide["m/z"] = list()
                elif line.startswith("TITLE"):
                    peptide["Title"] = line.split("TITLE=")[1]
                elif line.startswith("PEPMASS"):
                    pepmass = float(line.split("PEPMASS=")[1].split(" ")[0])
                    peptide["Pepmass"] = float(pepmass)
                elif line.startswith("CHARGE"):
                    charge = line.split("CHARGE=")[1][0]
                    sign = line.split("CHARGE=")[1][1]
                    if sign == "-":
                        peptide["Charge"] = int(charge) * (-1)
                    else:
                        peptide["Charge"] = int(charge)
                elif line.startswith("RTINSECONDS"):
                    peptide["Rtinseconds"] = int(line.split("RTINSECONDS=")[1])
                elif line.startswith("SCANS"):
                    peptide["Scans"] = int(line.split("SCANS=")[1])
                elif line.startswith("END IONS"):
                    peptide["Pepmass"] = peptide["Pepmass"] * peptide["Charge"]
                    peptide["PepmassUncharged"] = peptide["Pepmass"] - 1.007 * abs(peptide["Charge"])
                    peptide["PepmassSingleCharged"] = peptide["PepmassUncharged"] + 1.007
                    self.__peptides["Peptides"].append(peptide.copy())
                elif line.startswith("MASS"):
                    self.__peptides["Name"] = line
                else:
                    if (line != ""):
                        mz, relAbundance = line.split(" ")
                        peptide["m/z"].append(list([float(mz),float(relAbundance)]))

    def getPeptides(self):
        return self.__peptides
