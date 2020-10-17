"""
@topic: Detecting amount of AA in Helices, Betasheets and Loops
@author: Philippe Sanio
@version: 1.0.0
"""

from pymol import cmd

aaList = ["ALA", "ARG", "ASN", "ASP", "CYS", "GLU", "GLN", "GLY", "HIS", "ILE", "LEU", "LYS", "MET", "PHE", "PRO",
          "SER", "THR", "TRP", "TYR", "VAL"]
aaDict = dict()  # example: {"ALA":[cntHELIX, cntBETASHEET, cntLOOPS]}


def countStructuuresForAllAA():
    """
    Detects how often an amino acid has been detected in a helix, betasheet or a loop
    and displays how often it has been counted.
    """
    for aa in aaList:
        aaDict[aa] = countStructuresForAA(aa)

    for aa in aaDict:
        print(
            aa + " => " + str(aaDict[aa])
            + " \t\t=>\t HELIX: " + str(round(aaDict[aa][0] / sum(aaDict[aa]), 3)) + " %"
            + " \t BETA: " + str(round(aaDict[aa][1] / sum(aaDict[aa]), 3)) + " %"
            + " \t LOOP: " + str(round(aaDict[aa][2] / sum(aaDict[aa]), 3)) + " %"
        )


def countStructuresForAA(aa="ALA", elemName="all"):
    """
    Detects for a single amino acid how
    :param aa: 3 letter code amino acid
    :param elemName: all
    :return: list of counts in [helices, betasheet, loops, total amount]
    """
    structurList = [0, 0, 0]  # [H,S,L]
    objList = cmd.get_names("objects")
    obj = objList[0]  # elements (right in the protein browser in pymol)
    models = cmd.get_model(obj + " and resn " + aa + " and name CA")

    for i in models.atom:
        if i.ss == "H":
            structurList[0] += 1  # counting the Helix
        elif i.ss == "S":
            structurList[1] += 1
        else:
            structurList[2] += 1
    return structurList

#make function usable in PyMol
cmd.extend("countStructuresForAA", countStructuresForAA)
cmd.extend("countStructuresForAllAA", countStructuuresForAllAA)
