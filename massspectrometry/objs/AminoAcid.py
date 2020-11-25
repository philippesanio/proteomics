class AminoAcid:

    __AA = list()

    def __init__(self,aaCSV):
        self.__AA = self.__readCSVFile(aaCSV)

    def __readCSVFile(self, csvFile):
        aa = list()
        featureNames = list()
        with open(csvFile) as reader:
            cnt = 0
            for line in reader:
                lPre = line.replace("\n","")
                l = lPre.split(";")
                if cnt == 0: #get feature names (Single Letter Code, Three Letter Code, ...
                    featureNames = l
                else:
                    if len(featureNames)>0:
                        aaDict = dict()
                        for i in range(len(featureNames)):
                            aaDict[featureNames[i]] = l[i]
                        aa.append(aaDict.copy())
                cnt += 1
        return aa

    def getAminoAcids(self):
        return self.__AA