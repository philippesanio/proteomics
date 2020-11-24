class SequenzDB:

    __db = list()

    def __init__(self, databaseFile):
        self.__db = self.__readDataBase(databaseFile)

    def __readDataBase(self, databaseFile):
        db = list()
        with open(databaseFile) as reader:
            for line in reader:
                d = dict()
                s = line.replace("\n","").split(',')
                #mass, seq = s[0],s[1]
                d["PepmassUncharged"] = s[0]
                d["Sequence"] = s[1]
                db.append(d)
        return db


    def getDataBase(self):
        return self.__db

