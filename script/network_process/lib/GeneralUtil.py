class GeneralUtil:

    @staticmethod
    def readGeneIDDatabase(fileName):
        data = {}
        csv1 = open(fileName)
        for line in csv1:
            if "ensembl" in line:
                continue

            array=line.strip().split("\t")
            if array[0] not in data:
                data[array[0]] = {}
            data[array[0]] = {"ensembl":array[0],"symbol":array[1],"entrez":array[2]}
            try:
                data[array[0]]["uniprot"]=array[3]
            except:
                data[array[0]]["uniprot"]=""
        return data