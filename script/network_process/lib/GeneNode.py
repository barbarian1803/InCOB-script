class GeneNode:
    geneID = ""
    outWardEdge = set()
    inWardEdge = set()
    nodeType = ""
    logFC = ""
    score = ""
    connectedNodes = set()

    def __init__(self,geneID):
        self.geneID = geneID
        self.outWardEdge = set()
        self.inWardEdge = set()
        self.connectedNodes = set()
        self.logFC = 0
        self.score = 0
        
        if "ENSG" in geneID:
            self.nodeType = "gene"
        else:
            self.nodeType = "protein"

    def printNode(self):
        print "GeneID : "+self.geneID
        print "Type : "+self.nodeType
        print "Score : "+str(self.score)
        print "OutwardEdge:"
        print self.printEdge(self.outWardEdge)
        print "InwardEdge:"
        print self.printEdge(self.inWardEdge)

    def addInWardNode(self,node):
        self.inWardEdge.add(node)

    def addOutWardNode(self,node):
        self.outWardEdge.add(node)	

    def printEdge(self,edge):
        for n in edge:
            print n

    def printEdgeRel(self,toPrint):
        for n in self.outWardEdge:
            if self.nodeType=="gene" and "ENSG" not in n:  # if source is gene and target is protein
                toPrint.add(self.geneID+"\t"+n+"\t"+"to_protein"+"\n")
            else:
                toPrint.add(self.geneID+"\t"+n+"\t"+"bind_to"+"\n")
                
        for n in self.inWardEdge:
            if self.nodeType=="protein" and "ENSG" in n: #if source is gene and target is protein
                toPrint.add(n+"\t"+self.geneID+"\t"+"to_protein"+"\n")
            else:
                toPrint.add(n+"\t"+self.geneID+"\t"+"bind_to"+"\n")
                    
    def isEndNode(self):
        outward = [x for x in self.outWardEdge if x is not None]
        return (len(outward)==0)
        
    def isStartNode(self):
        inward = [x for x in self.inWardEdge if x is not None]
        return (len(inward)==0)
    
    def getInWardNode(self):
        return [x for x in self.inWardEdge if x is not None]
    
    def getOutWardNode(self):
        return [x for x in self.outWardEdge if x is not None]
    
    def removeInWardNode(self,nodes):
        self.inWardEdge.remove(nodes)

    def removeOutWardNode(self,nodes):
        self.outWardEdge.remove(nodes)
