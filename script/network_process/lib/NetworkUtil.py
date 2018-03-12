from GeneNode import GeneNode
from sets import Set

class NetworkUtil:

    @staticmethod   ## helper method to add node to the network from edge list element
    def createNode(network,element):
        source = element[0]
        target = element[1]

        if source not in network:
            network[source] = GeneNode(source)

        if target not in network:
            network[target] = GeneNode(target)

        network[source].addOutWardNode(target)
        network[target].addInWardNode(source)

        return network

    @staticmethod   ## build network data structure using edge list csv file
    def buildNetwork(fileName,split):
        network = {}
        edge_list = open(fileName,"r")

        for l in edge_list:
            #skip header line
            if "Source" in l:
                continue

            array = l.strip().split(split)
            network = NetworkUtil.createNode(network,array)

        return network

    @staticmethod
    def getEndNodesOfNetwork(network):
        endNodeSet = Set()
        for node in network:
            if network[node].isEndNode():
                endNodeSet.add(node)
        return endNodeSet

    @staticmethod
    def getStartNodesOfNetwork(network):
        startNodeSet = Set()
        for node in network:
            if network[node].isStartNode():
                startNodeSet.add(node)
        return startNodeSet

    @staticmethod
    def load_fimo_network_result(fname):
        checked_string = []
        table = []
        fobj = open(fname)
        for line in fobj:
            line = line.rstrip()
            if line[0] == "#":
                continue
            if line in checked_string:
                continue
            checked_string.append(line)
            split = line.split("\t")
            table.append({"source":split[1],"target":split[0],"type":"regulate"})
        return table

    @staticmethod
    def load_transfac_network(fname):
        table = []
        fobj = open(fname)
        for line in fobj:
            line = line.rstrip()
            if line[0] == "#":
                continue
            split = line.split("\t")
            to = split[0].replace("$","_")
            tfs = split[3].split("|")
            for tf in tfs:
                table.append({"source":tf,"target":to,"type":"protein"})

        return table

    @staticmethod   ## input file is fimo network and transfac gene-protein network and create merged network file
    def network_builder(fimo_network, transfac_network,output):
        fimo = NetworkUtil.load_fimo_network_result(fimo_network)
        transfac = NetworkUtil.load_transfac_network(transfac_network)

        network = fimo+transfac

        fobj = open(output,"w")
        fobj.write("Source\tTarget\tType\n")

        for n in network:
            fobj.write(n["source"].replace("V_","")+"\t"+n["target"].replace("V_","")+"\t"+n["type"]+"\n")


    @staticmethod   ##recursive method that is called for longest chain calculation
    def calculateLongestChain(network,node,level,visited):

        if network[node].isEndNode() or node in visited:
            return level+1


        visited.add(node)
        maxVal = -999
        for n in network[node].getOutWardNode():
            result = NetworkUtil.calculateLongestChain(network,n,level+1,visited)
            if result>maxVal:
                maxVal = result

        return maxVal

    @staticmethod   ## calculate longest chain from every starting nodes in the network
    def longestNodesChain(network):
        startNodes = NetworkUtil.getStartNodesOfNetwork(network)
        output = {}
        for n in startNodes:
            visited = Set()
            level = NetworkUtil.calculateLongestChain(network,n,0,visited)
            output[n] = level

        return output

    @staticmethod   ##print network into edge list file, you can skip node if you want to exclude some node
    def printNetworkTofile(network,skipNode,fName):
        fName.write("Source\tTarget\tType\n")
        toPrint = Set()
        for n in network:
            network[n].printEdgeRel(toPrint)

        for n in toPrint:
            split = n.split("\t")
            if split[0] in skipNode or split[1] in skipNode:
                continue
            fName.write(n)

    @staticmethod   ## load node value from inputted file
    def loadNodeValue(fname):
        nodeValue = {}
        for line in fname:
            if "Gene" in line:
                continue

            array = line.split("\t")

            nodeValue[array[0]] = array[1]

        return nodeValue

    @staticmethod   ## set node value and store it in the score attribute
    def setNodeValue(network,valueList):
        for n in network:
            if n in valueList:
                network[n].score = float(valueList[n])
            else:
                network[n].score = float(0)
        return network

    @staticmethod
    def removeNode(network,node):   ##remove node from the network
        for n in network[node].getInWardNode():
            network[n].removeOutWardNode(node)
        for n in network[node].getOutWardNode():
            network[n].removeInWardNode(node)
        network.pop(node)
        return network

    @staticmethod   ##remove edge that connect exactly same source and target
    def removeSamePath(network):
        for n in network:
            if network[n].nodeType=="protein":
                continue

            connectedNodes = Set()
            outWardEdge = network[n].getOutWardNode()

            for edge in outWardEdge:

                if len(network[edge].getInWardNode())>1:
                    continue

                edge_outward = Set(network[edge].getOutWardNode())
                intersect = connectedNodes.intersection(edge_outward)
                for i in intersect:
                    network[edge].removeOutWardNode(i)
                    network[i].removeInWardNode(edge)
                connectedNodes = connectedNodes.union(edge_outward)
        return network

    @staticmethod ##remove protein nodes that are in the starting node or end node
    def filterProteinNode(network):
        start = NetworkUtil.getStartNodesOfNetwork(network)
        for n in start:
            if network[n].nodeType=="protein":
                NetworkUtil.removeNode(network,n)
        end = NetworkUtil.getEndNodesOfNetwork(network)
        for n in end:
            if network[n].nodeType=="protein":
                NetworkUtil.removeNode(network,n)

        return network

    @staticmethod
    def filterNetwork(network):
        output = NetworkUtil.longestNodesChain(network)
        while 1 in output.values() or 2 in output.values() or 3 in output.values():
            for n in output:
                if output[n]<=3:
                    network = NetworkUtil.removeNode(network,n)
            network =  NetworkUtil.filterProteinNode(network)
            output = NetworkUtil.longestNodesChain(network)

        network = NetworkUtil.removeSamePath(network)
        return network

    @staticmethod
    def calculateDepth(network,node,score,visited):
        if network[node].score<score:
            network[node].score = score
        if node in visited:
            return

        visited.add(node)

        for n in network[node].getOutWardNode():
            if n in visited:
                continue
            NetworkUtil.calculateDepth(network,n,score+1,visited)
        return

    @staticmethod
    def printChain(network,node,toPrint):
        for n in network[node].getInWardNode():
            if network[n].score==network[node].score-1:
                toPrint.add(n+"\t"+node+"\n")
                NetworkUtil.printChain(network,n,toPrint)

    @staticmethod
    def printLongestChain(network,startnode,fout):
        for n in network:
            network[n].score=0

        visited = Set()
        NetworkUtil.calculateDepth(network,startnode,1,visited)

        maxScore = -1
        node = Set()
        for n in network:
            if maxScore<network[n].score:
                maxScore = network[n].score
                node = Set()
                node.add(n)
            elif maxScore==network[n].score:
                node.add(n)
        toPrint = Set()
        while node:
            NetworkUtil.printChain(network,node.pop(),toPrint)

        for n in toPrint:
            fout.write(n)

    @staticmethod
    def assignLogFC(network,logFC_file_name):
        logFCData = {}
        logFC_file = file(logFC_file_name)
        for line in logFC_file:
        #skip header line
            if "Gene" in line:
                continue

            values = line.split("\t")

            if values[0] in network:
                network[values[0]].logFC = float(values[1])

        return network

    @staticmethod
    def removeNonDiffNode(network,threshold):
        toRemove = Set()
        for n in network:
            if network[n].nodeType=="protein":
                continue
            node = network[n]
            if abs(node.logFC) <= threshold:
                toRemove.add(n)

        for n in toRemove:
            network = NetworkUtil.removeNode(network,n)

        return network

    @staticmethod
    def getEdgeType(network,nodeSrc):
        if network[nodeSrc].nodeType=="protein":
            return "regulate"
        else:
            return "protein"

    # @staticmethod
    # def printNodeSequence(netw):
