import codecs
import csv
import os
import statsmodels.stats.multitest as smm
import pandas as pd
import sys
import zipfile
from operator import itemgetter
from operator import add
from multiprocessing import Pool

def splitCSV(FILENAME,TMP_FOLDER,IS_ZIP):
    ENCODING = "utf-8"
    if IS_ZIP=="True":
        file_dir = "/".join(FILENAME.split("/")[:-1])
        print "Zip file, extracting"
        zip = zipfile.ZipFile(FILENAME)
        zip.extractall(file_dir)
        FILENAME = FILENAME.replace(".zip",".txt")

    if not os.path.isdir(TMP_FOLDER):
        os.makedirs(TMP_FOLDER)
        
    with codecs.open(FILENAME, "r", ENCODING) as fp:
        reader = csv.reader(fp,delimiter='\t')
        ## read rest of file
        line = 0
        for row in reader:
            if "#" in "-".join(row):
                continue
            line+=1
            if line%10000000==0:
                print "Current line: "+str(line)
            ## Print each row
            row[0] = row[0].replace("V_","")
            row[2] = row[2].split(":")[0]
            out = open(TMP_FOLDER+"/"+row[2]+".csv","a")
            out.write('\t'.join(row)+"\n")
            out.close()
    
    if IS_ZIP=="True":
        os.remove(FILENAME)

def preprocess_table(fileObj):
    targetGene = {}
    for line in fileObj:
        if line[0]=="#":
            continue
        
        lineEl = parseLine(line)
        if lineEl["target_name"] not in targetGene:
            targetGene[lineEl["target_name"]] = {}
            
        if lineEl["matrix_name"] not in targetGene[lineEl["target_name"]]:
            targetGene[lineEl["target_name"]][lineEl["matrix_name"]] = []
        
        #if lineEl not in targetGene[lineEl["target_name"]][lineEl["matrix_name"]]:
        targetGene[lineEl["target_name"]][lineEl["matrix_name"]].append(lineEl)
        
    retval = sortByPval(targetGene)
    
    return retval
      
def parseLine(line):
    retval = {}
    element = line.split("\t")
    
    retval["matrix_name"] = element[0]
    retval["target_name"] = element[2].split(":")[0]
    retval["start"] = int(element[3])
    retval["stop"] = int(element[4])
    if element[5]=="+":
        retval["strand"] = 1
    else:
        retval["strand"] = -1
    retval["p_val"] = float(element[7])
    
    return retval
    
def sortByPval(table):
    output = {}
    for target in table.keys():
        output[target] = {}
        for matrix in table[target].keys():
            output[target][matrix] = sorted(table[target][matrix], key=itemgetter('p_val'))
    return output
    
def calculateQVal(tableInput):
    for target in tableInput.keys():
        for matrix in tableInput[target].keys():
            matchPValuesList = tableInput[target][matrix]
            Q = len(matchPValuesList)
            pvals = []
            for i in xrange(0,len(matchPValuesList)):
                pval = matchPValuesList[i]["p_val"]
                pvals.append(pval)
                rank = i+1
                matchPValuesList[i]["rank"] = rank
                matchPValuesList[i]["Q"] = Q
            
            qvals = smm.multipletests(pvals,method='fdr_bh')[1]
            
            for i in xrange(0,len(matchPValuesList)):
                matchPValuesList[i]["q_val"] = qvals[i]
            
            tableInput[target][matrix] = matchPValuesList
    return tableInput
    
def qval_filter(tableInput,qv_thresh):
    output = {}
    for target in tableInput.keys():
        ## bonferoni correction, qvalue will be divided by the number of matched matrix for a given promoter region
        bonferroni = len(tableInput[target].keys())
        #print "no mat: "+str(bonferroni)
        output[target] = {}
        for matrix in tableInput[target].keys():
            output[target][matrix] = []
            matchList = tableInput[target][matrix]
            for el in matchList:
                if el["q_val"] <= (qv_thresh/bonferroni):
                    output[target][matrix].append(el)
    return output

def printTable(tableInput,fileName):
    fobj = open(fileName,"w")
    header = ["## Target","Matrix name","start","stop","strand","pvalue","rank","Q","qvalue"]
    header = "\t".join(header)+"\n"
    fobj.write(header)
    for target in tableInput.keys():
        for matrix in tableInput[target].keys():
            matchList = tableInput[target][matrix]
            for el in matchList:
                data = [el["target_name"],el["matrix_name"],str(el["start"]),str(el["stop"]),str(el["strand"]),str(el["p_val"]),str(el["rank"]),str(el["Q"]),str(el["q_val"])]
                data = "\t".join(data)+"\n"
                fobj.write(data)
                
def printNetwork(tableInput,fileName):
    fobj = open(fileName,"w")
    fobj.write("## Target\tMatrix name\n")
    pairs = []
    for target in tableInput.keys():
        for matrix in tableInput[target].keys():
            matchList = tableInput[target][matrix]
            for el in matchList:
                el["matrix_name"] = el["matrix_name"].replace("V_","")
                if el["target_name"]+"-"+el["matrix_name"] not in pairs:
                    fobj.write(el["target_name"]+"\t"+el["matrix_name"]+"\n")
                    pairs.append(el["target_name"]+"-"+el["matrix_name"])
                
def startQvalCalcPerFile(params):
    params = params.split(":")
    TMP_FOLDER = params[0]
    OUTPUT = params[1]
    THRESHOLD = float(params[2])
    finput = params[3]
    
    fileObj = open(TMP_FOLDER+"/"+finput)
    table = preprocess_table(fileObj)
    fileObj.close()
    table_qval = calculateQVal(table)
    qv_thresh = THRESHOLD
    filtered_table = qval_filter(table_qval,qv_thresh)
    printNetwork(filtered_table,OUTPUT+"_network/"+finput)
    printTable(filtered_table,OUTPUT+"/"+finput)

def startQvalCalc(TMP_FOLDER,OUTPUT,THRESHOLD):
    listFile = os.listdir(TMP_FOLDER)
    
    listParams = map(lambda x: ":".join([TMP_FOLDER,OUTPUT,str(THRESHOLD),x]),listFile)

    if not os.path.isdir(OUTPUT):
        os.makedirs(OUTPUT)
        os.makedirs(OUTPUT+"_network")

    p = Pool()
    dfs = p.map(startQvalCalcPerFile, listParams)

def mergeNetwork(NET_DIR,OUTPUT):
    listFile = os.listdir(NET_DIR)

    network_df = pd.DataFrame(columns=["Target_transcript","Target_gene","Source"])


    for n in listFile:
        fname = NET_DIR+"/"+n
        df = pd.read_csv(fname,sep="\t",header=0,index_col=False)
        df.columns = ["Target","Source"]
        df["Target_transcript"] = df["Target"].apply(lambda x: x.split("-")[-1])
        df["Target_gene"] = df["Target"].apply(lambda x: x.split("-")[0])
        df = df[["Target_transcript","Target_gene","Source"]]
        network_df = pd.concat([network_df,df])

    network_df.to_csv(OUTPUT.replace(".txt","_per_transcript.txt"),sep="\t",index=False)
    network_df = network_df[["Target_gene","Source"]]
    network_df = network_df.drop_duplicates()
    network_df.to_csv(OUTPUT,sep="\t",index=False)

def mergeFullNetwork(NETWORK_FILE,ANNOT_FILE,OUTPUT):
    network = pd.read_csv(NETWORK_FILE,sep="\t",header=0,index_col=None)
    network.columns = ["Target","Source"]
    network["Type"] = ["bind_to"]*len(network)
    annot_network = pd.read_csv(ANNOT_FILE,sep="\t",header=None,index_col=None)
    annot_network.columns = ["Target","Source","Type"]

    complete_network = pd.concat([network,annot_network])
    complete_network = complete_network[["Source","Target","Type"]]
    complete_network.to_csv(OUTPUT,sep="\t",index=False)


if __name__ == "__main__":
    FILENAME = sys.argv[1] # fimo result file name
    IS_ZIP = sys.argv[2] # True or False
    QVAL_THRESHOLD = sys.argv[3]
    ANNOT_FILE = sys.argv[4]
    OUTPUT = sys.argv[5]
    
    file_dir = "/".join(FILENAME.split("/")[:-1])
    file_name = FILENAME.split("/")[-1]
    file_name = ".".join(file_name.split(".")[:-1])
    
    tmp_folder = file_dir+"/tmp_"+file_name
    
    print "Split FIMO result file into smaller file..."
    splitCSV(FILENAME,tmp_folder,IS_ZIP)
        
    qval_folder = file_dir+"/qval_"+file_name
    
    print "Start calculating q-value for fimo result..."
    startQvalCalc(tmp_folder,qval_folder,QVAL_THRESHOLD)
    
    print "Generate network file from fimo result..."
    raw_network = file_dir+"/raw_network_"+file_name+".txt"
    mergeNetwork(qval_folder+"_network",raw_network)
    mergeFullNetwork(raw_network,ANNOT_FILE,OUTPUT)
    
    print "Removing all temporary file"
    os.system("rm -rf "+tmp_folder)
    os.system("rm -rf "+qval_folder)
    os.system("rm -rf "+qval_folder+"_network")
    os.system("rm -f "+raw_network)
    os.system("rm -f "+raw_network.replace(".txt","_per_transcript.txt"))
    
    print "Finish!"
