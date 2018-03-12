## adding differential expression data from cancer dataset to network.
## accept 2 parameters, input folder and aoutput folder
## cancer data is hard coded

import sys
import os
import shutil
import numpy as np
import pandas as pd
from lib.Network import *
from lib.GeneralUtil import *
from lib.johnson import *

NETWORK_INPUT_FOLDER = sys.argv[1]
NETWORK_OUTPUT_FOLDER = sys.argv[2]

network_files = os.listdir(NETWORK_INPUT_FOLDER)
if not os.path.isdir(NETWORK_OUTPUT_FOLDER):
    os.mkdir(NETWORK_OUTPUT_FOLDER)

def createNetwork(network_file,filename,colname,output):
    logfc = pd.read_csv(filename,sep="\t",header=0,index_col=0)
    logfc = logfc[~logfc.index.duplicated(keep='first')]
    raw_network = Network(network_file,"\t")

    for n in raw_network.getNodes():
        try:
            raw_network.getNode(n).logFC = logfc.loc[n,colname]
        except:
            continue
    # logFC file has been filtered by p-value to be less than <0.05, then node without logFC value (0 logFC) should be removed
    raw_network.removeNonDiffNode() 
    raw_network.filterNetwork()
    raw_network.printNetworkTofile(set(),file(output,"w"))
    return raw_network 

cancers = {}
cancers["bdc"] = "../../data/cancer_data/BDC.csv"
cancers["crc"] = "../../data/cancer_data/CRC.csv"
cancers["hcc"] = "../../data/cancer_data/HCC.csv"
cancers["luad"] = "../../data/cancer_data/LUAD.csv"
cancers["ossrc"] = "../../data/cancer_data/OSSRC.csv"

for network_file in network_files:
    for cancer in cancers:
        network_type = network_file.split(".")[0]
        output = NETWORK_OUTPUT_FOLDER+"/"+network_type+"/output_"+cancer+".csv"
        if not os.path.isdir(NETWORK_OUTPUT_FOLDER+"/"+network_type):
            os.mkdir(NETWORK_OUTPUT_FOLDER+"/"+network_type)
        createNetwork(NETWORK_INPUT_FOLDER+"/"+network_file,cancers[cancer],"log2FoldChange",output)

