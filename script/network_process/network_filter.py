## Script filter network generated from FIMO
## accept a parameter that is the folder of the network files

import sys
import os
import shutil
sys.path.append("..")
import numpy as np
import pandas as pd
from lib.Network import *
from lib.GeneralUtil import *
from lib.johnson import * 
 
 
folder = sys.argv[1]
files = os.listdir(folder)
for f in files:
    raw_network = Network(folder+"/"+f,"\t")
    raw_network.filterNetwork()
    fout = file(folder+"/filtered_"+f,"w")
    raw_network.printNetworkTofile(set(),fout)
    fout.close()
