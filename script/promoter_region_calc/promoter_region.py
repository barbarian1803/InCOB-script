import GTF
import pprint
import sys
from operator import itemgetter
import os

def ReadFilterGene(file_gene_filter):
    gene_filter = {}

    file_object = open(file_gene_filter,"r")
    for line in file_object:
        line = line.rstrip()
        if line in gene_filter:
            gene_filter[line] = gene_filter[line]+1
        else:
            gene_filter[line] = 1
    return list(gene_filter.keys())

def CalculateAllPromoterRegions(result,promoter_length,gene_filter_list):
    promoter = {}
    for index, row in result.iterrows():
        if (row["gene_id"] in gene_filter_list) and (row["feature"]=="exon") and (row["exon_number"]=="1"):
            region = CalculatePromoterRegion(row,promoter_length)
            if row["gene_id"] not in promoter:
                promoter[row["gene_id"]] = [region]
            else:
                promoter[row["gene_id"]].append(region)

    return promoter

def CalculateAllPromoterRegions2(result,promoter_length,gene_filter_list):
    promoter = {}
    for index, row in result.iterrows():
        if (row["gene_id"] in gene_filter_list) and (row["feature"]=="exon") and (row["exon_number"]=="1"):
            region = CalculatePromoterRegion(row,promoter_length)
            key = row["gene_id"]+"-"+row["transcript_id"]
            if key not in promoter:
                promoter[key] = [region]
            else:
                promoter[key].append(region)

    return promoter

def CalculatePromoterRegion(row,promoter_length):
    chrom = row["seqname"]
    gene = row["gene_id"]
    strand = row["strand"]
    start = int(row["start"])
    end = int(row["end"])
    promoter_start = -1
    promoter_end = -1
    if strand=="+":
        promoter_start = start-promoter_length-1
        promoter_end = start-1
    else:
        promoter_start = end+1
        promoter_end = end+1+promoter_length

    if promoter_start<0:
        promoter_start = 0
    return {"chrom":chrom,"promoter_start":promoter_start,"promoter_end":promoter_end,"strand":strand}

def printToFile(promoter,filename):
    file_object = open(filename,"w")
    for name in promoter.keys():
        rows = promoter[name]
        for row in rows:
            file_object.write(row["chrom"]+"\t"+str(row["promoter_start"])+"\t"+str(row["promoter_end"])+"\t"+name+"\t"+"1\t"+row["strand"]+"\n")

def sortRegions(regions):
    return sorted(regions, key=itemgetter('promoter_start'))


def filterOverlaps(promoter):
    for name in promoter.keys():
        regions = sortRegions(promoter[name])
        toExclude = []
        for i in xrange(len(regions),1,-1):
            idx1 = i-1
            idx2 = i-2

            reg1 = regions[idx1]
            reg2 = regions[idx2]

            print "CalculateOverlap"
            print reg1
            print reg2

            lenreg1 = reg1["promoter_end"]-reg1["promoter_start"]+1
            lenreg2 = reg2["promoter_end"]-reg2["promoter_start"]+1
            distance_reg1_reg2 = reg1["promoter_end"]-reg2["promoter_start"]+1

            if distance_reg1_reg2 < lenreg1+lenreg2:    #overlap
                regions[idx2]["promoter_start"] = min(reg1["promoter_start"],reg2["promoter_start"])
                regions[idx2]["promoter_end"] = max(reg1["promoter_end"],reg2["promoter_end"])
                toExclude.append(idx1)

        regions = [elem for i, elem in enumerate(regions) if i not in toExclude]
        promoter[name] = regions
    return promoter

param_1= sys.argv[1] # gtf file
param_2= sys.argv[2] # promoter region length
param_3= sys.argv[3] # list gene to calculate
param_4= sys.argv[4] # output file
param_5= sys.argv[5] # per_transcript or per_gene
param_6= sys.argv[6] # human genome reference fasta


promoter_length = int(param_2)
file_gene_filter = param_3
file_output = param_4


print "Read GTF file into memory"
result = GTF.dataframe(param_1)

print "Read gene list"
gene_filter_list = ReadFilterGene(file_gene_filter)

print "Calculate promoter region to bed format"
if param_5=="per_gene":
    promoter = CalculateAllPromoterRegions(result,promoter_length,gene_filter_list)
    promoter = filterOverlaps(promoter)
else:
    promoter = CalculateAllPromoterRegions2(result,promoter_length,gene_filter_list)

print "Writing temporary bed file"
temp_file = file_output+".tmp.txt"
printToFile(promoter,temp_file)

print "Calculating promoter region finished, extracting fasta sequence"
os.system("bedtools getfasta -fi "+param_6+" -bed "+temp_file+" -s -name >> "+file_output)

print "Deleting temp file"
os.system("rm -f "+temp_file)

print "Finish..."
