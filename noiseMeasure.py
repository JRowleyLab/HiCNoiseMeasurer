#!/usr/bin/env python
import numpy as np
import argparse
import subprocess
import math
import statsmodels
from statsmodels import tsa
from statsmodels.tsa import stattools
import warnings
warnings.filterwarnings("ignore")

parser = argparse.ArgumentParser(description='Tool used to either downsample or upsample Hi-C reads\n\n')

parser.add_argument("-i", "--hicfile", dest="hic_file", help="Input .hic file(s), can be comma separated list\n\n", metavar="FILE", required=True)
parser.add_argument('--res', dest='res', metavar='str', help="Resolution(s) to process the Hi-C file, can be comma separated list\n\n", required=True)
parser.add_argument('-c', dest='chrom', default = "1", help="Chromosome to examine\n\n")
parser.add_argument('-j', dest="juicertools", help="Path to juicer tools\n\n", required=True)
parser.add_argument('-o', dest="outputfile", default=0, help="Path to output file. If none specified, will output to terminal\n\n")
args = parser.parse_args()

inputfiles = args.hic_file
myfiles = inputfiles.split(",")
myfiles = [str(x) for x in myfiles]

inputres = args.res
myres = inputres.split(",")
myres = [str(x) for x in myres]

fullnoiselist=[[]]*len(myres)
mychr = str(args.chrom)
startbin = 1000000
endbin = 10000000

 
for r in range(0,len(myres)):
    
    res=float(myres[r]) 
    binnum = int(float(endbin-startbin)/res)
    #endbin = startbin+(binnum*res)
    noiselist=[]
    intres=int(res)
    location=":".join([str(mychr),str(int(startbin)),str(int(endbin))])
    startoffset = int(startbin/res)

    for interactionfile in myfiles:
        print("processing " + str(interactionfile) + " at " + str(intres) + "  resolution")
        expecteds=[]
        with subprocess.Popen(['java', '-jar', args.juicertools, 'dump', 'expected', 'NONE', interactionfile, mychr, 'BP', str(intres)], stdout=subprocess.PIPE, stderr=subprocess.PIPE) as p:
            output, errors = p.communicate()
        lines = output.decode('utf-8').splitlines()
        for line in lines:
            myE = line.split()

            if myE[0] == "WARN":
                continue
            if myE[0] == "INFO":
                continue
            expecteds.append(float(myE[0]))
            #print("expecteds from juicer")
        matsize=int((endbin-startbin)/intres+1)
        dense=np.zeros((matsize,matsize))
        with subprocess.Popen(['java', '-jar', args.juicertools, 'dump', 'observed', 'NONE', interactionfile, location, location, 'BP', str(intres)], stdout=subprocess.PIPE, stderr=subprocess.PIPE) as p:
            output, errors = p.communicate()
        resultstmp = output.decode('utf-8').splitlines()

        for result in resultstmp:
            re = result.split()
            if re[0] == "INFO":
                continue
            if re[0] == "WARN":
                continue


            left = math.floor(int(re[0])/res) - startoffset
            right = math.floor(int(re[1])/res) - startoffset

            dist = int((int(re[1])-int(re[0]))/res)
            score = (float(re[2])+1)/(float(expecteds[dist])+1)
            dense[left,right] = score
            dense[right,left] = score
        noise=[]
        for i in range(1,len(dense)):
            try:
                noise.append(statsmodels.tsa.stattools.acf(dense[i-1:i,:], nlags=1)[1])
            except:
                continue
        noiselist.append(1/np.nanmean(noise))
    fullnoiselist[r] = noiselist
if args.outputfile == 0:
    for j in range(0,len(myres)):
        for h in range(0,len(myfiles)):
            
            print(str(myfiles[h]) + " at " + str(myres[j]) + " resolution, the noise is: " + str(fullnoiselist[j][h]))
else:
    with open(args.outputfile, "w") as ofile:
        for j in range(0,len(myres)):
            for h in range(0,len(myfiles)):

                linetowrite=(str(myfiles[h]) + "\t" + str(myres[j]) + "\t" + str(fullnoiselist[j][h]) + "\n")
                ofile.write(linetowrite)
