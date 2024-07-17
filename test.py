import gzip
import sys




"""


path = sys.argv[1]
with gzip.open(path, 'rt') as fp:
    for line in fp:
        if line[0] == "a" or line[0] == "c" or line[0] == "t" or line[0] == "g":
        	print(line)
        else: continue






linecountex = 0
linecountin = 0
path = sys.argv[1]
klen = 3
kex = {}
kin = {}

validex = ["A", "C", "T", "G"]
validin = ["a", "c", "t", "g"]
for i in range(klen):

	




with gzip.open(path, 'rt') as fp:
    for line in fp:
        invalidex = False
        invalidin = False
        if line[0] == '#': continue
        if line[0] == '>': continue
        for i in range (len(line)-klen+1):
        	kmer = line[i:i+klen]
        	for j in range(len(kmer)):
        		if kmer[j] not in validex:
        			invalidex = True
        		if kmer[j] not in validin:
        			invalidin = True
        	if invalidex == False:
        		if kmer not in kex:
        			kex[kmer] = 0
        		kex[kmer] += 1
        		linecountex += 1
        	if invalidin == False:
        		if kmer not in kin:
        			kin[kmer] = 0
        		kin[kmer] += 1
        		linecountin += 1
"""