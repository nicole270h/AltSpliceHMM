import gzip
import sys
import math


def f(path, klen):
	kex = {}
	kin = {}
	exoncount = 0
	introncount = 0
		
	with gzip.open(path, 'rt') as fp:
		for line in fp:
			invalidex = False
			invalidin = False
			if line[0] == '#' or line[0] == '>': continue
			if line[0].isupper():
			 exoncount += len(line)
			else:
				introncount += len(line)
			for i in range (len(line)-klen+1):
				kmer = line[i:i+klen]
				if "n" in kmer or "\n" in kmer: continue
				if kmer.isupper():
					exoncount += 1
					kmer = kmer.lower()
					if kmer not in kex:
						kex[kmer] = 0
					kex[kmer] += 1
					if kmer not in kin:
						kin[kmer] = 0
				else:
					if kmer not in kin:
						kin[kmer] = 0
					kin[kmer] += 1
					if kmer not in kex:
						kex[kmer] = 0
	
	klistex = list(kex.keys())
	knumex = list(kex.values())
	klistin = list(kin.keys())
	knumin = list(kin.values())
	logmin = 0
	logmax = 0
	for i in range(len(klistex)):
		logodds = math.log2((knumex[i]*(introncount-klen+1-knumin[i]))/(knumin[i]*(exoncount-klen+1-knumex[i])))
		if logodds > logmax:
			logmax = logodds
			logmaxdex = klistex[i]
		elif logodds < logmin:
			logmin = logodds
			logmindex = klistex[i]
		#print(str(klistex[i]),str(knumex[i]), str(knumin[i]), logodds, sep="\t")
	#print("logmin: "+str(logmin)+" at "+str(logmindex))
	#print("logmax: "+str(logmax)+" at "+str(logmaxdex))
	return kex, kin

