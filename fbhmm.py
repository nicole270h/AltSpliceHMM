import argparse
import gzip
import json
import math
import sys


def log2(a):
	if a != 0:
		return math.log2(a)
	elif a == 0:
		return -99.0

def logsum(logp, logq): #given log p and log q, return log(p+q)
		logr = logp+log2(1+2**(logq-logp))
		return logr


def read_fasta(filename):
	"""iteratively read records from a FASTA file"""
	if   filename == '-':          fp = sys.stdin
	elif filename.endswith('.gz'): fp = gzip.open(filename, 'rt')
	else:                          fp = open(filename)
	name = None
	seqs = []
	while True:
		line = fp.readline()
		if line == '': break
		line = line.rstrip()
		if line.startswith('>'):
			if len(seqs) > 0:
				yield(name, ''.join(seqs))
				name = line[1:]
				seqs = []
			else:
				name = line[1:]
		else:
			seqs.append(line)

	yield(name, ''.join(seqs))
	fp.close()


#in log and new
def fwdhmmlog(sampobs):
	fwd = {}
	for state in hmm['states']:
			fwd[state] = [-99.0]*emm
	fwd["exon1"] = [0.0]*emm
	#^ created dict, every state has its own list which will have values appended to items
	for i in range(emm, len(sampobs)): #length of sequence minus emm length, because exon2 is where it ends
		for nextstate in fwd: #l
			totalsum = -99
			for nowstate in fwd: #k
				totalsum = logsum(totalsum, fwd[nowstate][i-1]+log2(hmm['transitions'][nowstate][nextstate])) #"i" is wrong here
			#now multiply totalsum with emission prob from x = hmm['emissions'][nextstate]
			if	 nextstate == 'exon1' or nextstate == 'exon2':	window = sampobs[i-emm:i+1]
			elif nextstate == 'intron':							window = sampobs[i-imm:i+1]
			else:												window = sampobs[i]
			forwardval = log2(hmm['emissions'][nextstate][window])+totalsum
			fwd[nextstate].append(forwardval)
	return fwd


def bwdemission(state, kmer):
	total = 0
	rest = kmer[1:]
	for nt in'ACGT':
		total += hmm['emissions'][state][nt+''.join(rest)]
	ans = hmm['emissions'][state][kmer]/total
	return ans
"""
def bwdhmm(sampobs):
	bwd = {}
	for state in hmm['states']:
		bwd[state] = [0.0]*emm
	bwd["exon2"] = [1.0]*emm
	for i in range(len(sampobs)-emm-1, -1, -1): #moving backwards, so step is -1
		for prevstate in bwd: #k
			totalsum = 0
			for nowstate in bwd: #l
				if 		prevstate == 'exon1' or prevstate == 'exon2':	window = sampobs[i:i+emm+1]
				elif 	prevstate == 'intron':						window = sampobs[i:i+imm+1]
				else: 												window = sampobs[i]
				eprob = bwdemission(prevstate,window)
				totalsum += hmm['transitions'][prevstate][nowstate]*eprob*bwd[nowstate][len(sampobs)-i-2]
			bwd[prevstate].append(totalsum)
	for state in bwd:
		bwd[state].reverse()
	return bwd
"""

def bwdhmmlog(sampobs):
	bwd = {}
	for state in hmm['states']:
		bwd[state] = [-99.0]*emm
	bwd["exon2"] = [0.0]*emm
	for i in range(len(sampobs)-emm-1, -1, -1): #moving backwards, so step is -1
		for prevstate in bwd: #k
			totalsum = -99
			for nowstate in bwd: #l
				if 		prevstate == 'exon1' or prevstate == 'exon2':	window = sampobs[i:i+emm+1]
				elif 	prevstate == 'intron':						window = sampobs[i:i+imm+1]
				else: 												window = sampobs[i]
				eprob = bwdemission(prevstate,window)
				totalsum = logsum(totalsum, log2(hmm['transitions'][prevstate][nowstate])+log2(eprob)+bwd[nowstate][len(sampobs)-i-2])
			bwd[prevstate].append(totalsum)
	for state in bwd:
		bwd[state].reverse()
	return bwd		


def chart(x): #formats the fwd or bwd list so it prints legibly
	for key in x:
		print(key, end="\t")
		for i in range(len(x[key])):
			print(f"{x[key][i]:.2f}", end="\t")
		print("")
	return ""

## CLI ##

parser = argparse.ArgumentParser(description='Forward/Backward Decoder')
parser.add_argument('hmm', help='hmm model file')
parser.add_argument('fasta', help='fasta file')
parser.add_argument('--local', action='store_true',
	help='allow local initializations')
parser.add_argument('--debug', action='store_true',
	help='show debugging output')
arg = parser.parse_args()

# run with python3 fbhmm.py dm.hmm

hmm = json.load(open(arg.hmm))
emm = hmm['states']['exon1']
imm = hmm['states']['intron']

for defline, seq in read_fasta(arg.fasta):
	useq = seq.upper()
	fwd = fwdhmmlog(useq)
	bwd = bwdhmmlog(useq)
	print(f'>{defline}')
	cols = ['pos', 'nt', 'label']
	cols.extend(list(fwd.keys()))
	print('\t'.join(cols))
	for i in range(len(seq)):
		label = 'exon' if seq[i].isupper() else 'intron'
		print(i, seq[i], label, sep='\t', end='')
		for state in fwd:
			fb = fwd[state][i] + bwd[state][i]
			print(f'\t{fwd[state][i]:.1f}', end='')
		print()

"""
print("", end = "\t")
for nt in sample:
	print(nt, end = "\t")
print("")

print("Forward:")
print(chart(fwdhmmlog(sample)))
print("Backward:")
print(chart(bwdhmmlog(sample)))


forward = fwdhmmlog(sample)
backward = bwdhmmlog(sample)
fb = {}
for key in forward:
	fb[key] = []
	for i in range(len(forward[key])):
		fb[key].append(forward[key][i]+backward[key][i])
print(chart(fb))
"""


				

	









#sample = "CTACTATTCGACATTTTCATGCGTCTCAATCTTCCGGACTgtgagtgtccctgattgaaattctcttcaattaacattgaacaattatcttactcagCTGCACCGCAACCGAATGGAGACAACGAATTGTCGCCTAA"
#sample = sample.upper()

"""
# not in log
def fwdhmm(sampobs):
	fwd = {}
	for state in hmm['states']:
			fwd[state] = [0.0]*emm
	fwd["exon1"] = [1.0]*emm
	#^ created dict, every state has its own list which will have values appended to items
	for i in range(emm, len(sampobs)): #length of sequence minus emm length, because exon2 is where it ends
		for nextstate in fwd: #l
			totalsum = 0
			for nowstate in fwd: #k
				totalsum += fwd[nowstate][i-1]*hmm['transitions'][nowstate][nextstate] #"i" is wrong here
			#now multiply totalsum with emission prob from x = hmm['emissions'][nextstate]
			#forwardval = hmm['emissions'][nextstate][sampobs[i]]*totalsum #something is wrong here
			if	 nextstate == 'exon1' or nextstate == 'exon2':	window = sampobs[i-emm:i+1]
			elif nextstate == 'intron':							window = sampobs[i-imm:i+1]
			else:												window = sampobs[i]
			forwardval = hmm['emissions'][nextstate][window]*totalsum
			fwd[nextstate].append(forwardval)
	return fwd
"""

