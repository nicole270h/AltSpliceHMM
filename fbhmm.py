import argparse
import gzip
import json
import math
import statistics
import sys
import numpy as np
import matplotlib.pyplot as plt

#call with python3 fbhmm.py dm2.hmm (fasta file)

def log2(a):
	if a != 0:
		return math.log2(a)
	elif a == 0:
		return -9999.0

def logsum(logp, logq): #given log p and log q, return log(p+q)
	if abs(logp-logq) > 30: return max(logp, logq)
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
		elif 'n' in line or 'N' in line:
			if len(seqs) > 0:
				seqs.pop()
		else:
			seqs.append(line)

	yield(name, ''.join(seqs))
	fp.close()


#in log and new
def fwdhmmlog(sampobs):
	fwd = {}
	for state in hmm['states']:
			fwd[state] = [-9999.0]*emm
	fwd["exon1"] = [0.0]*emm
	#^ created dict, every state has its own list which will have values appended to items
	for i in range(emm, len(sampobs)): #length of sequence minus emm length, because exon2 is where it ends
		for nextstate in fwd: #l
			totalsum = -9999
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

def bwdhmmlog(sampobs):
	bwd = {}
	for state in hmm['states']:
		bwd[state] = [-9999.0]*emm
	bwd["exon2"] = [0.0]*emm
	for i in range(len(sampobs)-emm-1, -1, -1): #moving backwards, so step is -1
		for prevstate in bwd: #k
			totalsum = -9999
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

def scatterplot(x, y, line, time):
	plt.ylim(-0.05, 1.05)
	plt.scatter(x, y, s=1)
	plt.axvspan(0,40, alpha=0.2, color='slategray')
	plt.axvspan(len(x)-40, len(x), alpha=0.2, color='slategray')
	plt.axvspan(40,45, alpha=0.1, color='lightseagreen')
	plt.axvspan(40,55, alpha=0.1, color='lightseagreen')
	plt.axvspan(len(x)-46, len(x)-40, alpha=0.1, color='lightseagreen')
	plt.axvspan(len(x)-55,len(x)-40, alpha=0.1, color='lightseagreen')
	plt.title(line)
	if time == True:
		plt.show(block = False)
		plt.pause(1)
		plt.close()
	else:
		plt.show()


## CLI ##

parser = argparse.ArgumentParser(description='Forward/Backward Decoder')
parser.add_argument('hmm', help='hmm model file')
parser.add_argument('fasta', help='fasta file')
parser.add_argument('--scatterplot', action='store_true',
	help='make graphs')
parser.add_argument('--textout', action='store_true',
	help='send text to stdout')
parser.add_argument('--debug', action='store_true',
	help='experimental')
arg = parser.parse_args()

# run with python3 fbhmm.py dm.hmm

hmm = json.load(open(arg.hmm))
emm = hmm['states']['exon1']
imm = hmm['states']['intron']


for defline, seq in read_fasta(arg.fasta):
	useq = seq.upper()
	fwd = fwdhmmlog(useq)
	bwd = bwdhmmlog(useq)
	xcoord = []
	ycoord = []
	for i in range(len(seq)):
		probsum = -9999
		#intronprob = fwd['intron'][i] + bwd['intron'][i]
		for state in fwd:
			probsum = logsum(probsum, fwd[state][i] + bwd[state][i])
		xcoord.append(i)
		ycoord.append(2**(fwd['intron'][i] + bwd['intron'][i] - probsum)) #change here
	if arg.debug:
		iscores = []
		for i in range(55, len(seq)-55):
			fbi = fwd['intron'][i] + bwd['intron'][i]
			fbe = max((fwd['exon1'][i] + bwd['exon1'][i], fwd['exon2'][i] + bwd['exon2'][i]))
			iscores.append(fbi-fbe)
		print(defline, len(seq), statistics.mean(iscores), statistics.stdev(iscores))
		#for i in range(len(iscores)):
		#	print(i, iscores[i])



	if arg.textout:
		for i in range(len(seq)):
			label = 'exon' if seq[i].isupper() else 'intron'
			print(i, seq[i], label, sep='\t', end='')
			for state in fwd:
				fb = fwd[state][i] + bwd[state][i]
				print(f'\t{fb:.1f}', end='')
			print("")
	xcoord = np.array(xcoord)
	ycoord = np.array(ycoord)
	if arg.scatterplot: scatterplot(xcoord, ycoord, defline, False)



"""
for defline, seq in read_fasta(arg.fasta):
	useq = seq.upper()
	fwd = fwdhmmlog(useq)
	bwd = bwdhmmlog(useq)
	print(f'>{defline}')
	cols = ['pos', 'nt', 'label']
	cols.extend(list(fwd.keys()))
	print('\t'.join(cols))
	print("")
	for i in range(len(seq)):
		label = 'exon' if seq[i].isupper() else 'intron'
		print(i, seq[i], label, sep='\t', end='')
		for state in fwd:
			fb = fwd[state][i] + bwd[state][i]
			print(f'\t{fb:.1f}', end='')
		print("")



"""


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

