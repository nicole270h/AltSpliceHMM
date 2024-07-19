import argparse
import gzip
import json
import mcb185

parser = argparse.ArgumentParser(description='Splice model trainer')
parser.add_argument('file', help='fasta file of intron')
parser.add_argument('--don', type=int, default=5,
	help='donor length [%(default)i')
parser.add_argument('--acc', type=int, default=6,
	help='acceptor length [%(default)i]')
arg = parser.parse_args()

# counts
don = [ {'A':0, 'C':0, 'G':0, 'T':0} for i in range(arg.don)]
acc = [ {'A':0, 'C':0, 'G':0, 'T':0} for i in range(arg.acc)]
total = 0
for defline, seq in mcb185.read_fasta(arg.file):
	if 'N' in seq: continue               # some introns have Ns
	if not seq.startswith('GT'): continue # canonical only
	if not seq.endswith('AG'): continue   # canonical only
	dseq = seq[:arg.don]
	aseq = seq[-arg.acc:]
	for i, nt in enumerate(dseq): don[i][nt] += 1
	for i, nt in enumerate(aseq): acc[i][nt] += 1
	total += 1

# probs
for i in range(arg.don):
	for nt in don[i]:
		don[i][nt] /= total
for i in range(arg.acc):
	for nt in acc[i]:
		acc[i][nt] /= total

# out
splice = {'donor': don, 'acceptor': acc}
print(json.dumps(splice, indent=4))
