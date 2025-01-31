import argparse
import gzip
import mcb185
import json

parser = argparse.ArgumentParser(description='Markov model trainer')
parser.add_argument('file', help='fasta file')
parser.add_argument('order', type=int)
arg = parser.parse_args()

# counts
model = {}
for defline, seq in mcb185.read_fasta(arg.file):
	if 'N' in seq: continue
	for i in range(len(seq) - arg.order):
		ctx = seq[i:i+arg.order]
		nt = seq[i+arg.order]
		if ctx not in model: model[ctx] = {'A':0, 'C':0, 'G':0, 'T':0}
		model[ctx][nt] += 1

# probs
for ctx in model:
	total = sum(model[ctx].values())
	for nt in model[ctx]: model[ctx][nt] /= total

# out
#for ctx in sorted(model):
#	for nt in sorted(model[ctx]):
#		print(ctx, nt, model[ctx][nt])


for ctx in sorted(model):
	for nt in sorted(model[ctx]):
		out = {"ctx": ctx, "nt" : nt, "prob" : model[ctx][nt]}
		print(out)

