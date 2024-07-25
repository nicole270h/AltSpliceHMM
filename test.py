import gzip
import sys
import numpy as np
import matplotlib.pyplot as plt


file = sys.argv[1]

def read_fasta(filename):
	fp = gzip.open(filename, 'rt')
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

for name, seq in read_fasta(file):
	print(name)
	print(seq)
	print("")
"""
x = [2, 10, 6]
x = np.array(x)
print(x)
"""
