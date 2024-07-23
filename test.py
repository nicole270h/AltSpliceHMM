import gzip
import sys
"""
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

#FBgn0003559.FBtr0306545.exon9-exon10
yep = "CTACTATTCGACATTTTCATGCGTCTCAATCTTCCGGACTgtgagtgtccctgattgaaattctcttcaattaacattgaacaattatcttactcagCTGCACCGCAACCGAATGGAGACAACGAATTGTCGCCTAA"

for l in yep:
	print(l, end="\t")