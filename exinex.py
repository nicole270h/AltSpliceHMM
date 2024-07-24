import gzip
import sys

from grimoire.genome import Reader

efh = open('exons.fa', 'w')
ifh = open('introns.fa', 'w')
eie = open('eie.fa', 'w')

for chrom in Reader(sys.argv[1], sys.argv[2]):
	for gene in chrom.ftable.build_genes():
		if len(gene.transcripts()) == 0: continue
		tx = gene.transcripts()[0]
		if len(tx.introns) == 0: continue
		for i in range(1, len(tx.exons)):
			exon1 = tx.exons[i-1]
			exon2 = tx.exons[i]
			intron = tx.introns[i-1]
			if exon1.length < 40: continue
			if exon2.length < 40: continue
			if intron.length < 50: continue
			if intron.length > 500: continue
			
			efh.write(f'>{gene.id}.{tx.id}.exon{i}\n{exon1.seq_str()}\n')
			ifh.write(f'>{gene.id}.{tx.id}.intron{i}\n{intron.seq_str()}\n')
			
			if exon1.strand == '+':
				eie.write(f'>{gene.id}.{tx.id}.exon{i}-exon{i+1}\n')
				eie.write(f'{exon1.seq_str()[-40:]}\n')
				eie.write(f'{intron.seq_str().lower()}\n')
				eie.write(f'{exon2.seq_str()[:40]}\n')
			else:
				eie.write(f'>{gene.id}.{tx.id}.exon{i+1}-exon{i}\n')
				eie.write(f'{exon2.seq_str()[-40:]}\n')
				eie.write(f'{intron.seq_str().lower()}\n')
				eie.write(f'{exon1.seq_str()[:40]}\n')