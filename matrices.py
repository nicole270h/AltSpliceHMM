import kmeroddsfunc
import gzip
import sys

filepath = sys.argv[1]
klen = int(sys.argv[2]) #can only be a integer between 1 and 4 because matrix()

kex, kin = kmeroddsfunc.f(filepath, klen)


def matrix(n): #n-dimensional matrix for any number n between 1 and 4
	if n == 1:
		return[0 for i in range(4)]
	elif n == 2:
		return[[0 for i in range(4)] for j in range(4)]
	elif n == 3:
		return[[[0 for i in range(4)] for j in range(4)] for k in range(4)]
	elif n == 4:
		return[[[[0 for i in range(4)] for j in range(4)] for k in range(4)] for l in range(4)]
	"""
	matrix = 0
	for i in range(n):
		matrix = [matrix for j in range(4)]
	return matrix
	"""
def c(x): #converts nt to index value to search in matrix
	if x == "a":
		return 0
	elif x == "c":
		return 1
	elif x == "g":
		return 2
	elif x == "t":
		return 3


#initializing the matrix with 0
kmer3ex = matrix(klen)
kmer3in = matrix(klen)

for key in kex:
	for i in range(len(key)):
		print(key[i])
	print("-------")


