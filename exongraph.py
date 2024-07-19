import gzip
import sys
import math

path = sys.argv[1]

introncount = 0
intronlist = []
cumulative = []

"""
def last(x, xlist): #in a sorted list, returns last occurence of a number x or less
	revlist = xlist[::-1]
	for i in range(len(xlist)):
		if revlist[i] <= x:
			return len(revlist)-i-1
		else: continue
	return 0

testlist = [1,1,3,4,4,5,5,5,6,8]
testcumulative = []
totalc = 0
for i in range(testlist[-1]):
	counts = testlist.count(i+1)
	totalcounts += counts
	testcumulative.append(totalc)
print(testcumulative)
"""
totalc = 0
with gzip.open(path, 'rt') as fp:
    for line in fp:
        if line[0] == '#' or line[0] == '>': continue
        if line[0].islower():
        	introncount += 1
        	intronlist.append(len(line))
    intronlist.sort()
    for i in range(intronlist[-1]):
    	c = intronlist.count(i+1)
    	totalc += c
    	#cumulative.append(totalc)
    	print(i+1, totalc/35110, sep=",")
    #print(cumulative)
