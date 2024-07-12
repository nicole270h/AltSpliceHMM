import math


sampobs = "6611"

def max2(a,b):
	if a >= b:
		return a
	if b > a:
		return b

def logsum(logp, logq): #given log p and log q, return log(p+q)
	logr = logp+math.log2(1+2**(logq-logp))
	return logr


def a(prev, current):
	if prev == 0:
		if current == 0:
			return 1 #??????? otherwise my equation initialization is weird
		if current == "f":
			return 2/3 #initial probability of fair state
		if current == "l":
			return 1/3 #initial probability of loaded state
	if prev == "f":
		if current == "f":
			return 0.95
		if current == "l":
			return 0.05
	if prev == "l":
		if current == "f":
			return 0.1
		if current == "l":
			return 0.9
#^ matrix for initialization & transition probabilities between fair and loaded

def e(l, xi):
	if l == "f":
		return 1/6
	if l == "l":
		if xi == 1 or xi == 2 or xi == 3 or xi == 4 or xi == 5:
			return 1/10
		if xi == 6:
			return 1/2
#^ emission matrix for probability that you see xi given state l

# in log
def fwd():
	sumf = 0
	suml = 0
	fwdflist = []
	fwdllist = []
	for i in range(len(sampobs)):
		fwdf = math.log2(e("f", int(sampobs[i])))+sumf
		fwdl = math.log2(e("l", int(sampobs[i])))+suml
		sumf = logsum(fwdf+math.log2(a("f", "f")), fwdl+math.log2(a("l", "f")))
		suml = logsum(fwdf+math.log2(a("f", "l")), fwdl+math.log2(a("l", "l")))
		fwdflist.append(fwdf)
		fwdllist.append(fwdl)
	return fwdflist, fwdllist


def bwd():
	prevf = 0
	prevl = 0
	bwdflist = []
	bwdllist = []
	for i in range(len(sampobs)):
		previ = int(sampobs[-1-i])
		subsumf1 = math.log2(a("f", "f"))+math.log2(e("f", previ))+prevf
		subsumf2 = math.log2(a("f", "l"))+math.log2(e("l", previ))+prevl
		bwdf = logsum(subsumf1, subsumf2)
		subsuml1 = math.log2(a("l", "f"))+math.log2(e("f", previ))+prevf
		subsuml2 = math.log2(a("l", "l"))+math.log2(e("l", previ))+prevl
		bwdl = logsum(subsuml1, subsuml2)
		prevf = bwdf
		prevl = bwdl
		bwdflist.append(bwdf)
		bwdllist.append(bwdl)
	bwdflist.reverse()
	bwdllist.reverse()
	return bwdflist, bwdllist
	
forwardfair, forwardload = fwd()
backwardfair, backwardload = bwd()

for i in range(len(forwardfair)):
	avgfair = 2**(logsum(forwardfair[i], backwardfair[i])-1)
	avgload = 2**(logsum(forwardload[i], backwardload[i])-1)
	print(avgfair, avgload)
