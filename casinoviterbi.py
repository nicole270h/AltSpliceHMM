import math

sampobs = "1111111611111666661"

test = "315116246446644245311321631164152133625144543631656626566666"
sampstates = "FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFLLLLLLLLLLLLLLL"
#^ to test/compare with


def max2(a,b):
	if a >= b:
		return a
	if b > a:
		return b


k = ["f", "l"] #possible states

def a(prev, current):
	if prev == 0:
		if current == 0:
			return 1 #??????? otherwise my equation initialization is weird
		if current == "f":
			return 2/3 #initial probability of fair state
		if current == "l":
			return 1/3 #initial probability of fair state
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

"""
tempvit = 1
tempval = 0
for i in range(len(sampobs)):
	viterbif = e("f", int(sampobs[i]))*tempvit*a(tempval, "f")
	viterbil = e("l", int(sampobs[i]))*tempvit*a(tempval, "l")
	tempvit = max2(viterbif, viterbil)
	if tempvit == viterbif:
		tempval = "f"
	elif tempvit == viterbil:
		tempval = "l"
	print(i, sampobs[i], tempvit, tempval)
"""

"""
# this is the log version because underflow
tempvit = 1
tempval = 0
for i in range(len(sampobs)):
	viterbif = math.log(e("f", int(sampobs[i])))+tempvit+math.log(a(tempval, "f"))
	viterbil = math.log(e("l", int(sampobs[i])))+tempvit+math.log(a(tempval, "l"))
	tempvit = max2(viterbif, viterbil)
	if tempvit == viterbif:
		tempval = "f"
	elif tempvit == viterbil:
		tempval = "l"
	print(i, sampobs[i], tempvit, viterbif, viterbil, tempval)
"""

tempvit = 1
tempval = 0
for i in range(len(sampobs)):
	viterbif = e("f", int(sampobs[i]))*tempvit*a(tempval, "f")
	viterbil = e("l", int(sampobs[i]))*tempvit*a(tempval, "l")
	tempvit = max2(viterbif, viterbil)
	if tempvit == viterbif:
		tempval = "f"
	elif tempvit == viterbil:
		tempval = "l"
	print(i, sampobs[i], tempvit, viterbif, viterbil, tempval)