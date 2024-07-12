import math
import random

#function that creates a sample sequence of states/observations so I can compare
def hmm(length):
	states = ""
	observations = ""
	val = random.choice("ffl")
	for i in range(length):
		states += val
		if val == "f":
			val = random.choice("fffffffffffffffffffl")
			obs = random.choice("123456")
		if val == "l":
			val = random.choice("lllllllllf")
			obs = random.choice("1234566666")
		observations += obs
	return states, observations


states, observations = hmm(200)
print(observations)
print(states)




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


# in log because underflow
sampobs = observations
tempvitf = 0
tempvitl = 0
tempval = 0
nottempval = 0
total = []
for i in range(len(sampobs)):
	maxkf = max2(tempvitf+math.log2(a(tempval, "f")), tempvitl+math.log2(a(nottempval, "f")))
	maxkl = max2(tempvitf+math.log2(a(tempval, "l")), tempvitl+math.log2(a(nottempval, "l")))
	viterbif = math.log2(e("f", int(sampobs[i])))+maxkf
	viterbil = math.log2(e("l", int(sampobs[i])))+maxkl
	tempvitf = viterbif
	tempvitl = viterbil
	if viterbif >= viterbil:
		tempval = "f"
		nottempval = "l"
	else:
		tempval = "l"
		nottempval = "f"
	total.append(tempval)
finalstring = "".join(total)
print(finalstring)



#below is to see how accurately stuff matches
check = ""
for i in range(len(states)):
	if finalstring[i] != states[i]:
		check += "!"
	else:
		check += "."
print(check)
accuracy = check.count(".")
print(str((accuracy/len(states))*100)+"% accurate")
	
