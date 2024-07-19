import json
import sys

#run with python3 reader.py splice.model
#splice.model was taken from the output of python3 splicetrainer.py exons.fa.gz

#fp = open(sys.argv[1])
#file = json.load(fp)
#print(file)
#produces python dictionary with donor and acceptor site probabilities

#correct vers
fp = open(sys.argv[1])
sm = json.load(fp)
print(json.dumps(sm, indent=4))

