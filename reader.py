import json
import sys

fp = open(sys.argv[1])
sm = json.load(fp)
print(json.dumps(sm, indent=4))
