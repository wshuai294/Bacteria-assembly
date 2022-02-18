import sys

raw = sys.argv[1]
out = sys.argv[2]
mean_depth = sys.argv[3]
min_ratio = sys.argv[4]
min_depth = float(mean_depth) * float(min_ratio) 


f = open(out, 'w')
for line in open(raw):
    line = line.strip()
    array = line.split()
    dp = int(array[4])
    if dp < min_depth:
        continue
    print (line, file = f)

f.close()

