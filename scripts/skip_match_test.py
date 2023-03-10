import sys

graph_file = sys.argv[1]
path_file = sys.argv[2]

out = open(path_file, "w")
for line in open(graph_file):
    array = line.strip().split()
    if array[0] == "SEG":
        print (array[1] + "+", file = out)
out.close()
