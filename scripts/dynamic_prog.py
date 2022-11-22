import copy
import sys






# graph = {"1+":["5+"], "2+":["5+"], "5+":["3+", "2+"]}
# score = {"1+":{"5+":-1}, "2+":{"5+":-1}, "5+":{"3+":-1, "2+":-1}}
# copy_number = {"1+":1, "2+":1, "3+": 1,  "4+": 1, "5+":2}

# graph = {"1+":["2+"], "2+":["3+"], "3+":["4+", "5+"], "4+":["3+"]}
# # score = {"1+":{"5+":-1}, "2+":{"5+":-1}, "5+":{"3+":-1, "2+":-1}}
# copy_number = {"1+":1, "2+":1, "3+": 2,  "4+": 1, "5+":1}

def read_graph():
    g = open(graph_file)
    graph = {}
    copy_number = {}
    for line in g:  
        array = line.strip().split()
        if array[0] == "SEG":
            copy_number[array[1]+"+"] = int(array[2])
            copy_number[array[1]+"-"] = int(array[2])
        else:
            node1 = array[1] + array[2]
            node2 = array[3] + array[4]
            if node1 not in graph:
                graph[node1] = [node2]
            else:
                graph[node1] += [node2]
    g.close()
    return graph, copy_number

def DFS(node, path, paths, copy_number):
    # print (paths, path, node)
    if copy_number[node] == 0:
        paths.append(path)
        # path.append(node)
        return 0

    path.append(node)
    # copy_number[node] -= 1
    if node[:-1]+"+" in copy_number:
        copy_number[node[:-1]+"+"] -= 1
    if node[:-1]+"-" in copy_number:
        copy_number[node[:-1]+"-"] -= 1

    if node not in graph:
        # print (path, paths)
        paths.append(path)
        # cnv = copy.deepcopy(copy_number)
        return 0

    for near in graph[node]:
        my_path = copy.deepcopy(path)
        my_copy_number = copy.deepcopy(copy_number)
        DFS(near, my_path, paths, my_copy_number)
        # path.pop()
        # path = path[:-1]

def find_start(graph):
    in_num = {}
    for node in graph:
        for near in graph[node]:
            if near not in in_num:
                in_num[near] = 1
            else:
                in_num[near] = 1
    for node in graph:
        if node not in in_num and copy_number[node] > 0:
            return node
    return list(graph.keys())[0]


def update_graph(graph, paths):
    new_graph = {}
    record = {}
    for path in paths:
        for p in range(1, len(path)):
            from_node = path[p-1] 
            to_node = path[p]
            if from_node not in record:
                record[from_node] = [to_node]
            else:
                record[from_node] += [to_node]
    for node in graph:
        if copy_number[node] == 0:
            continue
        elif node not in record:
            new_graph[node] = graph[node]
    
        else:
            new_next = []
            for near in graph[node]:
                if near not in record[node] and copy_number[near] > 0:
                    new_next.append(near)
            if len(new_next) > 0:
                new_graph[node] = new_next
    return new_graph

def get_longest(paths):
    longest_path = ''
    max_len = 0
    for p in paths:
        if len(p) >= max_len:
            max_len = len(p)
            longest_path = p
    return longest_path

def remove_repeat_segs(paths):
    flag = True
    if len(paths) > 1:
        while flag:
            flag = False
            for i in range(len(paths) - 1):
                if len(paths[i]) < len(paths[i+1]):
                    a = paths[i]
                    paths[i] = paths[i+1]
                    paths[i+1] = a
    f = open(out_file, 'w')
    for p in paths:
        # print ("len", len(p))
        for seg in p:
            if origin_copy_number[seg] > 0:
                if seg[:-1]+"+" in origin_copy_number:
                    origin_copy_number[seg[:-1]+"+"] -= 1
                if seg[:-1]+"-" in origin_copy_number:
                    origin_copy_number[seg[:-1]+"-"] -= 1
                print (seg, end = "\t", file = f)
                print (seg, end = "\t")
        print ('', file = f)
        print ('')
    f.close()



graph_file = sys.argv[1]
out_file = sys.argv[2]
final_paths = []
graph, copy_number = read_graph()
origin_copy_number = copy.deepcopy(copy_number)
for p in graph:
    print (p, len(graph[p]), copy_number[p])


while len(graph) > 0:
    
    path = []
    paths = []
    
    node = find_start(graph)
    # print (path, )
    print ("iteration*********", node, copy_number[node], graph[node])
    DFS(node, path, paths, copy_number)
    longest_pa = get_longest(paths)
    final_paths.append(longest_pa)
    graph = update_graph(graph, final_paths)
    # print (paths)
    # print ("paths:")
    for p in paths:
        print (len(p), p[:15])
# print (final_paths)
    # print (len(final_paths), len(graph))
print ("finish---------")
# print (final_paths)
remove_repeat_segs(final_paths)

