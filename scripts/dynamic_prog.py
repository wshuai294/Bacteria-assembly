import copy
import sys




graph_file = sys.argv[1]
out_file = sys.argv[2]

# graph = {"1+":["5+"], "2+":["5+"], "5+":["3+", "2+"]}
# score = {"1+":{"5+":-1}, "2+":{"5+":-1}, "5+":{"3+":-1, "2+":-1}}
# copy_number = {"1+":1, "2+":1, "3+": 1,  "4+": 1, "5+":2}

graph = {"1+":["2+"], "2+":["3+"], "3+":["4+", "5+"], "4+":["6+"], "5+":["8+"]}
# score = {"1+":{"5+":-1}, "2+":{"5+":-1}, "5+":{"3+":-1, "2+":-1}}
copy_number = {"1+":1, "2+":1, "3+": 3,  "4+": 1, "5+":1, "6+":1, "8+":1}

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
    if copy_number[node] == 0:
        path.append(node)
        return 0
    # print (paths, node)
    path.append(node)
    copy_number[node] -= 1
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
        if node not in in_num:
            return node
    return list(graph.keys())[0]

def given_graph(graph):
    # print (graph)
    begin_node = find_start(graph)
    record = {}
    for node in graph:
        record[node] = -1
        for near in graph[node]:
            record[near] = -1
    seg_list = []
    return begin_node, record, seg_list

def update_graph(graph, paths):
    new_graph = {}
    record = {}
    for path in paths:
        for p in range(1, len(path)):
            from_node = path[p-1] 
            to_node = path[p]
    #         record[p] = 1
    # for node in graph:
    #     if node in record:
    #         continue
    #     if copy_number[node] == 0:
    #         continue
    #     new_graph[node] = []
    #     for near in graph[node]:
    #         if near in record:
    #             continue
    #         if copy_number[near] == 0:
    #             continue
    #         new_graph[node].append(near)
    # return new_graph

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

# node = "1+"
final_paths = []
# graph, copy_number = read_graph()
origin_copy_number = copy.deepcopy(copy_number)


while len(graph) > 0:
    print ("iteration*********")
    path = []
    paths = []
    
    node, record, seg_list = given_graph(graph)
    # print (path, )
    DFS(node, path, paths, copy_number)
    longest_pa = get_longest(paths)
    final_paths.append(longest_pa)
    graph = update_graph(graph, final_paths)
    # print (paths)
# print (final_paths)
    print (final_paths, graph, origin_copy_number)
print ("finish---------")
remove_repeat_segs(final_paths)

