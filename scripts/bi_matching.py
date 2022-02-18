


graph_file = "/mnt/d/breakpoints/assembly/simulation/test/test.graph.txt"
figure_file = "/mnt/d/breakpoints/assembly/simulation/test/test.graph.pdf"
out_file = "/mnt/d/breakpoints/assembly/simulation/test/test.solve.path.txt"

def read_graph():
    g = open(graph_file)
    graph = {}
    score = {}
    for line in g:  
        array = line.strip().split()
        node1 = array[0] + array[1]
        node2 = array[2] + array[3]
        if node1 not in graph:
            graph[node1] = [node2]
            score[node1] = {node2:-1}
        else:
            graph[node1] += [node2]
            score[node1][node2] = -1
    g.close()
    return graph, score

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
        for p in path:
            record[p] = 1
    for node in graph:
        if node in record:
            continue
        new_graph[node] = []
        for near in graph[node]:
            if near in record:
                continue
            new_graph[node].append(near)
    return new_graph

def BFS(graph, node, seg_list, record):
    seg_list.append(node)
    next_node = ''
    next_node_num = 0
    if node in graph:
        for near in graph[node]:
            if near == node or record[near] != -1:
                continue
            next_node_num += 1
            next_node = near
    if next_node_num == 1:
        record[node] = next_node
        BFS(graph, next_node, seg_list, record)
    else:
        # print (seg_list)
        path.append(seg_list)
        record[node] = 0
        graph = update_graph(graph, record)
        if len(graph) > 0:
            begin_node, record, seg_list = given_graph(graph)
            BFS(graph, begin_node, seg_list, record)
        else:
            return 1

def remove_repeat_segs(paths):
    print ("000")
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
                origin_copy_number[seg[:-1]+"+"] -= 1
                origin_copy_number[seg[:-1]+"-"] -= 1
                print (seg, end = "\t", file = f)
                print (seg, end = "\t")
        print ('', file = f)
        print ('')
    f.close()

    # print (path)

def simu_copy():
    copy_number = {}
    for key in graph:
        copy_number[key] = 1
        if key[:-1] == "NODE_6_length_604_cov_756.573055:1-604":
            copy_number[key] = 2
        for near in graph[key]:
            copy_number[near] = 1
            if near[:-1] == "NODE_6_length_604_cov_756.573055:1-604":
                copy_number[near] = 2
    return copy_number

def dynamic_programming(node, node_num):
    # print (node)
    if copy_number[node] == 0:
        return node_num
    copy_number[node] -= 1
    if node not in graph or len(graph[node]) == 0:      
        return node_num+1
    max_num = 0
    for near in graph[node]:
        test_node_num = dynamic_programming(near, node_num)
        # print (near, test_node_num)
        if test_node_num >= max_num:
            max_num = test_node_num
            next_node = near
    node_num = max_num + 1
    score[node][next_node] = node_num
    return node_num

def trace_back(begin_node, score):
    path = [begin_node]
    while True:
        if begin_node not in score:
            break
        next_node = ''
        next_num = -1
        for near in score[begin_node]:
            if score[begin_node][near] >= next_num:
                next_num = score[begin_node][near]
                next_node = near
        path.append(next_node)
        begin_node = next_node
    print (path)
    return path

    


graph, score = read_graph()
copy_number = simu_copy()
origin_copy_number = copy_number.copy()
paths = []

while len(graph) > 0:
    begin_node, record, seg_list = given_graph(graph)
    start_node = begin_node
    # print ("begin", begin_node, "-------------", graph)

    node_num = 0
    node_num = dynamic_programming(begin_node, node_num)
    path = trace_back(begin_node, score)
    
    paths.append(path)
    graph = update_graph(graph, paths)

print ("finish---------")
print (origin_copy_number)
remove_repeat_segs(paths)









# graph = read_graph()
# print (graph)
# copy_number = simu_copy()
# print (copy_number)
# begin_node, record, seg_list = given_graph(graph)
# # print ("begin", begin_node )
# path = []
# BFS(graph, begin_node, seg_list, record)
# # print (path)
# remove_repeat_segs(path)


















# M=[]
# class DFS_hungary():

#     def __init__(self, nx, ny, edge, cx, cy, visited):
#         self.nx, self.ny=nx, ny
#         self.edge = edge
#         self.cx, self.cy=cx,cy
#         self.visited=visited

#     def max_match(self):
#         res=0
#         for i in self.nx:
#             if self.cx[i]==-1:
#                 for key in self.ny:         # 将visited置0表示未访问过
#                     self.visited[key]=0
#                 res+=self.path(i)
#         print (cx)
#         return res

#     def path(self, u):
#         for v in self.ny:
#             if self.edge[u][v] and (not self.visited[v]):
#                 self.visited[v]=1
#                 if self.cy[v]==-1:
#                     self.cx[u] = v
#                     self.cy[v] = u
#                     M.append((u,v))
#                     return 1
#                 else:
#                     M.remove((self.cy[v], v))
#                     if self.path(self.cy[v]):
#                         self.cx[u] = v
#                         self.cy[v] = u
#                         M.append((u, v))
#                         return 1
#         return 0

# if __name__ == '__main__':
#     nx, ny = ['A', 'B', 'C', 'D'], ['E', 'F', 'G', 'H']
#     edge = {'A':{'E': 1, 'F': 0, 'G': 1, 'H':0}, 'B':{'E': 0, 'F': 1, 'G': 0, 'H':1}, 'C':{'E': 1, 'F': 0, 'G': 0, 'H':1}, 'D':{'E': 0, 'F': 0, 'G': 1, 'H':0}} # 1 表示可以匹配， 0 表示不能匹配
#     cx, cy = {'A':-1,'B':-1,'C':-1,'D':-1}, {'E':-1,'F':-1,'G':-1,'H':-1}
#     visited = {'E': 0, 'F': 0, 'G': 0,'H':0}

#     print (DFS_hungary(nx, ny, edge, cx, cy, visited).max_match())