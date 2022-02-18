


graph_file = "/mnt/d/breakpoints/assembly/simulation/test/test.graph.txt"
figure_file = "/mnt/d/breakpoints/assembly/simulation/test/test.graph.pdf"
out_file = "/mnt/d/breakpoints/assembly/simulation/test/test.solve.path.txt"

def read_graph():
    g = open(graph_file)
    graph = {}
    for line in g:  
        array = line.strip().split()
        node1 = array[0] + array[1]
        node2 = array[2] + array[3]
        if node1 not in graph:
            graph[node1] = [node2]
        else:
            graph[node1] += [node2]
    g.close()
    return graph


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

def update_graph(graph, record):
    new_graph = {}
    for node in graph:
        if record[node] != -1:
            continue
        new_graph[node] = []
        for near in graph[node]:
            if record[near] != -1:
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
        pathes.append(seg_list)
        record[node] = 0
        graph = update_graph(graph, record)
        if len(graph) > 0:
            begin_node, record, seg_list = given_graph(graph)
            BFS(graph, begin_node, seg_list, record)
        else:
            return 1

def remove_repeat_segs(pathes):
    flag = True
    if len(pathes) > 1:
        while flag:
            flag = False
            for i in range(len(pathes) - 1):
                if len(pathes[i]) < len(pathes[i+1]):
                    a = pathes[i]
                    pathes[i] = pathes[i+1]
                    pathes[i+1] = a
    f = open(out_file, 'w')
    record_exist = {}
    for path in pathes:
        for seg in path:
            if seg[:-1] not in record_exist:
                print (seg, end = "\t", file = f)
                record_exist[seg[:-1]] = 1
        print ('', file = f)
    f.close()

    # print (pathes)






graph = read_graph()
begin_node, record, seg_list = given_graph(graph)
# print ("begin", begin_node )
pathes = []
BFS(graph, begin_node, seg_list, record)
# print (pathes)
remove_repeat_segs(pathes)


















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