import sys
from collections import defaultdict

def reverse_dir(dir):
    if dir == "+":
        return "-"
    else:
        return "+"

def build_graph(file_name):
    graph = defaultdict(list)
    with open(file_name, 'r') as file:
        for line in file:
            columns = line.strip().split()
            if columns[0].startswith("JUNC"):
                node1 = columns[1] + columns[2]
                node2 = columns[3] + columns[4]
                graph[node1].append(node2)
                graph[columns[3] + reverse_dir(columns[4])].append(columns[1] + reverse_dir(columns[2]))
                # graph[node2].append(node1)  # Uncomment this line if the graph is undirected
    return graph

def dfs(graph, start, end, visited=None):
    if visited is None:
        visited = set()
    if start not in graph or end not in graph:
        return None
    path = [start]
    visited.add(start)

    if start == end:
        return path

    for node in graph[start]:
        if node not in visited:
            sub_path = dfs(graph, node, end, visited)
            if sub_path is not None:
                return path + sub_path
    return None

def find_paths(file_a, file_b):
    graph_a = build_graph(file_a)
    graph_b = build_graph(file_b)
    paths = set()
    # Combine the graphs
    combined_graph = defaultdict(list, graph_a)
    for node, adjacents in graph_b.items():
        for adj in adjacents:
            new_path = dfs(graph_a, node, adj)
            if new_path != None and len(new_path) > 2:
                paths.add("\t".join(new_path))
    for item in paths:
        print(item)
if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python script.py <file_a> <file_b>")
    else:
        find_paths(sys.argv[1], sys.argv[2])
