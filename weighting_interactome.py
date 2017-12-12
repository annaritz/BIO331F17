
##Program for reducing the size of  the interactome in order to weight it appropriately

# /Users/kathythompson/Documents/school_stuff/college/bio_331/group_331/BIO331F17/weighting_interactome.py

def main():
    interactome = "interactome-flybase.txt" #interactome file
    positives = "positive-ids.txt" #Positive node file
    edges, nodes = read_edge_file(interactome) 
    positives = read_id_file(positives,nodes)
    edges, nodes = read_edge_file(interactome) 
    adj_list = make_adj_list(edges, nodes)
    nodes, edges = remove_by_dist(adj_list, positives) #only removing the nodes with distance 2 away, just to reduce the size of the file
    write_to_file(edges, "interactome_flybase_intscoreinput.txt")
    return

#Input: Text file containing edges in the interactome
#Output: Set of edges and set of nodes in the whole interactome
def read_edge_file(filename): ##taken from L.T.'s code and then edited by K.T (labtime), commented by EEK
    nodes = set()
    edges = set()
    with open (filename, 'r') as f:
        s = f.readline() #takes away header
        for line in f:
            k = line.strip().split() #this is only for the toy dataset. split on tabs otherwise
            k = k[0:2]
            k.append(1) ## EEK, adds edge weight 1 to every edge, used for calculating get_adj_list_with_weights
            edges.add(tuple(k))
            nodes.add(k[0])
            nodes.add(k[1])
    print('Number of edges: '+str(len(edges)))
    print('Number of nodes: '+str(len(nodes)))
    return edges,nodes
# Input: text file containing positive nodes, set of all nodes
# Output: set of all positive nodes in the interactome
def read_id_file(filename,nodes): #K.T (labtime)
    positives = set()
    with open (filename, 'r') as f:
        for line in f:
            k = line.strip().split()
            if k[0] in nodes:
                positives.add(k[0])
    print('Number of positive nodes: '+str(len(positives)))
    return positives
#Input: adjacency list, starting node, and connected component set (visited)
#Output: D - dictionary of number of visits per node, visited- the total number of nodes we were able to reach with BFS
def BFS(adj_list, s, visited): #K.T(labtime) ##from HW3.py
    LARGE_NUM = 100000000000
    D = {n:LARGE_NUM for n in adj_list}
    D[s] = 0
    q = [s]
    while len(q) != 0:
        w = q.pop(0)
        visited.add(w)
        for neighbor in adj_list[w]:
            if D[neighbor] == 100000000000:
                D[neighbor] = D[w]+1
                q.append(neighbor)
    return D,visited

##function returns an unweighted adjacency list
#Input: set of edges and set of nodes
# Output: adjacency list dictionary with nodes as keys and neighbor lists as values
def make_adj_list(edges,nodes): #K.T(labtime), but copied from Lab6 (anna)
    adj_list = {n:[] for n in nodes}  ## another way to initialize dictonaries
    for e in edges:
        adj_list[e[0]].append(e[1]) #use adjacency list
        adj_list[e[1]].append(e[0])
    return adj_list


#Runs BFS with every known positive node as a source node,
# adds if a node is within or equal to 4 units away
#Input: adjacency list, set of positives
#Output: set of nodes and set of edges containing nodes 4 or fewer paths from a positive node
def remove_by_dist(adj_list,positives): #K.T, with debugging done by the entire group
    print("Running remove_by_dist")##EEK
    nodes = set()
    visited = set()
    for p in positives:
        if p in adj_list: 
            test_distance, visited = BFS(adj_list, p , visited)
        for node in test_distance:
            if test_distance[node] <= 3:
                nodes.add(node)
    print('Number of processed nodes: ',len(nodes))
    #         edge.append(1) ## doing this to make Steiner work
    edges = set()
    seen = set() #seen keeps track of redundant nodes
    for v in adj_list:
        if v in nodes:
            for u in adj_list[v]:
                if u in nodes and u not in seen:
                    edges.add(tuple([v,u,1]))
        seen.add(v)
    print('Length of processed edges',len(edges))
    print("Done with remove_by_dist!")
    return nodes, edges
#Input: an edge list of the interactome, and then 
def write_to_file(edgels, filename):
    with open(filename, "w") as k:
        for edge in edgels:
            k.write(edge[0]+"\t"+edge[1]+"\n")
    return

if __name__ == '__main__':
    main()