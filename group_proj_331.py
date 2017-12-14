from __future__ import print_function # (needed for python2 vs. python3)
from graphspace_python.api.client import GraphSpace
from graphspace_python.graphs.classes.gsgraph import GSGraph
from datetime import datetime

def main(): # EEK added comments to this
    #Input files

    # interactome = "toy_dataset.txt" #interactome file (added by KT)
    # positives = "toy_pos_set.txt" #Positive node file (added by KT)
    print('Start: ' + str(datetime.now()))

    interactome = "interactome-flybase.txt" #interactome file
    positives = "positive-ids.txt" #Positive node file


    #read interactome and positive node files
    edges, nodes = read_edge_file(interactome)
    positives = read_id_file(positives,nodes)

    #Make adjacency list from nodes and edges in interactome
    adj_list = make_adj_list(edges, nodes)

    #removes nodes more that 4 nodes away from any positive node, reassigns nodes and edges
    nodes, edges = remove_by_dist(adj_list, positives)

    print('Done with Pre-Processing: ' + str(datetime.now()))
    # #generates a steiner tree, and set of non terminal nodes, and adj_list
    
    #steiner_tree,nonterminal_ST_nodes,steiner_adj_list = SteinerApprox(nodes,edges,positives)

     # returns steiner tree nodes(from steiner edges out) as list of nodes
    #all_nodes = steiner_edges_out(steiner_tree,'tree_edges')

    # steiner non-positive terminals (list of nodes)
    #steiner_nodes_out(all_nodes, nonterminal_ST_nodes, 'tree_nodes')
    
    #steiner_adj_list_file(steiner_adj_list, 'steiner_adj_list')

    steiner_adj_list = adj_list_read('steiner_adj_list.txt')
    print('Done with Steiner Tree: ' + str(datetime.now()))
    # # runs BFS on the processed nodes, adj_list from the steiner tree, and positive set
    bfs_dict = bfs_rank(nodes,steiner_adj_list,positives)
    # BFS rank
    BFS_rank_out(bfs_dict,'BFS_rank')

    print('Done with BFS Rank: ' + str(datetime.now()))

    #Computes shortest paths given a node and adjacency list
    #pos_node_dict, SP_nonterminal_nodes = shortest_paths(nodes, edges, positives)
   
    # new_shortest_paths input (dictionary with key = non pos node, value = upstream pos node)
    #shortest_paths_out(pos_node_dict, 'new_shortest_paths')

    #print('Done with Shortest Paths: ' + str(datetime.now()))
    # #Reassigns nodes and edges to be a subgraph 
    # nodes,edges = select_subgraph_to_post(edges,nonterminal_ST_nodes,positives,steiner_tree,bfs_dict)

    
    
    # #Posts subgraph to GraphSpace
    # title = 'Interactome draft'+str(datetime.now())
    # post_graph(nodes,edges,nonterminal_ST_nodes,positives,steiner_tree,bfs_dict,title)
   

    print('Program Complete: ' + str(datetime.now()))
    





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
    #print('Number of edges: '+str(len(edges)))
    #print('Number of nodes: '+str(len(nodes)))
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
    #print('Number of positive nodes: '+str(len(positives)))
    return positives

##Network Pre-processing
#Update edges given a set of nodes
#removes all edges/nodes from the graph that do not include nodes in given set
#Input: visited - a connected component, set of edges
#Output: nes edge set that only containes edges in the connected component
def update_edges(visited,edges): #KT
    removing_edges = set()
    for edge in edges:
        if edge[0] not in visited:
            removing_edges.add(edge)
        if edge[1] not in visited:
            removing_edges.add(edge)
    edges = edges - removing_edges
    return edges


##Not used
##Checks to see if the positives and negatives are in the graph
#Input: graph(nodes and edges) and sets of positive and negative nodes
#Output: modifies pos and neg by removing nodes not in graph
def check_pos_negs(graph, positives, negatives): # K.T.
    for a in positives:
        if a not in nodes:
            trash = pos.pop(a)
    for a in negatives:
        if a not in nodes:
            trash.neg.pop(a)


##make sure the graph is connected, if not, takes the largest component by running BFS
#Input: adjacency list and list of nodes
#Output: set of nodes in the connected component
def check_connected(adj_list, nodes): #K.T(labtime)
    visited = set()
    for node in nodes:
        if node not in visited:
            distances,visited = BFS(adj_list, node, visited)
    return visited


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
    #print("Running remove_by_dist")##EEK
    nodes = set()
    visited = set()
    for p in positives:
        if p in adj_list: 
            test_distance, visited = BFS(adj_list, p , visited)
        for node in test_distance:
            if test_distance[node] <= 4:
                nodes.add(node)
    #print('Number of processed nodes: ',len(nodes))
    #         edge.append(1) ## doing this to make Steiner work
    edges = set()
    seen = set() #seen keeps track of redundant nodes
    for v in adj_list:
        if v in nodes:
            for u in adj_list[v]:
                if u in nodes and u not in seen:
                    edges.add(tuple([v,u,1]))
        seen.add(v)
    #print('Length of processed edges',len(edges))
    #print("Done with remove_by_dist!")
    return nodes, edges


#Input: adjacency list
#Output: edge list (set)
def adj_to_edge(adj_list): ##labtime, possibly EEK
    edges = set()
    for a in adj_list:
        for n in adj_list[a]:
            edge = [a,n]
            if [n,a] not in edges:
                edges.add(edge)
    return edges

# Input: set of nodes, list of edges, and a set of terminals as inputs
# Output: the metric closure, which is composed of terminals for nodes, and weighted, minimum shortest distances as edges, and adj_list.
# Wyatt modified the dijkstra & we_adj_list functions, unnesting the latter so it is only run once.
def get_metric_closure(nodes,edges,terminals):
    mc_edges = [] #construct a list of metric closure edges
    pi_dict = {}
    distance_dict = {}
    adj_list = get_adj_list_with_weights(edges) #adj_list for edges of G
# The rest of the function builds a list of edges for the metric closure, using two for loops, such that every terminal node gets connected to every other.
    for v in terminals:
        D,pi = dijkstra(nodes,adj_list,v)
        pi_dict[v] = pi
        distance_dict[v] = D
        for u in terminals:
#this if Statement just checks to exclude redundancy and self-loops.
            mc_edge = [v,u,D[u]]
            same_edge = [u,v,D[u]]
            if mc_edge not in mc_edges and same_edge not in mc_edges and v != u:
                mc_edges.append(mc_edge) # and adds it to the MC edges list.
    return mc_edges,adj_list,pi_dict,distance_dict


## Function uses a dictionary pi (see dijkstra's algorithm) implicitly including starting node 's', and an ending node as the second argument.
#Input: 
#Output:
def get_path(pi,node):
    path = [node] ## path starts with the ending node (& works backwords)
##So long as there is a previous path, pi[path[0]] does not return None.  In that case, the loop ends.
    while pi[path[0]]:
## + is used to append previous node to the start of the list, so it will become path[0] on next iteration.
        path = [pi[path[0]]] + path
    return path



## Make an adjacency list that contains the weights of each edge.(Anna)
## e.g., for edge (u,v), you can access the weight of that edge
## with adj_list[u][v] OR adj_list[v][u]
## INPUT: 3-element list of edges [node1,node2,weight]
## OUTPUT: dictionary of dictionaries
def get_adj_list_with_weights(edges):
    adj_list = {}
    for u,v,w in edges: ## another way to specify elements of key
        ## We want to add the key-value pair (v,w) to adj_list[u].
        ## First see if u is a key in adj_list.
        if u not in adj_list:
            adj_list[u] = {}  ## add the key (value is a DICTIONARY)
        ## Add the key-value pair (v,w) to adj_list[u]
        adj_list[u][v] = w

        ## We want to add the key-value pair (u,w) to adj_list[v].
        ## First see if v is a key in adj_list.
        if v not in adj_list:
            adj_list[v] = {}  ## add the key (value is a DICTIONARY)
        ## Add the key-value pair (u,w) to adj_list[v]
        adj_list[v][u] = w

    return adj_list


## Code reused from Lab6, which was developed collaboratively in class. 
#Input: nodes and edges of graph G
#Output: the minimum spanning tree, which is a list of edges
def kruskal(nodes,edges):
    T = [] #spanning tree started as a list
    C = set() # set of connected components initialized
    for node in nodes:
        C.add(frozenset([node]))
    edges = sorted(edges, key=lambda x:x[2]) #in-line function selects the index to for the sorted function to look at
    for edge in edges:
    	if acyclic(edge[0],edge[1],C): #renamed based off T property
    		T.append(edge) #build spanning tree
    		update_c(edge[0],edge[1],C)
    return T
## code developed collaboratively in class.  Modifications noted in #
def acyclic(node1,node2,C):
	for item in C: ##for each item in a connected component
## check if new edge (node1,node2) creates a cycle
		if (node1 in item) & (node2 in item):
			return False # If so, the addition will not be acyclic
	return True

## Run Dijkstra's in the weighted, undirected graph. (from Anna)
## INPUT: set of nodes, 3-element list of edges [node1,node2,weight], source s
## OUTPUT: Dictionary of distances (D), Dictionary of predecessors (pi)
def dijkstra(nodes,adj_list,s):
    #print("Running Dijkstra's") #EEK
    ## Build adjacency list that contains the weights of the edge.
    ## e.g., for edge (u,v), you can access the weight of that edge
    ## with adj_list[u][v] OR adj_list[v][u]

    LARGE_NUM = 1000000 ## like "infinity" here.

    ## initialize distances dictionary D.
    D = {n:LARGE_NUM for n in nodes}

    ## initialize predecessor dictionary pi.
    pi = {n:None for n in nodes}
    ## set distance to s to be 0
    D[s] = 0

    ## Queue is a dictionary (slow implementation)
    ## This could be sped up with a proper priority queue,
    ## but is fine for this homework.
    ## The queue values start as the distances for each node.
    Q = {n:D[n] for n in nodes}

    while len(Q) > 0: ## While we haven't visited all the nodes...
        ## Find the node with the minimum weight.
        w = None
        for n in Q: ## for every node in the Queue...
            if w == None or Q[n] < Q[w]: ## if we haven't set w yet or n is better...
                w = n ## set w to be this node.

        ## remove w from queue
        del Q[w]

        ## Iterate through the neighbors of w
        for x in adj_list[w]:
            # print('x ', x)
            # print('w ', w)
            # print('Dx ', D[x])
            # print('Dw ', D[w])
            # print('adj ', adj_list[w][x])
            ## If the current distance to x is larger than coming from w, update
            if D[x] > D[w] + adj_list[w][x]:
                D[x] = D[w] + adj_list[w][x] ## update the distance
                pi[x] = w ## update the predecessor (we came from w)
                Q[x] = D[x] ## update the entry in the queue

    return D,pi



#Input: a list of nodes, edges, and terminal nodes L
#Output: the Steiner Tree of the graph as a set of edges and a list of Steiner Tree nodes
def SteinerApprox(nodes,edges,terminals): ##Miriam
    #print("Beginning Steiner Approximation") ##EEK
    # Following solves for weighted edges of the metric closure.  The adj_list is not dependent on a start node, so it is run once and passed throughout the algorithm.
    mc_edges,steiner_adj_list,pi_dict,distance_dict = get_metric_closure(nodes,edges,terminals)
    ## Following function reused from Lab6.  It returns the minimum spanning tree for the metric closure of G.
    Tmc = kruskal(terminals,mc_edges)
    # T will build the full Steiner tree as a list of edges.
    T = set()
    for edge in Tmc: #for each edge in the metric closure
    # dijkstra's is rerun to solve for pi, so previous paths can be reconstructed from 's' (edge[0]) to end (edge[1]).
        D,pi = dijkstra(nodes, steiner_adj_list,edge[0])
        P = get_path(pi,edge[1]) # Reconstructs subpath from 's' to end.
        for i in range(len(P)): # for each node in subpath P
            if i <= len(P)-2: # Up until the second to last index
                if tuple([P[i],P[i+1]]) not in T and tuple([P[i+1],P[i]]) not in T:
                    T.add(tuple([P[i],P[i+1]])) # Add the edge to T
    #print('steiner tree: '+str(T))
    nonterminal_ST_nodes = set()
    for i in T:
        if i[0] not in terminals:
            nonterminal_ST_nodes.add(i[0])
        if i[1] not in terminals:
            nonterminal_ST_nodes.add(i[1])
    #print('nonterminal_ST_nodes: '+str(nonterminal_ST_nodes))
    return T, nonterminal_ST_nodes, steiner_adj_list


## Function updates the connected component based on an input of previous connected components and two nodes (the latter are the new edge in the min spanning tree).  Returns updated connected component.
#Input: previous connected components, 2 nodes
#Output: the new connected component
def update_c(node1,node2,C):
## c1 and c2 will include a minimum of one new node
    c1 = []
    c2 = []
## for each connected components
    for item in C:
        if node1 in item: ## if a node is part of that connected component
            c1 = item ## store that connected component as c1
            continue
        if node2 in item: ##if a second node is part of that connected component
            c2 = item ## store that connected component as c1
    C.remove(c1) ## remove both stored connected
    C.remove(c2)
## and add the union of the two stored cc's into the set of cc sets, C.
    C.add(frozenset(c1).union(frozenset(c2)))
    return



#Input: path as a list of nodes, adjacency list
#Output: the path as a list of weighted edges
def path_to_edges(path, adjacency): ##Miriam
    path_edges = []
    for i in range(len(path)-1):
        node1 = path[i] #first node in the edge - source
        node2 = path[i+1] #second node in the edge
        weight = adjacency[node1][node2] #look up the weight of the path between node1 and node2
        path_edges.append([node1,node2,weight])
    return path_edges

#Input: list of nodes, an adj_ls, and list of terminal nodes
#Output: a BFS dictionary of distances from positive nodes calculated by dijkstra's
def bfs_rank(nodes,adj_list,terminals): ##Wyatt
    #print('Running BFS rank')
    bfs_dict = {}
    for node in nodes:
        if node not in terminals:
            bfs_dict[node]=0
    for t in terminals:
        D, pi = dijkstra(nodes,adj_list,t)
        for key in D:
            if key not in terminals:
                bfs_dict[key] += (1.0/D[key])
    normalize_bfs_rank(bfs_dict)
    #print('BFS ranking completed:'+str(bfs_dict))
    return bfs_dict


#Input: takes in a dictionary of nodes with ranks from bfs_rank
#Output: a normalized dictionary using the maximum rank
def normalize_bfs_rank(rank):##Wyatt
    mxm = 0
    for node in rank:
        if rank[node] > mxm:
            mxm = rank[node]
    for node in rank:
        rank[node] = rank[node]/mxm #normalizes the values
    return 

#Input: edges, non terminal nodes that were included in the steiner tree, 
#a list of positive nodes, the list of edges from the steiner tree, 
#and the bfs dictionary
#Output: a new graph integrating this information
def select_subgraph_to_post(edges,nonterminal_ST_nodes,positives,steiner_tree,bfs_dict):##Wyatt
    bfs_list = []
    subedges = set()
    for node in bfs_dict:
        bfs_list.append([node,bfs_dict[node]])
    subnodes = sorted(bfs_list, key=lambda x:x[1],reverse = True)
    subnodes = [item[0] for item in subnodes]
    if len(subnodes) >= 100:
        subnodes = subnodes[0:100]
    print('100 highest BFS rank: '+str(subnodes)) ##KT:sorts out the 100 best BFS ranked nodes
    subnodes = subnodes + list(positives)
    subnodes = set(subnodes).union(nonterminal_ST_nodes)
    for edge in update_edges(subnodes,edges):
        subedges.add(tuple(sorted(edge)[1:])) #creates a new connected component with these nodes
    subedges.union(steiner_tree) #adds the steiner tree
    print('subedges: '+str(subedges))
    return subnodes,subedges

## Template code provided by Anna in Lab4 pdf.
#Input: floats
#Output: a hexadecimal color combining the floats
def rgb_to_hex(red,green,blue):
    maxHexValue = 255
    r = int(red*maxHexValue)
    g = int(green*maxHexValue)
    b = int(blue*maxHexValue)
    RR = format(r, '02x')
    GG = format(g, '02x')
    BB = format(b, '02x')
    return '#'+RR+GG+BB


#Input: takes in a list of nodes, edges, nonterminal_ST_nodes, terminals, the edges of the steiner_tree, and the BFS rank, and a title
#Output: an uploaded graph
def post_graph(nodes,edges,nonterminal_ST_nodes,terminals,steiner_tree,BFS_rank,title): ##Collaborative
    ## connect to GraphSpace
    USERNAME = 'wgormley@reed.edu'
    PASSWORD = 'side_flop'
    if USERNAME == 'FILL IN':
        sys.exit('ERROR: add your username and password in the post_graph() function.  Exiting.')
    graphspace = GraphSpace(USERNAME,PASSWORD)

    # create Graph instance, set title and tags.
    G = GSGraph()
    m = 2
    G.set_name(title)
    G.set_tags(['Hw5'])
    for n in nodes:
        if n in nonterminal_ST_nodes:
            color = rgb_to_hex(0,1,0)
        if n in BFS_rank:
            color = rgb_to_hex(1*BFS_rank[n],0,0)
        if n in nonterminal_ST_nodes and n in BFS_rank:
            color = rgb_to_hex(1*BFS_rank[n],1,0)

        popup=None
        if n in terminals:
            color='#0C7999'
            popup = 'terminal node'
        G.add_node(n,label=n,popup=popup)
        G.add_node_style(n,color=color,shape='ellipse',height=30,width=30)
    for e in edges:
        G.add_edge(e[0],e[1])
        G.add_edge_style(e[0],e[1],edge_style='dotted',color='#B1B1B1')
    for e in steiner_tree:
        G.add_edge_style(e[0],e[1],color='#000000',width=2)
        G.add_edge_style(e[1],e[0],color='#000000',width=2)
    G.set_data(data={'Regulators':'Blue','top 100 BFS Rank':'Red','nonterminal_ST_nodes+':'green','Spanning Tree Edge':'Black','ST Node and BFS ranked':'Yellow (R+G)'})
    try:
        graph = graphspace.update_graph(G)
        print('updated graph with title',title)
    except:
        graph = graphspace.post_graph(G)
    print('posted graph with title',title)
    return


##New Formulation code (KT)
#This will compute shortest paths to a particular node
#keeps track of the positive it is going from in a dictionary--{key is non-positive node: value is upstream pos nodes it came from}
def shortest_paths(nodes,edges,terminals): ##KT, with help from Anna
    print("Beginning shortest paths call") 
    pos_node_dict={} #will keep track of what positive comes with each node
    # T will build the full Steiner tree as a list of edges.
    adj_list = get_adj_list_with_weights(edges) #adj_list for edges of G
    T = set()
    SP_nonterminal_nodes = set()
    D,pi = dijkstra(nodes, adj_list, "FBgn0265434") #We can make one single call because this is an undirected graph
    for node in terminals: #for each node in the node
    # dijkstra's is rerun to solve for pi, so previous paths can be reconstructed from 's' (edge[0]) to end (edge[1]).
        P = get_path(pi,node) # Reconstructs subpath from 's' to end.
        for i in range(len(P)): # for each node in subpath P
            if i <= len(P)-2: # Up until the second to last index
                if tuple([P[i],P[i+1]]) not in T and tuple([P[i+1],P[i]]) not in T:
                    T.add(tuple([P[i],P[i+1]])) # Add the edge to T
                    for i in T:
                        if i[0] not in terminals:
                            SP_nonterminal_nodes.add(i[0]) #adds the node if it hasn't seen in before
                            if i[0] not in pos_node_dict:
                                pos_node_dict[i[0]]=set() #adds the key to the dictionary if it doesn't exist yet
                            pos_node_dict[i[0]].add(node) #adds the positive we were using to the set of positives that reach this node
                        if i[1] not in terminals:
                            SP_nonterminal_nodes.add(i[1]) #again, adds the node if it hasn't seen it before
                            if i[1] not in pos_node_dict:
                                pos_node_dict[i[1]]=set() #adds the key to the dictionary if it doesn't exist yet
                            pos_node_dict[i[1]].add(node) #adds the positive we were using to the set of positives that reach this node
    print("pos_node_dict",pos_node_dict)
    print('nonterminal_ST_nodes: '+str(SP_nonterminal_nodes))

    return pos_node_dict, SP_nonterminal_nodes


'''
Elaine's output functions!
'''

##Input is a set of tuples (edges)
## Build list of nodes in the tree for use in nodes_out
## Output is two columns, one per node in edge
def steiner_edges_out(tree_edge_set, filename): # Steiner
    out_file = open(str(filename)+'.txt','w')
    all_nodes = set()
    out_file.write('Edge1'+'\t'+'Edge2'+'\n')
    for m in tree_edge_set:
        for i in range(len(m)):
            if m[i] not in all_nodes:
                all_nodes.add(m[i])
            if i == 0:
                out_file.write(str(m[i]) + '\t')
            else:
                out_file.write(str(m[i]) + '\n')
    out_file.close()
    print(all_nodes)
    return all_nodes
'''
This works now!
'''


##Input all_nodes is set of nodes from egdes_out and input non_pos_nodes is SET of non positive nodes
#output is two columns, one is the node and the other is whether it is a positive node (N/Y) 
def steiner_nodes_out(all_nodes, non_term_nodes, filename): # Steiner
    out_file = open(str(filename)+'.txt','w')
    out_file.write('Node'+'\t'+'Terminal(Y/N)'+'\n')
    for node in all_nodes:
        out_file.write(str(node) + '\t')
        if node in non_term_nodes:
            out_file.write('N' + '\n')
        else:
            out_file.write('Y' + '\n')
    out_file.close()
'''
This works too!
'''
        
#Input BFS_rank_list is a list of two item lists [[node,float],[node1, float1] ]
#Output is two columns, one is the node and the other is the BFS rank
def BFS_rank_out(BFS_rank_list, filename): # BFS rank 
    out_file = open(str(filename)+'.txt','w')
    out_file.write('Node'+'\t'+'BFS_Rank'+'\n')
    for m in BFS_rank_list:
        out_file.write(str(m) + '\t' + str(FS_rank_list[m]) + '\n')
    out_file.close()

'''
This works!
'''

## Input dict is a dictionary with key = non pos node, value = upstream pos node
## output is two columns, with one as non_pos_node and the other as upstream pos node
def shortest_paths_out(dict, filename): # new shortest paths
    out_file = open(str(filename)+'.txt','w')
    out_file.write('Node'+'\t'+'Up_pos_nodes'+'\n')
    for key in dict:
        out_file.write(str(key) + '\t')
        out_file.write(str(dict[key]) + '\n')
    out_file.close()

'''
This works!
'''

def steiner_adj_list_file(adj_list, filename):
    out_file = open(str(filename)+'.txt','w')
    out_file.write(str(adj_list))
    out_file.close()

def adj_list_read(filename):
    with open (filename, 'r') as f:
        for line in f:
            k = line 
    return k


if __name__ == '__main__':
    main()
