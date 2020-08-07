# Computational Investigation of Fog Signaling Pathway
#
# What it does: This code will generate a list of potential regulators from a protein-protein interactome and a file of positive regulators
# containing interactions between nodes using pre-processing techniques, Steiner Tree Approximations, shortestpaths-ranking, and a shortest paths algorithm. 
# 
# Input: a plain text file, formatted into at least 2 columns to indicate node to node interaction, a text file of regulators
# Output: Text file of the edge list of the Steiner Tree and a file containing all nodes within it , text files of potential regulators ranked by the shortestpaths-ranking
#         and text files of potential regulators calculated by the shortest-paths algorithm
# 2017 Authors: Miriam Bern, Wyatt Gormley, Elaine Kushkowski, Kathy Thompson, Logan Tibbetts, and Anna Ritz
# Runtime note: This code takes at least 3 hours to run.
#
# Updated Jan 2020 by Anna Ritz

from __future__ import print_function 
from datetime import datetime


### ORIGINAL (uncompressed) INTERACTOME
# INTERACTOME_FILE = "inputs/interactome-flybase.txt" #interactome file
# POSITIVES_FILE = "inputs/positive-ids.txt" #Positive node file
# OUTPUT_PREFIX="outputs/"
# NODE_COLS = [0,1]

### COLLAPSED INTERACTOME
INTERACTOME_FILE = "inputs/interactome-flybase-collapsed-evidence.txt" #interactome file
POSITIVES_FILE = "inputs/positive-ids.txt" #Positive node file
OUTPUT_PREFIX="outputs/collapsed_"
NODE_COLS = [3,4]

def main(): # EEK, KT added comments to this
    #Input files

    print('Start: ' + str(datetime.now()))    
    #read interactome and positive node files
    edges, nodes = read_edge_file(INTERACTOME_FILE)
    positives = read_id_file(POSITIVES_FILE,nodes)
    #positives = set([p for p in positives][:20])

    #Make adjacency list from nodes and edges in interactome
    adj_list = make_adj_list(edges, nodes)
    print('Done reading files and making adjacency list: ' + str(datetime.now()))

    #removes nodes more that 4 nodes away from any positive node, reassigns nodes and edges
    #nodes, edges = remove_by_dist(adj_list, positives)

    print('Done with Pre-Processing: ' + str(datetime.now()))

    # #generates a steiner tree, and set of non terminal nodes, and adj_list
    steiner_tree,nonterminal_ST_nodes,steiner_adj_list, pi_dict, distance_dict = SteinerApprox(nodes,edges,positives)

    # returns steiner tree nodes(from steiner edges out) as list of nodes
    all_nodes = steiner_edges_out(steiner_tree,OUTPUT_PREFIX+'tree_edges')
    steiner_nodes_out(all_nodes, nonterminal_ST_nodes, OUTPUT_PREFIX+'tree_nodes')
    print('Done with Steiner Tree: ' + str(datetime.now()))
    
    # shortestpaths rank
    shortestpaths_rank_dict,shortestpaths_rank_list = shortestpaths_rank(nodes, steiner_adj_list, positives, pi_dict, distance_dict)
    shortestpaths_rank_out(shortestpaths_rank_list,OUTPUT_PREFIX+'shortestpaths_rank')
    print('Done with shortestpaths Rank: ' + str(datetime.now()))

    #Computes shortest paths given a node and adjacency list
    pos_node_dict, SP_nonterminal_nodes = shortest_paths(nodes, edges, positives, pi_dict)
    shortest_paths_out(pos_node_dict, OUTPUT_PREFIX+'new_shortest_paths')
    print('Done with Shortest Paths: ' + str(datetime.now()))

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
            row = line.strip().split("\t")
            k = [row[NODE_COLS[0]],row[NODE_COLS[1]],1] ## EEK, adds edge weight 1 to every edge, used for calculating get_adj_list_with_weights 
            edges.add(tuple(k))
            nodes.add(k[0])
            nodes.add(k[1])
    print('%d edges and %d nodes' % (len(edges),len(nodes)))
    return edges,nodes

# Input: text file containing positive nodes, set of all nodes in the interactome
# Output: set of all positive nodes in the interactome
def read_id_file(filename,nodes): #K.T (labtime)
    positives = set()
    with open (filename, 'r') as f:
        for line in f:
            k = line.strip().split()
            if k[0] in nodes:
                positives.add(k[0])
    return positives

##Network Pre-Processing

#Update edges given a set of nodes
#removes all edges/nodes from the graph that do not include nodes in given set
#Input: visited - a connected component, set of edges
#Output: new edge set that only containes edges in the connected component
def update_edges(visited,edges): #KT
    removing_edges = set()
    for edge in edges:
        if edge[0] not in visited:
            removing_edges.add(edge)
        if edge[1] not in visited:
            removing_edges.add(edge)
    edges = edges - removing_edges
    return edges


##make sure the graph is connected, if not, takes the largest component by running BFS
#Input: adjacency list and list of nodes
#Output: set of nodes in the connected component
def check_connected(adj_list, nodes): #K.T(labtime)
    visited = set()
    for node in nodes:
        if node not in visited:
            distances,visited = BFS(adj_list, node, visited) #runs BFS on each node, and checks if we can reach it with breadth first search
    return visited


#Input: adjacency list, starting node, and connected component set (visited)
#Output: D - dictionary of number of visits per node, visited- the total number of nodes we were able to reach with BFS
def BFS(adj_list, s, visited): #K.T(labtime) ##from HW3.py
    LARGE_NUM = 100000000000
    D = {n:LARGE_NUM for n in adj_list} # assigns everything that we are considering a distance of "infinity"
    D[s] = 0 #initializes start node's distance to be 0
    pi = {n:None for n in adj_list}

    q = [s] #puts the start node in the queue to "search" for neighbors
    while len(q) != 0: #while we can reach something that we haven't seen
        w = q.pop(0) # remove the current node from the queue and find neighbors
        visited.add(w) # add that to "visited"
        for neighbor in adj_list[w]: #search for neighbors of that node
            if D[neighbor] == LARGE_NUM: # if we haven't seen it before
                D[neighbor] = D[w]+1 # reassign distance to be the distance from the node being considered, plus 1 (to account for it being a neighbor)
                pi[neighbor]= w
                q.append(neighbor) # now go through the neighbors of the neighbor
    return D,pi,visited 


##function returns an unweighted adjacency list
#Input: set of edges and set of nodes
# Output: adjacency list dictionary with nodes as keys and neighbor lists as values
def make_adj_list(edges,nodes): #K.T(labtime), but copied from Lab6 (anna)
    adj_list = {n:set() for n in nodes}  ## another way to initialize dictonaries
    for e in edges:
        adj_list[e[0]].add(e[1]) 
        adj_list[e[1]].add(e[0])
    return adj_list


#Runs BFS with every known positive node as a source node,
# adds if a node is within or equal to 4 units away
#Input: adjacency list, set of positives, max distance (default 4)
#Output: set of nodes and set of edges containing nodes 4 or fewer paths from a positive node
def remove_by_dist(adj_list,positives,max_dist=4): #K.T, with debugging done by all
    print("Running remove_by_dist")##EEK
    print('%d positives' % (len(positives)))

    ## make super-source node and add it to all positives.
    ss = 'supersource'
    adj_list[ss] = set()
    for p in positives:
        if p in adj_list:
            adj_list[p].add(ss)
            adj_list[ss].add(p)

    nodes = set()
    test_distance, test_pi, visited = BFS(adj_list, ss , set()) # gets the distance from each node in the adjacency list to the considered positive
    for node in test_distance: # for each node in the distance dictionary from that positive
        if node != ss and test_distance[node] <= max_dist+1: # ss adds an extra edge to distance.
            nodes.add(node)
    
    edges = set() # initializes new edge set
    seen = set() #seen keeps track of redundant nodes
    for v in adj_list: # for each node in the adjacency list
        if v != ss and v in nodes: # if it is in the new set of nodes
            for u in adj_list[v]: # for each neighbor of that node
                if u != ss and u in nodes and u not in seen: # if that neighbor is in the new set of nodes and it hasn't been seen
                    edges.add(tuple([v,u,1])) # add it to the new edge set
        seen.add(v) # then add the node we considered, because we have seen it now
    print('%d edges and %d nodes after removing nodes >%d from positives' % (len(edges),len(nodes),max_dist))
    return nodes, edges

#Input: adjacency list
#Output: edge list (set)
def adj_to_edge(adj_list): ##labtime, EEK
    edges = set()
    for a in adj_list: # for each node
        for n in adj_list[a]: #goes through all nodes that are neighbors of the top-level node 
            edge = [a,n] #creates an edge to show they are neighbors
            if [n,a] not in edges: #checks for duplicates of the edge created
                edges.add(edge) # adds it to the edge set
    return edges

# Input: set of nodes, list of edges, and a set of terminals as inputs
# Output: the metric closure, which is composed of terminals for nodes, and weighted, minimum shortest distances as edges, and adj_list.
# Wyatt modified the shortestpaths & we_adj_list functions, unnesting the latter so it is only run once.
def get_metric_closure(nodes,edges,terminals):
    mc_edges = [] #construct a list of metric closure edges
    pi_dict = {}
    distance_dict = {}
    adj_list = make_adj_list(edges, nodes)
    # The rest of the function builds a list of edges for the metric closure, using two for loops, such that every terminal node gets connected to every other.
    terminal_list = list(terminals)
    for i in range(len(terminal_list)):
        print('-->terminal #%d of %d' % (i,len(terminal_list)))
        v = terminal_list[i]
        D,pi,visited = BFS(adj_list, v, set())
        pi_dict[v] = pi
        distance_dict[v] = D
        for j in range(i+1,len(terminal_list)):
            u = terminal_list[j]            
            mc_edge = [v,u,D[u]]
            same_edge = [u,v,D[u]]
            if mc_edge not in mc_edges and same_edge not in mc_edges:
                mc_edges.append(mc_edge) # and adds it to the MC edges list.
    out = open(OUTPUT_PREFIX+'all_terminals_shortest_paths.txt','w')
    for i in range(len(terminal_list)):
        for j in range(i+1,len(terminal_list)):
            pi = pi_dict[terminal_list[i]]
            P = get_path(pi,terminal_list[j])
            out.write('\t'.join(P)+'\n')
    out.close()
    print('wrote to %sall_terminals_shortest_paths.txt' % (OUTPUT_PREFIX))
    return mc_edges,adj_list,pi_dict,distance_dict


## Function uses a dictionary pi (see shortestpaths's algorithm) implicitly including starting node 's', and an ending node as the second argument.
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
## Input: 3-element list of edges [node1,node2,weight]
## Output: dictionary of dictionaries
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

    return adj_list #AR

## Code reused from Lab6, which was developed collaboratively in class.
#Input: nodes and edges of graph G
#Output: the minimum spanning tree, which is a list of edges
def kruskal(nodes,edges):
    T = [] #spanning tree started as a list
    C = set() # set of connected components initialized
    for node in nodes:
        C.add(frozenset([node]))
    edges = sorted(edges, key=lambda x:x[2]) #in-line function selects the index for the sorted function to look at
    for edge in edges:
    	if acyclic(edge[0],edge[1],C): #checks if the tree created is acyclic
    		T.append(edge) #build spanning tree
    		update_c(edge[0],edge[1],C)
    return T
## code developed collaboratively in class.  Modifications noted in #

#Checks to see if a tree is acyclic, given a new edge and connected component
#Input: 2 nodes, and a connected component
#Output: boolean value evaluating if the new edge would introduce a cycle into the component
def acyclic(node1,node2,C):
	for item in C: ##for each item in a connected component
## check if new edge (node1,node2) creates a cycle
		if (node1 in item) & (node2 in item):
			return False # If so, the addition will not be acyclic
	return True

## Run shortestpaths's in the weighted, undirected graph. 
## INPUT: set of nodes, 3-element list of edges [node1,node2,weight], source s, set of targets (default is None)
## OUTPUT: Dictionary of distances (D), Dictionary of predecessors (pi). If target_set is specified, function will
## terminate when all nodes in target_set have been reached.
def shortestpaths(nodes,adj_list,s,target_set=None):
    #print("Running shortestpaths's") #EEK
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
    if target_set != None:
        seen_targets = set()

    while len(Q) > 0: ## While we haven't visited all the nodes...
        ## Find the node with the minimum weight.
        w = min(Q.keys(), key=Q.get)

        ## if this is a target, add it to the seen_targets. 
        if target_set != None and w in target_set:
            seen_targets.add(w)
        ## if this is the last target to find; return.
        if target_set != None and len(seen_targets) == len(target_set):
            #print('returning early b/c %d targets found' % (len(target_set)))
            return D,pi

        ## remove w from queue
        del Q[w]

        ## Iterate through the neighbors of w
        for x in adj_list[w]:
            ## If the current distance to x is larger than coming from w, update
            if D[x] > D[w] + adj_list[w][x]:
                D[x] = D[w] + adj_list[w][x] ## update the distance
                pi[x] = w ## update the predecessor (we came from w)
                Q[x] = D[x] ## update the entry in the queue

    return D,pi 



#Input: a list of nodes, edges, and terminal nodes L
#Output: the Steiner Tree of the graph as a set of edges and a list of Steiner Tree nodes
def SteinerApprox(nodes,edges,terminals): ##MB
    print("Beginning Steiner Approximation") ##EEK
    # Following solves for weighted edges of the metric closure.  The adj_list is not dependent on a start node, so it is run once and passed throughout the algorithm.
    mc_edges,steiner_adj_list,pi_dict,distance_dict = get_metric_closure(nodes,edges,terminals)
    print("Got metric Closure")
    ## Following function reused from Lab6.  It returns the minimum spanning tree for the metric closure of G.
    Tmc = kruskal(terminals,mc_edges)
    # T will build the full Steiner tree as a list of edges.
    print('MST')
    T = set()
    for edge in Tmc: #for each edge in the metric closure
        T = expand_edge(edge,pi_dict,T)
    print('expanded tree')
    nonterminal_ST_nodes = set()
    terminal_ST_nodes = set()
    for i in T:
        if i[0] not in terminals:
            nonterminal_ST_nodes.add(i[0])# each node part of the steiner tree and that is not a terminal node is added
        else:
            terminal_ST_nodes.add(i[0])

        if i[1] not in terminals:
            nonterminal_ST_nodes.add(i[1])
        else:
            terminal_ST_nodes.add(i[1])
    print('done')
    return T, nonterminal_ST_nodes, steiner_adj_list, pi_dict, distance_dict

def expand_edge(edge,pi_dict,T):
    tree_nodes = set([t[0] for t in T]).union([t[1] for t in T])
    path = get_path(pi_dict[edge[0]],edge[1]) # Reconstructs subpath from 's' to end.

    # get p_i and p_j as first and last verts 
    # already in T.
    p_i = None
    p_j = None
    for n in path:
        # if we haven't set p_i and n is in the nodeset, 
        # set p_i to this node.
        if not p_i and n in tree_nodes: 
            p_i = n
        # if n is in the nodeset, set p_j to this node.
        # note that this will get overwritten multiple times, 
        # but it will end with the correct node.
        if n in tree_nodes:
            p_j = n

    if p_i == None or p_i == p_j: 
        ## if there are 0 or 1 vertices, update the entire path.
        print('updating entire path')
        for i in range(len(path)-1):
            u = path[i]
            v = path[i+1]
            T.add(tuple([u,v])) # Add the edge to T

    elif p_i == edge[0] and p_j == edge[1] and not is_connected(T,edge[0],edge[1]):
        print('updating entire path -- trees are not connected')
        ## if u and v are the first and last nodes
        ## in the graph, BUT u and v aren't connected, then they
        ## are part of two different trees. Connect them with
        ## the entire path.
        for i in range(len(path)-1):
            u = path[i]
            v = path[i+1]
            T.add(tuple([u,v])) # Add the edge to T
    else:
        # there are 2 vertices so the path must be added in pieces.  
        # go through each edge in the path, and add it if we haven't
        # reached p_i OR we have reached p_j.
        print('adding in pieces with p_i=%s and p_j=%s' % (p_i,p_j))
        skip = False
        for i in range(len(path)-1):
            if path[i] == p_i: # we have reached p_i. Skip edges.
                skip = True
            if path[i] == p_j: # we have reached p_j, Stop skipping edges.
                skip = False
            
            if not skip:
                u = path[i]
                v = path[i+1]
                T.add(tuple([u,v])) # Add the edge to T
    return T

def is_connected(edges,u,v):
    ## This is a function for the EGFR case to merge trees.
    ## returns True if u and v are connected.
    conn_nodes = set([u]) # set of connected nodes in this "layer" of neighbors.
    old_neighbors = set() # all neighbors we've seen.
    while len(conn_nodes) > 0: # while there are more connected nodes to explore...
        # get new neighbors we haven't seen yet. If v is one of them, return true.
        new_neighbors = get_new_neighbors(edges,conn_nodes,old_neighbors)
        if v in new_neighbors:
            return True
        # update old_neighbors and conn_nodes.
        old_neighbors.update(conn_nodes)
        conn_nodes = new_neighbors
    # if we've reached here, v wasn't in any neighbor set. REturn false.
    return False

def get_new_neighbors(edges,node_set,old_neighbors):
    ## utility function to loop through edges and get all neighbors
    ## of the node set that haven't already been seen. 
    ## this could be made faster with a different data structure.
    neighbors = set()
    for u,v in edges:
        if u in node_set and v not in old_neighbors:
            neighbors.add(v)
        if v in node_set and u not in old_neighbors:
            neighbors.add(u)
    return neighbors

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
def path_to_edges(path, adjacency): ##MB
    path_edges = []
    for i in range(len(path)-1):
        node1 = path[i] #first node in the edge - source
        node2 = path[i+1] #second node in the edge
        weight = adjacency[node1][node2] #look up the weight of the path between node1 and node2
        path_edges.append([node1,node2,weight])
    return path_edges

#shortestpaths's Ranking
#Input: list of nodes, an adj_ls, list of terminal nodes,
# pi dict (dict of pi dicts), and distance dict (dict of D dicts)(distance dictionary from metric_closure function)
#Output: a shortestpaths's ranked dictionary proportional to distances from positives
def shortestpaths_rank(nodes,adj_list,terminals,pi_dict,distance_dict): ##WG
    print('Running shortestpaths rank')
    shortestpaths_rank_dict = {} # initializes the dictionary
    for node in nodes: # for each node in the node list 
        if node not in terminals:
            shortestpaths_rank_dict[node]=0 # if not a positive, assign a rank of 0
    for t in terminals: # for every node in the positives
        D = distance_dict[t] #references the distance from a particular positive to all other nodes in the steiner tree
        for key in D: # for each node distance from the positive
            if key not in terminals: # if the node is not a positive itself
                shortestpaths_rank_dict[key] += (1.0/D[key]) #assigns it a rank based on the distance to nodes
    shortestpaths_rank_list = normalize_shortestpaths_rank(shortestpaths_rank_dict) #creates a list of two-element lists containing [node, rank]
    print('shortestpaths ranking completed:')
    return shortestpaths_rank_dict,shortestpaths_rank_list


#Input: takes in a dictionary of nodes with score from shortestpaths_rank
#Output: an ordered list of ranked nodes (highest to lowest), with
#normalized (according to the maximum) scores
def normalize_shortestpaths_rank(rank):##WG
    mxm = 0 # Used to find maximum value
    for node in rank: # for each ranked node
        if rank[node] > mxm: # If its score is higher than current maximum
            mxm = rank[node] # Make maximum equal to that score
    for node in rank: # For every ranked node
        rank[node] = rank[node]/mxm #normalizes it by the maximum

    shortestpaths_rank_list = [] # Build a list to sort
    for node in rank: # For each ranked node
        shortestpaths_rank_list.append([node,rank[node]]) # Add [node,score] to the list
    shortestpaths_rank_list = sorted(shortestpaths_rank_list, key=lambda x:x[1],reverse = True) #return a sorted list according to rank values
    return shortestpaths_rank_list

##New Formulation code (KT)
#This will compute shortest paths to a particular node
#keeps track of the positive it is going from in a dictionary--{key is non-positive node: value is upstream pos nodes it came from}
def shortest_paths(nodes,edges,terminals,pi_dict): ##KT, with help from Anna
    ## variables for proteins
    ZIP = "FBgn0265434"
    SQH = "FBgn0003514"

    print("Beginning shortest paths call")
    out = open(OUTPUT_PREFIX+'shortest_paths_from_sqh.txt','w')
    pos_node_dict={} #will keep track of what positive comes with each node
    # T will build the full Steiner tree as a list of edges.
    adj_list = get_adj_list_with_weights(edges) #adj_list for edges of G
    T = set()
    SP_nonterminal_nodes = set()
    for node in terminals: #for each node in the node
    # shortestpaths's is rerun to solve for pi, so previous paths can be reconstructed from 's' (edge[0]) to end (edge[1]).
        P = get_path(pi_dict[SQH],node) # Reconstructs subpath from 's' to end.
        out.write('\t'.join(P)+'\n')
        for n in P: # for each node in path P 
            if n not in terminals:
                SP_nonterminal_nodes.add(n) #adds the node if it hasn't seen in before
                if n not in pos_node_dict:
                    pos_node_dict[n]=set() #adds the key to the dictionary if it doesn't exist yet
                pos_node_dict[n].add(node) #adds the positive we were using to the set of positives that reach this node
           
    print("%d in pos_node_dict" % len(pos_node_dict))
    print('%d nonterminal_ST_nodes: ' % len(SP_nonterminal_nodes))
    out.close()
    print("wrote to %sshortest_paths_from_sqh.txt" % (OUTPUT_PREFIX))
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
    out_file.write('#Edge1\tEdge2\n')
    for m in tree_edge_set:
        for i in range(len(m)):
            all_nodes.add(m[0])
            all_nodes.add(m[1])
            out_file.write('%s\t%s\n' % (m[0],m[1]))
    out_file.close()
    #print(all_nodes)
    return all_nodes


##Input all_nodes is set of nodes from egdes_out and input non_pos_nodes is SET of non positive nodes
#output is two columns, one is the node and the other is whether it is a positive node (N/Y)
def steiner_nodes_out(all_nodes, non_term_nodes, filename): # Steiner
    out_file = open(str(filename)+'.txt','w')
    out_file.write('#Node\tTerminal(Y/N)\n')
    for node in all_nodes:
        out_file.write(str(node) + '\t')
        if node in non_term_nodes:
            out_file.write('N' + '\n')
        else:
            out_file.write('Y' + '\n')
    out_file.close()

#Input BFS_rank_list is a list of two item lists [[node,float],[node1, float1] ]
#Output is two columns, one is the node and the other is the BFS rank
def shortestpaths_rank_out(BFS_rank_list, filename): # BFS rank
    out_file = open(str(filename)+'.txt','w')
    out_file.write('#Node\tBFS_Rank\n')
    for m in BFS_rank_list:
        for i in range(len(m)):
            if i == 0:
                out_file.write(str(m[i]) + '\t')
            else:
                out_file.write(str(m[i]) + '\n')

    out_file.close()


## Input dict is a dictionary with key = non pos node, value = upstream pos node
## output is two columns, with one as non_pos_node and the other as upstream pos node
def shortest_paths_out(up_dict, filename): # new shortest paths
    out_file = open(str(filename)+'.txt','w')
    out_file.write('Node\t#Up_pos_nodes\tUp_pos_nodes\n')
    for key in up_dict:
        out_file.write('%s\t%d\t%s\n' % (key,len(up_dict[key]),';'.join(up_dict[key])))
    out_file.close()



if __name__ == '__main__':
    main()
