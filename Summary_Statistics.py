## Summary Statistics Modification of Hw3
## Wyatt Gormley
## New coded is meant to be commented with a single #, triple quotes are excluded.
## Remove by distance currently uses BFS, may need to switch to dijkstra's,
## Summary stats 1a and 1b include degree distribution and average AND per k,
## Which can be useful in the report.  Summary stat 2 reviews the quanitity of
## nodes and edges within path 0<r<7 from a positive, whih maxes out at 17715
## when all nodes of the largest connected omponent are included/
## Import Statements
from __future__ import print_function # (needed for python2 vs. python3)
import matplotlib.pyplot as plt
from copy import deepcopy
import numpy as np
from datetime import datetime
import math
print(str(datetime.now()))

def main():
    edges,nodes = read_edge_file("interactome-flybase.txt") # get graph
    positives = read_id_file("positive-ids.txt",nodes) # get pos set
    negatives = read_id_file("negative-ids.txt",nodes) # get neg set
    adj_list = get_adj_list_with_weights(edges) # get adj_list w/ weights
    # Sumary Stats 1a
    AND_dict = get_AND(adj_list) # construct dict {node1 : Ave. Nv dv,node2...}
    # Summary Stats 1b
    # get 'degree hist' dict of:
    # { k=1 : average AND for k=1, k=2 : ave AND for k=2,...}
    # and k_count dict of {k: No. of ks in set, k+1: No. of Ks, ...}
    # for the positive, negative, & full interactome sets
    pos_degree_hist,pos_k_count = get_histogram(AND_dict,adj_list,positives)
    neg_degree_hist,neg_k_count = get_histogram(AND_dict,adj_list,negatives)
    interactome_degree_hist,k_count = get_histogram(AND_dict,adj_list,nodes)
    # Plot these two data sets as scatterplots
    create_K_count_plot(k_count,pos_k_count,neg_k_count) #plots histogram
    create_K_AND_plot(interactome_degree_hist, pos_degree_hist,neg_degree_hist)#plots
    """create_plot(d,'path_hist') #histogram of shortest_path lengths created"""

    # The following keeps track of how many nodes (dict value) are some minimum
    # distance (numbers which become dictionary keys) from any positibe node.  Records some
    separation_nodes_dict,separation_edges_dict = create_separation_dict(adj_list,positives)
    create_separation_nodes_plot(separation_nodes_dict)
    create_separation_edges_plot(separation_edges_dict)
    return
#end of main function
## reads .txt file as a set of edge tupples (u, v, 1) and a set of nodes.
def read_edge_file(filename): ##taken from L.T.'s code and then edited by K.T (labtime)
    nodes = set()
    edges = set()
    with open (filename, 'r') as f:
        s = f.readline() #takes away header
        for line in f:
            k = line.strip().split() #this is only for the toy dataset. split on tabs otherwise
            k = k[0:2]
            k.append(1) ## EEK, adds edge weight 1 to every edge, used for calculating get_adj_list_with_weights
            if k[0] != k[1]:
                edges.add(tuple(k))
                nodes.add(k[0])
                nodes.add(k[1])
    print(s)
    print('len edges: '+str(len(edges)))
    print('len nodes: '+str(len(nodes)))
    return edges,nodes
## reads positive or negative .txt file, checking that each node is contained in
## the main graph, and then builds a set of read nodes to return.
def read_id_file(filename,nodes): #K.T (labtime)
    terminals = set()
    with open (filename, 'r') as f:
        for line in f:
            k = line.strip().split()
            if k[0] in nodes:
                terminals.add(k[0])
    print('len terminals: '+str(len(terminals)))
    return terminals
# Make an adjacency list that contains the weights of each edge.(Anna)
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
        adj_list[u][v] = w ## Add the key-value pair (v,w)
        ## We want to add the key-value pair (u,w) to adj_list[v].
        ## First see if v is a key in adj_list.
        if v not in adj_list:
            adj_list[v] = {}  ## add the key (value is a DICTIONARY)
        ## Add the key-value pair (u,w) to adj_list[v]
        adj_list[v][u] = w
    return adj_list
## Uses an adjacency list to return a dictionary {node: AND, ...} where AND
## indicates 'average neighbor degree.'
def get_AND(adj_list):
    print('Starting Summary Statistics 1a: Average Neighbor Degree of Nodes')
    AND_dict = {} ##initialize an ave. Nv degree {node:AND}
    for key in adj_list: ##for each node
        AND = 0 ##initialize ave. Nv degree to 0
        dv = len(adj_list[key]) ##store degree of a node
        for neighbor in adj_list[key]: ##for each neighboring node
            AND = AND + len(adj_list[neighbor]) ##increment
        if dv > 0: ##prevents dividing by 0
            AND = float(AND)/float(dv) ##complete AND computation for node.
        AND_dict[key] = AND
    print('AND_dict len: '+str(len(AND_dict)))##produces a histogram dict {k:average AND}
    return AND_dict
# Code rewritten from HW3.  Function uses adj_list, previously constructed Ave
# Nv dv dict, and a node set (terminals) to analyze, producing a histogram in
# the form { k : ave AND, ... }.  An assumption was made to remove 0 as a key to
# prevent log(0) error in plot, maybe this should get handled elsewhere?
def get_histogram(AND_dict,adj_list,terminals):
    print('Starting Summary Statistics 1b: K-Count and K-AND')
    histogram = {} # create output dictionary {k: ave AND}
    k_count = {} # records the number of k occurences
    for node in terminals: # for each node in the ave Nv degree dictionary
        k = len(adj_list[node]) # k equals the number of its neighbors, or dv
        if k not in histogram.keys(): # If such value is not yet in AND_dict
            histogram[k] = AND_dict[node] # create pair {k:AND}
            k_count[k] = 1 # and initialize k_count[k] to one
        else:
            histogram[k] += AND_dict[node] # continue calc ave AND
            k_count[k] +=1 # and continue recording k_count
    for k in histogram:# for each degree
        histogram[k] = float(histogram[k])/float(k_count[k]) # take the average
    if 0 in histogram.keys(): #if 0 is a k value in histogram
        del histogram[0] #delete it, as it will produce a log[0] error elsewhere
    print('lengths: '+str(len(histogram)))
    return histogram, k_count
#This is similar to the get_histogram() function, changes will be indicated
#Histogram does not compare all nodes to all other nodes, but each node (loop 1)
#to each other unpaired node (loop 2 and unpaired_nodes set)
"""def get_SPL(nodes,adj_list):
    histogram = {} # initialize a histogram dictionary of {SPL : count, ...}
    unpaired_nodes = set(nodes) # initizlize set to keep track o
    for node in nodes: # for each terminal or starting node
        distance = shortest_path(nodes,node,adj_list) #(equivalent to the 'value')
        unpaired_nodes.discard(node)                #in the get_histogram function
        for node in unpaired_nodes: #loop looks at unpaired nodes to generate dict
            if distance[node] not in histogram:
                histogram[distance[node]] = 1
            else:
                histogram[distance[node]] += 1
    if (len(nodes)+1) in histogram.keys(): #this removes the key for unconnected nodes
        del histogram[len(nodes)+1]     #which get initialized to len(nodes)+1 in distances function
    return histogram
"""
# Inputs adj_list and positive set
# Outputs dictionaries two dictionaries, the first of form {0:~104, 1:3,613,...}
# indiating there are 104 positives 0 distance away from the positives, 3,613
# a distance of 1 or less from any positive, which will continue up until key=7.
# Remove by distane is used so as to return data on edges as well.
def create_separation_dict(adj_list,positives):
    print('Starting Summary Statistics 2:\nSeparation')
    separation_nodes_dict = {0:len(positives)}
    separation_edges_dict = {} #Excludes 0:0 so log transformed plot does nor
    # raise an issue
    nodes = set() # sets needs to start off at 0 and enables the while loop
    removal_distance = 1
    while len(nodes) <= 17714: #17715 observed largest connected component
        nodes, edges = remove_by_dist(adj_list,positives, removal_distance)
        separation_nodes_dict[deepcopy(removal_distance)] = deepcopy(len(nodes))
        separation_edges_dict[deepcopy(removal_distance)] = deepcopy(len(edges))
        removal_distance += 1
    return separation_nodes_dict,separation_edges_dict

# Application may be a stretch, as x and y values are usually related as a series
# of n terms, i.e. (x1 y1), (x2, y2), . . . (xn, yn), while here x and y terms
# follow the log transformation of the k: ave-AND data set.  In this analysis,
# negatice slopes indicate disasortative networks, positive slopes indicated
# asortative networks. x_ and y_ indicate averages.
def get_P_correlation_coefficient(k_aveAND):
    sigma_x = float(sum(k_aveAND.keys())) # take the sum of the keys
    x_ = sigma_x/len(k_aveAND) # and divide by the No. of keys to get x average
    sigma_y = float(sum(k_aveAND.values())) # take the sum of the values
    y_ = sigma_y/len(k_aveAND) # and divide by the No. of values to get y ave.
    # pcc = num/(math.sqrt(denom_x)*math.sqrt(denom_y))
    numerator = 0 # = sigma (xi - x_average)*(yi - y_average)
    denominator_x = 0 # = sigma (xi - x_average)^2
    denominator_y = 0 # = sigma (yi - y_average)^2
    for key in k_aveAND: # perform sumation process
        xi_x_ = float(key-x_) # = (xi - x_average)
        yi_y_ = float(k_aveAND[key]-y_) # = (yi - y_average)
        numerator += (xi_x_ * yi_y_)
        denominator_x += (xi_x_ ** 2)
        denominator_y += (yi_y_ ** 2)
    pcc = numerator/(math.sqrt(denominator_x)*math.sqrt(denominator_y))
    return pcc
## Borrowed from group assignment
# Removal distance taken as an input argument
## Runs BFS with every known positive node as a source node,
# 'visited' had been part of the code shared with BFS, here it is removed.
def remove_by_dist(adj_list,positives,removal_distance): #K.T
    print("Running remove_by_dist")##EEK
    nodes = set() # Build new seet of include nodes
    for p in positives: # for each positive regulator
        test_distance = BFS(adj_list, p)# get a dict of dist to all other nodes
        for node in test_distance:## for each node
        ## if it is less than the removal distance to any positive regulator
            if test_distance[node] <= removal_distance:
                nodes.add(node) ##add that node to the new nodes set
    print('Length of processed nodes',len(nodes))
    edges = set()  ## Build new set of included edges
    seen = set() ## seen keeps track of redundant nodes
    for v in adj_list:  ##for each node in the adj_list (all in graph G)
        if v in nodes: ## verify it is in the new nodes set
            for u in adj_list[v]: ## for each of v's neighbors, u
            ## if current u was not prior v and u is in new nodes
                if u in nodes and u not in seen:
                    edges.add(tuple([v,u,1])) ##rebuild an edge in new edge set
        seen.add(v)  ## record that we've seen v to optimize our algorithm
    print('Length of processed edges',len(edges))
    print("Done with remove_by_dist!")
    return nodes, edges

##Used in remove by distance, perhaps this should get replaced by dijkstra's
def BFS(adj_list, s): ##K.T(labtime)
    LARGE_NUM = 100000000000 ## Initialize ~infinity
    D = {n:LARGE_NUM for n in adj_list} ## build a distance dict for each node
    D[s] = 0  ## Initilize starting node as 0 distance from itself
    q = [s]  ## add the start to the front of the queue!
    while len(q) != 0: ## While there is a queue
        w = q.pop(0) ## grab the first one off the stack
        for neighbor in adj_list[w]: ## and for each of its neighbors
            if D[neighbor] == 100000000000:  ## if their distance is ~infinity
                D[neighbor] = D[w]+1 ## make it just one more than current space
                q.append(neighbor) ## and add this neighbor to the queue
    return D  ## return the dctionary

def create_K_count_plot(k_count,pos_k_count,neg_k_count):    #histogram of degrees created
    print('Plotting Summary Statistics 1bI: K-Count')
    fig = plt.figure(figsize=(12,6))

    x1 = k_count.keys()
    y1 = [k_count[xval] for xval in x1]
    logx1 = [math.log(a) for a in x1]
    logy1 = [math.log(b) for b in y1]
    line1 = 'k:.'
    x2 = pos_k_count.keys()
    y2 = [pos_k_count[xval] for xval in x2]
    logx2 = [math.log(a) for a in x2]
    logy2 = [math.log(b) for b in y2]
    line2 = 'c:.'
    x3 = neg_k_count.keys()
    y3 = [neg_k_count[xval] for xval in x3]
    logx3 = [math.log(a) for a in x3]
    logy3 = [math.log(b) for b in y3]
    line3 = 'r:.'

    plt.subplot(1,2,1)  #subplot 1 uses values from each data set and a line code
    plt.xlim([0,100])
    plt.plot(x1,y1,line1)
    plt.plot(x2,y2,line2)
    plt.plot(x3,y3,line3)

    plt.xlabel('k')
    plt.ylabel('Number of Nodes')
    plt.title('Degree Distribution')

    plt.subplot(1,2,2)  #generates subplot 2.
    plt.plot(logx1,logy1,line1)
    plt.plot(logx2,logy2,line2)
    plt.plot(logx3,logy3,line3)

    plt.xlabel('log(x)')
    plt.ylabel('log(y)')
    plt.title('Degree Distribution (log)')

    plt.tight_layout()

    plt.savefig('Degree Distribution'+str(datetime.now())+'.png')
    print('wrote to '+'Degree Distribution'+'.png')
    return

#I had trouble with running a loop using matplotlib, likely user error.
#uses 'h_entry' to allow for generation of each histogram, degree, AND, & SPL.
def create_K_AND_plot(degree_histogram, pos_degree_hist,neg_degree_hist):
    print('Plotting Summary Statistics 1bII: K-average AND')
    fig = plt.figure(figsize=(12,6))

    x1 = degree_histogram.keys()
    y1 = [degree_histogram[xval] for xval in x1]
    logx1 = [math.log(a) for a in x1]
    logy1 = [math.log(b) for b in y1]
    line1 = 'k.'

    (m,b) = np.polyfit(logx1,logy1,1)
    pcc1 = get_P_correlation_coefficient(degree_histogram)
    print('Pearson correlation coefficient for full interactome, log transformed: '+str(round(pcc1,4)))
    print('linear fit of the log transformed full interactome set')
    print('y = '+str(round(m, 3))+'x + '+str(round(b,3)))
    lr1 = np.polyval([m,b],x1)
    line1lr = 'k.'

    x2 = pos_degree_hist.keys()
    y2 = [pos_degree_hist[xval] for xval in x2]
    logx2 = [math.log(a) for a in x2]
    logy2 = [math.log(b) for b in y2]
    line2 = 'c.'

    (n,d) = np.polyfit(logx2,logy2,1)
    pcc2 = get_P_correlation_coefficient(pos_degree_hist)
    print('Pearson correlation coefficient for positive set, log transformed: '+str(round(pcc2,4)))
    print('linear fit of the log transformed positive set')
    print('y = '+str(round(n, 3))+'x + '+str(round(d,3)))
    lr2 = np.polyval([n,d],x2)
    line2lr = 'c.'

    x3 = neg_degree_hist.keys()
    y3 = [neg_degree_hist[xval] for xval in x3]
    logx3 = [math.log(a) for a in x3]
    logy3 = [math.log(b) for b in y3]
    line3 = 'r.'

    (o,p) = np.polyfit(logx3,logy3,1)
    pcc3 = get_P_correlation_coefficient(neg_degree_hist)
    print('Pearson correlation coefficient for negative set, log transformed: '+str(round(pcc3,4)))
    print('linear fit of the log transformed negative set')
    print('y = '+str(round(o, 3))+'x + '+str(round(p,3)))
    lr3 = np.polyval([n,d],x3)
    line3lr = 'r--'

    plt.subplot(1,2,1)  #subplot 1 uses values from each data set and a line code
    plt.xlim([0,150])
    plt.plot(x1,y1,line1)
    plt.plot(x2,y2,line2)
    plt.plot(x3,y3,line3)

    plt.xlabel('k')
    plt.ylabel('average AND')
    plt.title('Average AND for k')

    plt.subplot(1,2,2)  #generates subplot 2.
    plt.plot(logx1,logy1,line1)
    plt.plot(logx2,logy2,line2)
    plt.plot(logx3,logy3,line3)
    #plt.plot(logx1, lr1, line1lr)
    #plt.plot(logx2, lr2, line2lr)
    #plt.plot(logx3, lr3, line3lr)

    plt.xlabel('log(x)')
    plt.ylabel('log(y)')
    plt.title('Average AND for k (log)')

    plt.tight_layout()

    plt.savefig('Average AND plot'+str(datetime.now())+'.png')
    print('wrote to '+'Average AND Histogram'+'.png')
    return

def create_separation_nodes_plot(separation_nodes_dict):
    print('Plotting Summary Statistics 2: Nodes returned using remove_by_dist')
    fig = plt.figure(figsize=(12,6))

    x1 = separation_nodes_dict.keys()
    y1 = [separation_nodes_dict[xval] for xval in x1]
    #logx1 = [math.log(a) for a in x1]
    logy1 = [math.log(b) for b in y1]
    line1 = 'k:.'

    plt.subplot(1,2,1)  #subplot 1 uses values from each data set and a line code
    plt.xlim([0,6])
    plt.plot(x1,y1,line1)

    plt.xlabel('Nodes of x Pathlength to Positives')
    plt.ylabel('Number of Nodes')
    plt.title('Separation Distribution')

    plt.subplot(1,2,2)  #generates subplot 2.
    plt.plot(x1,logy1,line1)

    plt.xlabel('Nodes of Pathlength x to Positives')
    plt.ylabel('log(y)')
    plt.title('Separation Distribution (log)')

    plt.tight_layout()

    plt.savefig('Node Count by Some Distance of Separation from Positives '+str(datetime.now())+'.png')
    print('wrote to '+'Node Count by Some Distance of Separation from Positives'+'.png')

def create_separation_edges_plot(separation_edges_dict):
    print('Plotting Summary Statistics 2: Edges returned using remove_by_dist')
    fig = plt.figure(figsize=(12,6))

    x2 = separation_edges_dict.keys()
    y2 = [separation_edges_dict[xval] for xval in x2]
    logx2 = [math.log(a) for a in x2]
    logy2 = [math.log(b) for b in y2]
    line2 = 'k:.'

    plt.subplot(1,2,1)  #subplot 1 uses values from each data set and a line code
    plt.xlim([0,6])
    plt.plot(x2,y2,line2)

    plt.xlabel('Distance of Separation from Positives')
    plt.ylabel('Number of Edges')
    plt.title('Separation Distribution')

    plt.subplot(1,2,2)  #generates subplot 2.
    plt.plot(x2,logy2,line2)

    plt.xlabel('log(x)')
    plt.ylabel('log(y)')
    plt.title('Separation Distribution (log)')

    plt.tight_layout()

    plt.savefig('Edge Count by Some Distance of Separation from Positives'+str(datetime.now())+'.png')
    print('wrote to '+'Edge Count by Some Distance of Separation from Positives'+'.png')

    return


if __name__ == '__main__':
    main()
