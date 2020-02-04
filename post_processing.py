import matplotlib.pyplot as plt
from matplotlib_venn import venn3
from graphspace_python.api.client import GraphSpace
from graphspace_python.graphs.classes.gsgraph import GSGraph
import sys

'''
Use of these functions assumes that the following files 
have been generated and are located in the same directory:

'tree_edges.txt'
'tree_nodes.txt'
'new_shortest_paths.txt'
'Dijkstra_rank.txt'

This code also requires:
'nodes-flybase.txt' for conversion from FlyBase IDs to common names.
'
'''

### ORIGINAL (uncompressed) INTERACTOME
# INTERACTOME_FILE = 'inputs/interactome-flybase.txt'
# POSITIVE_FILE = "inputs/positive-ids.txt" 
# OUTPUT_PREFIX = 'outputs/'

#### COLLASPED INTERACTOME
INTERACTOME_FILE = "inputs/interactome-flybase-collapsed-evidence.txt" #interactome file
POSITIVE_FILE = "inputs/positive-ids.txt" #Positive node file
OUTPUT_PREFIX="outputs/collapsed_"

TREE_NODE_FILE = OUTPUT_PREFIX+'tree_nodes.txt'
TREE_EDGE_FILE = OUTPUT_PREFIX+'tree_edges.txt'
SHORTEST_PATHS_FILE = OUTPUT_PREFIX+'new_shortest_paths.txt'
PATHS_RANK_FILE = OUTPUT_PREFIX+'shortestpaths_rank.txt'
SHORTEST_PATHS_FROM_SQH = OUTPUT_PREFIX+'shortest_paths_from_sqh.txt'
ALL_PAIRS_FILE = OUTPUT_PREFIX+'all_terminals_shortest_paths.txt'

def main():

	# creates a discionary of FB identifiers and common names
	common_names_dict = read_common_names('inputs/nodes-flybase.txt')

	# reads in output files
	steiner_terminal_nodes,steiner_non_terminal_nodes = read_tree_nodes(TREE_NODE_FILE)
	
	shortest_paths,sp_dict = read_new_shortest_paths(SHORTEST_PATHS_FILE)
	
	dijkstra,dij_dict = read_dijkstra_rank(PATHS_RANK_FILE)

	# converts node sets for each single output into common names
	steiner_nt_common = FB_to_common(steiner_non_terminal_nodes,common_names_dict)
	new_shortest_common = FB_to_common(shortest_paths, common_names_dict)
	dijkstra_common = FB_to_common(dijkstra, common_names_dict)

	make_venn(steiner_nt_common,dijkstra_common,new_shortest_common,'Nodes in Steiner Tree\napproximation','Nodes on paths to\nmany positives','Nodes on paths from Sqh to positives','outputs/venn.png')

	## connect to GraphSpace
	graphspace = GraphSpace('aritz@reed.edu', 'platypus')
	post_graphspace_graphs(common_names_dict,sp_dict,dij_dict,graphspace)

	print_to_screen(steiner_nt_common,new_shortest_common,dijkstra_common)

	return

def print_to_screen(steiner_nt_common,new_shortest_common,dijkstra_common):
	print('-'*100)

	# alphabetizes the common names for single outputs
	steiner_alpha = alphabetize(steiner_nt_common)
	dijkstra_alpha = alphabetize(dijkstra_common)
	shortest_alpha = alphabetize(new_shortest_common)

	# gets intersections of single outputs
	steiner_shortest_set = in_both_sets(steiner_nt_common,new_shortest_common)
	dijkstra_shortest_set = in_both_sets(dijkstra_common, new_shortest_common)
	steiner_dijkstra_set =  in_both_sets(dijkstra_common, steiner_nt_common)

	# alphabetizes the intersections 
	steiner_shortest_alpha = alphabetize(steiner_shortest_set)
	dijkstra_shortest_alpha = alphabetize(dijkstra_shortest_set)
	steiner_dijkstra_alpha = alphabetize(steiner_dijkstra_set)

	# gets intersection of all sets
	all_set = in_both_sets(steiner_shortest_set, dijkstra_common)

	# alphabetizes intersection of all sets
	all_set_alpha = alphabetize(all_set)

	#gets all proteins seen in any output
	union_of_steiner_dijkstra = steiner_nt_common.union(dijkstra_common)
	union_of_all = union_of_steiner_dijkstra.union(new_shortest_common) 
	union_alpha = alphabetize(union_of_all)
	

	# Prints alphabetized outputs
	print('\nPRINTING OUTPUTS')
	print()
	print('Steiner Candidates:', len(steiner_alpha))
	print('Dijkstra Candidates:', len(dijkstra_alpha))
	print('Shortest Paths to NM-II Candidates:', len(shortest_alpha))
	print()
	print('Steiner and Shortest Paths Candidates:', ','.join(steiner_shortest_alpha), len(steiner_shortest_alpha))
	print('Dijkstra and Shortest Paths : ',','.join(dijkstra_shortest_alpha),len(dijkstra_shortest_alpha))
	print('Steiner and Dijkstra Candidates: ', ','.join(steiner_dijkstra_alpha), len(steiner_dijkstra_alpha))
	print()
	print('Candidates present in all outputs:', ','.join(all_set_alpha), len(all_set_alpha))
	print()
	print('Union of all candidates: ', ','.join(union_alpha), len(union_alpha))

## Input file: 'nodes-flybase.txt'
## Outputs a dictionary with flybase IDs as keys and common names as values
def read_common_names(filename):
    name_dict = {}
    with open (filename, 'r') as f:
        for line in f:
            k = line.strip().split()
            if k[2] is not k[0]:
                name_dict[k[0]] = k[2]
    #print('Name dictionary:' + str(name_dict))
    return name_dict

## Input file: 'tree_nodes.txt'
## Outputs sets of terminal and non terminal nodes from the steiner tree 
def read_tree_nodes(filename):
	terminal_nodes = set()
	non_terminal_nodes = set()
	with open (filename, 'r') as f:
		s = f.readline() #takes away header 
		for line in f:
			k = line.strip().split()
			#print(k)
			if k[1] == 'Y':
				terminal_nodes.add(k[0])
			else:
				non_terminal_nodes.add(k[0])
	#print('Terminal Nodes:', terminal_nodes)
	#print('Non-terminal Nodes:', non_terminal_nodes)
	return terminal_nodes, non_terminal_nodes

def read_tree(filename):
	nodes = set()
	edges = set()
	with open(filename) as fin:
		fin.readline() # takes away header
		for line in fin:
			row = line.strip().split()
			nodes.add(row[0])
			nodes.add(row[1])
			edges.add((row[0],row[1]))
	return nodes,edges

def read_paths(filename):
	edges = set()
	nodes = set()
	with open(filename) as fin:
		for line in fin:
			#no header
			row = line.strip().split()
			nodes.update(set(row))
			for i in range(len(row)-1):
				u = row[i]
				v = row[i+1]
				if (u,v) not in edges and (v,u) not in edges:
					edges.add((u,v))
	return nodes,edges

## Input file: 'Dijkstra_rank.txt'
## Outputs set of nodes in Dijkstra
def read_dijkstra_rank(filename):
	Dijkstra_nodes = set()
	dijkstra_dict = {}
	with open (filename, 'r') as f:
		s = f.readline() #takes away header 
		for line in f:
			k = line.strip().split()
			if float(k[1]) > 0.7:
				Dijkstra_nodes.add(k[0])
				dijkstra_dict[k[0]] = float(k[1])
	#print('Dijkstra Nodes:', Dijkstra_nodes)
	return Dijkstra_nodes,dijkstra_dict


## Input file: 'new_shortest_paths.txt'
## Outputs set of nodes that are keys in new_shortest_paths
def read_new_shortest_paths(filename):
	shortest_paths_nodes = set()
	shortest_paths_dict = {}
	with open (filename, 'r') as f:
		s = f.readline() #takes away header 
		for line in f:
			k = line.strip().split()
			shortest_paths_nodes.add(k[0])
			shortest_paths_dict[k[0]] = int(k[1])
	#print('Shortest Paths Nodes:', shortest_paths_nodes)
	return shortest_paths_nodes,shortest_paths_dict

## converts FlyBase identifiers to common names
## leaces FlyBase id if there is no common name
def FB_to_common(node_set, common_names_dict):
	new_set = set()
	for node in node_set:
		if node in common_names_dict:
			new_set.add(common_names_dict[node])
		else:
			new_set.add(node)
	#print(new_set)
	return new_set

def read_graph(filename):
	edges = set()
	nodes = set()
	with open(filename) as fin:
		fin.readline() # ignore header
		for line in fin:
			row = line.strip().split()
			nodes.update(set([row[0],row[1]]))
			if row[0] != row[1]: # ignore self-loops for viz.
				edges.add((row[0],row[1]))
	return nodes,edges

## returns intersection of two sets (common elements)
def in_both_sets(set1, set2):
	both_set = set1.intersection(set2)

	return both_set

def alphabetize(set):
	sorted_list = []
	for item in set:
		sorted_list.append(item)
	sorted_list.sort()

	return sorted_list


def make_venn(list1,list2,list3,name1,name2,name3,figname):
	plt.figure(figsize=(6.5,5))
	venn3([list1,list2,list3],(name1,name2,name3))
	plt.title('Candidate proteins from each method',fontweight='bold',fontsize=14)
	plt.tight_layout()
	plt.savefig(figname)
	print("wrote to %s" % figname)
	return

def post_graphspace_graphs(common_names_dict,sqh_dict,rank_dict,graphspace): # these are hard-coded in.
	positives = set([s.strip() for s in open(POSITIVE_FILE).readlines()])
	positives = FB_to_common(positives,common_names_dict)
	print('%d positives' % (len(positives)))

	## TODO update popup text.
	post_steiner_graph(positives,common_names_dict,graphspace)

	post_sqh_graph(positives,common_names_dict,sqh_dict,graphspace)

	post_rank_graph(positives,common_names_dict,rank_dict,graphspace)
	
	return

def post_steiner_graph(positives,common_names_dict,graphspace):
	## get Steiner Tree
	tree_nodes,tree_edges = read_tree(TREE_EDGE_FILE)
	tree_nodes = FB_to_common(tree_nodes,common_names_dict)
	mapped_edges = set()
	for e in tree_edges:
		mapped_edges.add((common_names_dict.get(e[0],e[0]),common_names_dict.get(e[1],e[1])))
	tree_edges = mapped_edges
	print('%d nodes and %d edges from steiner tree'  % (len(tree_nodes),len(tree_edges)))
	post_graph(tree_nodes,tree_edges,positives,graphspace,'Steiner Tree Approximation')
	return 

def post_sqh_graph(positives,common_names_dict,sqh_dict,graphspace):
	## get paths from Sqh
	sqh_nodes,sqh_edges = read_paths(SHORTEST_PATHS_FROM_SQH)
	sqh_nodes = FB_to_common(sqh_nodes,common_names_dict)
	mapped_edges = set()
	for e in sqh_edges:
		mapped_edges.add((common_names_dict.get(e[0],e[0]),common_names_dict.get(e[1],e[1])))
	sqh_edges = mapped_edges
	mapped_sqh_dict = {}
	for key,value in sqh_dict.items():
		mapped_sqh_dict[common_names_dict.get(key,key)] = value
	sqh_dict =mapped_sqh_dict
	print('%d nodes and %d edges from sqh'  % (len(sqh_nodes),len(sqh_edges)))
	post_graph(sqh_nodes,sqh_edges,positives,graphspace,'Paths From Sqh',ranks=sqh_dict)
	return

def post_rank_graph(positives,common_names_dict,rank_dict,graphspace):
	all_nodes,all_edges = read_graph(INTERACTOME_FILE)
	#all_pairs_nodes,all_pairs_edges = read_paths(ALL_PAIRS_FILE)
	print('%d nodes and %d edges from interactome'  % (len(all_nodes),len(all_edges)))	
	all_nodes = FB_to_common(all_nodes,common_names_dict)
	mapped_edges = set()
	for e in all_edges:
		mapped_edges.add((common_names_dict.get(e[0],e[0]),common_names_dict.get(e[1],e[1])))
	all_edges = mapped_edges
	
	## map ranking dictionary
	mapped_rank_dict = {}
	for key,value in rank_dict.items():
		mapped_rank_dict[common_names_dict.get(key,key)] = value
	rank_dict =mapped_rank_dict
	sorted_items = sorted(rank_dict.items(), key=lambda x: x[1], reverse=True)
	#sorted_items = sorted_items[:10] # only take first 10 

	rank_nodes = positives.union(set([i[0] for i in sorted_items]))
	rank_edges = set()
	for e in all_edges:
		if e[0] in rank_nodes and e[1] in rank_nodes and ((e[0] not in positives and e[1] in positives) or (e[1] not in positives and e[0] in positives)):
			rank_edges.add(e)

	new_rank_nodes = set()
	for n in rank_nodes:
		s = sum([1 for e in rank_edges if e[0]==n or e[1]==n ])
		if s == 0:
			print(n,rank_dict.get(n,n),'has no edges!!  Ignoring if positive.')
			if n not in positives:
				new_rank_nodes.add(n)
		else:
			new_rank_nodes.add(n)
	new_rank_edges = set()
	for e in rank_edges:
		if e[0] in new_rank_nodes and e[1] in new_rank_nodes:
			new_rank_edges.add(e)
	rank_nodes = new_rank_nodes
	rank_edges = new_rank_edges

	print('%d nodes and %d edges from ranked nodes plus positives'  % (len(rank_nodes),len(rank_edges)))			
	post_graph(rank_nodes,rank_edges,positives,graphspace,'Shortest Paths Rank',ranks=rank_dict, weight_edges = True)
	return

def post_graph(nodes,edges,positives,graphspace,graph_name,ranks=None,weight_edges=False):
	sizes = {}
	if ranks != None:
		max_size = max(ranks.values())
		min_size = min(ranks.values())
		for key in ranks:
			sizes[key] = (ranks[key]-min_size)/(max_size-min_size)

	G = GSGraph()
	G.set_name(graph_name)
	G.set_tags(['NMII'])
	G.set_data(data={'description': 'Nodes are proteins and edges are protein interactions. Gray nodes are positives.'})

	for n in nodes:
		scale_factor = 1
		rank_num = 'unranked'
		norm_num = 'unranked'
		if ranks != None:
			scale_factor = 1+sizes.get(n,0)
			rank_num = str(ranks.get(n,'unranked'))
			norm_num = str(sizes.get(n,'unranked'))
		if n in positives:
			color = '#AAAAAA'
		else:
			color = '#AA2255'
		G.add_node(n,label=n,popup='scale factor = %f<br>rank = %s<br>normalized rank = %s' % (scale_factor,rank_num,norm_num))
		G.add_node_style(n,color=color,shape='ellipse',height=50*scale_factor,width=50*scale_factor)

	for e in edges:
		scale_factor = 1
		color = 'k'
		max_norm_rank = -1
		if ranks != None and weight_edges:
			max_norm_rank = max(sizes.get(e[0],0),sizes.get(e[1],0))
			scale_factor = max((1-max_norm_rank)**2,0.3)
			print(e[0],e[1],[sizes.get(e[0],0),sizes.get(e[1],0)],scale_factor)
			color_max = min(max_norm_rank,0.5)
			color = rgb_to_hex(color_max,color_max,color_max)
			popup_str = 'linewidth = %f <br> max_norm_rank = %f<br> line color = (%f,%f,%f)' % (scale_factor,max_norm_rank,color_max,color_max,color_max)
		else:
			popup_str = 'linewidth = %f' % (scale_factor)
		G.add_edge(e[0],e[1],popup=popup_str)
		G.add_edge_style(e[0],e[1],color=color,width=scale_factor)

	# post graph
	post(G,graphspace)
	print('Done posting %s' % (graph_name))
	return

def rgb_to_hex(r,g,b):
    return "#{:02x}{:02x}{:02x}".format(int(r*255),int(g*255),int(b*255))

def post(G,gs):
	try:
		graph = gs.update_graph(G)
	except:
		graph = gs.post_graph(G)

	return graph

if __name__ == '__main__':
	main()