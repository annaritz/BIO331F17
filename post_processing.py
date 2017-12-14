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
def main():

	# creates a discionary of FB identifiers and common names
	common_names_dict = read_common_names('nodes-flybase.txt')

	# reads in output files
	steiner_terminal_nodes,steiner_non_terminal_nodes = read_tree_nodes('tree_nodes.txt')
	shortest_paths = read_new_shortest_paths('new_shortest_paths.txt')
	dijkstra = read_dijkstra_rank('Dijkstra_rank.txt')

	# converts node sets for each single output into common names
	steiner_nt_common = FB_to_common(steiner_non_terminal_nodes,common_names_dict)
	new_shortest_common = FB_to_common(shortest_paths, common_names_dict)
	dijkstra_common = FB_to_common(dijkstra, common_names_dict)

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
	print('Steiner Candidates:', steiner_alpha, len(steiner_alpha))
	print('Dijkstra Candidates: ',dijkstra_alpha, len(dijkstra_alpha))
	print('Shortest Paths to NM-II Candidates:', shortest_alpha, len(shortest_alpha))
	print('Steiner and Shortest Paths Candidates:', steiner_shortest_alpha, len(steiner_shortest_alpha))
	print('Dijkstra and Shortest Paths : ',dijkstra_shortest_alpha)
	print('Steiner and Dijkstra Candidates: ', steiner_dijkstra_alpha, len(steiner_dijkstra_alpha))

	print('Candidates present in all outputs:', all_set_alpha, len(all_set_alpha))

	print('Union of all candidates: ', union_alpha, len(union_alpha))


	return


## Input file: 'nodes-flybase.txt'
## Outputs a dictionary with flybase IDs as keys and common names as values
def read_common_names(filename):
    name_dict = {}
    with open (filename, 'r') as f:
        for line in f:
            k = line.strip().split()
            if k[2] is not k[0]:
                name_dict[k[0]] = k[2]
    print('Name dictionary:' + str(name_dict))
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
			print(k)
			if k[1] == 'Y':
				terminal_nodes.add(k[0])
			else:
				non_terminal_nodes.add(k[0])
	#print('Terminal Nodes:', terminal_nodes)
	#print('Non-terminal Nodes:', non_terminal_nodes)
	return terminal_nodes, non_terminal_nodes


## Input file: 'Dijkstra_rank.txt'
## Outputs set of nodes in Dijkstra
def read_dijkstra_rank(filename):
	Dijkstra_nodes = set()
	with open (filename, 'r') as f:
		s = f.readline() #takes away header 
		for line in f:
			k = line.strip().split()
			if float(k[1]) > 0.7:
				Dijkstra_nodes.add(k[0])
	#print('Dijkstra Nodes:', Dijkstra_nodes)
	return Dijkstra_nodes


## Input file: 'new_shortest_paths.txt'
## Outputs set of nodes that are keys in new_shortest_paths
def read_new_shortest_paths(filename):
	shortest_paths_nodes = set()
	with open (filename, 'r') as f:
		s = f.readline() #takes away header 
		for line in f:
			k = line.strip().split()
			shortest_paths_nodes.add(k[0])
	#print('Shortest Paths Nodes:', shortest_paths_nodes)
	return shortest_paths_nodes

## converts FlyBase identifiers to common names
## leaces FlyBase id if there is no common name
def FB_to_common(node_set, common_names_dict):
	new_set = set()
	for node in node_set:
		if node in common_names_dict:
			new_set.add(common_names_dict[node])
		else:
			new_set.add(node)
	print(new_set)
	return new_set


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





if __name__ == '__main__':
	main()