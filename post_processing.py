'''
Use of these functions assumes that the following files 
have been generated and are located in the same directory:

'tree_edges.txt'
'tree_nodes.txt'
'new_shortest_paths.txt'
'BFS_rank.txt'

This code also requires:
'nodes-flybase.txt' for conversion from FlyBase IDs to common names.
'
'''




def main():
	common_names_dict = read_common_names('nodes-flybase.txt')
	steiner_terminal_nodes,steiner_non_terminal_nodes = read_tree_nodes('tree_nodes.txt')
	#read_BFS_rank('BFS_rank.txt')
	shortest_paths = read_new_shortest_paths('new_shortest_paths.txt')
	steiner_nt_common = FB_to_common(steiner_non_terminal_nodes,common_names_dict)
	new_shortest_common = FB_to_common(shortest_paths, common_names_dict)

	both_set = in_both_sets(steiner_nt_common,new_shortest_common)

	print('Both:', both_set)
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
	print('Terminal Nodes:', terminal_nodes)
	print('Non-terminal Nodes:', non_terminal_nodes)
	return terminal_nodes, non_terminal_nodes



def read_BFS_rank(filename):
	BFS_nodes = set()
	with open (filename, 'r') as f:
		s = f.readline() #takes away header 
		for line in f:
			k = line.strip().split()
			print(k)
			BFS_nodes.add(k[0])
	print('BFS Nodes:', BFS_nodes)
	return BFS_nodes

def read_new_shortest_paths(filename):
	shortest_paths_nodes = set()
	with open (filename, 'r') as f:
		s = f.readline() #takes away header 
		for line in f:
			k = line.strip().split()
			shortest_paths_nodes.add(k[0])
	print('Shortest Paths Nodes:', shortest_paths_nodes)
	return shortest_paths_nodes

## converst 
def FB_to_common(node_set, common_names_dict):
	new_set = set()
	for node in node_set:
		if node in common_names_dict:
			new_set.add(common_names_dict[node])
		else:
			new_set.add(node)
	print(new_set)
	return new_set

def in_both_sets(set1, set2):
	both_set = set()
	if len(set1)>len(set2):
		for node in set1:
			if node in set2:
				both_set.add(node)
	else:
		for node in set2:
			if node in set1:
				both_set.add(node)
	return both_set




if __name__ == '__main__':
	main()