# BIO331F17
Group project to identify potential regulators of NMII and Fog pathway members in a fly interactome.  

## Authors

### 2020 Update
- **Anna Ritz**, professor

### 2017 Final Project
- **Miriam Bern**, student 
- **Wyatt Gormley**, student
- **Elaine Kushkowski**, student
- ** Kathy Thompson**, student
- **Logan Tibbetts**, student

## Instructions 

The `group_proj_331.py` script will generate a list of potential regulators from a protein-protein interactome and a file of positive regulator containing interactions between nodes using pre-processing techniques, Steiner Tree Approximations, Dijkstra-ranking, and a shortest paths algorithm. 
- *Input:* a plain text file, formatted into at least 2 columns to indicate node to node interaction, a text file of regulators.
- *Output:* text file of the edge list of the Steiner Tree and a file containing all nodes within it, text files of potential regulators ranked by the Dijkstra-ranking, and text files of potential regulators calculated by the shortest-paths algorithm.
- *Output note:* "tree_nodes.txt" and "tree_edges.txt" may be different with each run because we are using a Steiner tree approximation.
- *Runtime note:* This code takes at least 3 hours to run on the original interactome.


Make sure that all files to be analyzed are in the same file folder as `group_proj_331.py`. When running, change the interactome variable in the `main()` function to be whatever text file interactome you are analyzing, and change the positives variable in the `main()` function to be a text file of positive regulators for that genome.

The purpose of `post_processing.py` is to produce output files from input obtained from the main file. It was also used for testing output functions separately from the main code. 
