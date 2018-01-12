import glob

#Note: DroID-flybase.txt was taken from https://github.com/annaritz/fly-interactome/tree/master/databases/DroID

def main():
    #edit_file function takes the DroID portion of FlyBase interactome and removes the extra info
    #the output file just contains the DroID edges/interactions
    edit_file("DroID-flybase.txt", "DroID-flybase-ppis.txt")

    #edge_list function creates a list of edges from the file containing DroID FlyBase interactions
    edges = edge_list("DroID-flybase-ppis.txt")

    #node_list function creates a file of all the DroID FlyBase nodes 
    nodes = node_list(edges, "DroID-flybase-nodes.txt")

    #DroID_files is a list of all the DroID database files containing all weighted interactions between
    #DroID FlyBase nodes and any other protein in the database (therefore they contain some interactions
    #that are not in the FlyBase interactome)
    DroID_files = ["DroID_Database_PPIs/DroID1.txt",
    "DroID_Database_PPIs/DroID2.txt",
    "DroID_Database_PPIs/DroID3.txt",
    "DroID_Database_PPIs/DroID4.txt",
    "DroID_Database_PPIs/DroID5.txt",
    "DroID_Database_PPIs/DroID6.txt",
    "DroID_Database_PPIs/DroID7.txt",
    "DroID_Database_PPIs/DroID8.txt",
    "DroID_Database_PPIs/DroID9.txt",
    "DroID_Database_PPIs/DroID10.txt",
    "DroID_Database_PPIs/DroID11.txt",
    "DroID_Database_PPIs/DroID12.txt",
    "DroID_Database_PPIs/DroID13.txt",
    "DroID_Database_PPIs/DroID14.txt",
    "DroID_Database_PPIs/DroID15.txt",
    "DroID_Database_PPIs/DroID16.txt",
    "DroID_Database_PPIs/DroID17.txt",
    "DroID_Database_PPIs/DroID18.txt",
    "DroID_Database_PPIs/DroID19.txt",
    "DroID_Database_PPIs/DroID20.txt",
    "DroID_Database_PPIs/DroID21.txt",
    ]

    #This for loop goes through the database files and cleans them up by removing extra info
    for i in range(len(DroID_files)):
        infile = "DroID_Database_PPIs/DroID" + str(i+1) + ".txt"
        outfile = "DroID-weighted-interactions" + str(i+1) + ".txt"
        edit_database_interactions(infile,outfile)

    #The following line grabs all the edited (i.e. extra info removed) DroID database files
    interaction_files = glob.glob("DroID-weighted/*") 

    #This line combines the DroID database files into one file
    all_DroID = combine_files(interaction_files, "all-DroID-weighted.txt")

    #This line takes the file of all DroID database files and converts it into a dictionary
    #Key = edge, value = weight of the edge
    Dro_dict = convert_to_dict("all-DroID-weighted.txt")

    #This line takes the DroID FlyBase interactions and checks the dictionary to see if each interaction
    #(i.e. [a,b] and [b,a] is a key in the dictionary. If it is, it takes this edge and its weight
    #and puts it into a new file. This new file contains all the weighted DroID FlyBase interactions 
    #and serves as the weighted interactome. 
    weight_interactome(Dro_dict, "DroID-flybase-ppis.txt", "flybase-weighted-interactome.txt")

    return


#Takes as input "DroID-flybase.txt" which contains all the DroID interactions in the FlyBase interactome
#Outputs a two-column file "DroID-flybase-ppis.txt" of the interactions without the extra information
def edit_file(infile, outfile):
    with open(infile, 'r') as f:
        with open(outfile, 'w') as ids:
            for line in f:
                k = line.strip().split("\t")
                if k == ['#FBID1', 'FBID2', 'PubMedIDs', 'Evidence']:
                    pass
                else:
                    k = k[0:2]
                    pr1, pr2 = k[0], k[1]
                    ids.write(pr1 + "\t" + pr2 + "\n")
    return

#Takes as input "DroID-flybase-ppis.txt" and returns a list of edges 
def edge_list(infile):
    edges = []
    with open(infile, 'r') as ppi:
        for line in ppi:
            e = line.strip().split("\t")
            edges.append(e)
    return edges

#Takes as input the list of edges and an outfile name and produces a file that contains all
#the nodes in the DroID portion of the FlyBase interactome. These nodes are used to search the DroID
#database to find all scored interaction for each node/protein.
def node_list(edges, outfile):
    nodes = []
    with open(outfile, "w") as out:
        for e in edges:
            if e[0] not in nodes:
                nodes.append(e[0])
                out.write(e[0] + "," + "\n")
            if e[1] not in nodes:
                nodes.append(e[1])
                out.write(e[1] + "," + "\n")
    #print("# Nodes: ", len(nodes))
    return nodes

#Takes as input the weighted DroID interactions downloaded from the DroID database website and 
#outputs a three-column file with the interaction and its confidence score (removes extra info)
def edit_database_interactions(infile, outfile):
    with open(infile, 'r') as f:
        with open(outfile, 'w') as ids:
            for line in f:
                k = line.strip().replace('"', '').split(",")
                if k == ['Gen1', 'Gene2', 'Symbol1', 'Symbol2', 'Table Names', 'Weighted Correlation', 'Confidence']:
                    pass
                else:
                    pr1, pr2, conf = k[0], k[1], k[6]
                    ids.write(pr1 + "\t" + pr2 + "\t" + conf + "\n")
    return

#Takes all the three-column DroID interaction files and combines them into one file
def combine_files(files, outfile):
    with open(outfile, "wb") as o:
        for f in files:
            #print("Reading file", f)
            with open(f, "rb") as wppi:
                o.write(wppi.read())
    return

#Takes as input the file of all weighted interactions found in the DroID database and
#outputs a dictionary where the keys are a tuple representing the edges and the values are the 
#confidence scores
def convert_to_dict(all_weighted):
    weighted_DroID_dict = {} #initialize empty dictionary
    with open(all_weighted, "r") as weighted:
        for line in weighted:
            l = line.strip().split("\t")
            intn = (l[0], l[1]) #creates a tuple of the interaction
            score = l[2]
            weighted_DroID_dict[intn] = score #enters it into the dictionary
    return weighted_DroID_dict


#Takes as input the dictionary of weighted interactions and the file of DroID FlyBase interactions.
#Checks if the flybase interaction is a key in the dictionary. If it is, it adds it to a new file 
#and updates the interaction with its score. 
#Outputs this file of weighted FlyBase interactions

def weight_interactome(Dro_dict, FB_interactions, outfile):
    with open(outfile, "w") as weighted:
        with open(FB_interactions, "r") as FI: #reads the DroID FlyBase interactions file
            for line in FI:
                l = line.strip().split("\t")
                prot1 = l[0]
                prot2 = l[1]
                edge1 = (prot1, prot2)
                edge2 = (prot2, prot1)
                if edge1 in Dro_dict:
                    weighted.write(prot1 + "\t" + prot2 + "\t" + str(Dro_dict[edge1]) + "\n")
                if edge2 in Dro_dict:
                    weighted.write(prot2 + "\t" + prot1 + "\t" + str(Dro_dict[edge2]) + "\n")
    return


if __name__ == "__main__":
    main()