'''
Import both the imputgmt and STRING files
Generate a dict of dicts wherein each key is a gene and each value is a dict with connected genes and their scores
Create a dict of the loci and each gene at each loci

Overall must do a t-test (or a stat test) on the number of edges (edge density = number of edges/number of nodes)
between the control (random generate graphs which random loci) and the subnetworks that were generated with the FA loci
'''

####Gets the dict of loci and all the genes within each
f = open("C:\\Users\\ascol\\Downloads\\7711_Day3_Assignment2\\Assignment2\\Input.gmt.txt", "r")
loci_gene_dict = {}
loci = 0
for line in f:
    line = line.strip("\n") #takes out the \n at the end of each line
    line_list = line.split("\t") #creates list divided by tabs
    loci_gene_dict[loci] = line_list[2:] #removes the first two entries which are the FA locus number and the tagSNP
    loci += 1 # keeps track of the loci number used for the keys
f.close()

####I need to create a dictionary of dicts wherein each entry is a gene and the genes that
f = open("C:\\Users\\ascol\\Downloads\\7711_Day3_Assignment2\\Assignment2\\STRING.txt", "r")
gene_interaction_dict = {}
connected_gene_weight = {}
for line in f:
    line = line.strip("\n") #takes the \n off of each line
    line = line.split("\t") #takes the string and converts it into a list of three elements
    if gene_interaction_dict.get(line[0]) is None: #checks to make sure that the key is exists
        connected_gene_weight = {} #if it doesnt this generates a new empty connected_gene weight dictionary
    connected_gene_weight[line[1]] = line[2] #adds to the connected_gene_weight dictionary of the current key
    gene_interaction_dict[line[0]] = connected_gene_weight #updates the gene_interaction_dict which is a dictionary of nodes (keys) with values of a dictionaary contraining the connected nodes and the weight of the connection
f.close()

####need to create a dctionary of nodes and their degrees (number of connections)
num_connections_dict = {}
for key in gene_interaction_dict:
    num_connections_dict[key] = len(gene_interaction_dict[key]) #for each gene (key) counts how many connections exist (value)

####need to bin by similar key values


'''print(len(gene_interaction_dict))
first_value = list(gene_interaction_dict.values())[0]
second_value = list(gene_interaction_dict.values())[1]
print('First Value: ', first_value, "\n\n", "Second value: ", second_value)'''






