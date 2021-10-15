'''
Import both the imputgmt and STRING files
Generate a dict of dicts wherein each key is a gene and each value is a dict with connected genes and their scores
Create a dict of the loci and each gene at each loci

Overall must do a t-test (or a stat test) on the number of edges (edge density = number of edges/number of nodes)
between the control (random generate graphs which random loci) and the subnetworks that were generated with the FA loci
'''
import random
import scipy.stats as stats
import argparse as arg

parser = arg.ArgumentParser(description="Create two networks and determine whether their edge density is significantly different") #starting to make arguments
parser.add_argument("--GI", "-gene_interaction_file", type=str, help="gene interaction file path (default STRING.txt)", default="STRING.txt")
parser.add_argument("--ef", "-experimental_file", type=str, help="experimental loci file path (default Input.gmt.txt)", default="Input.gmt.txt")
parser.add_argument("--nb", "-n_bins", type=int, help="number of bins (default: 128)", default=128)
parser.add_argument("--nt", "-n_trials", type=int, help="number of trials (default: 1000)", default=1000)

args = parser.parse_args()
print(args)

####Gets the dict of loci and all the genes within each
def dict_loci(f_experimental): #input .gmt file format path
    f = open(f_experimental,"r")
    loci_gene_dict = {} #empty dict of loci (key) and list of genes (value) is initialized
    loci = 0 #number of current loci initialized
    for line in f: #for each line in the file
        line = line.strip("\n") #takes out the \n at the end of each line
        line_list = line.split("\t") #creates list divided by tabs
        loci_gene_dict[loci] = line_list[2:] #removes the first two entries which are the FA locus number and the tagSNP
        loci += 1 #keeps track of the loci number used for the keys
    f.close()
    return loci_gene_dict #returns a dict of loci (key) and list of genes at loci (values)

####Creates a dictionary of dicts wherein each key is a gene and value is dict of the genes (keys) that connect to it and the weight of their functional similarity (value)
def dict_gene_interactions (f_cofunction_net): #input STRING.txt database file path
    f = open(f_cofunction_net,"r")
    gene_interaction_dict = {} #empty dict of dicts is intializied
    connected_gene_weight = {} #empty dict of connected gene (key) and weight (value) of interaction is initialized
    for line in f: #for every line in the file
        line = line.strip("\n") #takes the \n off of each line
        line = line.split("\t") #takes the string and converts it into a list of three elements
        if gene_interaction_dict.get(line[0]) is None: #checks to make sure that the key exists
            connected_gene_weight = {} #if it doesnt this generates a new empty connected_gene weight dictionary
        connected_gene_weight[line[1]] = line[2] #adds to the connected_gene_weight dictionary of the current key
        gene_interaction_dict[line[0]] = connected_gene_weight #updates the gene_interaction_dict which is a dictionary of nodes (keys) with values of a dictionaary contraining the connected nodes and the weight of the connection
    f.close()
    return gene_interaction_dict #returns a dict of gene (key) with values of a dict containing gene (key) and edge weight (value)

####Creats a dictionary of nodes and their degrees (number of connections), creates a sorted list of them by value from low to high and then creates a dictionary of bins and the genes in each bin
####Generates a graph using the 12 loci and measure the edge density (in this case the number of edges) and a graph of 12 nodes from the cofunction netwoek wherein the nodes have a similar degree
def listofedgecounts (f_experimental,GI_file,n_bins,n_trials): #input .gmt file path, STRING.txt database file path, number of bins, number of trials
    loci_gene_dict = dict_loci(f_experimental)
    gene_interaction_dict = dict_gene_interactions(GI_file)
    num_connections_dict = {} #dict of gene (key) and the number of nodes it connects to (value) is initialized
    for key in gene_interaction_dict:
        num_connections_dict[key] = len(gene_interaction_dict[key]) #for each gene (key) counts how many connections exist (value)
    sorted_num_connections_list = sorted(num_connections_dict.items(), key=lambda x: x[1]) #generates a sorted list of gene,degree pairs
    gene_bin_dict = {} #dict of genes (key) and bin number (value) initialized
    bin_genesinbin_dict = {} #dict of bin number (key) and list of genes in bin (value) initialized
    for bin in range(0,n_bins): #creating 128 bins as this was the number of bins in the paper
        list_of_genes = [] #empty list of genes in a bin
        if bin != n_bins: #for every bin except for the last one
            new_list_gene_degree_pairs = sorted_num_connections_list[bin*(len(sorted_num_connections_list)//(n_bins-1)):(bin+1)*(len(sorted_num_connections_list)//(n_bins-1))] #creates respective bin list of gene, degree pairs, since we want 128 bins we divide by 127 and the last bin will be the remainder
        else:
            new_list_gene_degree_pairs = sorted_num_connections_list[-(len(sorted_num_connections_list)%(n_bins-1)):] #creates a list of gene, degree pairs for the final bin which is the remainder of the total genes divided by 127 (our chosen number of bins)
        for pair in new_list_gene_degree_pairs:
            gene_bin_dict[pair[0]] = bin #adds a gene (key) and bin number (value) to dict
            list_of_genes.append(pair[0]) #adds the gene to a list of genes in one bin
        bin_genesinbin_dict[bin] = list_of_genes #creates bin (key) and list of genes at that bin (value) and adds to dictionary

    edge_counts_experimental = [] #empty list of edge densities from n trials taken from a random set of genes from the FA loci is intialized
    edge_counts_CFN = [] #empty list of edge densities from n trials taken from random CFN genes is intialized
    for num in range(0,n_trials): # n trials
        sum_counts = 0  #initialized sum of the number of the connections that each node has every trial
        CFN_sum_counts = 0 #initialized sum of the number of the connections that each node has every trial from the CFN
        for key in loci_gene_dict: #iterated over ever loci in FA file (0-11)
            rand_gene = random.choice(loci_gene_dict[key]) #selects a random gene is slected from the current loci
            if rand_gene in num_connections_dict: #checks to make sure that the gene is in the dict of connections taken from the STRING database
                CFN_bin = gene_bin_dict[rand_gene] #gets bin number from the randomly chosen gene in the experimental loci
                CFN_rand_gene = random.choice(bin_genesinbin_dict[int(CFN_bin)]) #a gene is randomly chosen from the same bin that the rand_gene was in
                CFN_sum_counts += num_connections_dict[CFN_rand_gene] #takes the number of connected genes and adds it to CFN_sum_counts in order to get total edge count
                sum_counts += num_connections_dict[rand_gene] #takes the number of connected genes and adds it to sum_counts in order to get the total number of edges in a network of 12 genes 1 from each of the 12 loci
        edge_counts_CFN.append(CFN_sum_counts) #after sum is calculated from the 12 randomly chosen genes and appended to the CFN_sum_counts list
        edge_counts_experimental.append(sum_counts) #after the sum is calculated it is appended to the edge_count list
    return [edge_counts_CFN,edge_counts_experimental] #returns a list of lists containing edgecounts for each trial from the experimental and Cofunction networks

####Stubbed GA, this function will eventually optimize the edge count list
def stub_GA(experimental_and_CFN_edgecount_list): #inputs list of of lists with edgecounts for each trial
    return experimental_and_CFN_edgecount_list #returns input

####Performs a statistical test on the two lists of edge counts
def rank_based_test(f_experimental,GI_file,n_bins,n_trials): #input .gmt file path, STRING.txt database file path, number of bins, number of trials
    list_of_ec = listofedgecounts(f_experimental, GI_file, n_bins, n_trials) #generates list of edge counts for both networks
    GA_edgecounts = stub_GA(list_of_ec) #stubbed genetic algorithm function
    pval = stats.wilcoxon(GA_edgecounts[0],GA_edgecounts[1])[1] #gets the p-value using the wilcoxon rank sum tset that returns a tuple where the index 1 contains the p-value
    print("CFN edge counts: ", GA_edgecounts[0])
    print("FA Network edge counts: ", GA_edgecounts[1])
    print("Rank-based test p-value: ", pval)

#rank_based_test(args.ef, args.GI, args.nb, args.nt)


rank_based_test("C:\\Users\\ascol\\Downloads\\7711_Day3_Assignment2\\Assignment2\\Input.gmt.txt",
                "C:\\Users\\ascol\\Downloads\\7711_Day3_Assignment2\\Assignment2\\STRING.txt" ,128,1000 )




