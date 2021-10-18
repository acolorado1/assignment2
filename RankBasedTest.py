'''
RankBasedTest uses multiple functions wherein a STRING.txt database file and .gmt formated file are input and two lists
of edge node densities are generated. A p-value is calculated using a Wilcoxon rank sum test to determine if the two
lists are significantly different.
'''

import random
import scipy.stats as stats
import argparse as arg

'''
Added four arguments to parser: 
1. --GI parameter that takes gene interaction file path as a string 
2. --ef parameter that takes loci and genes at loci file path as a string
3. --nb parameter that takes the number of bins as an integer
4. --nt parameter that takes the number of trials as an integer
'''

parser = arg.ArgumentParser(description="Create two networks and determines whether their edge density is significantly "
                                        "different")
parser.add_argument("--GI", "-gene_interaction_file", type=str, help="gene interaction file path (default STRING.txt)",
                    default="STRING.txt")
parser.add_argument("--ef", "-experimental_file", type=str, help="experimental loci file path (default Input.gmt.txt)",
                    default="Input.gmt.txt")
parser.add_argument("--nb", "-n_bins", type=int, help="number of bins (default: 128)", default=128)
parser.add_argument("--nt", "-n_trials", type=int, help="number of trials (default: 1000)", default=1000)

args = parser.parse_args()

'''
Parameter f_experimental takes .gmt formatted file containing loci and genes at that loci.
Function opens file and reads each line creating a dictionary where keys are the loci number and the values are a list
of genes at that loci. 
Returns loci_gene_dict dictionary. 
'''
def dict_loci(f_experimental):
    f = open(f_experimental,"r")
    loci_gene_dict = {}
    loci = 0
    for line in f:
        line = line.strip("\n")
        line_list = line.split("\t")
        loci_gene_dict[loci] = line_list[2:]
        loci += 1
    f.close()
    return loci_gene_dict

'''
Parameter f_cofunction_net takes STRING.txt database file containing gene interactions. 
dict_genes_interactions opens files and reads each line creating a dictionary of dictionaries where the keys are genes 
and the values are dictionaries containing gene names (keys) and the strength of their interaction (values). 
Returns this gene_interactions_dict dictionary.
'''
def dict_gene_interactions (f_cofunction_net):
    f = open(f_cofunction_net,"r")
    gene_interaction_dict = {}
    connected_gene_weight = {}
    for line in f:
        line = line.strip("\n")
        line = line.split("\t")
        if gene_interaction_dict.get(line[0]) is None:
            connected_gene_weight = {}
        connected_gene_weight[line[1]] = line[2]
        gene_interaction_dict[line[0]] = connected_gene_weight
    f.close()
    return gene_interaction_dict

'''
Parameter f_experiemental takes .gmt formatted file as a string 
Parameter GI_file takes STRING.txt database file as a string 
Parameter n_bins takes an integer representing the number of bins 
Parameter n_trials takes an integer representing the number of trials that will be run 
listofedgecounts calls functions dict_loci to create a dictionary of loci (keys) and genes at the loci (values) and 
dict_gene_interactions to create a dictionary of gene interactions. Then creates num_connections_dict dictionary where
keys are gene names and values are the number of genes they interact with. num_connections_dict is sorted from lowest 
to highest value and converted into a list of tuples with gene and number of interaction as key value pairs. Binning is 
done by size wherein the total number is divided by n_bins and the remainder is calculated. gene_bin_dict is initialized
to create a dictionary of genes (keys) and their corresponding bin number (value). bin_genesinbin_dict is initialized to 
create a dictionary of bins (keys) and a list of genes in each bin (values). 
edge_counts_experimental list is initialized where each item is the edge density of the network calculated for each n_trial
based a randomly chosen gene for each loci in the loci_genes_dict. edge_counts_CFN list is initialized where each item 
is the edge density of the full network of genes calculated for each n_trial based on a randomly chosen gene in the same 
bin as gene chosen at each loci. 
Returns a list of the two lists (edge_counts_experimental and edge_counts_CFN)
'''
def listofedgecounts (f_experimental,GI_file,n_bins,n_trials):
    loci_gene_dict = dict_loci(f_experimental)
    gene_interaction_dict = dict_gene_interactions(GI_file)
    num_connections_dict = {}
    for key in gene_interaction_dict:
        num_connections_dict[key] = len(gene_interaction_dict[key])
    sorted_num_connections_list = sorted(num_connections_dict.items(), key=lambda x: x[1])
    gene_bin_dict = {}
    bin_genesinbin_dict = {}
    for bin in range(0,n_bins):
        list_of_genes = []
        if bin != n_bins:
            new_list_gene_degree_pairs = sorted_num_connections_list[bin*(len(sorted_num_connections_list)//(n_bins-1)):
                                                        (bin+1)*(len(sorted_num_connections_list)//(n_bins-1))]
        else:
            new_list_gene_degree_pairs = sorted_num_connections_list[-(len(sorted_num_connections_list)%(n_bins-1)):]
        for pair in new_list_gene_degree_pairs:
            gene_bin_dict[pair[0]] = bin
            list_of_genes.append(pair[0])
        bin_genesinbin_dict[bin] = list_of_genes

    edge_counts_experimental = []
    edge_counts_CFN = []
    for num in range(0,n_trials):
        sum_counts = 0
        CFN_sum_counts = 0
        for key in loci_gene_dict:
            rand_gene = random.choice(loci_gene_dict[key])
            if rand_gene in num_connections_dict:
                CFN_bin = gene_bin_dict[rand_gene]
                CFN_rand_gene = random.choice(bin_genesinbin_dict[int(CFN_bin)])
                CFN_sum_counts += num_connections_dict[CFN_rand_gene]
                sum_counts += num_connections_dict[rand_gene]
        edge_counts_CFN.append(CFN_sum_counts)
        edge_counts_experimental.append(sum_counts)
    return [edge_counts_CFN,edge_counts_experimental]

'''
Parameter experimental_and_CFN_edgecount_list is a list of two lists containing edge counts 
Stubbed function returns the input 
'''
def stub_GA(experimental_and_CFN_edgecount_list):
    return experimental_and_CFN_edgecount_list

'''
Parameter f_experiemental takes .gmt formatted file as a string 
Parameter GI_file takes STRING.txt database file as a string 
Parameter n_bins takes an integer representing the number of bins 
Parameter n_trials takes an integer representing the number of trials that will be run 
rank_based_test calls of listofedgecounts function to get a list of two lists containing the calculated edgecounts, calls 
the stubbed genetic algorithm on the lists of edgecounts and run a wilcoxon ranked sum test on the two lists 
Returns nothing, but prints the list of edge counts from the subnetwork of the experimental loci, the list of edge
counts from the randomly generated network, and the p-value that was calculated.
'''
def rank_based_test(f_experimental,GI_file,n_bins,n_trials):
    list_of_ec = listofedgecounts(f_experimental, GI_file, n_bins, n_trials)
    GA_edgecounts = stub_GA(list_of_ec)
    pval = stats.wilcoxon(GA_edgecounts[0],GA_edgecounts[1])[1]
    print("CFN edge counts: ", GA_edgecounts[0])
    print("FA Network edge counts: ", GA_edgecounts[1])
    print("Rank-based test p-value: ", pval)

rank_based_test(args.ef, args.GI, args.nb, args.nt)