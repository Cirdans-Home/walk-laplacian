# This function generates a tree with a given number of nodes and a given branching factor using
# the NetworkX library. The script has to be executed through the MATLAB engine API in order to
# be able to call the function from MATLAB and return the adjacency matrix of the generated tree
# as a MATLAB sparse matrix.

import networkx as nx
import numpy as np

def generate_tree(n_nodes, seed):
    G = nx.generators.random_unlabeled_tree(n_nodes,seed=seed)
    # Convert the graph to a sparse adjacency matrix
    A = nx.adjacency_matrix(G).toarray()
    pos = nx.spring_layout(G)
    return A, pos

A, pos = generate_tree(n_nodes, seed)