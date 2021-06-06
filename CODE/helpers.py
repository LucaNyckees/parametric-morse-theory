import networkx as nx
import numpy as np
from numpy import linalg as LA
import scipy
from simplicial import *
import gudhi as gd
import json
import sys
from scipy.stats import bernoulli

def stochastic_2block_model(n,p,q):
    
    if (n % 2) != 0:
        print("Please take an even integer n.")
        return 0
    
    G = nx.Graph()
    G.add_nodes_from(range(n))
    for i in range(int(n/2)):
        for j in range(int(n/2)):
            if i!=j and bernoulli.rvs(p):
                G.add_edge(i,j)
    for i in range(int(n/2),n):
        for j in range(int(n/2),n):
            if i!=j and bernoulli.rvs(p):
                G.add_edge(i,j)
    for i in range(int(n/2)):
        for j in range(int(n/2),n):
            if i!=j and bernoulli.rvs(q):
                G.add_edge(i,j)
    return G

def stochastic_4block_model(n,p,q):
    
    if (n % 4) != 0:
        print("Please take n as a multiple of 4.")
        return 0
    
    G = nx.Graph()
    G.add_nodes_from(range(n))
    for i in range(int(n/4)):
        for j in range(int(n/4)):
            if i<j and bernoulli.rvs(p):
                G.add_edge(i,j)
    for i in range(int(n/4),int(n/2)):
        for j in range(int(n/4),int(n/2)):
            if i<j and bernoulli.rvs(p):
                G.add_edge(i,j)
    for i in range(int(n/2),int(3*n/4)):
        for j in range(int(n/2),int(3*n/4)):
            if i<j and bernoulli.rvs(p):
                G.add_edge(i,j)
    for i in range(int(3*n/4),int(n)):
        for j in range(int(3*n/4),int(n)):
            if i<j and bernoulli.rvs(p):
                G.add_edge(i,j)
    for i in range(n):
        for j in range(n):
            if i<j and bernoulli.rvs(q):
                G.add_edge(i,j)
    return G



def build_series(G, functions, start, time_step, time_count):
    
    """
    Args:
        G : networkX.Graph(), a graph
        functions : list of functions (one function for each node of G)
        T : list of time slices (for which we sample the functions)
        
    Returns:
        This function returns a sequence of node labelings on K, where K is 
        the clique complex built on the graph G (there is a node labeling 
        for each time slice in T), put in the form of an array np.ndarray() 
        of size |G| x |T|.
    """
    T = []
    for i in range(int(time_count)):
        T.append(start)
        start += time_step
    
    n = G.number_of_nodes()
    time_count = len(T)
    
    assert n == len(functions), "There has to be a function for each node."
    
    K = gd.SimplexTree()
    for e in G.edges:
        K.insert(e)
    for n in G.nodes:
        K.insert([n])
    K.expansion(3)
    
    complex_series = []
    complex_series += time_count * [K]
    
    for j in range(time_count):
        
        nodes = list(complex_series[j].get_skeleton(0))
        
        for i in range(n):
            
            complex_series[j].assign_filtration(nodes[i][0],functions[i](T[j])) 

    
    return complex_series


def get_skel(K,p):
    
    skel = list(K.get_skeleton(p))
    skel_new = []
    for simplex in skel:
        if len(simplex[0]) == p+1:
            skel_new.append(simplex)
    return skel_new
    


def build_morse_function(K,d,g,noise):
    
    """
    K : simplicial complex (gudhi.SimplexTree)
    d : dimension of K 
    g : labels on nodes given by a dictionary
    """
    
    
    Flag = {}
    f = {} 
    
    for p in range(d):
        for simplex in list(K.get_skeleton(p)):
            Flag[str(simplex[0])]=0


    for simplex in list(K.get_skeleton(0)):
        f[str(simplex[0])] = g[str(simplex[0])]
        
        
        
    for p in range(1,d+1):
        
        for simplex in get_skel(K,p):
            Faces = [face[0] for face in list(K.get_boundaries(simplex[0]))]
            Faces = [str(face) for face in Faces]
            Faces = sorted(Faces, key=lambda x: f[x], reverse=True) 
            
            gamma_0 = Faces[0]
            gamma_1 = Faces[1]
            
            
            if (Flag[Faces[0]]==0 and f[gamma_0]>f[gamma_1]):
                f[str(simplex[0])]=(f[gamma_0]+f[gamma_1]) / 2
                
                Flag[gamma_0]=1
            
            else:
                epsilon = np.random.uniform(low=0, high=noise)
                f[str(simplex[0])] = f[gamma_0] + epsilon
    
    return [f, Flag]


def magic_function(G, functions, start, time_step, time_count, noise): 
    
    """
    Returns a list of pairs, of the form [... [S,f_i] ...] where f_i is a 
    discrete Morse function that is built on the simplicial complex S.
    """
    
    series = build_series(G, functions, start, time_step, time_count)
    
    list_g = []
    list_f = []
    
    for i in range(len(series)):
        
        g = {}
        
        for j in range(G.number_of_nodes()):
            
            simplex = get_skel(series[i],0)[j]
            g[str(simplex[0])] = simplex[1] 
        
        list_g.append(g)
    
    for i in range(len(series)):
        
        f = build_morse_function(series[i], series[i].dimension(), list_g[i],noise)[0]
    
        list_f.append(f)
    
    list_pairs = [[series[i],list_f[i]] for i in range(len(series))]
    
    #print("We have a sequence of {} complexes.\n".format(len(series)))
        
    return list_pairs

