import networkx as nx
import numpy as np
from numpy import linalg as LA
import scipy
from simplicial import *
import gudhi as gd
import json
import sys


def build_series(G):
    
    """
    This function returns a series of simplicial complexes that are built as clique 
    complexes on a series of labeled graphs (labels based on spectral theory).
    
    Arguments:
        G : unlabeled, undirected graph
    """
    
    L_sparse = nx.linalg.laplacianmatrix.laplacian_matrix(G)
    L = scipy.sparse.csr_matrix.toarray(L_sparse)
    
    # use scipy to compute first k eigen pairs (in case graph too big)
    # think about ordering the eigen pairs (order on eigenvalues)
    # add function : distr of eigenvalues (at what point are they insignificant?)
    
    w, v = LA.eig(L)
    
    graph_series = []
    complex_series = []
    
    w_res = []
    indices = []
    count = -1
    for i in w:
        count += 1
        if i not in w_res:
            w_res.append(i)
            indices.append(count)
            
    # taking out the multiplicities
    
    v_res = [v[i] for i in indices] 
    
    # sorting the eigenvalues in increasing order
    
    pairs = [[w_res[i],v_res[i]] for i in range(len(v_res))]
    
    pairs = sorted(pairs,key=lambda pairs: pairs[0])
    
    w_res = [pairs[i][0] for i in range(len(pairs))]
    v_res = [pairs[i][1] for i in range(len(pairs))]

        
    for i in indices:
        st = gd.SimplexTree()
        for e in G.edges:
            st.insert(e)
        for n in G.nodes:
            st.insert([n])
        st.expansion(3)
        complex_series.append(st)
    
    for i in range(len(indices)):
        nodes = list(complex_series[i].get_skeleton(0))

        for j in range(len(nodes)):
            complex_series[i].assign_filtration(nodes[j][0],v_res[i][j])
            
            
    
    return complex_series

def get_skel(K,p):
    
    skel = list(K.get_skeleton(p))
    skel_new = []
    for simplex in skel:
        if len(simplex[0]) == p+1:
            skel_new.append(simplex)
    return skel_new
    


def build_morse_function(K,d,g):
    
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
                epsilon = np.random.uniform(low=0.0, high=0.5)
                f[str(simplex[0])] = f[gamma_0] + epsilon
    
    return [f, Flag]


def magic_function(G): 
    
    """
    Returns a list of pairs, of the form [... [S,f_i] ...] where f_i is a 
    discrete Morse function that is built on the simplicial complex S.
    """
    
    series = build_series(G)
    
    list_g = []
    list_f = []
    
    for i in range(len(series)):
        
        g = {}
        
        for j in range(G.number_of_nodes()):
            
            simplex = get_skel(series[i],0)[j]
            g[str(simplex[0])] = simplex[1] 
        
        list_g.append(g)
    
    for i in range(len(series)):
        
        f = build_morse_function(series[i], series[i].dimension(), list_g[i])[0]
    
        list_f.append(f)
    
    list_pairs = [[series[i],list_f[i]] for i in range(len(series))]
    
    print("We have a sequence of {} complexes.\n".format(len(series)))
        
    return list_pairs
