import networkx as nx
import numpy as np
from numpy import linalg as LA
import scipy
from simplicial import *
import gudhi as gd
import json
import sys
from DMT import *
from helpers import *

def are_connected(K,V,W,s1,s2):
    
    """
    Arguments:
        K : gudhi.SimplexTree
        V : list of pairs of cells
        W : list of pairs of cells 
        s1 : cell critical for V
        s2 : cell critical for W
    Returns:
        True if s1 and s2 are connected,
        False otherwise (according to def.3.1. of PMT).
    """
                  
    if len(s1)!=len(s2) or s1 == s2:
        
        return False
    
    else:
        
        k = len(s1)-1  
        
        k_simplices = [simplex[0] for simplex in get_skel(K,k)]

        for simplex in k_simplices:
            
            #print(simplex)
                        
            for path1 in V_paths(K,V,s1,simplex):
                
                for path2 in V_paths(K,W,simplex,s2):
                    
                    if path1 == [s1] or path2 == [s2]:
                        
                        #print("Path from {} to {} via the subpaths {} and {}.".format(
                            s1,s2,path1,path2))
                        
                        
                        return True
                    
                    
                    elif len(path1[1])==k and len(path2[1])==k+2:
                    
                        #print("WE have a path from {} to {} via the subpaths {} and {}.".format(
                            s1,s2,path1,path2))
                        
                        return True
           
    return False

        
def connecting_critical_cells(G,C=[],V=[]):
    
    """
    + COMPLETE DMT ANALYSIS OF G
    + Returns the list of pairs of critical cells for f_i and f_j that are connected,
    for all indices i and j, as 4-tuples (alpha,beta,i,j) indicating the time slices.
    If no C and V are given as input, it deals with G via the spectral theory principle.
    Else, it considers a time-series with Morse theory characteristics given by C and V.
    """
    
    all_connected_pairs = []
    
    if C!=[] and V!=[]:
        
        for i in range(len(C)):
            print("At time slice {}, K has critical cells {}.\n".format(i,C[i]))
            print("At time slice {}, K has gradient vector field {}.\n".format(i,V[i]))
            print("----------------------------------------------------\n")
        
        K = gd.SimplexTree()
        for e in G.edges:
            K.insert(e)
        for n in G.nodes:
            K.insert([n])
        
        index = len(C)
        
    else:
        
        list_dmts = magic_function(G)

        C = []
        V = []

        for i in range(len(list_dmts)):
            K = list_dmts[i][0]
            f = list_dmts[i][1]
            C.append(critical_cells(K,f))
            print("At time slice {}, K has critical cells {}.\n".format(i,C[i]))
            V.append(gradient(K, f))
            print("At time slice {}, K has gradient vector field {}.\n".format(i,V[i]))
            print("----------------------------------------------------\n")
            
        index = len(list_dmts)
        
    print(K.dimension())
        
        
    for i in range(index):
        
        for j in range(index):  
        
            for k in range(K.dimension()+1):
                
                if i!=j:
                
                    k_crit_i = [cell for cell in C[i] if len(cell)==k+1]
                    k_crit_j = [cell for cell in C[j] if len(cell)==k+1]
                
                    for s1 in k_crit_i:
                    
                        for s2 in k_crit_j:
                            
                            if are_connected(K,V[i],V[j],s1,s2):
                            
                                all_connected_pairs.append([s1,s2,i,j])
                            
    print("All connected pairs of critical cells across time slices are : {}.\n".format(all_connected_pairs))
    
    return None
