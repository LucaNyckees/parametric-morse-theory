import networkx as nx
import numpy as np
from numpy import linalg as LA
import scipy
from simplicial import *
import gudhi as gd
import json
import sys
from helpers import *

def check_Morse(K,f):
    
    """
    Returns a boolean value determining if f is a discrete Morse function on K.
    """
    
    simplices = [simplex[0] for simplex in list(K.get_simplices())]
    
    for simplex in simplices:
        faces_1 = K.get_boundaries(simplex)
        cofaces_1 = K.get_cofaces(simplex,1)
        count_f = 0
        for face in faces_1:
            if f[str(face[0])] > f[str(simplex)]:
                count_f +=1
        count_c = 0
        for coface in cofaces_1:
            if f[str(coface[0])] < f[str(simplex)]:
                count_c +=1
        if count_c > 1 or count_f > 1:
            return False
    return True

def critical_cells(K,f):
    
    """
    Returns the set of critical cells of (K,f) (in order of increasing dim).
    
    Arguments:
        K : gudhi.SimplexTree (simplicial complex)
        f : dict (dictionary with labels on simplices, f is a DMT on K)
    Return type : list
    """
    
    C = []
    
    simplices = [simplex[0] for simplex in list(K.get_simplices())]
    
    for simplex in simplices:
        faces_1 = K.get_boundaries(simplex)
        cofaces_1 = K.get_cofaces(simplex,1)
        count_f = 0
        for face in faces_1:
            if f[str(face[0])] > f[str(simplex)]:
                count_f +=1
        count_c = 0
        for coface in cofaces_1:
            if f[str(coface[0])] < f[str(simplex)]:
                count_c +=1
        if count_c == 0 and count_f == 0:
            C.append(simplex)
            
        C = sorted(C,key=len)
    
    return C


def gradient(K,f):
    
    """
    Returns the gradient vector field of f (as a list of pairs of simplices).
    """
    
    V = []
    
    simplices = [simplex[0] for simplex in list(K.get_simplices())]
    
    for simplex in simplices:
        cofaces_1 = K.get_cofaces(simplex,1)
        for coface in cofaces_1: 
            if f[str(coface[0])] < f[str(simplex)] and [simplex,coface[0]] not in V:
                V.append([simplex,coface[0]])
    
    return V


def Hasse_diagram(K,V):
    
    H = nx.DiGraph()
    
    simplices = [simplex[0] for simplex in list(K.get_simplices())]
    
    for cell in simplices:
        H.add_node(str(cell))
    
    for cell in simplices:
        cofaces_1 = [coface[0] for coface in  K.get_cofaces(cell,1)]
        H.add_edges_from([(str(coface),str(cell)) for coface in cofaces_1])
        
    for edge in V:
        H.remove_edge(str(edge[1]),str(edge[0]))
        H.add_edge(str(edge[0]),str(edge[1]))
    
    return H


def V_paths(K,V,s1,s2):
    
    H = Hasse_diagram(K, V)
    
    all_paths = list(nx.all_simple_paths(H,str(s1),str(s2)))
    
    all_paths = [[json.loads(cell) for cell in path] for path in all_paths]
    
    if s1 == s2:
        all_paths.append([s1])
    
    return all_paths
