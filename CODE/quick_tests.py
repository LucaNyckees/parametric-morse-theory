import networkx as nx
import numpy as np
from numpy import linalg as LA
import scipy
from simplicial import *
import gudhi as gd
import json
import sys
import matplotlib.pyplot as plt
import numpy as np

from helpers import *
from PMT import *
from DMT import *


G3 = nx.Graph()
for i in range(1,8):
    G3.add_node(i)
    G3.add_edge(i,i+1)
G3.add_edge(1,8)

C = []

C.append([[5],[1,2]])
C.append([[7],[4,5]])
C.append([[5],[7],[5,6],[1,2]])
C.append([[1],[5],[3,4],[7,8]])
C.append([[8],[5,6]])
C.append([[4],[7,8]])
C.append([[1],[6],[2,3],[7,8]])
C.append([[5],[4,5]])

V = []

V.append([[[1],[1,8]],[[8],[7,8]],[[7],[6,7]],[[6],[5,6]],[[4],[4,5]],[[3],[3,4]],[[2],[2,3]]])
V.append([[[1],[1,8]],[[8],[7,8]],[[6],[6,7]],[[5],[5,6]],[[4],[3,4]],[[3],[2,3]],[[2],[1,2]]])
V.append([[[1],[1,8]],[[8],[7,8]],[[6],[6,7]],[[4],[4,5]],[[3],[3,4]],[[2],[2,3]]])
V.append([[[8],[1,8]],[[7],[6,7]],[[6],[5,6]],[[4],[4,5]],[[3],[2,3]],[[2],[1,2]]])
V.append([[[1],[1,8]],[[7],[7,8]],[[6],[6,7]],[[5],[4,5]],[[4],[3,4]],[[3],[2,3]],[[2],[1,2]]])
V.append([[[8],[1,8]],[[7],[6,7]],[[6],[5,6]],[[5],[4,5]],[[3],[3,4]],[[2],[2,3]],[[1],[1,2]]])
V.append([[[8],[1,8]],[[7],[6,7]],[[5],[5,6]],[[4],[4,5]],[[3],[3,4]],[[2],[1,2]]])
V.append([[[1],[1,8]],[[8],[7,8]],[[7],[6,7]],[[6],[5,6]],[[4],[3,4]],[[3],[2,3]],[[2],[1,2]]])


"""BD_coordinates = life_coordinates(G3,C,V)

print(BD_coordinates)

print("----------------------------------------------------\n")

plot_persistence_diagram(BD_coordinates)

#plot_BD_barcode(BD_coordinates)
"""

"""
K = gd.SimplexTree()
for e in G3.edges:
    K.insert(e)
for n in G3.nodes:
    K.insert([n])
all_connected_pairs = connecting_critical_cells(G3,C,V)
coordinates = life_coordinates_of_cell(K,[5],C,V,1,all_connected_pairs)

print(coordinates)

print(are_connected(K,V[6],V[7],[6],[5]))
print(are_connected(K,V[7],V[6],[5],[6]))
"""

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

G2 = stochastic_2block_model(10,0.7,0.05)


#G3 = nx.erdos_renyi_graph(10,0.5)

BD_coordinates = life_coordinates(G2)
plot_PD(BD_coordinates)





    







