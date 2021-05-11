import networkx as nx
import numpy as np
from numpy import linalg as LA
import scipy
from simplicial import *
import gudhi as gd
import matplotlib.pyplot as plt
import json
import sys
import numpy as np
from DMT import *
from helpers import *

plt.rcParams["figure.figsize"] = (10,10)

def right_dimension(path,k):
    
    for i in range(len(path)):
        
        if i % 2 == 0 and len(path[i]) != k+1:
            
            return False
        
        if i % 2 == 1 and len(path[i]) != k:
            
            return False
        
    return True
     
        

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
                  
    if len(s1)!=len(s2):
        
        return False
    
    if s1 == s2:
        
        return True
    
    else:
        
        k = len(s1)-1  
        
        k_simplices = [simplex[0] for simplex in get_skel(K,k)]

        for simplex in k_simplices:
            
            #print(simplex)
                        
            for path1 in V_paths(K,V,s1,simplex):
                
                for path2 in V_paths(K,W,simplex,s2):
                    
                    if path1 == [s1] or path2 == [s2]:
                        
                        #print("Path from {} to {} via the subpaths {} and {}.".format(
                            #s1,s2,path1,path2))
                        
                        
                        return True
                    
                    
                    elif right_dimension(path1, k) and right_dimension(path2,k+1):
                    
                        #print("WE have a path from {} to {} via the subpaths {} and {}.".format(
                            #s1,s2,path1,path2))
                        
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
    
    return all_connected_pairs


def life_coordinates_of_cell_at_i(K,s,C,V,i,all_connected_pairs):
    
    """
    Args
    - s is a simplex that is critical for gradient field V[i] (i.e. s in C[i])
    - all_connected_pairs is as returned by connecting_critical_cells()
    
    Determines whether s dies/is born at time i
    """    
    
    dies_at_i = True
    born_at_i = True
    
    if i != len(C)-1:
    
        for k in range(K.dimension()+1):
            
            k_crit_i_plus_1 = [cell for cell in C[i+1] if len(cell)==k+1]
            
            for s1 in k_crit_i_plus_1:
                
                if are_connected(K,V[i],V[i+1],s,s1)==True:
                    
                    if are_connected(K,V[i+1],V[i],s1,s)==True:
                        
                        dies_at_i = False
                        
                        break
                    
    elif i != len(C)-1 and s in C[i+1]:
            
        dies_at_ie = False
                    
    elif i == len(C)-1:
        
        dies_at_i = False
                        
    if i != 0:
        
        for k in range(K.dimension()+1):
         
            k_crit_i_minus_1 = [cell for cell in C[i-1] if len(cell)==k+1]
                
            for s1 in k_crit_i_minus_1:
                    
                if are_connected(K,V[i-1],V[i],s1,s)==True:
                        
                    if are_connected(K,V[i],V[i-1],s,s1)==True:
                            
                        born_at_i = False
                        
                        break
                    
    elif i != 0 and s in C[i-1]:
            
        born_at_i = False
                
    elif i==0:
        
        born_at_i = False
        
    
    #print("Death : {}\n Birth : {}\n for i = {}\n".format(dies_at_i,born_at_i,i))
    return (dies_at_i,born_at_i)



def life_coordinates_of_cell(K,s,C,V,i,all_connected_pairs):
    
    """
    Args
    - s is a simplex that is critical for gradient field V[i] (i.e. s in C[i])
    - all_connected_pairs is as returned by connecting_critical_cells()
    
    Return the time of birth and death of cell s
    """    
    
    time_slices = [i for i in range(len(C)) if s in C[i]]
    
    deaths = []
    births = []
    
    for i in time_slices:
        
        if life_coordinates_of_cell_at_i(K,s,C,V,i,all_connected_pairs)[0]:
            
            deaths.append(i)
            
        if life_coordinates_of_cell_at_i(K,s,C,V,i,all_connected_pairs)[1]:   
            
            births.append(i)
            
    if deaths == []:
        death = max(time_slices)
    else:
        death = max(deaths)
        
    if births == []:
        birth = min(time_slices)
    else:
        birth = min(births)
            
    return (birth,death)
    
    



def life_coordinates(G,C=[],V=[]):
    
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
    
    print("----------------------------------------------------\n")
    
    coordinates = []
    
    for i in range(index): 
    
        for k in range(K.dimension()+1):
        
            k_crit_i = [cell for cell in C[i] if len(cell)==k+1]
    
            for s in k_crit_i:
                
                l = life_coordinates_of_cell(K,s,C,V,i,all_connected_pairs)
                
                b = l[0]
                d = l[1]
                
                # 3-tuples (simplex, birth, death)
        
                coordinates.append([s,b,d])
    
    
    return coordinates


def plot_BD_barcode(coordinates):
    
    
    cells = [str(trio[0]) for trio in coordinates]
    births = [str(trio[1]) for trio in coordinates]
    
    plt.style.use('ggplot')
    
    plt.barh(cells,births)
    plt.title('Life coordinates of critical cells')
    plt.ylabel('Critical cells')
    plt.xlabel('Life time')
    plt.show()

def plot_PD(coordinates):

    maximum = max([trio[2] for trio in coordinates])
    
    coord_0 = [trio for trio in coordinates if len(trio[0])==1]
    coord_1 = [trio for trio in coordinates if len(trio[0])==2]
    coord_2 = [trio for trio in coordinates if len(trio[0])==3]
    coord_3 = [trio for trio in coordinates if len(trio[0])==4]
    
    births_0 = [trio[1] for trio in coord_0]
    deaths_0 = [trio[2] for trio in coord_0]
    births_1 = [trio[1] for trio in coord_1]
    deaths_1 = [trio[2] for trio in coord_1]
    births_2 = [trio[1] for trio in coord_2]
    deaths_2 = [trio[2] for trio in coord_2]
    births_3 = [trio[1] for trio in coord_3]
    deaths_3 = [trio[2] for trio in coord_3]
    

    N = len(coord_0)
    life_0 = [(trio[1],trio[2]) for trio in coord_0]
    z = range(maximum+1)

    multiplicities_0 = []

    for i in range(N):

        mult = life_0.count(life_0[i])
        multiplicities_0.append(mult)
        multiplicities_0[i]=multiplicities_0[i]**2+20
        
    colors_0 = [0.2 for i in range(N)]
        
    M = len(coord_1)
    life_1 = [(trio[1],trio[2]) for trio in coord_1]

    multiplicities_1 = []

    for i in range(M):

        mult = life_1.count(life_1[i])
        multiplicities_1.append(mult)
        multiplicities_1[i]=multiplicities_1[i]**2+20

    colors_1 = [0.2 for i in range(M)]
    
    O = len(coord_2)
    life_2 = [(trio[1],trio[2]) for trio in coord_2]

    multiplicities_2 = []

    for i in range(O):

        mult = life_2.count(life_2[i])
        multiplicities_2.append(mult)
        multiplicities_2[i]=multiplicities_2[i]**2+20

    colors_2 = [0.2 for i in range(O)]
    
    P = len(coord_3)
    life_3 = [(trio[1],trio[2]) for trio in coord_3]

    multiplicities_3 = []

    for i in range(P):

        mult = life_3.count(life_3[i])
        multiplicities_3.append(mult)
        multiplicities_3[i]=multiplicities_3[i]**2+20

    colors_3 = [0.2 for i in range(P)]

    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, sharex=True, sharey=True)

    ax1.scatter(births_0, deaths_0, s=multiplicities_0, c=colors_0, alpha=0.5)
    ax1.plot(z,z,'g')
    ax1.set_title("0-cells")
    ax2.scatter(births_1, deaths_1, s=multiplicities_1, c=colors_1, alpha=0.5)
    ax2.plot(z,z,'g')
    ax2.set_title("1-cells")
    ax3.scatter(births_2, deaths_2, s=multiplicities_2, c=colors_2, alpha=0.5)
    ax3.plot(z,z,'g')
    ax3.set_title("2-cells")
    ax4.scatter(births_3, deaths_3, s=multiplicities_3, c=colors_3, alpha=0.5)
    ax4.plot(z,z,'g')
    ax4.set_title("3-cells")
    fig.suptitle("Life coordinates of critical cells")
    
    for i in range(N):
        
        ax1.plot([births_0[i],births_0[i]],[births_0[i],deaths_0[i]],'r--')
        
    for i in range(M):
        
        ax2.plot([births_1[i],births_1[i]],[births_1[i],deaths_1[i]],'r--')
        
    for i in range(O):
        
        ax3.plot([births_2[i],births_2[i]],[births_2[i],deaths_2[i]],'r--')
        
    for i in range(P):
        
        ax4.plot([births_3[i],births_3[i]],[births_3[i],deaths_3[i]],'r--')
        

    plt.show()










