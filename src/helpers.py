import networkx as nx
import numpy as np
import gudhi as gd
from scipy.stats import bernoulli
from gudhi.simplex_tree import SimplexTree


def add_stochastic_community(interval1: tuple[int, int], interval2: tuple[int, int], proba: float, G: nx.Graph) -> None:
    """
    Adds to the graph G a community (or 'block') built by forming edges (i,j) is in interval1\times interval2.

    Args:
        - G: graph on which to add the block
        - interval1: range of node indices
        - interval2: range of node indices
        - proba: edge probability with which the block is built

    Note:
        Only edges are added with this procedure, not nodes.
    """
    for i in range(int(interval1[0]), int(interval1[1])):
        for j in range(int(interval2[0]), int(interval2[1])):
            if i < j and bernoulli.rvs(proba):
                G.add_edge(i, j)


def stochastic_block_model(n: int, p: float, q: float, k: int) -> nx.Graph:
    """
    Returns the planted partition model (https://en.wikipedia.org/wiki/Stochastic_block_model)
    with parameters n, p, q and k defined as follows.

    Args:
        - n : number of nodes in the graph
        - p : within-communities edge probability
        - q : between-communities edge probability
        - k : number of communities (blocks)
    """
    if (n % k) != 0 or not (0 <= p <= 1 and 0 <= q <= 1):
        raise ValueError("n must be a multiple of k, and p and q must be in [0,1].")

    G = nx.Graph()
    G.add_nodes_from(range(n))
    for i in range(k):
        add_stochastic_community(interval1=(i, i + 1), interval2=(i, i + 1), proba=p, G=G)

    for c1 in range(k):
        for c2 in range(k):
            if c1 < c2:
                add_stochastic_community(interval1=(c1 * n / k, (c1 + 1) * n / k),
                                         interval2=(c2 * n / k, (c2 + 1) * n / k),
                                         proba=q,
                                         G=G)
    return G


def build_complex_series(G: nx.Graph, functions: list, start: int, time_step: float, count: int) -> list:
    """
    Args:
        - G : networkX.Graph(), a graph
        - functions : list of functions (one function for each node of G)
        - (start, time_step, count) : define the list of time slices (for which we sample the functions)

    Returns:
        This function returns a sequence of node labelings on K, where K is
        the clique complex built on the graph G (there is a node labeling
        for each time slice in T), put in the form of an array np.ndarray()
        of size |G| x |T|.
    """
    T = list(map(lambda i: start + i * time_step, range(count)))
    n = G.number_of_nodes()
    assert n == len(functions), "There has to be a function for each node."

    K = gd.SimplexTree()
    for e in G.edges:
        K.insert(e)
    for n in G.nodes:
        K.insert([n])
    K.expansion(3)

    complex_series = count * [K]

    for j in range(count):
        nodes = list(complex_series[j].get_skeleton(0))

        for i in range(n):
            complex_series[j].assign_filtration(nodes[i][0], functions[i](T[j]))

    return complex_series


def get_skeleton(K: SimplexTree, p: int) -> list:
    skel = list(K.get_skeleton(p))
    skel_new = [s for s in skel if len(s[0]) == p + 1]
    return skel_new


def build_morse_function(K: SimplexTree, d: int, g: dict, noise: float) -> dict:
    """
    Args:
        - K : simplicial complex (gudhi.SimplexTree)
        - d : dimension of K
        - g : labels on nodes given by a dictionary
        - noise: noise term to ensure injectivity of resulting morse function

    Returns:
        A dictionary containing both the flag complex of K and the function f we built on it.
    """

    flag = {}
    f = {}

    for p in range(d):
        for simplex in list(K.get_skeleton(p)):
            flag[str(simplex[0])] = 0

    for simplex in list(K.get_skeleton(0)):
        f[str(simplex[0])] = g[str(simplex[0])]

    for p in range(1, d + 1):
        for simplex in get_skeleton(K, p):
            faces = sorted([str(face[0]) for face in list(K.get_boundaries(simplex[0]))],
                           key=lambda x: f[x], reverse=True)
            gamma_0 = faces[0]
            gamma_1 = faces[1]

            if flag[faces[0]] == 0 and f[gamma_0] > f[gamma_1]:
                f[str(simplex[0])] = (f[gamma_0] + f[gamma_1]) / 2
                flag[gamma_0] = 1
            else:
                epsilon = np.random.uniform(low=0, high=noise)
                f[str(simplex[0])] = f[gamma_0] + epsilon

    return {'flag': flag, 'f': f}


def build_function_series(G: nx.Graph,
                          functions: list,
                          start: int,
                          time_step: float,
                          count: int,
                          noise: float) -> list[dict]:
    """
    Returns a list of dictionaries of the form {'complex': K_i, 'f': f_i} where f_i is a discrete Morse function that
    is built on the simplicial complex K_i, where the latter results from applying function functions[i] on the graph G.

    Args:
        - G : networkX.Graph(), a graph
        - functions : list of functions (one function for each node of G)
        - (start, time_step, count) : define the list of time slices (for which we sample the functions)
        - noise: noise term to ensure injectivity of resulting morse functions
    """

    series = build_complex_series(G, functions, start, time_step, count)

    r = range(len(series))

    def build_vertex_function(i) -> dict:
        g = {}
        for j in range(G.number_of_nodes()):
            simplex = get_skeleton(series[i], 0)[j]
            g[str(simplex[0])] = simplex[1]
        return g

    list_g = list(map(lambda i: build_vertex_function(i), r))
    list_f = list(map(lambda i: build_morse_function(K=series[i],
                                                     d=series[i].dimension(),
                                                     g=list_g[i],
                                                     noise=noise)['f'], r))
    list_complex_functions = [{'complex': series[i], 'f': list_f[i]} for i in r]

    return list_complex_functions
