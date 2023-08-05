import networkx as nx
import numpy as np
import gudhi as gd
from scipy.stats import bernoulli


def add_stochastic_community(interval1: tuple[int, int], interval2: tuple[int, int], proba: float, G: nx.Graph) -> None:
    for i in range(interval1[0], interval1[1]):
        for j in range(interval2[0], interval2[1]):
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
            complex_series[j].assign_filtration(nodes[i][0], functions[i](T[j]))

    return complex_series


def get_skel(K, p):
    skel = list(K.get_skeleton(p))
    skel_new = []
    for simplex in skel:
        if len(simplex[0]) == p + 1:
            skel_new.append(simplex)
    return skel_new


def build_morse_function(K, d, g, noise):
    """
    K : simplicial complex (gudhi.SimplexTree)
    d : dimension of K
    g : labels on nodes given by a dictionary
    """

    Flag = {}
    f = {}

    for p in range(d):
        for simplex in list(K.get_skeleton(p)):
            Flag[str(simplex[0])] = 0

    for simplex in list(K.get_skeleton(0)):
        f[str(simplex[0])] = g[str(simplex[0])]

    for p in range(1, d + 1):
        for simplex in get_skel(K, p):
            Faces = [face[0] for face in list(K.get_boundaries(simplex[0]))]
            Faces = [str(face) for face in Faces]
            Faces = sorted(Faces, key=lambda x: f[x], reverse=True)

            gamma_0 = Faces[0]
            gamma_1 = Faces[1]

            if Flag[Faces[0]] == 0 and f[gamma_0] > f[gamma_1]:
                f[str(simplex[0])] = (f[gamma_0] + f[gamma_1]) / 2

                Flag[gamma_0] = 1

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
            simplex = get_skel(series[i], 0)[j]
            g[str(simplex[0])] = simplex[1]

        list_g.append(g)

    for i in range(len(series)):
        f = build_morse_function(series[i], series[i].dimension(), list_g[i], noise)[0]

        list_f.append(f)

    list_pairs = [[series[i], list_f[i]] for i in range(len(series))]

    # print("We have a sequence of {} complexes.\n".format(len(series)))

    return list_pairs
