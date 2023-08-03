import networkx as nx
from simplicial import *
import matplotlib.pyplot as plt
from DMT import V_paths, critical_cells, gradient
from helpers import get_skel, magic_function

plt.rcParams["figure.figsize"] = (10, 10)


def right_dimension(path, k):
    for i in range(len(path)):
        if i % 2 == 0 and len(path[i]) != k + 1:
            return False

        if i % 2 == 1 and len(path[i]) != k:
            return False

    return True


def are_connected(K, V, W, s1, s2):
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

    if len(s1) != len(s2):
        return False

    if s1 == s2:
        return True

    else:
        k = len(s1) - 1

        k_simplices = [simplex[0] for simplex in get_skel(K, k)]

        for simplex in k_simplices:
            # print(simplex)

            for path1 in V_paths(K, V, s1, simplex):
                for path2 in V_paths(K, W, simplex, s2):
                    if path1 == [s1] or path2 == [s2]:
                        # print("Path from {} to {} via the subpaths {} and {}.".format(
                        # s1,s2,path1,path2))

                        return True

                    elif right_dimension(path1, k) and right_dimension(path2, k + 1):
                        # print("WE have a path from {} to {} via the subpaths {} and {}.".format(
                        # s1,s2,path1,path2))

                        return True

    return False


def parametric(K, V, C, drawing=False):
    "Returns a list of birth-death coordinates of all critical cells appearing along the sequence of DVFs on K."

    temp_critical_cells = []
    critical_cells = []
    nb_slices = len(C)

    coordinates = []
    full_coordinates = []

    for i in range(nb_slices):
        for s in C[i]:
            temp_critical_cells.append([s, i])
            critical_cells.append([s, i])

    graph = abstract_diagram(K, V, C)

    for indexed_cell in critical_cells:
        if indexed_cell in temp_critical_cells:
            for node_ in graph.nodes:
                if graph.nodes[node_]["index"] == indexed_cell:
                    node = node_
                    break

            magic_path = diagram_path(graph, node, nb_slices)
            time = indexed_cell[1]

            if len(magic_path) == 1:
                coordinates.append([time, time])
                full_coordinates.append([indexed_cell[0], time, time])
                cell = graph.nodes[node_]["index"]
                if cell in temp_critical_cells:
                    temp_critical_cells.remove(cell)

            else:
                for cell in magic_path:
                    index_cell_ = graph.nodes[cell]["index"]
                    if index_cell_ in temp_critical_cells:
                        temp_critical_cells.remove(index_cell_)

                coordinates.append([time, time + len(magic_path) - 1])
                full_coordinates.append(
                    [indexed_cell[0], time, time + len(magic_path) - 1]
                )

    if drawing:
        nx.draw(graph)

    return [coordinates, full_coordinates]


def per_time_slice(graph, i):
    slice_nodes = [
        node for node in list(graph.nodes) if graph.nodes[node]["index"][1] == i
    ]

    return slice_nodes


def abstract_diagram(K, V, C):
    critical_cells = []

    for i in range(len(C)):
        for s in C[i]:
            critical_cells.append([s, i])

    N = len(critical_cells)

    graph = nx.DiGraph()
    node = 1
    for cell in critical_cells:
        graph.add_node(node, index=cell)
        node += 1

    for node1 in range(1, N + 1):
        for node2 in range(1, N + 1):
            i = graph.nodes[node1]["index"][1]
            j = graph.nodes[node2]["index"][1]

            if j - i == 1:
                s1 = graph.nodes[node1]["index"][0]
                s2 = graph.nodes[node2]["index"][0]

                if are_connected(K, V[i], V[j], s1, s2) and are_connected(
                    K, V[j], V[i], s2, s1
                ):
                    graph.add_edge(node1, node2)

    return graph


def diagram_path(graph, node, nb_slices):
    path = [node]
    temp_node = node

    for i in range(nb_slices):
        slice_nodes = per_time_slice(graph, i + 1)

        for next_node in slice_nodes:
            if (temp_node, next_node) in graph.edges:
                path.append(next_node)
                temp_node = next_node
                break

    return path


def pipeline(G, functions=[], start=0, time_step=0.05, time_count=12, noise=0.05):
    list_dmts = magic_function(G, functions, start, time_step, time_count, noise)

    C = []
    V = []

    for i in range(len(list_dmts)):
        K = list_dmts[i][0]
        f = list_dmts[i][1]
        C.append(critical_cells(K, f))
        V.append(gradient(K, f))

    return parametric(K, V, C)


def plot_PD(coordinates):
    maximum = max([trio[2] for trio in coordinates])

    coord_0 = [trio for trio in coordinates if len(trio[0]) == 1]
    coord_1 = [trio for trio in coordinates if len(trio[0]) == 2]
    coord_2 = [trio for trio in coordinates if len(trio[0]) == 3]
    coord_3 = [trio for trio in coordinates if len(trio[0]) == 4]

    births_0 = [trio[1] for trio in coord_0]
    deaths_0 = [trio[2] for trio in coord_0]
    births_1 = [trio[1] for trio in coord_1]
    deaths_1 = [trio[2] for trio in coord_1]
    births_2 = [trio[1] for trio in coord_2]
    deaths_2 = [trio[2] for trio in coord_2]
    births_3 = [trio[1] for trio in coord_3]
    deaths_3 = [trio[2] for trio in coord_3]

    z = range(maximum + 1)

    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, sharex=True, sharey=True)

    ax1.scatter(births_0, deaths_0, color="skyblue", alpha=0.5)
    ax1.plot(z, z, color="orange")
    ax1.set_title("0-cells")
    ax2.scatter(births_1, deaths_1, color="skyblue", alpha=0.5)
    ax2.plot(z, z, color="orange")
    ax2.set_title("1-cells")
    ax3.scatter(births_2, deaths_2, color="skyblue", alpha=0.5)
    ax3.plot(z, z, color="orange")
    ax3.set_title("2-cells")
    ax4.scatter(births_3, deaths_3, color="skyblue", alpha=0.5)
    ax4.plot(z, z, color="orange")
    ax4.set_title("3-cells")
    fig.suptitle("Life coordinates of critical cells")

    for i in range(len(coord_0)):
        ax1.plot([births_0[i], births_0[i]], [births_0[i], deaths_0[i]], "m--")

    for i in range(len(coord_1)):
        ax2.plot([births_1[i], births_1[i]], [births_1[i], deaths_1[i]], "m--")

    for i in range(len(coord_2)):
        ax3.plot([births_2[i], births_2[i]], [births_2[i], deaths_2[i]], "m--")

    for i in range(len(coord_3)):
        ax4.plot([births_3[i], births_3[i]], [births_3[i], deaths_3[i]], "m--")

    plt.show()
