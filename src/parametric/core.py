import networkx as nx
from gudhi.simplex_tree import SimplexTree
from itertools import product

from discrete.core import v_paths, critical_cells, gradient
from helpers import get_skeleton, build_function_series


def has_correct_dims(path: list, k: int) -> bool:
    """
    Returns True if path contains simplices of alternating dimensions k + 1 and k, else False.
    """

    even_indexed_checks = [len(s) == k + 1 for i, s in enumerate(path) if i % 2 == 0]
    odd_indexed_checks = [len(s) == k for i, s in enumerate(path) if i % 2 == 1]
    return all(even_indexed_checks + odd_indexed_checks)


def are_connected(K: SimplexTree, V: list, W: list, s1: list, s2: list) -> bool:
    """
    Returns True if s1 and s2 are connected, False otherwise (according to def.3.1. of PMT).

    Args:
        K : simplicial complex
        V : list of pairs of cells
        W : list of pairs of cells
        s1 : cell critical for V
        s2 : cell critical for W
    """

    if len(s1) != len(s2):
        return False

    if s1 == s2:
        return True

    else:
        k = len(s1) - 1
        k_simplices = [s[0] for s in get_skeleton(K, k)]
        for s in k_simplices:
            v_paths1 = v_paths(K, V, s1, s)
            v_paths2 = v_paths(K, W, s, s2)
            if any([has_correct_dims(p1, k) and has_correct_dims(p2, k + 1) for p1, p2 in product(v_paths1, v_paths2)]):
                return True
            if any([p1 == [s1] or p2 == [s2] for p1, p2 in product(v_paths1, v_paths2)]):
                return True

    return False


def parametric_coordinates(K: SimplexTree, V: list, C: list, drawing=False) -> list:
    """
    Returns a list of birth-death coordinates of all critical cells appearing in the sequence of DVFs on K.
    """

    temp_critical_cells = [[s, i] for i, s_ in enumerate(C) for s in s_]
    critical_cells = temp_critical_cells.copy()

    coordinates = []
    full_coordinates = []

    G = abstract_diagram(K, V, C)

    for indexed_cell in critical_cells:
        if indexed_cell in temp_critical_cells:
            for node_ in G.nodes:
                if G.nodes[node_]["index"] == indexed_cell:
                    node = node_
                    break

            magic_path = diagram_path(G, node, len(C))
            time = indexed_cell[1]

            if len(magic_path) == 1:
                coordinates.append([time, time])
                full_coordinates.append([indexed_cell[0], time, time])
                cell = G.nodes[node_]["index"]
                if cell in temp_critical_cells:
                    temp_critical_cells.remove(cell)

            else:
                for cell in magic_path:
                    index_cell_ = G.nodes[cell]["index"]
                    if index_cell_ in temp_critical_cells:
                        temp_critical_cells.remove(index_cell_)

                coordinates.append([time, time + len(magic_path) - 1])
                full_coordinates.append([indexed_cell[0], time, time + len(magic_path) - 1])

    if drawing:
        nx.draw(G)

    return [coordinates, full_coordinates]


def per_time_slice(G: nx.Graph, i: int) -> list:
    slice_nodes = [v for v in list(G.nodes) if G.nodes[v]["index"][1] == i]
    return slice_nodes


def abstract_diagram(K: SimplexTree, V: list, C: list) -> nx.Graph:
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


def diagram_path(graph: nx.Graph, node: int, nb_slices: int) -> list:
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


def parametric_pipeline(
        G: nx.Graph,
        functions: list = [],
        start: int = 0,
        time_step: float = 0.05,
        count: int = 12,
        noise: float = 0.05) -> list:
    list_dmts = build_function_series(G, functions, start, time_step, count, noise)

    C = list(map(lambda dmt: critical_cells(K=dmt['complex'], f=dmt['f']), list_dmts))
    V = list(map(lambda dmt: gradient(K=dmt['complex'], f=dmt['f']), list_dmts))

    return parametric_coordinates(K=list_dmts[-1]['complex'], V=V, C=C)
