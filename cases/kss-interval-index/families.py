"""Graph families from arXiv:2610.39966."""
import networkx as nx


def dumbbell(n):
    """D_n (n >= 8): two disjoint C4 joined by a path of length n-7."""
    G = nx.Graph()
    G.add_edges_from([(0, 1), (1, 2), (2, 3), (3, 0)])            # first C4, attach at 0
    off = 4
    G.add_edges_from([(off, off + 1), (off + 1, off + 2), (off + 2, off + 3), (off + 3, off)])  # second C4, attach at 4
    L = n - 7  # path length (edges) between 0 and 4 ; internal vertices: L-1
    path = [0] + list(range(8, 8 + L - 1)) + [4]
    nx.add_path(G, path)
    assert G.number_of_nodes() == n, (n, G.number_of_nodes())
    return G


def two_triangles(n):
    G = nx.Graph()
    G.add_edges_from([(0, 1), (1, 2), (2, 0), (3, 4), (4, 5), (5, 3)])
    L = n - 5
    path = [0] + list(range(6, 6 + L - 1)) + [3]
    nx.add_path(G, path)
    assert G.number_of_nodes() == n
    return G


