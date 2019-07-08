from ._misc import powerset

__all__ = ['spanning_trees']


def spanning_trees(graph, verts_in=None, edge_list=None, cur_vert=0):
    """
    Yields all spanning trees of a given graph.

    This function returns a generator that yields spanning trees,
    as a list of edges of the graph.

    This function should be called only with the input graph as
    its argument.  The other arguments are used for its recursion.

    Parameters
    ----------
    graph : NetworkX graph
        This graph should be a simple graph.  The correctness of this
        function has not been tested with other graphs.

    verts_in : list, optional (default=None)
        The current list of vertices to have their neighbors checked.

    edge_list : list, optional (default=None)
        The current list of edges in the (potential) spanning tree.

    cur_vert : integer, optional (default=None)
        The index in `verts_in` that will be checked.  This value
        increases by 1 with each iteration.

    Yields
    ------
    tree : list
        A list of edges in the spanning tree as pairs of vertices.
    """
    if edge_list is None:
        degree_one = [n for n, d in graph.degree() if d == 1]
        init_vert = []
        init_edge = []
        start_index = 0
        if degree_one:
            for node in degree_one:
                if node not in init_vert:
                    neighbor = list(graph[node])[0]
                    init_edge.append((node, neighbor))
                    init_vert.append(neighbor)
            # is it necessary to sort init_verts?
            init_vert = degree_one + init_vert[:1]
            start_index = len(degree_one)
        else:
            init_vert.append(list(graph)[0])
        for tree in spanning_trees(graph, init_vert, init_edge, start_index):
            yield tree
    else:
        if len(edge_list) == len(graph)-1:
            yield edge_list  # is it necessary to sort edgelist?
        else:
            if cur_vert < len(verts_in):
                neighbors = [n for n in graph[verts_in[cur_vert]] if n not in verts_in]
                for subset in powerset(neighbors):
                    new_edges = [(verts_in[cur_vert], n) for n in subset]
                    for tree in spanning_trees(graph, verts_in + list(subset), edge_list + new_edges, cur_vert + 1):
                        yield tree


# use nx.connected_components() instead.
# def components(graph):
#     """
#     Returns representatives of the components of a graph as a list.
#
#     Every vertex in the graph will be connected to exactly one
#     element of the returned list.  As such, the length of this
#     list is precisely the number of components of the graph.
#
#     Parameters
#     ----------
#     graph : NetworkX graph
#         This graph should be a simple graph.  The correctness of this
#         function has not been tested with other graphs.
#
#     Returns
#     -------
#     seeds : list
#          A list of representatives of the components of the graph.
#          This is chosen arbitrarily.
#     """
#     import networkx as nx
#     graph = nx.Graph(graph)
#     visited = []
#     seeds = []
#     for e in graph.nodes():
#         if e in visited:
#             continue
#         seeds.append(e)
#         stack = [e]
#         while stack:
#             cur_node = stack.pop()
#             visited.append(cur_node)
#             neighbors = list(graph[cur_node])
#             stack.extend([x for x in neighbors if x not in visited])
#     return seeds
