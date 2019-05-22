from ._misc import powerset

__all__ = ['spanning_trees', 'components']


def spanning_trees(graph, verts_in=None, edge_list=None, cur_vert=0):
    """Generate and yield all spanning trees of a given graph."""
    if edge_list is None:
        degree_one = [n for n, d in graph.degree if d == 1]
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


def components(graph):
    import networkx as nx
    graph = nx.Graph(graph)
    visited = []
    seeds = []
    for e in graph.nodes():
        if e in visited:
            continue
        seeds.append(e)
        stack = [e]
        while stack:
            cur_node = stack.pop()
            visited.append(cur_node)
            neighbors = list(graph[cur_node])
            stack.extend([x for x in neighbors if x not in visited])
    return seeds
