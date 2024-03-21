#!/usr/bin/env ipython

# depth_first_scc_recursive.py
#
# Copyright (c) 2020 Carlos Braga
# This program is free software; you can redistribute it and/or modify it under
# the terms of the MIT License. See accompanying LICENSE.md or
# https://opensource.org/licenses/MIT.

import sys
import networkx
import numpy
from matplotlib import pyplot
from ds.digraph import DiGraph
from ds.queue import Queue
from ds.stack import Stack
from depth_first_toposort_recursive import DepthFirstTopoSortRecursive

# -----------------------------------------------------------------------------
class DepthFirstSCC:
    """
    Compute the undirected connected components of a graph using breadth first
    search algorithm.
    """

    # ---- Special methods ----------------------------------------------------
    def __init__(self, graph):
        self.NEW = 0
        self.VISITED = 1

        # First pass - compute the topological order of the reverse graph
        toposort = DepthFirstTopoSortRecursive(DiGraph.reverse(graph))

        # Second pass - compute SCCs of graph in reverse topological order
        self._state = [self.NEW for _ in graph.vertices()]
        self._id = [numpy.inf for _ in graph.vertices()]
        self._count = 0
        for v in toposort:
            if not self.is_visited(v):
                self._dfs(graph, v)
                self._count += 1

    def __del__(self):
        pass

    # ---- API ----------------------------------------------------------------
    def is_visited(self, v):
        """ Has the specified vertex been visited? """
        return self._state[v] == self.VISITED

    def connected(self, v, w):
        """ Are v and w connected? """
        return self._id[v] == self._id[w]

    def id(self, v):
        """ Return the component id for the specified vertex """
        return self._id[v]

    def count(self):
        """ Return the number of strongly connnected components """
        return self._count

    def _dfs(self, graph, v):
        """ Compute the vertices reachable from a source vertex s """
        self._state[v] = self.VISITED
        self._id[v] = self._count
        for w in graph.adj(v):
            if self._state[w] == self.NEW:
                self._dfs(graph, w)


# -----------------------------------------------------------------------------
def draw(dfs, graph, cc):
    """ Draw the graph using networkx package """
    # Create the networkx graph
    G = networkx.DiGraph()
    for v in graph.vertices():
        if dfs.id(v) == cc:
            for w in graph.adj(v):
                G.add_edge(v, w)

    # Draw the networkx graph
    fig = pyplot.figure(figsize=(7,6))
    options = {
        "font_size": 10,
        "font_color": "#000000",
        "with_labels": False,
        "node_size": 10,
        "node_shape": "o",
        "node_color": "#ee0000",
        "edge_color": "#1f78b4",
        "linewidths": 1,
        "width": 1,
        "alpha" : 0.6,
    }
    pos = networkx.spring_layout(G, k=2.0, iterations=1000, seed=63)
    networkx.draw_networkx(G, pos, **options)
    pyplot.axis("off")
    pyplot.show()
    # pyplot.savefig('fig.pdf', bbox_inches='tight')


# -----------------------------------------------------------------------------
def main (argv):
    """ DepthFirstSCC test client """

    # Create a graph from stream
    with open('../data/testSCC.txt') as stream:
        graph = DiGraph.from_stream(stream)
    print(graph)

    # Compute the list of vertices reachable from 0
    dfs = DepthFirstSCC(graph)
    print("\n strong components ", dfs.count(), file=sys.stdout, end='\n')
    for scc in range(dfs.count()):
        for v in graph.vertices():
            if dfs.id(v) == scc:
                print(v, file=sys.stdout, end=' ')
        print(file=sys.stdout, end='\n')
    # for scc in range(dfs.count()):
    #     draw(dfs, graph, scc)

    # Compute the five largest strong components in the graph
    scc_sizes = [0 for _ in range(dfs.count())]
    for v in graph.vertices():
        scc_sizes[dfs.id(v)] += 1
    scc_sizes.sort(reverse=True)
    print("scc_sizes", scc_sizes)

    # Done
    return 0

if __name__ == '__main__':
    sys.exit(main(sys.argv))
