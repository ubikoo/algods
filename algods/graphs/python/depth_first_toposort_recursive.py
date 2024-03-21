#!/usr/bin/env ipython

# depth_first_toposort_recursive.py
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

# -----------------------------------------------------------------------------
class DepthFirstTopoSortRecursive:
    """
    Compute the undirected connected components of a graph using breadth first
    search algorithm.
    """

    # ---- Special methods ----------------------------------------------------
    def __init__(self, graph):
        self.NEW = 0
        self.VISITED = 1
        self._state = [self.NEW for _ in graph.vertices()]
        self._vertices = Stack()
        for v in graph.vertices():
            if self._state[v] == self.NEW:
                self._dfs(graph, v)

    def __del__(self):
        pass

    def __iter__(self):
        """ Iterate over the ordered vertices """
        for v in self._vertices:
            yield v;

    def __repr__(self):
        """ Return a string representation of the ordered vertex stack """
        return repr(self._vertices)

    # ---- API ----------------------------------------------------------------
    def is_visited(self, v):
        """ Has the specified vertex been visited? """
        return self._state[v] == self.VISITED

    def _dfs(self, graph, v):
        """ Compute the vertices reachable from a source vertex s """
        self._state[v] = self.VISITED
        for w in graph.adj(v):
            if self._state[w] == self.NEW:
                self._dfs(graph, w)
        self._vertices.push(v)

# -----------------------------------------------------------------------------
def draw(toposort, graph):
    """ Draw the graph using networkx package """
    # Create the networkx graph
    G = networkx.DiGraph()
    for v in toposort:
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
    pos = networkx.spring_layout(G, k=3.0, iterations=100, seed=63)
    networkx.draw_networkx(G, pos, **options)
    pyplot.axis("off")
    pyplot.show()
    # pyplot.savefig('fig.pdf', bbox_inches='tight')


# -----------------------------------------------------------------------------
def main (argv):
    """ DepthFirstTopoSortRecursive test client """

    # Create a graph from stream
    with open('../data/tinyDAG.txt') as stream:
        graph = DiGraph.from_stream(stream)
    print(graph)

    # Compute the list of vertices reachable from 0
    toposort = DepthFirstTopoSortRecursive(graph)
    print(toposort)
    # draw(toposort, graph)

    # Done
    return 0

if __name__ == '__main__':
    sys.exit(main(sys.argv))
