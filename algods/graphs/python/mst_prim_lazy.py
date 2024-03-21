#!/usr/bin/env ipython

# mst_prim_lazy.py
#
# Copyright (c) 2020 Carlos Braga
# This program is free software; you can redistribute it and/or modify it under
# the terms of the MIT License. See accompanying LICENSE.md or
# https://opensource.org/licenses/MIT.

import sys
import collections
import random
import numpy
import matplotlib
from matplotlib import pyplot
from mpl_toolkits.mplot3d import axes3d
import networkx

from ds.edge_weighted_graph import Edge, EdgeWeightedGraph
from ds.priority_queue import PriorityQueue
from ds.union_find import UnionFind

# -----------------------------------------------------------------------------
class MstPrimLazy:
    """
    MstPrimLazy implements the lazy version of Prims algorithm to compute the
    minimum spanning tree of an edge weighted graph.
    """

    # ---- MstPrimLazy special methods ----------------------------------------
    def __init__(self, graph):
        # Graph minimum spanning tree
        self._mst = []

        # Create a priority queue of crossing edges in ascending order
        # self._pq = PriorityQueue("min")
        self._pq = PriorityQueue(graph.n_edges(), (lambda x, y: x < y))

        # Array of vertex states
        self.NEW = 0
        self.VISITED = 1
        self._state = [self.NEW for _ in graph.vertices()]

        # Visit the first vertex assuming the graph is connected. Keep visiting
        # vertices while there are crossing edges between the MST and the graph
        self._visit(graph, 0)
        while not self._pq.is_empty():
            edge = self._pq.top()
            self._pq.pop()

            # If both vertices have been visited alreay, its not a crossing
            # edge anymore - discard it and continue
            v = edge.first()
            w = edge.second()
            if self._is_visited(v) and self._is_visited(w):
                continue

            # Add the edge to the minimum spanning tree
            self._mst.append(edge)

            # Either vertex v or vertex w is visited (but not both, see above)
            # Visit the unvisited vertex and add any crossing edges to the queue
            if not self._is_visited(v):
                self._visit(graph, v)
            if not self._is_visited(w):
                self._visit(graph, w)

    def __del__(self):
        pass

    def __iter__(self):
        for edge in self.edges():
            yield edge

    def __repr__(self):
        """ Return a string representation of this MstPrimLazy """
        s = []
        for edge in self.edges():
            s.append(
                "%d %d %lf\n" % (edge.first(), edge.second(), edge.weight()))
        return ''.join(s)

    def __len__(self):
        """ Return the number of keys in the queue """
        return len(self._mst)

    # ---- MstPrimLazy api ----------------------------------------------------
    def edges(self):
        """ Return the edges in the minimum spanning tree """
        return self._mst

    def weight(self):
        """ Return the weight of the minimum spanning tree """
        w = 0.0
        for edge in self.edges():
            w += edge.weight()
        return w

    def _is_visited(self, v):
        """ Has the specified vertex been visited? """
        return self._state[v] == self.VISITED

    def _visit(self, graph, v):
        """
        Visit the specified vertex in the graph.
        For every edge incident to the vertex, if the edge is a crossing edge,
        i.e., the other vertex has not been visited, add the edge to the queue.
        """
        self._state[v] = self.VISITED
        for edge in graph.adj(v):
            if not self._is_visited(edge.other(v)):
                self._pq.push(edge)


# -----------------------------------------------------------------------------
def draw(graph, mst):
    """ Draw the graph using networkx package """
    # Create the networkx graph
    G = networkx.Graph()
    for edge in graph.edges():
        G.add_edge(edge.first(), edge.second())

    # Create the minimum spanning tree
    T = networkx.Graph()
    for edge in mst:
        T.add_edge(edge.first(), edge.second())

    # Draw the networkx graph
    pyplot.figure(figsize=(7,6))
    Gopts = {
        "font_size": 10,
        "font_color": "#000000",
        "with_labels": False,
        "node_size": 10,
        "node_shape": "o",
        "node_color": "#eeeeee",
        "edge_color": "#eeeeee",
        "linewidths": 1,
        "width": 1,
        "alpha" : 0.2,
    }
    Topts = {
        "font_size": 10,
        "font_color": "#000000",
        "with_labels": False,
        "node_size": 10,
        "node_shape": "o",
        "node_color": "#ff0000",
        "edge_color": "#000000",
        "linewidths": 1,
        "width": 1,
        "alpha" : 0.6,
    }

    pos = networkx.spring_layout(T, k=2.0, iterations=10000, seed=63)
    networkx.draw_networkx(G, pos, **Gopts)
    networkx.draw_networkx(T, pos, **Topts)

    pyplot.axis("off")
    pyplot.show()
    # pyplot.savefig('fig.pdf', bbox_inches='tight')

# -----------------------------------------------------------------------------
def main (argv):
    """ MstPrimLazy test client """
    with open('../data/mediumEWG.txt') as stream:
        graph = EdgeWeightedGraph.from_stream(stream)
        mst = MstPrimLazy(graph)
    print(graph)
    print("mst\n", mst)
    print("mst weight", mst.weight())
    draw(graph, mst)

    # Done
    return 0

if __name__ == '__main__':
    sys.exit(main(sys.argv))
