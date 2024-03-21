#!/usr/bin/env ipython

# mst_kruskal.py
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
class MstKruskal:
    """
    MstKruskal implements the Kruskal algorithm to compute the minimum spanning
    tree of an edge weighted graph.
    """

    # ---- MstKruskal special methods ------------------------------------------
    def __init__(self, graph):
        """
        Create a priority queue of edges in ascending order with minimum weight
        edge on top priority.
        Create a union find structure of connected components of the minimum
        spanning tree.
        """
        # Graph minimum spanning tree
        self._mst = []

        # Minimum spanning tree components
        self._uf = UnionFind(graph.n_vertices())

        # Create a priority queue of edges in ascending order
        self._pq = PriorityQueue(graph.n_edges(), (lambda x, y: x < y))
        # print(graph.edges())
        for edge in graph.edges():
            self._pq.push(edge)
        # print(self._pq)

        # Add edges to the mst from the priority queue unless it creates a cycle
        while (not self._pq.is_empty()
               and len(self._mst) < graph.n_vertices() - 1):
            edge = self._pq.top()
            self._pq.pop()

            # If edge vertices belong to the same component, the edge creates
            # a cycle - discard it and continue
            v = edge.first()
            w = edge.second()
            if self._uf.connected(v, w):
                continue

            # Merge the components of v and w to the same tree
            self._uf.join(v, w)

            # Add the edge to the minimum spanning tree
            self._mst.append(edge)

    def __del__(self):
        pass

    def __iter__(self):
        for edge in self.edges():
            yield edge

    def __repr__(self):
        """ Return a string representation of this MstKruskal """
        s = []
        for edge in self.edges():
            s.append(
                "%d %d %lf\n" % (edge.first(), edge.second(), edge.weight()))
        return ''.join(s)

    def __len__(self):
        """ Return the number of keys in the queue """
        return len(self._mst)

    # ---- MstKruskal api -----------------------------------------------------
    def edges(self):
        """ Return the edges in the minimum spanning tree """
        return self._mst

    def weight(self):
        """ Return the weight of the minimum spanning tree """
        w = 0.0
        for edge in self.edges():
            w += edge.weight()
        return w


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
    """ MstKruskal test client """
    with open('../data/mediumEWG.txt') as stream:
        graph = EdgeWeightedGraph.from_stream(stream)
        mst = MstKruskal(graph)
    print(graph)
    print("mst\n", mst)
    print("mst weight", mst.weight())
    draw(graph, mst)

    # Done
    return 0

if __name__ == '__main__':
    sys.exit(main(sys.argv))
