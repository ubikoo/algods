#!/usr/bin/env ipython

# mst_prim_eager.py
#
# Copyright (c) 2020 Carlos Braga
# This program is free software; you can redistribute it and/or modify it under
# the terms of the MIT License. See accompanying LICENSE.md or
# https://opensource.org/licenses/MIT.

import select
import sys
import collections
import random
import numpy
import matplotlib
from matplotlib import pyplot
from mpl_toolkits.mplot3d import axes3d
import networkx

from ds.edge_weighted_graph import Edge, EdgeWeightedGraph
from ds.priority_index_queue import PriorityIndexQueue
from ds.union_find import UnionFind

# -----------------------------------------------------------------------------
class MstPrimEager:
    """
    MstPrimEager implements the eager version of Prims algorithm to compute the
    minimum spanning tree of an edge weighted graph.
    """

    # ---- MstPrimEager special methods ---------------------------------------
    def __init__(self, graph):
        # Vertex state
        self.NEW = 0
        self.VISITED = 1
        self._state = [self.NEW for _ in graph.vertices()]

        # Shortest edge from the tree to a vertex
        self._edge_to = [None for _ in graph.vertices()]

        # Shortest distance from the tree to a vertex
        # self._dist_to[v] = self._edge_to[v].weight()
        self._dist_to = [numpy.inf for _ in graph.vertices()]

        # Indexed priority queue of crossing edges in ascending order
        # self._pq = PriorityIndexQueue("min", graph.n_vertices())
        self._pq = PriorityIndexQueue(graph.n_vertices(), (lambda x, y: x < y))

        # Initialize pq with vertex 0, and corresponding weight 0.0
        self._dist_to[0] = 0.0
        self._pq.push(0, 0.0)
        while not self._pq.is_empty():
            v = self._pq.top_ix()
            self._pq.pop()
            self._visit(graph, v)

        # Compute the minimum spanning tree from the set of shortest edges
        # from the tree to every vertex in the graph. The first vertex has
        # no incident vertex, therefore it is ignored.
        self._mst = []
        for v in range(1, graph.n_vertices()):
            self._mst.append(self._edge_to[v])

    def __del__(self):
        pass

    def __iter__(self):
        for edge in self.edges():
            yield edge

    def __repr__(self):
        """ Return a string representation of this MstPrimEager """
        s = []
        for edge in self.edges():
            s.append("%d %d %lf\n" % (edge.first(), edge.second(), edge.weight()))
        return ''.join(s)

    def __len__(self):
        """ Return the number of keys in the queue """
        return len(self._mst)

    # ---- MstPrimEager api ---------------------------------------------------
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
        Visit the specified vertex in the graph. For every edge incident to
        vertex-v, if the edge is a crossing edge, i.e. the other vertex-w has
        not been visited, do:

        If the edge weight is smaller than the current weight distance to w,
        the edge is the new best connection of w to the tree T. Add it to
        the priority queue keyed by its distance or updated it if it already
        exists.
        """
        self._state[v] = self.VISITED
        for edge in graph.adj(v):
            w = edge.other(v)

            # Vertex w is already connected to the spanning tree
            if self._is_visited(w):
                continue

            # Edge is the new shortest connection from the tree to w
            if edge.weight() < self._dist_to[w]:
                self._edge_to[w] = edge
                self._dist_to[w] = edge.weight()

                # If priority queue already contains w, update its priority
                # with the new distance key. Otherwise, add w to the queue
                if self._pq.contains(w):
                    self._pq.change(w, self._dist_to[w])
                else:
                    self._pq.push(w, self._dist_to[w])


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
    """ MstPrimEager test client """
    with open('../data/mediumEWG.txt') as stream:
        graph = EdgeWeightedGraph.from_stream(stream)
        mst = MstPrimEager(graph)
    print(graph)
    print("mst\n", mst)
    print("mst weight", mst.weight())
    draw(graph, mst)

    # Done
    return 0

if __name__ == '__main__':
    sys.exit(main(sys.argv))
