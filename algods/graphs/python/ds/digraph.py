#!/usr/bin/env ipython

# digraph.py
#
# Copyright (c) 2020 Carlos Braga
# This program is free software; you can redistribute it and/or modify
# it under the terms of the MIT License.
# See accompanying LICENSE.md or https://opensource.org/licenses/MIT.

import sys
import numpy
import matplotlib
from matplotlib import pyplot
from mpl_toolkits.mplot3d import axes3d
import networkx

# -----------------------------------------------------------------------------
class DiGraph:
    """
    DiGraph represents a directed graph in adjacency list representation.
    """

    # ---- DiGraph special methods --------------------------------------------
    def __init__(self, n_vertices = 0):
        assert n_vertices >= 0, "invalid number of vertices"
        self._n_vertices = n_vertices
        self._n_edges = 0
        self._adj = []
        for _ in self.vertices():
            self._adj.append([])     # empty adjacency list for each vertex

    def __del__(self):
        pass

    def __repr__(self):
        """ Return a string representation of this graph """
        s = []
        s.append("%d vertices, %d edges\n" % (self._n_vertices, self._n_edges))
        for v in self.vertices():
            s.append("%d : " % (v))
            for w in self.adj(v):
                s.append("%d " % (w))
            s.append("\n")
        return ''.join(s)

    # ---- DiGraph API --------------------------------------------------------
    def is_valid(self, v):
        """ Is the specified vertex valid? """
        return v >= 0 and v < self._n_vertices

    def n_vertices(self):
        """ Return the number of vertices """
        return self._n_vertices

    def n_edges(self):
        """ Return the number of edges """
        return self._n_edges

    def degree(self, v):
        """ Return the outdegree of the specified vertex """
        assert self.is_valid(v), "invalid vertex"
        return len(self._adj[v])

    def vertices(self):
        """ Iterate over the vertices """
        for v in range(self._n_vertices):
            yield v

    def adj(self, v):
        """ Iterate over adjacency list of the specified vertex """
        assert self.is_valid(v), "invalid vertex"
        for w in reversed(self._adj[v]):
            yield w

    def add_edge(self, v, w):
        """ Add a new edge between the specified vertices. """
        assert self.is_valid(v) and self.is_valid(w), "invalid edge vertices"
        self._adj[v].append(w)
        self._n_edges += 1

    # ---- DiGraph factory methods --------------------------------------------
    @staticmethod
    def from_stream(stream):
        """ Factory method. Create a graph from an input stream. """
        n_vertices = int(stream.readline())
        assert n_vertices >= 0, "invalid number of vertices"

        n_edges = int(stream.readline())
        assert n_edges >= 0, "invalid number of edges"

        graph = DiGraph(n_vertices)
        for _ in range(n_edges):
            buffer = stream.readline().split()
            assert len(buffer) == 2, "invalid edge"

            v = int(buffer[0])
            w = int(buffer[1])
            # print("add", v, w, file=sys.stdout, end='\n');
            graph.add_edge(v, w)    # edge v -> w

        return graph

    @staticmethod
    def from_graph(other):
        """ Factory method. Create a graph from another graph. """
        # Create a graph with as many vertices as the other graph.
        graph = DiGraph(other.n_vertices())

        # Copy the adjacency lists from the other graph.
        graph._n_edges = other.n_edges()
        for v in other.vertices():
            graph._adj[v] = [w for w in other.adj(v)]
        return graph

    @staticmethod
    def reverse(other):
        """ Factory method. Return the reverse of another graph. """
        # Create a graph with as many vertices as the other graph.
        graph = DiGraph(other.n_vertices())

        # Reverse the adjacency lists of the other graph.
        for v in other.vertices():
            for w in other.adj(v):
                graph.add_edge(w, v)    # edge w -> v
        return graph


# -----------------------------------------------------------------------------
def draw(self):
    """ Draw the graph using networkx package """
    # Create the networkx graph
    G = networkx.DiGraph()
    for v in self.vertices():
        for w in self.adj(v):
            G.add_edge(v, w)

    # Draw the networkx graph
    pyplot.figure(figsize=(7,6))
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
    """ DiGraph test client """

    # Create a graph from stream
    with open('../../data/mediumDG.txt') as stream:
        graph1 = DiGraph.from_stream(stream)
    print(graph1)
    draw(graph1)

    # Create a graph from another graph
    graph2 = DiGraph.from_graph(graph1)
    print(graph2)
    draw(graph2)

    # Create a reverse graph from another graph
    graph3 = DiGraph.reverse(graph1)
    print(graph3)
    draw(graph3)

    # Done
    return 0

if __name__ == '__main__':
    sys.exit(main(sys.argv))
