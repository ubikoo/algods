#!/usr/bin/env ipython

# edge_weighted_graph.py
#
# Copyright (c) 2020 Carlos Braga
# This program is free software; you can redistribute it and/or modify it under
# the terms of the MIT License. See accompanying LICENSE.md or
# https://opensource.org/licenses/MIT.

import sys
import numpy
import matplotlib
from matplotlib import pyplot
from mpl_toolkits.mplot3d import axes3d
import networkx

# -----------------------------------------------------------------------------
class Edge:
    """
    Edge represents a weighted edge between two vertices.
    """

    # ---- Edge special methods -----------------------------------------------
    def __init__(self, v, w, weight):
        self._v = v
        self._w = w
        self._weight = weight

    def __del__(self):
        pass

    def __repr__(self):
        """ Return a string representation of this graph """
        return ("%d %d %f" % (self._v, self._w, self._weight))

    def __eq__(self, other):
        return (self._weight == other._weight)

    def __ne__(self, other):
        return (self._weight != other._weight)

    def __lt__(self, other):
        return (self._weight < other._weight)

    def __le__(self, other):
        return (self._weight <= other._weight)

    def __gt__(self, other):
        return (self._weight > other._weight)

    def __ge__(self, other):
        return (self._weight >= other._weight)

    # ---- Edge API -----------------------------------------------------------
    def first(self):
        """ Return the first vertex """
        return self._v

    def second(self):
        """ Return the second vertex """
        return self._w

    def other(self, v):
        """ Return the other vertex """
        if v == self._v :
            return self._w
        if v == self._w :
            return self._v
        raise ValueError('invalid vertex')

    def weight(self):
        """ Return the edge weight """
        return self._weight

# -----------------------------------------------------------------------------
class EdgeWeightedGraph:
    """
    EdgeWeightedGraph represents a undirected graph with weighted edges in
    adjacency list representation.
    """

    # ---- EdgeWeightedGraph special methods ----------------------------------
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
        s.append("%d vertices, %d edges\n" %
            (self._n_vertices, self._n_edges))
        for v in self.vertices():
            s.append("%d : " % (v))
            for edge in self.adj(v):
                s.append("%s, " % (repr(edge)))
            s.append("\n")
        return ''.join(s)

    # ---- EdgeWeightedGraph API ----------------------------------------------
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
        """ Return the degree of the specified vertex """
        assert self.is_valid(v), "invalid vertex"
        return len(self._adj[v])

    def vertices(self):
        """ Iterate over the graph vertices """
        for v in range(self._n_vertices):
            yield v

    def edges(self):
        """ Iterate over the graph edges """
        edges = []
        for v in self.vertices():
            for edge in self.adj(v):
                if v < edge.other(v):
                    edges.append(edge)
        return edges

    def adj(self, v):
        """ Iterate over the edges of the adjacency list in lifo order """
        assert self.is_valid(v), "invalid vertex"
        for edge in reversed(self._adj[v]):
            yield edge

    def add_edge(self, edge):
        """ Add a new edge between the specified vertices. """
        assert self.is_valid(edge.first()) and self.is_valid(edge.second()), \
            "invalid edge vertices"
        v = edge.first()
        w = edge.second()
        self._adj[v].append(edge)
        self._adj[w].append(edge)
        self._n_edges += 1

    # ---- EdgeWeightedGraph factory methods ----------------------------------
    @staticmethod
    def from_stream(stream):
        """ Factory method. Create a graph from an input stream. """
        n_vertices = int(stream.readline())
        assert n_vertices >= 0, "invalid number of vertices"

        n_edges = int(stream.readline())
        assert n_edges >= 0, "invalid number of edges"

        graph = EdgeWeightedGraph(n_vertices)
        for _ in range(n_edges):
            buffer = stream.readline().split()
            assert len(buffer) == 3, "invalid edge"

            v = int(buffer[0])
            w = int(buffer[1])
            weight = float(buffer[2])
            edge = Edge(v, w, weight)
            # print(
            #     "add",
            #     edge.first(),
            #     edge.second(),
            #     edge.weight(),
            #     file=sys.stdout, end='\n');
            graph.add_edge(edge)
        return graph

    @staticmethod
    def from_graph(other):
        """ Factory method. Create a graph from another graph. """
        # Create a graph with as many vertices as the other graph.
        graph = EdgeWeightedGraph(other.n_vertices())

        # Copy the adjacency lists from the other graph.
        graph._n_edges = other.n_edges()
        for v in other.vertices():
            graph._adj[v] = [edge for edge in other.adj(v)]
        return graph


# -----------------------------------------------------------------------------
def draw(graph):
    """ Draw the graph using networkx package """
    # Create the networkx graph
    G = networkx.Graph()
    for edge in graph.edges():
        G.add_edge(edge.first(), edge.second())

    # Draw the networkx graph
    pyplot.figure(figsize=(7,6))
    options = {
        "font_size": 10,
        "font_color": "#000000",
        "with_labels": False,
        "node_size": 10,
        "node_shape": "o",
        "node_color": "#dddddd",
        "edge_color": "#eeeeee",
        "linewidths": 1,
        "width": 1,
        "alpha" : 0.2,
    }
    pos = networkx.spring_layout(G, k=2.0, iterations=1000, seed=63)
    networkx.draw_networkx(G, pos, **options)
    pyplot.axis("off")
    pyplot.show()
    # pyplot.savefig('fig.pdf', bbox_inches='tight')


# -----------------------------------------------------------------------------
def main (argv):
    """ Graph test client """

    # Create a graph from stream
    with open('../../data/mediumEWG.txt') as stream:
        graph1 = EdgeWeightedGraph.from_stream(stream)
    print(graph1)
    draw(graph1)

    # Create a graph from another graph
    graph2 = EdgeWeightedGraph.from_graph(graph1)
    print(graph2)
    draw(graph2)

    # Done
    return 0

if __name__ == '__main__':
    sys.exit(main(sys.argv))
