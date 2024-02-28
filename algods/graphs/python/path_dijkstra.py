#!/usr/bin/env ipython

# path_dijkstra.py
#
# Copyright (c) 2020 Carlos Braga
# This program is free software; you can redistribute it and/or modify
# it under the terms of the MIT License.
# See accompanying LICENSE.md or https://opensource.org/licenses/MIT.

from multiprocessing import spawn
import sys
import collections
import random
import numpy
import matplotlib
from matplotlib import pyplot
from mpl_toolkits.mplot3d import axes3d
import networkx

from ds.edge_weighted_digraph import DiEdge, EdgeWeightedDiGraph
from ds.priority_queue import PriorityQueue
from ds.stack import Stack

# -----------------------------------------------------------------------------
class PathDijkstra:
    """
    PathDijkstra implements the lazy version of Dijkstra algorithm to compute
    single source shortest path algorithm in edge weighted directed graph.
    """
    class Pair:
        """ Class doc """
        def __init__(self, edge, weight):
            self.edge = edge
            self.weight = weight

        def __del__(self):
            pass

        def __repr__(self):
            """ Return a string representation of the pair """
            return ("%s %f" % (repr(self.edge), self.weight))

        def __eq__(self, other):
            return (self.weight == other.weight)

        def __ne__(self, other):
            return (self.weight != other.weight)

        def __lt__(self, other):
            return (self.weight < other.weight)

        def __le__(self, other):
            return (self.weight <= other.weight)

        def __gt__(self, other):
            return (self.weight > other.weight)

        def __ge__(self, other):
            return (self.weight >= other.weight)


    # ---- PathDijkstra special methods ---------------------------------------
    def __init__(self, graph, src):
        # Assert graph edges have non negative weights.
        for edge in graph.edges():
            assert edge.weight() >= 0, "edge has negative weight"

        # Array of vertex states
        self.NEW = 0
        self.VISITED = 1
        self._state = [self.NEW for _ in graph.vertices()]
        self._dist_to = [numpy.inf for _ in graph.vertices()]
        self._edge_to = [None for _ in graph.vertices()]

        # Priority queue of crossing edges in ascending order of path distance
        self._pq = PriorityQueue(graph.n_vertices(), (lambda x, y : x > y))

        # Relax the source vertex. Keep relaxing vertices while there are
        # crossing edges between the shortest path tree and the graph.
        # In the lazy approach, for a given crossing edge in the queue, vertex
        # w might have been relaxed already and in the tree, so we ignore it.
        self._dist_to[src] = 0.0
        self._relax(graph, src);

        # Keep relax
        while not self._pq.is_empty():
            pair = self._pq.top()
            self._pq.pop()

            edge = pair.edge
            v = edge.from_vertex()
            w = edge.to_vertex()

            if not self._is_visited(w):
                self._relax(graph, w)

    def __del__(self):
        pass

    def __iter__(self):
        for edge in self.edges():
            yield edge

    def __repr__(self):
        """ Return a string representation of this PathDijkstra """
        s = []
        for edge in self.edges():
            s.append(
                "%d %d %lf\n" % (edge.first(), edge.second(), edge.weight()))
        return ''.join(s)

    def __len__(self):
        """ Return the number of keys in the queue """
        return len(self._mst)

    # ---- PathDijkstra public api --------------------------------------------
    def dist_to(self, v):
        """
        Return the length of the shortest path from the source vertex to v
        """
        assert self._is_valid(v), "invalid vertex"
        return self._dist_to[v]

    def has_path_to(self, v):
        """ Is there a path from the source to the specified vertex? """
        assert self._is_valid(v), "invalid vertex"
        return self._dist_to[v] < float("inf")

    def path_to(self, v):
        """Returns a shortest path from the source vertex s to vertex v. """
        assert self._is_valid(v), "invalid vertex"
        if not self.has_path_to(v):
            return None
        path = Stack()
        e = self._edge_to[v]
        while e is not None:
            path.push(e)
            e = self._edge_to[e.from_vertex()]
        return path

    # ---- PathDijkstra private api -------------------------------------------
    def _is_valid(self, v):
        """ Is the speficified vertex a valid one? """
        return v >= 0 and v < len(self._dist_to)

    def _is_visited(self, v):
        """ Has the specified vertex been visited? """
        return self._state[v] == self.VISITED

    def _relax(self, graph, v):
        """
        Relax the specified vertex in the graph.
        """
        self._state[v] = self.VISITED
        for edge in graph.adj(v):
            w = edge.to_vertex()
            if self._dist_to[w] > self._dist_to[v] + edge.weight():
                self._dist_to[w] = self._dist_to[v] + edge.weight()
                self._edge_to[w] = edge
                self._pq.push(self.Pair(edge, self._dist_to[w]))


# -----------------------------------------------------------------------------
def draw(graph, sp):
    """ Draw the graph using networkx package """
    # Create the networkx graph
    G = networkx.DiGraph()
    for edge in graph.edges():
        G.add_edge(edge.from_vertex(), edge.to_vertex())

    # Create the minimum spanning tree
    T = networkx.DiGraph()
    for v in graph.vertices():
        if sp.has_path_to(v):
            for edge in sp.path_to(v):
                T.add_edge(edge.from_vertex(), edge.to_vertex())

    # Draw the networkx graph
    pyplot.figure(figsize=(7,6))
    Gopts = {
        "font_size": 10,
        "font_color": "#000000",
        "with_labels": False,
        "node_size": 50,
        "node_shape": "o",
        "node_color": "#dddddd",
        "edge_color": "#dddddd",
        "linewidths": 1,
        "width": 1,
        "alpha" : 0.2,
    }
    Topts = {
        "font_size": 10,
        "font_color": "#000000",
        "with_labels": True,
        "node_size": 100,
        "node_shape": "o",
        "node_color": "#ff0000",
        "edge_color": "#000000",
        "linewidths": 1,
        "width": 1,
        "alpha" : 0.6,
    }

    pos = networkx.spring_layout(T, k=2.0, iterations=10, seed=63)
    networkx.draw_networkx(G, pos, **Gopts)
    networkx.draw_networkx(T, pos, **Topts)

    pyplot.axis("off")
    pyplot.show()
    # pyplot.savefig('fig.pdf', bbox_inches='tight')

# -----------------------------------------------------------------------------
def main (argv):
    """ PathDijkstra test client """
    with open('../data/tinyEWD.txt') as stream:
        # Compute shortest paths
        graph = EdgeWeightedDiGraph.from_stream(stream)
        sp = PathDijkstra(graph, 0)

        # Print shortest paths
        for v in graph.vertices():
            if sp.has_path_to(v):
                print("%d to %d (%.2f):" % (0, v, sp.dist_to(v)), end=' ')
                for edge in sp.path_to(v):
                    print(edge, end=' ')
                print(end='\n')
            else:
                print("%d to %d no path", 0, v, end='\n')
        draw(graph, sp)
    # Done
    return 0

if __name__ == '__main__':
    sys.exit(main(sys.argv))
