#!/usr/bin/env ipython

# depth_first_paths.py
#
# Copyright (c) 2020 Carlos Braga
# This program is free software; you can redistribute it and/or modify
# it under the terms of the MIT License.
# See accompanying LICENSE.md or https://opensource.org/licenses/MIT.

import sys
import networkx
import numpy
from matplotlib import pyplot
from ds.graph import Graph
from ds.queue import Queue
from ds.stack import Stack

# -----------------------------------------------------------------------------
class DepthFirstPaths:
    """
    Compute the minimum distance paths between the set of reachable vertices and
    the source vertex using depth first search algorithm.
    """

    # ---- Special methods ----------------------------------------------------
    def __init__(self, graph, source):
        self.NEW = 0
        self.VISITED = 1
        self._source = source
        self._state = [self.NEW for _ in graph.vertices()]
        self._edge_to = [0 for _ in graph.vertices()]
        self._dfs(graph, source)

    def __del__(self):
        pass

    # ---- API ----------------------------------------------------------------
    def is_visited(self, v):
        """ Has the specified vertex been visited? """
        return self._state[v] == self.VISITED

    def path_to(self, v):
        """ Return the path from the source vertex to the specified vertex """
        if not self.is_visited(v):
            return None
        R = Stack()
        x = v
        while x != self._source:
            R.push(x)               # Add this vertex to the path
            x = self._edge_to[x]    # Move up the graph
        R.push(self._source)        # Add source vertex to complete the path
        return R

    def _dfs(self, graph, s):
        """ Compute the vertices reachable from a source vertex s """
        S = Stack()
        S.push(s)
        while not S.is_empty():
            v = S.pop()
            self._state[v] = self.VISITED
            for w in graph.adj(v):
                if self._state[w] == self.NEW:
                    self._edge_to[w] = v
                    S.push(w)


# -----------------------------------------------------------------------------
def draw(dfs, graph):
    """ Draw the graph using networkx package """
    # Create the networkx graph
    G = networkx.Graph()
    for v in graph.vertices():
        if dfs.is_visited(v):
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
    """ DepthFirstPaths test client """

    # Create a graph from stream
    with open('../data/tinyG.txt') as stream:
        graph = Graph.from_stream(stream)
    print(graph)

    # Compute the list of vertices reachable from v
    # for v in graph.vertices()
    #     dfs = DepthFirstPaths(graph, v)
    #     draw(dfs, graph)
    #     print("\nconnected to ", v, ": ", file=sys.stdout, end='\n')
    #     for w in graph.vertices():
    #         if dfs.is_visited(w):
    #             print(v, "->", w, ":", dfs.path_to(w), file=sys.stdout, end='')

    v = 12
    dfs = DepthFirstPaths(graph, v)
    draw(dfs, graph)
    print("\nconnected to ", v, ": ", file=sys.stdout, end='\n')
    for w in graph.vertices():
        if dfs.is_visited(w):
            print(v, "->", w, ":", dfs.path_to(w), file=sys.stdout, end='')

    # Done
    return 0

if __name__ == '__main__':
    sys.exit(main(sys.argv))
