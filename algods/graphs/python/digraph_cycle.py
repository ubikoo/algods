#!/usr/bin/env ipython

# depth_first_cycle.py
#
# Copyright (c) 2020 Carlos Braga
# This program is free software; you can redistribute it and/or modify
# it under the terms of the MIT License.
# See accompanying LICENSE.md or https://opensource.org/licenses/MIT.

import sys
import networkx
import numpy
from matplotlib import pyplot
from ds.digraph import DiGraph
from ds.stack import Stack

# -----------------------------------------------------------------------------
class DiGraphCycle:
    """
    Given a graph G(V,E) in adjacency list representation and a source vertex s,
    DiGraphCycle determines if G(V, E) has a cycle using depth first search.
    """

    # ---- Special methods ----------------------------------------------------
    def __init__(self, graph):
        self.NEW = 0
        self.ACTIVE = 1
        self.VISITED = 2

        self._state = [self.NEW for _ in graph.vertices()]
        self._edge_to = [-1 for _ in graph.vertices()]
        self._cycle = Stack()

        for s in graph.vertices():
            if self._state[s] == self.NEW:
                self._dfs(graph, s)

    def __del__(self):
        pass

    # ---- API ----------------------------------------------------------------
    def _dfs(self, graph, v):
        """ Run dfs from the specified vertex and check for cycles """
        # Mark the site search state as active
        self._state[v] = self.ACTIVE

        for w in graph.adj(v):
            # Short circuit the search if cycle already found
            if self.has_cycle():
                return

            # If the site is new continue the search. Otherwise if the site
            # state active, then we found a site in the path to the current
            # site that has not returned yet - a cycle.
            if self._state[w] == self.NEW:
                self._edge_to[w] = v
                self._dfs(graph, w)
            elif self._state[w] == self.ACTIVE:
                x = v
                while x != w:
                    self._cycle.push(x)
                    x = self._edge_to[x]
                self._cycle.push(w)
                self._cycle.push(v)

        # Mark site search state as visite.
        self._state[v] = self.VISITED

    def has_cycle(self):
        """ Does the graph have a cycle? """
        return not self._cycle.is_empty()

    def cycle(self):
        """ Return the cycle if found, return None otherwise """
        if not self.has_cycle():
            return None
        return self._cycle

# -----------------------------------------------------------------------------
def draw(dfs, graph):
    """ Draw the graph using networkx package """
    # Create the networkx graph
    if not dfs.has_cycle():
        return

    G = networkx.DiGraph()
    cycle = [v for v in dfs.cycle()]
    for i in range(1, len(cycle)):
        G.add_edge(cycle[i-1], cycle[i])

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
    """ DiGraphCycle test client """

    # Create a graph from stream
    with open('../data/tinyDG.txt') as stream:
        graph = DiGraph.from_stream(stream)
    print(graph)

    dfs = DiGraphCycle(graph)
    print("has_cycle", dfs.has_cycle(), file=sys.stdout, end='\n')
    if dfs.has_cycle():
        print([w for w in dfs.cycle()], file=sys.stdout, end='\n')
    draw(dfs, graph)

    # Done
    return 0

if __name__ == '__main__':
    sys.exit(main(sys.argv))
