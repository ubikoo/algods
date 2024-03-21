#!/usr/bin/env ipython

# union_find.py
#
# Copyright (c) 2020 Carlos Braga
# This program is free software; you can redistribute it and/or modify it under
# the terms of the MIT License. See accompanying LICENSE.md or
# https://opensource.org/licenses/MIT.

import sys
import collections

# -----------------------------------------------------------------------------
class UnionFind:
    """
    UnionFind implements a disjoint set data structure.
    It supports union and find operations on the sets, together with a count
    operation that returns the total number of components.
    This implementation uses the weighted-union by size with path compression
    during the find operation.

    UnionFind solves a dynamic connectivity problem. The input is taken as a
    sequence of pairs of integers, where each integer represents an object of
    some type.
    The pair p and q are to be interpreted as meaning p "is connected to" q.
    The term "is connected to" is assumed to be an equivalence relation:
        symmetric:  if p is connected to q, then q is connected to p.
        transitive: if p is connected to q and q is connected to r,
                    then p is connected to r.
        reflexive:  p is connected to p.

    The underlying data structure is a site indexed array which maintains the
    invariant that p and q are connected if and only if the root parent of p
    is equal to the root parent of q.
    The find operation returns the canonical element of the set containing a
    given element.
    The join operation merges the set containing element p with the set
    containing element q. Which root parent becomes child is determined by the
    rank (or size) of the sets.

    API UnionFind:
        UnionFind()     create an empty disjoint-set data structure.

        n_keys()        return the number of keys.
        n_sets()        return the number of sets.
        size(p)         return the size of the tree containing p
        contains(p)     does the union-find ensemble contain the p-key?
        connected(p, q) do p and q keys belong to the same component?
        find(p)         find the canonical component of element p
        join(p, q)      merge set containing p with set containing q

        std::map<self._sizet, std::vector<self._sizet>> sets{}
                        return a key-value map of all the components, where key
                        is the component parent index and value is a vector of
                        all the keys contained in the component.
    """

    # ---- UnionFind special methods ------------------------------------------
    def __init__(self, n_keys):
        """
        Create a collection of disjoint data sets with n keys. Each set
        contains initially a single key, it being its parent set.
        """
        assert n_keys > 0, "invalid number of keys"
        self._n_keys = n_keys
        self._n_sets = n_keys
        self._parent = [i for i in range(n_keys)] # each key is its parent set
        self._size = [1 for _ in range(n_keys)]   # with a single element

    def __del__(self):
        pass

    def __iter__(self):
        for k,v in self.sets():
            yield (k,v)

    def __repr__(self):
        """ Return a string representation of this UnionFind """
        s = []
        s.append("%d keys, %d sets\n" % (self._n_keys, self._n_sets))
        for k,v in self.sets().items():
            s.append("set id %d, size %d: " % (k, self._size[k]))
            for iv in v:
                s.append("%d " % iv)
            s.append("\n")
        return ''.join(s)

    def __len__(self):
        """ Return the number of keys in the queue """
        return self._n_keys

    # ---- UnionFind public api -----------------------------------------------
    def n_keys(self):
        """ Return the number of keys """
        return self._n_keys

    def n_sets(self):
        """ Return the number of sets """
        return self._n_sets

    def size(self, p):
        """ Return the size of the tree containing p """
        assert self.contains(p), "non existent key"
        return self._size[self.find(p)]

    def contains(self, p):
        """ Does the union-find ensemble contain the p-key? """
        return (p >= 0 and p < self._n_keys)

    def connected(self, p, q):
        """ Do p and q keys belong to the same component? """
        return (self.find(p) == self.find(q))

    def find(self, p):
        """
        Find the canonical component to which p-key belongs to.
        The find operation uses path compression - every element between the
        query value p and the root is set to to point to the root.
        Path compression can be implemented using two passes - one to find the
        root and the second to compress the path towards the root id.
        """
        assert self.contains(p), "non existent key"

        # Find the root element of p
        root = p
        while (root != self._parent[root]):
            root = self._parent[root]

        # Compress the path from p towards its root.
        key = p
        while (key != self._parent[key]):
            current = key
            key = self._parent[key]
            self._parent[current] = root

        return root

    def join(self, p, q):
        """
        Join the component of p-value and the component of q-value.
        If these are already connected, there is nothing to do.
        Otherwise, merge the canonical components of p and q.
        """

        # Get the components of p a q
        root_p = self.find(p)
        root_q = self.find(q)

        # Both p and q belong to the same component - nothing to do.
        if (root_p == root_q):
            return

        # Join the smaller set into the larger set to minimize tree depth after
        # the merge operation.
        # If the p-set is smaller than q-set, merge p-set into q-set and
        # increment q-size with p-size.
        # Otherwise, merge q-set into p-set and increment p-size.
        if (self._size[root_p] < self._size[root_q]):
            self._parent[root_p] = root_q
            self._size[root_q] += self._size[root_p]
        else:
            self._parent[root_q] = root_p
            self._size[root_p] += self._size[root_q]

        # For every merge of two components, the total number of component
        # sets is decreased by one.
        self._n_sets -= 1

    def sets(self):
        """
        Return a key-value map of all the components, where key is the parent
        index of the specified component set and value is a collection of all
        the keys belonging to that set.
        """
        # std::map<self._sizet, std::vector<self._sizet>> sets()
        components = collections.defaultdict(list)
        for p in range(self._n_keys):
            components[self.find(p)].append(p)
        return components


# -----------------------------------------------------------------------------
def main (argv):
    """ UnionFind test client """
    with open('../../data/tinyUF.txt') as stream:
        n_keys = int(stream.readline())
        assert n_keys >= 0, "invalid number of keys"

        uf = UnionFind(n_keys)

        buffer = stream.readline()
        while buffer:
            items = buffer.split()
            assert len(items) == 2, "invalid connection"

            p = int(items[0])
            q = int(items[1])
            print("join ", p, q, file=sys.stdout, end='\n')
            uf.join(p, q)

            buffer = stream.readline()
    print(uf)

    # Done
    return 0

if __name__ == '__main__':
    sys.exit(main(sys.argv))
