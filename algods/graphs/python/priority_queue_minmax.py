#!/usr/bin/env ipython

# priority_queue.py
#
# Copyright (c) 2020 Carlos Braga
# This program is free software; you can redistribute it and/or modify it under
# the terms of the MIT License. See accompanying LICENSE.md or
# https://opensource.org/licenses/MIT.

import sys
import collections
import random

# -----------------------------------------------------------------------------
class PriorityQueue:
    """
    PriorityQueue represents a heap-ordered complete binary tree.
    A binary tree is said to be heap-ordered if it maintains the invariant that
    for each node (a) with two children (a->b) and (a->c) the following is true:
         a (op) b = true
         a (op) c = true
    where (op) is a strict weak order binary comparison operator:
         x (op) x = false                           (irreflexivity)
         if x (op) y = true, then y (op) x = false  (asymmetry)
         if x (op) y = true, and
            y (op) z = true, then
            x (op) z = true                         (transitivity)

    The heap invariant is maintained when each node compares true against both
    of its children. The functions swim and sink restore the heap invariant by
    swapping the node with its parent or swapping the node with the extremum of
    its children.

    PriorityQueue API:
        is_empty()         is the heap empty?
        size()             number of keys in the heap
        clear()            remove all keys from the heap

        push(Key)          insert a key into the heap and return its index
        pop()              delete the topmost key from the heap
        top()              return the topmost key in the heap

        resize(capacity)   resize the heap underlying binary tree
        swap(i, j)         swap keys associated with index i and index j
        swim(k)            move the k-key up to restore the heap order
        sink(k)            move the k-key down to restore the heap order
    """

    # ---- PriorityQueue special methods --------------------------------------
    def __init__(self, order, capacity = 1):
        """
        Create a priority queue with a specified order and initial capacity
        """
        assert order == "max" or order == "min", "invalid priority queue order"
        assert capacity > 0, "invalid priority queue capacity"

        self._compare = self._less if order == "max" else self._greater
        self._capacity = capacity;
        self._count = 0
        self._pq = [None for _ in range(self._capacity + 1)]

    def __del__(self):
        pass

    def __iter__(self):
        ix = 1
        while ix <= self._count:
            key = self._pq[ix]
            yield key
            ix += 1

    def __repr__(self):
        """ Return a string representation of this PriorityQueue """
        s = []
        s.append("%d keys:" % (self._count))
        for key in self:
            s.append(" %s" % repr(key))
        s.append("\n")
        return ''.join(s)

    def __len__(self):
        """ Return the number of keys in the queue """
        return self.size()

    # ---- PriorityQueue public api -------------------------------------------
    def is_empty(self):
        """ Is the queue empty? """
        return self._count == 0

    def size(self):
        """ Return the size of the queue """
        return self._count

    def clear(self):
        """ Clear the queue """
        self._count = 0
        self._pq = [None for _ in range(self._capacity + 1)]

    def push(self, key):
        """ Insert an key into the queue """

        # Grow the priority queue if we reach its capacity
        if self._count == self._capacity - 1:
            self._resize(2 * self._capacity)

        # Add new key to the priority queue
        self._count += 1
        self._pq[self._count] = key
        self._swim(self._count)

    def pop(self):
        """ Remove the top key from the queue """
        assert not self.is_empty(), "priority queue underflow"

        # Remove the top key from the priority queue
        self._swap(1, self._count)  # swap the top-key with the bottom-key
        self._count -= 1            # discard the top-key at the bottom
        self._sink(1)               # sink the bottom-key from the top

        # Shrink the priority queue if its a quarter size of its capacity
        if (self._count > 0 and self._count == (self._capacity - 1) // 4):
            self._resize(self._capacity // 2)

    def top(self):
        """ Return the top key in the queue """
        assert not self.is_empty(), "priority queue underflow"
        return self._pq[1];

    # ---- PriorityQueue private api ------------------------------------------
    def _resize(self, capacity):
        """ Resize the priority queue with a new capacity """

        # Create a copy of the heap arrays with the new capacity.
        self._capacity = capacity
        copy = [None for _ in range(self._capacity)]

        # Copy the keys of the current heap.
        for i in range(1, self._count + 1):
            copy[i] = self._pq[i]
        self._pq = copy

    def _swap(self, i, j):
        """ Swap keys associated with index i and index j """
        self._pq[j], self._pq[i] = self._pq[i], self._pq[j]

    def _swim(self, k):
        """ Move the k-key up to restore the heap order """
        while k > 1 and self._compare(k // 2, k):
            self._swap(k // 2, k)
            k = k // 2

    def _sink(self, k):
        """ Move the k-key down to restore the heap order """
        while (2 * k <= self._count):
            j = 2 * k           # left child of parent k
            if (j < self._count and self._compare(j, j + 1)):
                j += 1          # take right child if compare left is true

            if self._compare(j, k):
                break           # compare parent k with child j is true
            self._swap(j, k)    # swap parent with its child
            k = j               # go down the childs branch

    def _less(self, i, j):
        """ Is i-key less than j-key? """
        return self._pq[i] < self._pq[j]

    def _greater(self, i, j):
        """ Is i-key greater than j-key? """
        return self._pq[i] > self._pq[j]


# -----------------------------------------------------------------------------
class Item:
    """ Item test class """
    def __init__(self, v):
        self.v = v

    def __del__(self):
        pass

    def __repr__(self):
        return "%d" % self.v

    def __eq__(self, other):
        return (self.v == other.v)

    def __ne__(self, other):
        return (self.v != other.v)

    def __lt__(self, other):
        return (self.v < other.v)

    def __le__(self, other):
        return (self.v <= other.v)

    def __gt__(self, other):
        return (self.v > other.v)

    def __ge__(self, other):
        return (self.v >= other.v)


# -----------------------------------------------------------------------------
def main (argv):
    """ PriorityQueue test client """
    n_items = 10

    # -------------------------------------------------------------------------
    # Test data
    data = [x for x in range(n_items)]
    random.shuffle(data)
    print("data", data, file=sys.stdout, end='\n')

    # Create a max PriorityQueue
    queue_max = PriorityQueue("max")

    # Insert 10 keys
    print("is_empty: ", queue_max.is_empty())
    for x in data:
        queue_max.push(x)
    print("is_empty: ", queue_max.is_empty())

    # Size and string representation
    print("len :", len(queue_max))
    print("repr:", repr(queue_max), file=sys.stdout, end='')

    # Iterate over the keys
    print("iterate :", file=sys.stdout, end=' ')
    for i in queue_max:
        print(i, file=sys.stdout, end=' ')
    print("")

    # Create a min PriorityQueue
    queue_min = PriorityQueue("min")

    # Delete the keys
    print("\nqueue max: ", queue_max, file=sys.stdout, end='\n')
    while not queue_max.is_empty():
        print("size ", queue_max.size(), file=sys.stdout, end=' ')
        print("top ", queue_max.top(), file=sys.stdout, end='\n')
        queue_min.push(queue_max.top())
        queue_max.pop()

    print("\nqueue min: ", queue_min, file=sys.stdout, end='\n')
    while not queue_min.is_empty():
        print("size ", queue_min.size(), file=sys.stdout, end=' ')
        print("top ", queue_min.top(), file=sys.stdout, end='\n')
        queue_min.pop()

    # -------------------------------------------------------------------------
    # Test items data
    items = [Item(11*i) for i in range(n_items)]
    print("items", items, file=sys.stdout, end='\n')

    # Create a max items PriorityQueue
    queue_items = PriorityQueue("max")
    for it in items:
        queue_items.push(it)

    print("\nqueue items: ", queue_items, file=sys.stdout, end='\n')
    while not queue_items.is_empty():
        print("size ", queue_items.size(), file=sys.stdout, end=' ')
        print("top ", queue_items.top(), file=sys.stdout, end='\n')
        queue_items.pop()

    # Done
    return 0

if __name__ == '__main__':
    sys.exit(main(sys.argv))
