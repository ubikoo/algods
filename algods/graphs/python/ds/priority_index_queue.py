#!/usr/bin/env ipython

# priority_queue_index.py
#
# Copyright (c) 2020 Carlos Braga
# This program is free software; you can redistribute it and/or modify
# it under the terms of the MIT License.
# See accompanying LICENSE.md or https://opensource.org/licenses/MIT.

import sys
import collections
import random

# -----------------------------------------------------------------------------
class PriorityIndexQueue:
    """
    PriorityIndexQueue represents an indexed heap-ordered complete binary tree.
    Each key is associated with an index between 0 and (max_keys - 1).
    The client uses the key index to change the key priority value or delete.
    The queue is implemented using a heap-ordered binary tree with an array to
    associate every key with the corresponding integers, in the given range.

    PriorityIndexQueue API:
        is_empty()          is the heap empty?
        size()              number of keys in the heap
        clear()             remove all keys from the heap

        is_valid(i)         is the specified index valid?
        contains(i)         is the specified index in the priority queue?

        push(i, Key)        insert a key into the queue at the specified index
        pop()               delete the topmost key in the queue

        top()               return the topmost key in the queue
        top_ix()            return the index of the topmost key in the queue
        key_of(i)           return the key associated the given index

        change(i, key)      change the priority key associated with index i
        delete(i)           remove key associated with index i from the queue

        swap(i, j)          swap keys associated with index i and index j
        swim(k)             move the k-key up to restore the heap order
        sink(k)             move the k-key down to restore the heap order
    """

    # ---- PriorityIndexQueue special methods --------------------------------------
    def __init__(self, capacity, compare):
        """
        Create a priority queue with a specified order and initial capacity
        """
        assert capacity > 0, "invalid priority queue capacity"
        assert compare != None, "invalid compare function"
        self._compare = compare
        self._capacity = capacity;
        self._count = 0
        self._keys = [None for _ in range(self._capacity + 1)]
        self._pq = [None for _ in range(self._capacity + 1)]
        self._qp = [-1 for _ in range(self._capacity + 1)]

    def __del__(self):
        pass

    def __iter__(self):
        ix = 0
        while ix < self._count:
            key = self.key_of(ix)
            yield key
            ix += 1

    def __repr__(self):
        """ Return a string representation of this PriorityIndexQueue """
        s = []
        s.append("%d keys:" % (self._count))
        for key in self:
            s.append(" %s" % repr(key))
        s.append("\n")
        return ''.join(s)

    def __len__(self):
        """ Return the number of keys in the queue """
        return self.size()

    # ---- PriorityIndexQueue public api -------------------------------------------
    def is_empty(self):
        """ Is the queue empty? """
        return self._count == 0

    def size(self):
        """ Return the size of the queue """
        return self._count

    def clear(self):
        """ Clear the queue """
        self._count = 0
        self._keys = [None for _ in range(self._capacity + 1)]
        self._pq = [None for _ in range(self._capacity + 1)]
        self._qp = [-1 for _ in range(self._capacity + 1)]

    def is_valid(self, i):
        """ Is the specified index valid? """
        return i >= 0 and i <= self._capacity

    def contains(self, i):
        """ Is the specified index in the priority queue? """
        assert self.is_valid(i), "index is outside priority queue range"
        return self._qp[i] != -1

    def push(self, i, key):
        """ Insert a key into the queue at the specified index  """
        assert not self.contains(i), "index is already in the priority queue"

        self._count += 1
        self._qp[i] = self._count
        self._pq[self._count] = i
        self._keys[i] = key
        self._swim(self._count)

    def pop(self):
        """ Remove the top key from the priority queue """
        assert not self.is_empty(), "priority queue underflow"
        _min = self._pq[1]              # store the minimum index

        self._swap(1, self._count)      # swap the top-key with the bottom-key
        self._count -= 1                # discard the top-key at the bottom
        self._sink(1)                   # sink the bottom-key from the top

        self._qp[_min] = -1             # remove the index from the queue
        self._keys[_min] = None         # remove the associated key
        self._pq[self._count + 1] = 0   # remove unused priority index

    def top(self):
        """ Return the topmost key in the queue """
        assert not self.is_empty(), "priority queue underflow"
        return self._keys[self._pq[1]];

    def top_ix(self):
        """ Return the index of the topmost key in the queu """
        assert not self.is_empty(), "priority queue underflow"
        return self._pq[1];

    def key_of(self, i):
        """ Return the key associated the given index """
        assert self.contains(i), "index is not in the priority queue"
        return self._keys[i];

    def change(self, i, key):
        """ Change the priority key associated with index i """
        assert self.contains(i), "index is not in the priority queue"
        self._keys[i] = key
        self._swim(self._qp[i])
        self._sink(self._qp[i])

    def delete(self, i):
        """ Remove key associated with index i from the queue """
        assert self.contains(i), "index is not in the priority queue"
        index = self._qp[i]
        self._swap(index, self._count)
        self._count -= 1
        self._swim(index)
        self._sink(index)
        self._keys[i] = None
        self._qp[i] = -1

    # ---- PriorityIndexQueue private api ------------------------------------------
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
        self._pq[i], self._pq[j] = self._pq[j], self._pq[i]
        self._qp[self._pq[i]] = i
        self._qp[self._pq[j]] = j

    def _swim(self, k):
        """ Move the k-key up to restore the heap order """
        while k > 1 and self._comp(k // 2, k):
            self._swap(k // 2, k)
            k = k // 2

    def _sink(self, k):
        """ Move the k-key down to restore the heap order """
        while (2 * k <= self._count):
            j = 2 * k           # left child of parent k
            if (j < self._count and self._comp(j, j + 1)):
                j += 1          # take right child if compare left is true

            if self._comp(j, k):
                break           # compare parent k with child j is true
            self._swap(j, k)    # swap parent with its child
            k = j               # go down the childs branch

    def _comp(self, i, j):
        """ Compare i-key with j-key """
        return self._compare(self._keys[self._pq[i]], self._keys[self._pq[j]])

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
    """ PriorityIndexQueue test client """
    n_items = 10

    # -------------------------------------------------------------------------
    # Test data
    data = [x for x in range(n_items)]
    random.shuffle(data)
    print("data", data, file=sys.stdout, end='\n')

    # Create a max PriorityIndexQueue
    queue_max = PriorityIndexQueue(n_items, (lambda x, y: x < y))

    # Insert 10 keys
    print("is_empty: ", queue_max.is_empty())
    for ix, key in enumerate(data):
        queue_max.push(ix, key)
    print("is_empty: ", queue_max.is_empty())

    # Size and string representation
    print("len :", len(queue_max))
    print("repr:", repr(queue_max), file=sys.stdout, end='')

    # Iterate over the keys
    print("iterate :", file=sys.stdout, end=' ')
    for i in queue_max:
        print(i, file=sys.stdout, end=' ')
    print("")

    # Create a min PriorityIndexQueue
    queue_min = PriorityIndexQueue(n_items, (lambda x, y: x > y))

    # Delete the keys
    print("\nqueue max: ", queue_max, file=sys.stdout, end='\n')
    while not queue_max.is_empty():
        print("size ", queue_max.size(), file=sys.stdout, end=' ')
        print("index ", queue_max.top_ix(), file=sys.stdout, end=' ')
        print("top ", queue_max.top(), file=sys.stdout, end='\n')
        ix = queue_max.top_ix()
        key = queue_max.top()
        queue_min.push(ix, key)
        queue_max.pop()

    print("\nqueue min: ", queue_min, file=sys.stdout, end='\n')
    while not queue_min.is_empty():
        print("size ", queue_min.size(), file=sys.stdout, end=' ')
        print("index ", queue_min.top_ix(), file=sys.stdout, end=' ')
        print("top ", queue_min.top(), file=sys.stdout, end='\n')
        queue_min.pop()

    # -------------------------------------------------------------------------
    # Test items data
    items = [Item(11*i) for i in range(n_items)]
    print("items", items, file=sys.stdout, end='\n')

    # Create a max items PriorityIndexQueue
    queue_items = PriorityIndexQueue(n_items, (lambda x, y: x < y))
    for ix, key in enumerate(items):
         queue_items.push(ix, key)

    print("\nqueue items: ", queue_items, file=sys.stdout, end='\n')
    while not queue_items.is_empty():
        print("size ", queue_items.size(), file=sys.stdout, end=' ')
        print("index ", queue_items.top_ix(), file=sys.stdout, end=' ')
        print("top ", queue_items.top(), file=sys.stdout, end='\n')
        queue_items.pop()

    # Done
    return 0

if __name__ == '__main__':
    sys.exit(main(sys.argv))
