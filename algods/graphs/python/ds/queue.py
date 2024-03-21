#!/usr/bin/env ipython

# queue.py
#
# Copyright (c) 2020 Carlos Braga
# This program is free software; you can redistribute it and/or modify it under
# the terms of the MIT License. See accompanying LICENSE.md or
# https://opensource.org/licenses/MIT.

import sys
import collections

# -----------------------------------------------------------------------------
class Queue:
    """
    Queue maintains a double ended queue of n items. An item is inserted at the
    end of the queue and removed from the front of the queue.
    """

    # ---- Queue special methods ----------------------------------------------
    def __init__(self):
        self._n_items = 0
        self._items = collections.deque()

    def __del__(self):
        pass

    def __iter__(self):
        ix = 0
        while ix < self._n_items:
            item = self._items[ix]
            yield item
            ix += 1

    def __repr__(self):
        """ Return a string representation of this Queue """
        s = []
        s.append("%d items " % (self._n_items))
        for item in self:
            s.append(" %s" % repr(item))
        s.append("\n")
        return ''.join(s)

    def __len__(self):
        """ Return the number of items in the queue """
        return self.size()

    # ---- Queue API ----------------------------------------------------------
    def is_empty(self):
        """ Is the queue empty? """
        return self._n_items == 0

    def size(self):
        """ Return the size of the queue """
        return self._n_items

    def clear(self):
        """ Clear the queue """
        self._n_items = 0
        self._items.clear()

    def push(self, item):
        """ Insert an item into the queue """
        self._items.append(item)
        self._n_items += 1

    def pop(self):
        """ Remove the top item from the queue """
        if self.is_empty():
            return None
        item = self._items.popleft()
        self._n_items -= 1
        return item

    def top(self):
        """ Return the top item in the queue """
        if self.is_empty():
            return None
        return self._items[0]


# -----------------------------------------------------------------------------
def main (argv):
    """ Queue test client """

    # Create an empty Queue
    queue = Queue()

    # Insert 10 items
    print("is_empty: ", queue.is_empty())
    for i in range(10):
        queue.push(i)
    print("is_empty: ", queue.is_empty())

    # Size and string representation
    print("len :", len(queue))
    print("repr:", repr(queue), file=sys.stdout, end='')

    # Iterate over the items
    print("iterate :", file=sys.stdout, end=' ')
    for i in queue:
        print(i, file=sys.stdout, end=' ')
    print("")

    # Delete the items
    while not queue.is_empty():
        print("top ", queue.top(), file=sys.stdout, end=' ')
        print("size ", queue.size(), file=sys.stdout, end=' ')
        print("pop ", queue.pop(), file=sys.stdout, end='\n')
    print("repr: ", repr(queue), file=sys.stdout, end=' ')

    # Done
    return 0

if __name__ == '__main__':
    sys.exit(main(sys.argv))
