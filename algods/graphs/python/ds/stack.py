#!/usr/bin/env ipython

# stack.py
#
# Copyright (c) 2020 Carlos Braga
# This program is free software; you can redistribute it and/or modify it under
# the terms of the MIT License. See accompanying LICENSE.md or
# https://opensource.org/licenses/MIT.

import sys
import collections

# -----------------------------------------------------------------------------
class Stack:
    """
    Stack maintains a double ended queue of n items. An item is inserted at the
    end of the queue and removed from the end of the queue.
    """

    # ---- Stack special methods ----------------------------------------------
    def __init__(self):
        self._n_items = 0
        self._items = collections.deque()

    def __del__(self):
        pass

    def __iter__(self):
        ix = self._n_items - 1
        while ix >= 0:
            item = self._items[ix]
            yield item
            ix -= 1

    def __repr__(self):
        """ Return a string representation of this Stack """
        s = []
        s.append("%d items " % (self._n_items))
        for item in self:
            s.append(" %s" % repr(item))
        s.append("\n")
        return ''.join(s)

    def __len__(self):
        """ Return the number of items in the stack """
        return self.size()

    # ---- Stack API ----------------------------------------------------------
    def is_empty(self):
        """ Is the stack empty? """
        return self._n_items == 0

    def size(self):
        """ Return the size of the stack """
        return self._n_items

    def clear(self):
        """ Clear the stack """
        self._n_items = 0
        self._items.clear()

    def push(self, item):
        """ Insert an item into the stack"""
        self._items.append(item)
        self._n_items += 1

    def pop(self):
        """ Remove the top item from the stack """
        if self.is_empty():
            return None
        item = self._items.pop()
        self._n_items -= 1
        return item

    def top(self):
        """ Return the top item in the stack """
        if self.is_empty():
            return None
        return self._items[-1]


# -----------------------------------------------------------------------------
def main (argv):
    """ Stack test client """

    # Create an empty Stack
    stack = Stack()

    # Insert 10 items
    print("is_empty: ", stack.is_empty())
    for i in range(10):
        stack.push(i)
    print("is_empty: ", stack.is_empty())

    # Size and string representation
    print("len :", len(stack))
    print("repr:", repr(stack), file=sys.stdout, end='')

    # Iterate over the items
    print("iterate :", file=sys.stdout, end=' ')
    for i in stack:
        print(i, file=sys.stdout, end=' ')
    print("")

    # Delete the items
    while not stack.is_empty():
        print("top ", stack.top(), file=sys.stdout, end=' ')
        print("size ", stack.size(), file=sys.stdout, end=' ')
        print("pop ", stack.pop(), file=sys.stdout, end='\n')
    print("repr: ", repr(stack), file=sys.stdout, end=' ')

    # Done
    return 0

if __name__ == '__main__':
    sys.exit(main(sys.argv))
