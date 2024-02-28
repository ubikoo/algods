//
// bstree.h
//
// Copyright (c) 2020 Carlos Braga
// This program is free software; you can redistribute it and/or modify it
// under the terms of the MIT License. See accompanying LICENSE.md or
// https://opensource.org/licenses/MIT.
//

#ifndef SEARCHING_BSTREE_H_
#define SEARCHING_BSTREE_H_

#include <string>
#include <sstream>
#include <vector>
#include <utility>
#include <functional>

#define BSTREE_T_DECL template<typename Key, typename Value, typename Compare>
#define BSTREE_C_DECL BSTree<Key, Value, Compare>

///
/// @brief BSTree represents a binary tree where each node has a comparable key
/// with an associated value and satisfies the constraint that the key in any
/// node (a) with children (a)->left and (a)->right, the following is invariant:
///     a.key (cmp) left.key = false
///     a.key (cmp) right.key = true
///
/// i.e., in a binary search tree the binary comparison operator is true for all
/// children right of the node and is false for all children left of the node.
///
/// The (cmp) operator is a strict weak order binary comparison operator:
///     x (cmp) x = false                                          (irreflexivity)
///     x (cmp) y = true => y (cmp) x = false                      (asymmetry)
///     x (cmp) y = true and y (cmp) z = true => x (cmp) z = true  (transitivity)
///
/// Each node contains a key, an associated value, a left link, a right link and
/// size of the subtree rooted at the node.
///
///  API BSTree<Key,Value,Compare>:
///      BSTree(compare)            create an empty tree with a comparator
///      ~BSTree()                  remove all key items and destroy the tree
///
///      size_t size() const        return the size of the tree
///      bool is_empty() const      is the binary search tree empty?
///      void clear()               clear all key items in the tree
///
///      void put(key, value)       insert a key-value item in the tree
///      const Value get(key)       get the value associated with the key
///      bool contains(key)         does the tree contain the key?
///
///      void erase_min()           remove the minimum key from the tree
///      void erase_max()           remove the maximum key from the tree
///      void erase(key)            remove the given key from the tree
///
///      const Key &min()           return the minimum key in the tree
///      const Key &max()           return the maximum key in the tree
///
///      const Key &floor(key)      return the largest item smaller than key
///      const Key &ceil(key)       return the smallest item larger then key
///      const Key &select(k)       return the(k+1)th smallest key in the tree
///      size_t rank(key)           return the number of keys smaller than key
///
///      vector<Key> keys()         return a vector with all the keys
///      vector<Key> keys(lo, hi)   return a vector with all keys in(lo,hi)
///
///      std::string graph()        return the graph of the tree
///      std::string graphviz()     return the graph of the tree
///
/// @see Algorithms 4th Edition by Robert Sedgewick, Kevin Wayne.
///
template<
    typename Key,
    typename Value,
    typename Compare = std::less<Key>>
class BSTree {
public:
    // Return the number of key-value pairs in the symbol table.
    size_t size() const {
        return size(root_);
    }

    // Does the symbol table contain any key?
    bool is_empty() const {
        return (size() == 0);
    }

    // Remove all key-value nodes from the symbol table.
    void clear() {
        clear(root_);
        root_ = nullptr;
    }

    // Inserts a given key-value pair into the symbol table. If the symbol table
    // already contains the specified key, overwrite the old value with the new
    // associated value.
    void put(const Key &key, const Value &val) {
        root_ = put(root_, key, val);
    }

    // Return the value associated with the given key.
    const Value &get(const Key &key) const {
        const Node *node = get(root_, key);
        assert(node != nullptr && "nonexistent key");
        return node->val_;
    }

    // Does the symbol table contain the given key?
    bool contains(const Key &key) const {
        const Node *node = get(root_, key);
        return (node != nullptr);
    }

    // Remove the smallest key and associated value from the symbol table.
    void erase_min() {
        assert(!is_empty() && "symbol table underflow");
        root_ = erase_min(root_);
    }

    // Remove the largest key and associated value from the symbol table.
    void erase_max() {
        assert(!is_empty() && "symbol table underflow");
        root_ = erase_max(root_);
    }

    // Remove the given key and associated value from the symbol table.
    void erase(const Key &key) {
        assert(!is_empty() && "symbol table underflow");
        if (!contains(key)) {
            return;
        }
        root_ = erase(root_, key);
    }

    // Return the minimum key from the symbol table.
    const Key &min() const {
        assert(!is_empty() && "undefined minimum on empty symbol table");
        const Node *node = min(root_);
        return node->key_;
    }

    // Return the maximum key from the symbol table.
    const Key &max() const {
        assert(!is_empty() && "undefined maximum on empty symbol table");
        const Node *node = max(root_);
        return node->key_;
    }

    // Return the largest key less than or equal to the given key.
    const Key &floor(const Key &key) const {
        assert(!is_empty() && "undefined floor on empty symbol table");
        const Node *node = floor(root_, key);
        assert(node != nullptr && "nonexistent, floor key is too small");
        return node->key_;
    }

    // Return the smallest key greater than or equal to the given key.
    const Key &ceil(const Key &key) const {
        assert(!is_empty() && "undefined floor on empty symbol table");
        const Node *node = ceil(root_, key);
        assert(node != nullptr && "nonexistent, ceil key is too large");
        return node->key_;
    }

    // Return the key in the symbol table with the given rank. The key has the
    // property that there are rank keys in the symbol table that are smaller:
    // the key is the (rank+1)th smallest key in the symbol table.
    const Key &select(const size_t k) const {
        assert(k >= 0 && k < size() && "invalid select argument");
        const Node *node = select(root_, k);
        return node->key_;
    }

    // Return the number of keys in the symbol table that are strictly less than
    // the given key.
    size_t rank(const Key &key) const {
        size_t k = rank(root_, key);
        assert(k != static_cast<size_t>(-1) && "nonexistent rank key argument");
        return k;
    }

    // Return an iterable vector of all keys in the symbol table.
    std::vector<Key> keys() const {
        return keys(min(), max());
    }

    // Return an iterable vector of all keys in the symbol table in the given
    // key range, lo <= key <= hi.
    std::vector<Key> keys(const Key &lo, const Key &hi) const {
        std::vector<Key> v;
        keys(root_, lo, hi, v);
        return v;
    }

    // Create a graph of the binary search tree.
    std::string graph() const {
        std::ostringstream ss;
        if (!is_empty()) {
            graph(root_, "", ss, false);
        }
        return ss.str();
    }

    // Create a graphviz graph of the red black tree.
    std::string graphviz() const {
        size_t nullcount = 0;
        std::ostringstream ss;
        ss << "digraph BSTree {\n";
        if (!is_empty()) {
            graphviz(root_->left_, root_->key_, nullcount, ss);
            graphviz(root_->right_, root_->key_, nullcount, ss);
        }
        ss << "}\n";
        return ss.str();
    }

    // BSTree constructor/destructor.
    explicit BSTree(const Compare &compare = Compare())
        : compare_(compare)
        , root_(nullptr) {}
    ~BSTree() { clear(); }

    // BSTree copy constructor/assignment.
    BSTree(const BSTree &other) {
        std::vector<Key> keys = other.keys();
        for (auto k : keys) {
            put(k, other.get(k));
        }
    }
    BSTree &operator=(const BSTree &other) {
        if (this == &other)
            return *this;
        std::vector<Key> keys = other.keys();
        for (auto k : keys) {
            put(k, other.get(k));
        }
        return *this;
    }

    // BSTree move constructor/assignment.
    BSTree(BSTree &&other) {
        compare_ = std::move(other.compare_);
        root_ = std::move(other.root_);
    }
    BSTree &operator=(BSTree &&other) {
        if (this == &other)
            return *this;
        compare_ = std::move(other.compare_);
        root_ = std::move(other.root_);
        return *this;
    }

private:
    // Node is a private nested class defining the nodes in the binary tree.
    // Each node contains a key, an associated value, a left link, a right link
    // and a node size.
    // The left link points to a tree with items whose keys obey the order:
    //      node.key(comp) item.key = false
    // The right link points to a tree with items whose keys obey the order:
    //      node.key(comp) item.key = true
    // The instance variable size gives the size of the of the subtree rooted at
    // this node.
    struct Node {
        // Node variables
        Key key_;           // comparable key
        Value val_;         // associated value
        size_t size_;       // #nodes in subtree rooted here
        Node *left_;        // link to left subtree
        Node *right_;       // link to right subtree

        // Node constructor/destructor
        Node(Key key, Value val, size_t size)
            : key_(key)
            , val_(val)
            , size_(size)
            , left_(nullptr)
            , right_(nullptr) {}
        ~Node() {}

        // Node copy constructor/assignment
        Node(const Node &other)
            : key_(other.key_)
            , val_(other.val_)
            , size_(other.size_)
            , left_(other.left_)
            , right_(other.right_) {}
        Node &operator=(const Node &other) {
            if (this == &other)
                return *this;
            key_ = other.key_;
            val_ = other.val_;
            size_ = other.size_;
            left_ = other.left_;
            right_ = other.right_;
            return *this;
        }

        // Node move constructor/assignment
        Node(Node &&other)
            : key_(std::move(other.key_))
            , val_(std::move(other.val_))
            , size_(std::move(other.size_))
            , left_(std::move(other.left_))
            , right_(std::move(other.right_)) {}
        Node &operator=(Node &&other) {
            if (this == &other)
                return *this;
            key_ = std::move(other.key_);
            val_ = std::move(other.val_);
            size_ = std::move(other.size_);
            left_ = std::move(other.left_);
            right_ = std::move(other.right_);
            return *this;
        }
    };

    // --- BSTree private API member functions --------------------------------
    size_t size(const Node *node) const;
    void clear(Node *node);

    Node *put(Node *node, const Key &key, const Value &val);
    const Node *get(const Node *node, const Key &key) const;

    Node *erase_min(Node *node);
    Node *erase_max(Node *node);
    Node *erase(Node *node, const Key &key);

    const Node *min(Node *node) const;
    const Node *max(Node *node) const;
    const Node *floor(const Node *node, const Key &key) const;
    const Node *ceil(const Node *node, const Key &key) const;
    const Node *select(const Node *node, const size_t k) const;
    size_t rank(const Node *node, const Key &key) const;

    void keys(
        const Node *node,
        const Key &lo,
        const Key &hi,
        std::vector<Key> &v) const;

    // ---- BSTree private debug functions ------------------------------------
    void graph(
        const Node *node,
        const std::string &prefix,
        std::ostringstream &ss,
        bool isleft) const;

    void graphviz(
        const Node *node,
        const Key &from,
        size_t &nullcount,
        std::ostringstream &ss) const;

    // ---- BSTree private member variables -----------------------------------
    Compare compare_;           // key comparator
    Node *root_;                // bst root link
};

/// ---- Binary-search tree query functions ------------------------------------
/// @brief Return the total number of key-value nodes in the tree.
///
BSTREE_T_DECL
inline size_t BSTREE_C_DECL::size(const BSTREE_C_DECL::Node *node) const
{
    if (node == nullptr) {
        return 0;
    }

    return node->size_;
}

///
/// @brief Remove all key-value nodes from the symbol table.
///
BSTREE_T_DECL
inline void BSTREE_C_DECL::clear(BSTREE_C_DECL::Node *node)
{
    if (node == nullptr) {
        return;
    }

    clear(node->left_);
    clear(node->right_);
    delete node;
}

/// ---- Binary-search tree accessor and mutator functions --------------------
/// @brief Inserts the given key-value pair into the symbol table. If the symbol
/// table already contains the specified key, overwrite the old value with the
/// new associated value.
///
BSTREE_T_DECL
inline typename BSTREE_C_DECL::Node *BSTREE_C_DECL::put(
    Node *node,
    const Key &key,
    const Value &val)
{
    // Fell out of the tree, create new node.
    if (node == nullptr) {
        return new Node(key, val, 1);
    }

    // If key is less than node->key, insert key on left subtree.
    // If node->key is less than key, insert key on right subtree.
    // Otherwise, key is equal to node->key, replace the value.
    if (compare_(key, node->key_)) {
        node->left_ = put(node->left_, key, val);
    } else if (compare_(node->key_, key)) {
        node->right_ = put(node->right_, key, val);
    } else {
        node->val_ = val;
    }

    node->size_ = 1 + size(node->left_) + size(node->right_);
    return node;
}

///
/// @brief Return the value associated with the given key.
///
BSTREE_T_DECL
inline const typename BSTREE_C_DECL::Node *BSTREE_C_DECL::get(
    const Node *node,
    const Key &key) const
{
    if (node == nullptr) {
        return nullptr;
    }

    if (compare_(key, node->key_)) {
        // key is less than node->key, get key from left subtree
        return get(node->left_, key);
    } else if (compare_(node->key_, key)) {
        // node->key is less than key, get key from right subtree
        return get(node->right_, key);
    } else {
        return node;
    }
}

/// ---- Binary-search tree deletion functions --------------------------------
/// @brief Remove the smallest key and associated value from the symbol table.
///
BSTREE_T_DECL
inline typename BSTREE_C_DECL::Node *BSTREE_C_DECL::erase_min(Node *node)
{
    if (node->left_ == nullptr) {
        Node *link = node->right_;
        delete node;
        return link;
    }

    node->left_ = erase_min(node->left_);
    node->size_ = 1 + size(node->left_) + size(node->right_);
    return node;
}

///
/// @brief Remove the largest key and associated value from the symbol table.
///
BSTREE_T_DECL
inline typename BSTREE_C_DECL::Node *BSTREE_C_DECL::erase_max(Node *node)
{
    if (node->right_== nullptr) {
        Node *link = node->left_;
        delete node;
        return link;
    }

    node->right_ = erase_max(node->right_);
    node->size_ = 1 + size(node->left_) + size(node->right_);
    return node;
}

///
/// @brief Remove the specified key and associated value from the symbol table.
/// 1) If node has no right child, just delete it and return its left child.
/// 2) If node has no left child instead, delete it and return its right child.
/// 3) If node is internal, use Hibbard's deletion - replace the node with
///    its successor. Because node contains a right child, its successor
///    is the smallest key in its right subtree.
///    The operation maintains the order of the tree because there are no
///    keys between node->key and its successor's key. So:
///      - save a link to the node to be deleted in link
///      - set node to point to its successor min(link->right_)
///      - relink node->right(pointing to a subtree with all keys larger than
///        node->key) to point to the subtree resultant from deleting the min
///        on the right subtree erase_min(link->right)
///      - relink node->left to the untouched left subtree of the saved link
///       (pointing to the subtree containing all keys smaller than node->key).
///
BSTREE_T_DECL
inline typename BSTREE_C_DECL::Node *BSTREE_C_DECL::erase(Node *node, const Key &key)
{
    // Fell out of the tree, return null node.
    if (node == nullptr) {
        return nullptr;
    }

    if (compare_(key, node->key_)) {
        // key is less than node->key, remove key from left subtree
        node->left_ = erase(node->left_, key);
    } else if (compare_(node->key_, key)) {
        // node->key is less than key, remove key from right subtree
        node->right_ = erase(node->right_, key);
    } else {
        // found node containing the key, delete it and return its successor
        if (node->right_ == nullptr) {
            Node *link = node->left_;
            delete node;
            return link;
        }

        if (node->left_ == nullptr) {
            Node *link = node->right_;
            delete node;
            return link;
        }

        const Node *link = min(node->right_);
        node->key_ = link->key_;
        node->val_ = link->val_;
        node->right_ = erase_min(node->right_);
    }

    node->size_ = 1 + size(node->left_) + size(node->right_);
    return node;
}

/// ---- Binary-search tree search functions ----------------------------------
/// @brief Return the minimum key from the symbol table.
///
BSTREE_T_DECL
inline const typename BSTREE_C_DECL::Node *BSTREE_C_DECL::min(Node *node) const
{
    if (node->left_ == nullptr) {
        return node;
    }

    return min(node->left_);
}

///
/// @brief Return the maximum key from the symbol table.
///
BSTREE_T_DECL
inline const typename BSTREE_C_DECL::Node *BSTREE_C_DECL::max(Node *node) const
{
    if (node->right_ == nullptr) {
        return node;
    }

    return max(node->right_);
}

///
/// @brief Return the largest key less than or equal to the given key.
///
BSTREE_T_DECL
inline const typename BSTREE_C_DECL::Node *BSTREE_C_DECL::floor(
    const Node *node,
    const Key &key) const
{
    if (node == nullptr) {
        return nullptr;
    }

    if (compare_(key, node->key_)) {
        // key is smaller than the node-key, floor key is on the left subtree.
        return floor(node->left_, key);
    } else if (compare_(node->key_, key)) {
        // node-key is smaller than key, try to find a floor key on the right
        // subtree, i.e., greater than node->key but smaller than key; if none
        // found, node->key is the floor.
        const Node *link = floor(node->right_, key);
        if (link != nullptr) {
            return link;
        }
    }

    return node;
}

///
/// @brief Return the smallest key greater than or equal to the given key.
///
BSTREE_T_DECL
inline const typename BSTREE_C_DECL::Node *BSTREE_C_DECL::ceil(
    const Node *node,
    const Key &key) const
{
    if (node == nullptr) {
        return nullptr;
    }

    if (compare_(key, node->key_)) {
        // key is smaller than node-key, try to find a ceiling key on the left
        // subtree, i.e., smaller than node->key but larger than key; if none
        // found, node->key is the ceiling.
        const Node *link = ceil(node->left_, key);
        if (link != nullptr) {
            return link;
        }
    } else if (compare_(node->key_, key)) {
        // node-key is smaller than key, ceiling key is on the right subtree.
        return ceil(node->right_, key);
    }

    return node;
}

///
/// @brief Return the key in the symbol table with the given rank.
/// The key has the property that there are rank keys in the symbol table that
/// are smaller - the key is the(rank+1)th smallest key in the symbol table.
///
BSTREE_T_DECL
inline const typename BSTREE_C_DECL::Node *BSTREE_C_DECL::select(
    const Node *node,
    const size_t k) const
{
    if (node == nullptr) {
        return nullptr;
    }

    size_t t = size(node->left_);
    if (k < t) {
        // rank is smaller than left subtree size, find key on left
        return select(node->left_, k);
    } else if (k > t) {
        // rank is larger than left subtree size - find key on right subtree by
        // appropriately subtracting left subtree size and this node contribution
        return select(node->right_, k - (t + 1));
    }

    return node;
}

///
/// Return the number of keys in the symbol table that are strictly less than the
/// given key.
///
BSTREE_T_DECL
inline size_t BSTREE_C_DECL::rank(const Node *node, const Key &key) const
{
    if (node == nullptr) {
        // fell out of the tree, return invalid rank.
        return -1;
    }

    if (compare_(key, node->key_)) {
        // key is less than node->key, rank-key is on the left subtree
        return rank(node->left_, key);
    } else if (compare_(node->key_, key)) {
        // node->key is less key, rank-key is on the right subtree.
        // the rank of this key is this node plus the size of the left subtree
        // plus the rank of the key on the right subtree.
        return (1 + size(node->left_) + rank(node->right_, key));
    } else {
        // node contains the key, return the size of its left subtree.
        return size(node->left_);
    }
}

///
/// @brief Return an iterable vector of all keys in the symbol table in the given
/// key range, lo <= key <= hi.
///
BSTREE_T_DECL
inline void BSTREE_C_DECL::keys(
    const Node *node,
    const Key &lo,
    const Key &hi,
    std::vector<Key> &v) const
{
    if (node == nullptr) {
        return;
    }

    // Search keys at the left subtree.
    keys(node->left_, lo, hi, v);

    // Insert key if within the range.
    if (!compare_(node->key_, lo) && !compare_(hi, node->key_)) {
        v.push_back(node->key_);
    }

    // Search keys on the right subtree.
    keys(node->right_, lo, hi, v);
}

/// ---- Binary-search tree debug functions -----------------------------------
/// @brief Create a graph of the subtree root at this node.
///
BSTREE_T_DECL
inline void BSTREE_C_DECL::graph(
    const Node *node,
    const std::string &prefix,
    std::ostringstream &ss,
    bool isleft) const
{
    if (node == nullptr) {
        return;
    }

    ss << prefix.c_str()
       << (isleft ? "└─<" : "└─>")
       << node->key_
       << "\n";
    graph(node->left_,  prefix + "   ", ss, true);
    graph(node->right_, prefix + "   ", ss, false);
}

///
/// @brief Create a graph of the subtree root at this node.
///
BSTREE_T_DECL
inline void BSTREE_C_DECL::graphviz(
    const Node *node,
    const Key &from,
    size_t &nullcount,
    std::ostringstream &ss) const
{
    if (node == nullptr) {
        std::string nullnode("null" + std::to_string(nullcount++));
        ss << "\t" << nullnode << " [shape=point];\n";
        ss << "\t" << from << " -> " << nullnode << "\n";
        return;
    }

    ss << "\t" << from << " -> " << node->key_ << ";\n";

    graphviz(node->left_, node->key_, nullcount, ss);
    graphviz(node->right_, node->key_, nullcount, ss);
}

#undef BSTREE_T_DECL
#undef BSTREE_C_DECL

#endif // SEARCHING_BSTREE_H_
