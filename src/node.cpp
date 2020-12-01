//
//  node.cpp
//  nori
//
//  Created by Nicolas Buisson on 04/03/2018.
//

#include <nori/node.h>
#include <stdio.h>

Node::Node() {
    for (int i = 0 ; i < 8 ; i++) {
        m_child[i] = nullptr;
    }
    m_triangles = nullptr;
}

Node::Node(std::vector<uint32_t>* triangles) {
    for (int i = 0 ; i < 8 ; i++) {
        m_child[i] = nullptr;
    }
    m_triangles = triangles;
}

Node::~Node() {
    // Call destruction for childs
    for (int i = 0 ; i < 8 ; i++) {
        delete m_child[i];
    }
    // Call destruction for stored triangles
    delete m_triangles;
}

Node* Node::getChild(int index) const {
    return m_child[index];
}

void Node::setChild(int index, Node* child) {
    m_child[index] = child;
}

std::vector<uint32_t>* Node::getTriangles() const {
    return m_triangles;
}

bool Node::isEmptyNode() const {
    // Is this a leaf node? Leaf nodes are not empty nodes!
    if (m_triangles != nullptr) {
        return false;
    }
    // Is there any child under this node? Nodes with childs are not empty nodes!
    for (int i = 0 ; i < 8 ; i++) {
        if (m_child[i] != nullptr) {
            return false;
        }
    }
    // Otherwise, this is an empty node
    return true;
}

bool Node::isLeafNode() const {
    if (m_triangles == nullptr) {
        return false;
    }
    return true;
}
