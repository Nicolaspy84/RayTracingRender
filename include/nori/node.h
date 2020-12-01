//
//  node.h
//  nori
//
//  Created by Nicolas Buisson on 04/03/2018.
//

#ifndef node_h
#define node_h

#include <vector>
#include <Eigen/Geometry>

class Node {
public:
    /// Default node constructor
    Node();
    /// Leaf node constructor (will setup the triangle indices vector list)
    Node(std::vector<uint32_t>* triangles);
    ~Node();
    
    Node* getChild(int index) const ;
    void setChild(int index, Node* child);
    
    /// Tests if the node is empty or not
    bool isEmptyNode() const ;
    /// Tests if the node is a leaf node or not
    bool isLeafNode() const ;
    
    std::vector<uint32_t>* getTriangles() const ;
    
private:
    /// Stores all child nodes
    Node* m_child[8];
    /// Stores potential triangles (if a node leaf)
    std::vector<uint32_t>* m_triangles;
};


#endif /* node_h */
