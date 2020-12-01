/*
 This file is part of Nori, a simple educational ray tracer
 
 Copyright (c) 2015 by Wenzel Jakob
 
 Nori is free software; you can redistribute it and/or modify
 it under the terms of the GNU General Public License Version 3
 as published by the Free Software Foundation.
 
 Nori is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 GNU General Public License for more details.
 
 You should have received a copy of the GNU General Public License
 along with this program. If not, see <http://www.gnu.org/licenses/>.
 */

#include <nori/customAccel.h>
#include <nori/timer.h>
#include <Eigen/Geometry>
#include <algorithm>
#include <array>

NORI_NAMESPACE_BEGIN

CustomAccel::~CustomAccel() {
    delete m_mainNode;
}

void CustomAccel::addMesh(Mesh *mesh) {
    if (m_mesh)
        throw NoriException("Accel: only a single mesh is supported!");
    m_mesh = mesh;
    m_bbox = m_mesh->getBoundingBox();
}

void CustomAccel::build() {
    
    std::cout << "Starting acceleration structure build (octree build)" << std::endl;
    
    Timer timer;
    
    /*
     * triangles stores the indices of all the triangles in the current mesh
     * No need to release memory in this function: it will be done when calling buildNode()
     */
    std::vector<uint32_t>* triangles = new std::vector<uint32_t>();
    
    for (uint32_t idx = 0; idx < m_mesh->getTriangleCount(); ++idx) {
        triangles->push_back(idx);
    }
    
    m_mainNode = buildNode(m_bbox, triangles, 0);
    
    std::cout << getStats() << "Acceleration structure build done (took " << timer.elapsedString() << ")" << std::endl;
}

Node* CustomAccel::buildNode(const BoundingBox3f &box, std::vector<uint32_t>* allTriangles, uint32_t numCalls) {
    
    // First check if there's any possible overlapping at all
    if (allTriangles->size() <= 0) {
        // Here we need to make sure that we release the memory of the empty vector
        delete allTriangles;
        return nullptr;
    }
    
    m_nbNodes++;
    
    // Returns a single leaf node with all the remaining triangles (stops there if we are too deep in the tree)
    if (allTriangles->size() <= 10 or numCalls >= 10) {
        // We don't need to worry about allTriangles anymore: the Node class will take care of its destruction
        Node *leafNode = new Node(allTriangles);
        m_nbStoredTriangles += allTriangles->size();
        m_nbLeafs++;
        // Allocating memory for 9 pointers (8 childs and a vector) as well as the remaining triangles
        m_memoryAllocated += sizeof(Node*) * 9 + sizeof(uint32_t) * allTriangles->size();
        return leafNode;
    }
    
    // Create a 8 vectors to store subnode triangles
    std::vector<uint32_t>* triangles[8];
    for (int i = 0 ; i < 8 ; i++) {
        // Memory allocation: make sure to delete this memory allocation later on
        triangles[i] = new std::vector<uint32_t>();
    }
    
    // Create 8 bounding boxes for the 8 subnodes (of the current node)
    BoundingBox3f innerBoxes[8];
    Vector3f corner1 = box.getCenter();
    for (int i = 0 ; i < 8 ; i++) {
        Vector3f corner2 = box.getCorner(i);
        Vector3f minCorner(std::min(corner1.x(), corner2.x()), std::min(corner1.y(), corner2.y()), std::min(corner1.z(), corner2.z()));
        Vector3f maxCorner(std::max(corner1.x(), corner2.x()), std::max(corner1.y(), corner2.y()), std::max(corner1.z(), corner2.z()));
        BoundingBox3f innerBox(minCorner, maxCorner);
        innerBoxes[i] = innerBox;
    }
    
    // Go through every triangle in the current node
    for (std::vector<uint32_t>::const_iterator it = allTriangles->begin() ; it != allTriangles->end() ; it++) {
        // For every subnode
        for (int i = 0 ; i < 8 ; i++) {
            // If the bounding box of the triangle is overlaping with the current subnode bounding box
            if (innerBoxes[i].overlaps(m_mesh->getBoundingBox(*it))) {
                // Add the triangle to its subnode list
                triangles[i]->push_back(*it);
            }
        }
    }
    
    // Here we don't need allTriangles anymore, so we should free the memory
    delete allTriangles;
    
    // Create a new node
    Node *node = new Node();
    for (int i = 0 ; i < 8 ; i++) {
        // Build subnodes
        node->setChild(i, buildNode(innerBoxes[i], triangles[i], numCalls + 1));
    }
    
    // Allocating memory for 9 pointers (8 childs and a vector)
    m_memoryAllocated += sizeof(Node*) * 9;
    
    return node;
}

// This function will explore the given node, setting the index of the closest triangle intersection with the ray (it will return if an intersection was found or not)
bool CustomAccel::exploreForIntersection(const BoundingBox3f &box, const Ray3f &ray_, const bool shadowRay, const Node* node, Intersection &its, uint32_t &index, int64_t &ignoreId) const {
    
    bool foundIntersection(false);
    
    // If this is an empty node
    if (node->isEmptyNode()) {
        // No intersection found so far
        foundIntersection = false;
    }
    
    // If this is a leaf node, the closest triangle has to be in there
    else if (node->isLeafNode()) {
        
        // Make a copy of the ray (we will need to update its '.maxt' value)
        Ray3f ray(ray_);
        
        // Test intersection for all triangles in the leaf
        for (std::vector<uint32_t>::const_iterator it = node->getTriangles()->begin() ; it != node->getTriangles()->end() ; it++) {
            float u, v, t;
            if (ignoreId != *it && m_mesh->rayIntersect(*it, ray, u, v, t)) {
                //An intersection was found!
                foundIntersection = true;
                //We can terminate immediately if this is a shadow ray query
                if (shadowRay)
                    break;
                ray.maxt = its.t = t;
                its.uv = Point2f(u, v);
                its.mesh = m_mesh;
                its.triId = *it;
                index = *it;
            }
        }
    }
    
    // Else, explore subnodes
    else {
        
        /*
         * Create 8 bounding boxes for the 8 subnodes, in the same way we create them during the build of the nodes
         * This will make sure they are in the same order as the child of the node when we call getChild(i)
         */
        BoundingBox3f innerBoxes[8];
        Vector3f corner1 = box.getCenter();
        for (int i = 0 ; i < 8 ; i++) {
            Vector3f corner2 = box.getCorner(i);
            Vector3f minCorner(std::min(corner1.x(), corner2.x()), std::min(corner1.y(), corner2.y()), std::min(corner1.z(), corner2.z()));
            Vector3f maxCorner(std::max(corner1.x(), corner2.x()), std::max(corner1.y(), corner2.y()), std::max(corner1.z(), corner2.z()));
            BoundingBox3f innerBox(minCorner, maxCorner);
            innerBoxes[i] = innerBox;
        }
        
        /*
         * Sort boxes by distance to ray origin and store the sorted indices in sortedIndices
         * This will basically sort subboxes by distance to the origin of the ray, so that we do test on the closest subboxes first
         * Note: if we find an intersection in a subnode, since they're sorted by distance to the origin of the ray,
         * we can stop looking in other nodes, since all of their triangles will be further away
         * We use the array sortedIndices to store the indices of the sorted list of subboxes by distance to the origin of the ray
         * We do this since innerBoxes[8] has to keep its order
         */
        std::array<uint, 8> sortedIndices = {{0, 1, 2, 3, 4, 5, 6, 7}};
        // Note: there's no need to sort if this is just a shadow query
        if (!shadowRay) {
            std::sort(sortedIndices.begin( ), sortedIndices.end( ), [&ray_, &innerBoxes](const uint indexA, const uint indexB)
                      {
                          return innerBoxes[indexA].squaredDistanceTo(ray_.o) < innerBoxes[indexB].squaredDistanceTo(ray_.o);
                      });
        }
        
        
        for (int unsorted_i = 0 ; unsorted_i < 8 ; unsorted_i++) {
            // We will use the sorted index instead of the "random" one
            uint i = sortedIndices[unsorted_i];
            // If the child exists and there's an intersection with the ray
            if (node->getChild(i) != nullptr && innerBoxes[i].rayIntersect(ray_)) {
                // Explore the subnode, if an intersection is found, we can break and terminate since there's no possible closer triangle
                if(exploreForIntersection(innerBoxes[i], ray_, shadowRay, node->getChild(i), its, index, ignoreId)) {
                    foundIntersection = true;
                    break;
                }
            }
        }
    }
    
    return foundIntersection;
}

bool CustomAccel::rayIntersect(const Ray3f &ray_, Intersection &its, bool shadowRay, int64_t ignoreTriangle) const {
    
    // Triangle index of the closest intersection
    uint32_t f = (uint32_t) -1;
    
    // Look for intersection
    bool foundIntersection = exploreForIntersection(m_bbox, ray_, shadowRay, m_mainNode, its, f, ignoreTriangle);
    
    if (foundIntersection && !shadowRay) {
        /* At this point, we now know that there is an intersection,
         and we know the triangle index of the closest such intersection.
         
         The following computes a number of additional properties which
         characterize the intersection (normals, texture coordinates, etc..)
         */
        
        /* Find the barycentric coordinates */
        Vector3f bary;
        bary << 1-its.uv.sum(), its.uv;
        
        /* References to all relevant mesh buffers */
        const Mesh *mesh   = its.mesh;
        const MatrixXf &V  = mesh->getVertexPositions();
        const MatrixXf &N  = mesh->getVertexNormals();
        const MatrixXf &UV = mesh->getVertexTexCoords();
        const MatrixXu &F  = mesh->getIndices();
        
        /* Vertex indices of the triangle */
        uint32_t idx0 = F(0, f), idx1 = F(1, f), idx2 = F(2, f);
        
        Point3f p0 = V.col(idx0), p1 = V.col(idx1), p2 = V.col(idx2);
        
        /* Compute the intersection positon accurately
         using barycentric coordinates */
        its.p = bary.x() * p0 + bary.y() * p1 + bary.z() * p2;
        
        /* Compute proper texture coordinates if provided by the mesh */
        if (UV.size() > 0)
            its.uv = bary.x() * UV.col(idx0) +
            bary.y() * UV.col(idx1) +
            bary.z() * UV.col(idx2);
        
        /* Compute the geometry frame */
        its.geoFrame = Frame((p1-p0).cross(p2-p0).normalized());
        
        if (N.size() > 0) {
            /* Compute the shading frame. Note that for simplicity,
             the current implementation doesn't attempt to provide
             tangents that are continuous across the surface. That
             means that this code will need to be modified to be able
             use anisotropic BRDFs, which need tangent continuity */
            
            its.shFrame = Frame(
                                (bary.x() * N.col(idx0) +
                                 bary.y() * N.col(idx1) +
                                 bary.z() * N.col(idx2)).normalized());
        } else {
            its.shFrame = its.geoFrame;
        }
    }
    
    return foundIntersection;
}

std::string CustomAccel::getStats() const {
    std::string output = "\nAcceleration structure stats : \n";
    output += "  -> Number of interior nodes : " + std::to_string(m_nbNodes);
    output += "\n  -> Number of leaf nodes : " + std::to_string(m_nbLeafs);
    output += "\n  -> Number of triangles stored : " + std::to_string(m_nbStoredTriangles);
    output += "\n  -> Average number of triangles per leaf node : " + std::to_string((float(m_nbStoredTriangles) / float(m_nbLeafs)));
    output += "\n  -> Size of a memory pointer : " + std::to_string(sizeof(Node*)) + " bytes";
    output += "\n  -> Size of a triangle id : " + std::to_string(sizeof(uint32_t)) + " bytes";
    output += "\n  -> Total memory allocated : " + std::to_string(m_memoryAllocated) + " bytes\n";
    return output;
}

NORI_NAMESPACE_END
