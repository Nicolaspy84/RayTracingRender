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

#pragma once

#include <nori/node.h>
#include <nori/mesh.h>
#include <vector>


NORI_NAMESPACE_BEGIN

/**
 * \brief Acceleration data structure for ray intersection queries
 *
 * The current implementation falls back to a brute force loop
 * through the geometry.
 */

class CustomAccel {
public:
    
    ~CustomAccel();
    
    /**
     * \brief Register a triangle mesh for inclusion in the acceleration
     * data structure
     *
     * This function can only be used before \ref build() is called
     */
    void addMesh(Mesh *mesh);
    
    /**
     * \brief Build the acceleration data structure as an octree
     */
    void build();
    Node* buildNode(const BoundingBox3f &box, std::vector<uint32_t>* allTriangles, uint32_t numCalls);
    
    /// Return an axis-aligned box that bounds the scene
    const BoundingBox3f &getBoundingBox() const { return m_bbox; }
    
    /**
     * \brief Intersect a ray against all triangles stored in the scene and
     * return detailed intersection information
     *
     * \param ray
     *    A 3-dimensional ray data structure with minimum/maximum extent
     *    information
     *
     * \param its
     *    A detailed intersection record, which will be filled by the
     *    intersection query
     *
     * \param shadowRay
     *    \c true if this is a shadow ray query, i.e. a query that only aims to
     *    find out whether the ray is blocked or not without returning detailed
     *    intersection information.
     *
     * \param ignoreTriangle
     *    Triangle id's to ignore (usually used to avoid casting a shadow on itself)
     *    Set it to -1 to try all triangles
     *
     * \return \c true if an intersection was found
     */
    bool rayIntersect(const Ray3f &ray, Intersection &its, bool shadowRay, int64_t ignoreTriangle) const;
    
    /**
     * \brief This function will explore the node provided, looking for an intersection with a triangle
     *
     * \param box
     *    The main bounding box
     *
     * \param ray
     *    A 3-dimensional ray data structure with minimum/maximum extent information
     *
     * \param shadowRay
     *    \c true if this is a shadow ray query (no detailed intersection information given here)
     *
     * \param node
     *    The main node
     *
     * \param its
     *    A detailed intersection record, which will be filled by the intersection query
     *
     * \param index
     *    Stores the index of the intersected triangle (if it exists)
     *
     * \return \c true if an intersection was found
     */
    bool exploreForIntersection(const BoundingBox3f &box, const Ray3f &ray, const bool shadowRay, const Node* node, Intersection &its, uint32_t &index, int64_t &ignoreId) const;
    
    /// This will just return a string with the build stats in it
    std::string getStats() const;
    
private:
    Mesh         *m_mesh = nullptr;        ///< Mesh (only a single one for now)
    Node         *m_mainNode = nullptr;    ///< Main node (octree)
    BoundingBox3f m_bbox;                  ///< Bounding box of the entire scene
    uint32_t      m_nbLeafs = 0;           ///< Total number of leafs after build() is called
    uint32_t      m_nbNodes = 0;           ///< Total number of nodes after build() is called
    uint32_t      m_nbStoredTriangles = 0; ///< Total number of triangles stored (including duplicates) after build() is called
    uint64_t      m_memoryAllocated = 0;
};

NORI_NAMESPACE_END
