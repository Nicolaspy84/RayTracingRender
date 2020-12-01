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

#include <nori/mesh.h>
#include <nori/bbox.h>
#include <nori/bsdf.h>
#include <nori/emitter.h>
#include <nori/warp.h>
#include <Eigen/Geometry>
#include <filesystem/resolver.h>

NORI_NAMESPACE_BEGIN

Mesh::Mesh(const PropertyList &props) {
    m_isVolume = props.getBoolean("volume", false);
    m_scattering = props.getFloat("scattering", 0.f);
    m_absorption = props.getFloat("absoption", 1.f);
    m_extinction = m_scattering + m_absorption;
    m_normalFile = getFileResolver()->resolve(props.getString("normalFile", "")).str();
    if (m_normalFile != getFileResolver()->resolve("").str())
        m_normalMap = new Bitmap(m_normalFile);
}

Mesh::~Mesh() {
    delete m_bsdf;
    delete m_emitter;
    delete m_dpdf;
}

void Mesh::activate() {
    if (!m_bsdf) {
        /* If no material was assigned, instantiate a diffuse BRDF */
        m_bsdf = static_cast<BSDF *>(
            NoriObjectFactory::createInstance("diffuse", PropertyList()));
    }
    
    if (m_emitter) {
        /* If an emitter has been assigned */
        createDPDF();
    }
    
    m_surfaceArea = 0.f;
    for (uint32_t idx = 0; idx <getTriangleCount(); ++idx) {
        m_surfaceArea += surfaceArea(idx);
    }
}

void Mesh::createDPDF() {
    // Allocate enough memory for the discrete PDF structure
    m_dpdf = new DiscretePDF(getTriangleCount());
    
    // For every triangle, add a new entry which corresponds to the surface area of the triangle
    for (uint32_t idx = 0; idx < getTriangleCount(); ++idx) {
        m_dpdf->append(surfaceArea(idx));
    }
    
    // Now we need to normalize, after this everything is setup
    m_dpdf->normalize();
}

Intersection Mesh::samplePosition(const Point2f &sample, float &pdf) const {
    
    // Create the intersection
    Intersection its;
    
    // Randomly select a triangle according to the specified distribution
    float x = sample.x();
    uint32_t idx = m_dpdf->sampleReuse(x, pdf);
    
    // Set pdf
    pdf /= surfaceArea(idx);
    
    // Compute the point
    float alpha = 1.f - sqrt(1.f - sample.y());
    float beta = x*(1.f - alpha);
    Point3f bary(alpha, beta, 1.f - alpha - beta);
    its.p = Point3f(m_V.col(m_F(0, idx))*bary.x() + m_V.col(m_F(1, idx))*bary.y() + m_V.col(m_F(2, idx))*bary.z());
    
    // Storing the triangle id could be useful for future use
    its.triId = idx;
    
    // Mesh pointer
    its.mesh = this;
    
    /* Vertex indices of the triangle */
    uint32_t idx0 = m_F(0, idx), idx1 = m_F(1, idx), idx2 = m_F(2, idx);
    
    Point3f p0 = m_V.col(idx0), p1 = m_V.col(idx1), p2 = m_V.col(idx2);
    
    /* Compute the intersection positon accurately
     using barycentric coordinates */
    its.p = bary.x() * p0 + bary.y() * p1 + bary.z() * p2;
    
    /* Compute proper texture coordinates if provided by the mesh */
    if (m_UV.size() > 0)
        its.uv = bary.x() * m_UV.col(idx0) +
        bary.y() * m_UV.col(idx1) +
        bary.z() * m_UV.col(idx2);
    
    /* Compute the geometry frame */
    its.geoFrame = Frame((p1-p0).cross(p2-p0).normalized());
    
    if (m_N.size() > 0) {
        /* Compute the shading frame. Note that for simplicity,
         the current implementation doesn't attempt to provide
         tangents that are continuous across the surface. That
         means that this code will need to be modified to be able
         use anisotropic BRDFs, which need tangent continuity */
        
        its.shFrame = Frame(
                            (bary.x() * m_N.col(idx0) +
                             bary.y() * m_N.col(idx1) +
                             bary.z() * m_N.col(idx2)).normalized());
    }
    else {
        its.shFrame = its.geoFrame;
    }
    
    // Return the final sampled intersection (ie point with data structure for the normal...)
    return its;
}

float Mesh::surfaceArea(uint32_t index) const {
    uint32_t i0 = m_F(0, index), i1 = m_F(1, index), i2 = m_F(2, index);

    const Point3f p0 = m_V.col(i0), p1 = m_V.col(i1), p2 = m_V.col(i2);

    return 0.5f * Vector3f((p1 - p0).cross(p2 - p0)).norm();
}

float Mesh::surfaceArea() const {
    return m_surfaceArea;
}

bool Mesh::rayIntersect(uint32_t index, const Ray3f &ray, float &u, float &v, float &t) const {
    uint32_t i0 = m_F(0, index), i1 = m_F(1, index), i2 = m_F(2, index);
    const Point3f p0 = m_V.col(i0), p1 = m_V.col(i1), p2 = m_V.col(i2);

    /* Find vectors for two edges sharing v[0] */
    Vector3f edge1 = p1 - p0, edge2 = p2 - p0;

    /* Begin calculating determinant - also used to calculate U parameter */
    Vector3f pvec = ray.d.cross(edge2);

    /* If determinant is near zero, ray lies in plane of triangle */
    float det = edge1.dot(pvec);

    if (det > -1e-8f && det < 1e-8f)
        return false;
    float inv_det = 1.0f / det;

    /* Calculate distance from v[0] to ray origin */
    Vector3f tvec = ray.o - p0;

    /* Calculate U parameter and test bounds */
    u = tvec.dot(pvec) * inv_det;
    if (u < 0.0 || u > 1.0)
        return false;

    /* Prepare to test V parameter */
    Vector3f qvec = tvec.cross(edge1);

    /* Calculate V parameter and test bounds */
    v = ray.d.dot(qvec) * inv_det;
    if (v < 0.0 || u + v > 1.0)
        return false;

    /* Ray intersects triangle -> compute t */
    t = edge2.dot(qvec) * inv_det;

    return t >= ray.mint && t <= ray.maxt;
}

const Vector3f Mesh::getNormal(Point2f uv) const {
    if (!normalMap())
        throw NoriException("Mesh: tried to access normal map when none are set up!");
    Color3f colNormal = m_normalMap->coeff((1.f - uv.y()) * m_normalMap->rows(), (1.f - uv.x()) * m_normalMap->cols());
    Vector3f normal = Vector3f(2.f * colNormal.x() - 1.f, 2.f * colNormal.y() - 1.f, 2.f * colNormal.z() - 1.f);
    return normal;
}

BoundingBox3f Mesh::getBoundingBox(uint32_t index) const {
    BoundingBox3f result(m_V.col(m_F(0, index)));
    result.expandBy(m_V.col(m_F(1, index)));
    result.expandBy(m_V.col(m_F(2, index)));
    return result;
}

Point3f Mesh::getCentroid(uint32_t index) const {
    return (1.0f / 3.0f) *
        (m_V.col(m_F(0, index)) +
         m_V.col(m_F(1, index)) +
         m_V.col(m_F(2, index)));
}

void Mesh::addChild(NoriObject *obj) {
    switch (obj->getClassType()) {
        case EBSDF:
            if (m_bsdf)
                throw NoriException(
                    "Mesh: tried to register multiple BSDF instances!");
            m_bsdf = static_cast<BSDF *>(obj);
            break;

        case EEmitter: {
                Emitter *emitter = static_cast<Emitter *>(obj);
                if (m_emitter)
                    throw NoriException(
                        "Mesh: tried to register multiple Emitter instances!");
                m_emitter = emitter;
            }
            break;

        default:
            throw NoriException("Mesh::addChild(<%s>) is not supported!",
                                classTypeName(obj->getClassType()));
    }
}

std::string Mesh::toString() const {
    return tfm::format(
        "Mesh[\n"
        "  name = \"%s\",\n"
        "  vertexCount = %i,\n"
        "  triangleCount = %i,\n"
        "  bsdf = %s,\n"
        "  emitter = %s\n"
        "]",
        m_name,
        m_V.cols(),
        m_F.cols(),
        m_bsdf ? indent(m_bsdf->toString()) : std::string("null"),
        m_emitter ? indent(m_emitter->toString()) : std::string("null")
    );
}

std::string Intersection::toString() const {
    if (!mesh)
        return "Intersection[invalid]";

    return tfm::format(
        "Intersection[\n"
        "  p = %s,\n"
        "  t = %f,\n"
        "  uv = %s,\n"
        "  shFrame = %s,\n"
        "  geoFrame = %s,\n"
        "  mesh = %s\n"
        "]",
        p.toString(),
        t,
        uv.toString(),
        indent(shFrame.toString()),
        indent(geoFrame.toString()),
        mesh ? mesh->toString() : std::string("null")
    );
}

NORI_NAMESPACE_END
