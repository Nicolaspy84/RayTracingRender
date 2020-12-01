//
//  simple.cpp
//  nori
//
//  Created by Nicolas Buisson on 13/03/2018.
//

#include <nori/integrator.h>
#include <nori/scene.h>
#include <Eigen/Core>

NORI_NAMESPACE_BEGIN

class SimpleIntegrator : public Integrator {
public:
    SimpleIntegrator(const PropertyList &props) {
        m_lightPosition = props.getPoint("position");
        m_lightEnergy = props.getColor("energy");
    }
    
    Color3f Li(const Scene *scene, Sampler *sampler, const Ray3f &ray) const {
        
        // Our first step is to find the point we are looking to render
        Intersection its;
        
        // If there's nothing, we render a black pixel
        if (!scene->rayIntersect(ray, its))
            return Color3f(0.0f);
        
        // Now we look if this point is visible
        
        // This stores the normalized vector from the point to the light
        Vector3f shadowDir = m_lightPosition - its.p;
        shadowDir.normalize();
        
        // This stores the angle between the normal of the point surface and the direction to the light
        float cosTheta = shadowDir.dot(its.shFrame.n);
        
        // This is a black point if the triangle is not visible
        if (cosTheta <= 0.f || scene->rayIntersect(Ray3f(its.p, shadowDir)))
            return Color3f(0.0f);
        
        // Now we now that this point is visible, so we render a color
        Color3f luminance = cosTheta/(4.f*M_PI*M_PI*(m_lightPosition-its.p).squaredNorm()) * m_lightEnergy;
        return luminance;
    }
    
    /// Return a human-readable description for debugging purposes
    std::string toString() const {
        return "SimpleIntegrator[]";
    }
    
private:
    // Stores the light's position
    Point3f m_lightPosition;
    // Stores the light's energy
    Color3f m_lightEnergy;
};

NORI_REGISTER_CLASS(SimpleIntegrator, "simple");
NORI_NAMESPACE_END
