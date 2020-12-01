//
//  ambocc.cpp
//  nori
//
//  Created by Nicolas Buisson on 13/03/2018.
//

#include <nori/integrator.h>
#include <nori/scene.h>
#include <Eigen/Core>
#include <nori/warp.h>

NORI_NAMESPACE_BEGIN

class AmbOccIntegrator : public Integrator {
public:
    AmbOccIntegrator(const PropertyList &props) {
        // Number of rays to compute the occlusion
        m_nbSamples = props.getInteger("nbSamples");
    }
    
    Color3f Li(const Scene *scene, Sampler *sampler, const Ray3f &ray) const {
        // Generate a warp class
        Warp warp;
        
        // Our first step is to find the point we are looking to render
        Intersection its;
        
        // If there's nothing, we render a black pixel
        if (!scene->rayIntersect(ray, its))
            return Color3f(0.0f);
        
        // Stores the number of occlusions found
        uint32_t occlusions = 0;
        
        // Now we compute the ambiant occlusion
        for (uint i = 0 ; i < m_nbSamples ; i++) {
            // Get a point on the hemisphere
            Vector3f point = warp.squareToCosineHemisphere(sampler->next2D());
            Frame normalFrame(its.shFrame.s, its.shFrame.t, its.shFrame.n);
            // Create the vector going from the intersection to the point on the hemisphere
            Vector3f rayDir = normalFrame.toWorld(point);
            // Create the occlusion test ray
            Ray3f occlusionRay(its.p, rayDir);
            // Test if this ray is occluded
            if (scene->rayIntersect(occlusionRay))
                occlusions += 1;
        }
        
        // Compute luminance (a simple average)
        float luminance = 1.f - float(occlusions)/float(m_nbSamples);
        Color3f color(luminance, luminance, luminance);
        return color;
    }
    
    /// Return a human-readable description for debugging purposes
    std::string toString() const {
        return "AmbiantOcclusionIntegrator[nbSamples=" + std::to_string(m_nbSamples) + "]";
    }
    
private:
    uint m_nbSamples;
};

NORI_REGISTER_CLASS(AmbOccIntegrator, "ao");
NORI_NAMESPACE_END
