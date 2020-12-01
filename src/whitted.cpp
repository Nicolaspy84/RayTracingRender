//
//  whitted.cpp
//  nori
//
//  Created by Nicolas Buisson on 24/03/2018.
//

#include <nori/integrator.h>
#include <nori/scene.h>
#include <nori/bsdf.h>
#include <nori/emitter.h>
#include <Eigen/Core>
#include <nori/dpdf.h>

NORI_NAMESPACE_BEGIN

class WhittedIntegrator : public Integrator {
public:
    WhittedIntegrator(const PropertyList &props) {
        // Nothing here yet
    }
    
    void preprocess(const Scene *scene) {
        
        m_emitters.clear();
        m_emitterDpdf.clear();
        
        // Preprocess for all lights (we generate a discrete pdf instance to randomly select lights)
        for (std::vector<Mesh*>::const_iterator it = scene->getMeshes().begin() ; it != scene->getMeshes().end() ; it++) {
            if((*it)->isEmitter()) {
                m_emitters.push_back((*it));
                m_emitterDpdf.append((*it)->surfaceArea());
            }
        }
        
        if (m_emitters.size() <= 0)
            throw NoriException("Integrator: critical error, no lights in scene");
        
        // We finally normalize the distribution
        m_emitterDpdf.normalize();
    }
    
    Color3f Li(const Scene *scene, Sampler *sampler, const Ray3f &ray) const {
        
        // The final radiance Li = Le + Lr
        Color3f radLi(0.f);
        Color3f radLe(0.f);
        Color3f radLr(0.f);
        
        // Our first step is to find the point we are looking to render
        Intersection its;
        
        // If there's nothing, we render a black pixel
        if (!scene->rayIntersect(ray, its))
            return Color3f(0.0f);
        
        // If the intersected triangle is an emitter, we have to account for Le
        if (its.mesh->isEmitter())
            radLe = its.mesh->getEmitter()->emittedRadiance(its.shFrame.toLocal(-ray.d));
        
        // If the mesh has a diffuse material
        if (its.mesh->getBSDF()->isDiffuse()) {
            // Sample for randomly choosing the light, the triangle on the light and then the point on the triangle
            Point2f sample = sampler->next2D();
            
            // Holds the pdf for choosing the light
            float lightChoicePdf(1.f);
            
            // We choose a light (reusing the y sample)
            Mesh* emitter = m_emitters[m_emitterDpdf.sampleReuse(sample.y(), lightChoicePdf)];
            
            // Now we sample a single point on the area light surface
            
            // This will store the probability for randomly choosing the point
            float lightPointPdf(1.f);
            
            // We sample a random point on its surface
            Intersection emitterIts = emitter->samplePosition(sample, lightPointPdf);
            
            // X->Y vector, and X-> normed vector
            Vector3f xy(emitterIts.p - its.p);
            Vector3f xyN = xy.normalized();
            
            // Distance along the ray
            emitterIts.t = xy.norm();
            
            // Compute the ray to test visibility
            Ray3f emitterRay(its.p, xy, Epsilon * (1.f + its.p.norm()), 1.f - (Epsilon * (1.f + its.p.norm())));
            
            // Now we test visibility
            if(!scene->rayIntersect(emitterRay)) {
                
                // Incident ray, in local coordinates
                Vector3f wi = its.shFrame.toLocal(-ray.d);
                // Outgoing ray, in local coodinates
                Vector3f wo = its.shFrame.toLocal(xyN);
                // Now we create a bsdf query object
                BSDFQueryRecord bsdfQuery(wi, wo, EMeasure::ESolidAngle);
                
                // Compute fr value
                Color3f fr = its.mesh->getBSDF()->eval(bsdfQuery);
                
                // Compute emmited light in the direction of the surface
                Color3f Le = emitterIts.mesh->getEmitter()->emittedRadiance(emitterIts.shFrame.toLocal(-xy));
                
                // Compute the geometric term
                float G = abs(its.shFrame.n.dot(xyN))*abs(emitterIts.shFrame.n.dot(-xyN))/xy.squaredNorm();
                
                // Now we can finally compute the value
                radLr = G * fr * Le / (lightPointPdf*lightChoicePdf);
            }
        }
        
        // Else, it's specular
        else {
            float randomizer = sampler->next1D();
            if (randomizer >= 0.95) {
                return Color3f(0.f);
            }
            
            // Create the BSDF query
            Vector3f wi = its.shFrame.toLocal(-ray.d);
            BSDFQueryRecord bsdfQuery(wi);
            
            // Generate next direction and next ray
            its.mesh->getBSDF()->sample(bsdfQuery, sampler->next2D());
            Ray3f nextRay(its.p, its.shFrame.toWorld(bsdfQuery.wo));
            
            // Compute the reflected/refracted value
            return (1.f/0.95f) * Li(scene, sampler, nextRay);
        }
        
        // We finally have the final value
        radLi = radLe + radLr;
        return radLi;
    }
    
    /// Return a human-readable description for debugging purposes
    std::string toString() const {
        return "WhittedIntegrator[]";
    }
    
private:
    std::vector<Mesh*> m_emitters;
    DiscretePDF m_emitterDpdf;
};

NORI_REGISTER_CLASS(WhittedIntegrator, "whitted");
NORI_NAMESPACE_END
