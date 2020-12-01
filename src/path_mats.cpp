//
//  path_mats.cpp
//  nori
//
//  Created by Nicolas Buisson on 12/04/2018.
//

#include <nori/integrator.h>
#include <nori/scene.h>
#include <nori/bsdf.h>
#include <nori/warp.h>
#include <nori/emitter.h>
#include <Eigen/Core>
#include <nori/dpdf.h>

NORI_NAMESPACE_BEGIN

class Path_mats : public Integrator {
public:
    Path_mats(const PropertyList &props) {
        m_probStopQuery = props.getFloat("probStopQuery", 0.05f);
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
        
        float q(m_probStopQuery);
        
        // Variables for the iterative loop
        // Note : final radiance is equal to add + mult
        Ray3f inRay(ray);
        Color3f add(0.f);
        Color3f mult(1.f);
        
        // Start the main loop
        while(true) {
            
            // The final radiance is L = Le + Li
            Color3f radLe(0.f);
            Color3f radLi(0.f);
            
            // Main intersection
            Intersection its;
            
            // If there's nothing, we break the loop (no further calculation needed)
            if (!scene->rayIntersect(inRay, its)) {
                mult = 0.f;
                break;
            }
            
            // If this is an emitter
            if (its.mesh->isEmitter())
                radLe = its.mesh->getEmitter()->emittedRadiance(its.shFrame.toLocal(-inRay.d));
            
            // Set the bsdf query
            BSDFQueryRecord bRec(its.shFrame.toLocal(-inRay.d.normalized()));
            
            // Indirect illumination
            radLi = its.mesh->getBSDF()->sample(bRec, sampler->next2D());
            
            // We update the next ray
            inRay = Ray3f(its.p, its.shFrame.toWorld(bRec.wo));
            
            add += mult * radLe;
            mult *= radLi / (1.f - q);
            
            // To avoid infinite loop
            if (sampler->next1D() <= q) {
                mult = 0.f;
                break;
            }
        }
        
        // We finally have the final value
        Color3f radL = add + mult;
        
        // If this is not a valid value, we prefer to abort for now
        if (!radL.isValid())
            throw NoriException("invalid radiance value\nadd=%s\nmult=%s", add.toString(), mult.toString());
        
        return radL;
    }
    
    /// Return a human-readable description for debugging purposes
    std::string toString() const {
        return "Path_mats_integrator[]";
    }
    
private:
    float m_probStopQuery;
    std::vector<Mesh*> m_emitters;
    DiscretePDF m_emitterDpdf;
};

NORI_REGISTER_CLASS(Path_mats, "path_mats");
NORI_NAMESPACE_END
