//
//  path_mis.cpp
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

class Path_mis : public Integrator {
public:
    Path_mis(const PropertyList &props) {
        m_probStopQuery = props.getFloat("probStopQuery", 0.05f);
    }
    
    ~Path_mis() {
        delete m_emitterDpdf;
    }
    
    void preprocess(const Scene *scene) {
        
        m_totalSurface = 0.f;
        m_triangleCount = 0;
        
        m_emitters.clear();
        if (!m_emitterDpdf)
            delete m_emitterDpdf;
        if (!m_triangleDpdf)
            delete m_triangleDpdf;
        
        // Preprocess for all lights
        for (std::vector<Mesh*>::const_iterator it = scene->getMeshes().begin() ; it != scene->getMeshes().end() ; it++) {
            if((*it)->isEmitter()) {
                m_emitters.push_back((*it));
                m_totalSurface += (*it)->surfaceArea();
                m_triangleCount += (*it)->getTriangleCount();
            }
        }
        
        if (m_emitters.size() <= 0)
            throw NoriException("Integrator: critical error, no lights in scene");
        
        // Generate new dpdfs
        m_emitterDpdf = new DiscretePDF(m_emitters.size());
        m_triangleDpdf = new DiscretePDF(m_triangleCount);
        for (std::vector<Mesh*>::const_iterator it = m_emitters.begin() ; it != m_emitters.end() ; it++) {
            
            // Add a new emitter
            m_emitterDpdf->append(1.f);
            
            // For every triangle, add a new entry which corresponds to the surface area of the triangle
            for (uint32_t idx = 0; idx < (*it)->getTriangleCount(); ++idx) {
                m_triangleDpdf->append((*it)->surfaceArea(idx));
            }
        }

        // We finally normalize the distributions
        m_triangleDpdf->normalize();
        m_emitterDpdf->normalize();
    }
    
    Color3f Li(const Scene *scene, Sampler *sampler, const Ray3f &ray) const {
        
        float q(m_probStopQuery);
        
        // Variables for the iterative loop
        // Note : final radiance is equal to add + mult
        Ray3f inRay(ray);
        bool prevSpec(true);
        Color3f add(0.f);
        Color3f mult(1.f);
        const Mesh* volume = nullptr;
        
        // Start the main loop
        while(true) {
            
            bool forceNextCast(false);
            
            // The final radiance is L = Le + Ld + Li
            Color3f radLe(0.f);
            Color3f radLd(0.f);
            Color3f radLi(0.f);
            
            // Main intersection
            Intersection its;
            
            // If there's nothing, we break the loop (no further calculation needed)
            if (!scene->rayIntersect(inRay, its)) {
                mult = 0.f;
                break;
            }
            
            // If we intersect a volume box and we are not in a volume
            if (its.mesh->isVolume() && volume == nullptr) {
                volume = its.mesh;
                inRay = Ray3f(its.p, inRay.d);
                forceNextCast = true;
            }
            
            // If we are already in a volume
            else if (volume != nullptr) {
                
                // We now have the maximum distance
                float tmax = its.t;
                
                // We can generate a random distance
                float t = volume->distance(sampler->next1D());
                
                if (t < tmax) {
                    // We update the next ray
                    inRay = Ray3f(inRay(t), Warp::squareToUniformSphere(sampler->next2D()));
                    
                    prevSpec = false;
                    
                    mult *= volume->albedo();
                    mult /= 1.f - q;
                }
                
                else {
                    
                    // Set the bsdf query
                    BSDFQueryRecord bRec(its.shFrame.toLocal(-inRay.d.normalized()));
                    bRec.uv = its.uv;
                    
                    // Sample direction
                    its.mesh->getBSDF()->sample(bRec, sampler->next2D());
                    
                    // We update the next ray
                    inRay = Ray3f(its.p, its.shFrame.toWorld(bRec.wo));
                    
                    if (!its.mesh->getBSDF()->isDiffuse()) {
                        prevSpec = true;
                        mult /= 1.f - q;
                    }
                    
                    else {
                        prevSpec = false;
                        mult *= its.mesh->getBSDF()->eval(bRec) / its.mesh->getBSDF()->pdf(bRec);
                        mult /= 1.f - q;
                    }
                    
                    // If the outgoing direction is outside of the volume
                    if (inRay.d.dot(its.geoFrame.n) >= 0.f)
                        volume = nullptr;
                }
            }
            
            // If this is a specular surface
            else if (!its.mesh->getBSDF()->isDiffuse()) {
                // Set the bsdf query
                BSDFQueryRecord bRec(its.shFrame.toLocal(-inRay.d.normalized()));
                bRec.uv = its.uv;
                
                // Sample direction
                its.mesh->getBSDF()->sample(bRec, sampler->next2D());
                
                // We update the next ray
                inRay = Ray3f(its.p, its.shFrame.toWorld(bRec.wo));
                
                prevSpec = true;
                
                mult /= 1.f - q;
                
            }
            
            // If this is a diffuse surface
            else {
                // Compute direct illumination from all lights
                if (its.mesh->getBSDF()->isDiffuse())
                    radLd = LdRandom(scene, sampler, inRay, its);
                
                // If this is an emitter
                if (its.mesh->isEmitter())
                    radLe = its.mesh->getEmitter()->emittedRadiance(its.shFrame.toLocal(-inRay.d));
                
                // Set the bsdf query
                BSDFQueryRecord bRec(its.shFrame.toLocal(-inRay.d.normalized()));
                bRec.uv = its.uv;
                
                // Indirect illumination
                radLi = its.mesh->getBSDF()->sample(bRec, sampler->next2D());
                
                // We update the next ray
                inRay = Ray3f(its.p, its.shFrame.toWorld(bRec.wo));
                
                // Test if the indirect new ray hits a light source
                Intersection lightIts;
                if (scene->rayIntersect(inRay, lightIts) && lightIts.mesh->isEmitter()) {
                    float pBRDF = its.mesh->getBSDF()->pdf(bRec);
                    float pLight = pow(lightIts.t, 2) / (lightIts.mesh->surfaceArea() * m_emitters.size() * lightIts.shFrame.cosTheta(lightIts.shFrame.toLocal(-its.shFrame.toWorld(bRec.wo))));
                    Color3f fr = its.mesh->getBSDF()->eval(bRec);
                    Color3f Le = lightIts.mesh->getEmitter()->emittedRadiance(lightIts.shFrame.toLocal(-inRay.d));
                    radLd += fr * Le * abs(its.shFrame.n.dot((lightIts.p - its.p).normalized())) / (pLight + pBRDF);
                }
                
                // Account for emitted radiance only if the previous mesh was specular
                if (prevSpec)
                    add += mult * (radLe + radLd);
                else
                    add += mult * radLd;
                mult *= radLi / (1.f - q);
                
                prevSpec = false;
            }
            
            // To avoid infinite loop
            if (sampler->next1D() <= q && !forceNextCast) {
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
    
    
    
    
    
    
    
    
    
    
    
    
    
    /*
    Color3f Li(const Scene *scene, Sampler *sampler, const Ray3f &ray) const {
        
        float q(m_probStopQuery);
        
        // Variables for the iterative loop
        // Note : final radiance is equal to add + mult
        Ray3f inRay(ray);
        bool prevSpec(true);
        Color3f add(0.f);
        Color3f mult(1.f);
        
        // Start the main loop
        while(true) {
            
            // The final radiance is L = Le + Ld + Li
            Color3f radLe(0.f);
            Color3f radLd(0.f);
            Color3f radLi(0.f);
            
            // Main intersection
            Intersection its;
            
            // If there's nothing, we break the loop (no further calculation needed)
            if (!scene->rayIntersect(inRay, its)) {
                mult = 0.f;
                break;
            }
            
            if (!its.mesh->getBSDF()->isDiffuse()) {
                // Set the bsdf query
                BSDFQueryRecord bRec(its.shFrame.toLocal(-inRay.d.normalized()));
                bRec.uv = its.uv;
                
                // Sample direction
                its.mesh->getBSDF()->sample(bRec, sampler->next2D());
                
                // We update the next ray
                inRay = Ray3f(its.p, its.shFrame.toWorld(bRec.wo));
                
                prevSpec = true;
                
                mult /= 1.f - q;
            }
            
            else {
                // Compute direct illumination from all lights
                if (its.mesh->getBSDF()->isDiffuse())
                    radLd = LdRandom(scene, sampler, inRay, its);
                
                // If this is an emitter
                if (its.mesh->isEmitter())
                    radLe = its.mesh->getEmitter()->emittedRadiance(its.shFrame.toLocal(-inRay.d));
                
                // Set the bsdf query
                BSDFQueryRecord bRec(its.shFrame.toLocal(-inRay.d.normalized()));
                bRec.uv = its.uv;
                
                // Indirect illumination
                radLi = its.mesh->getBSDF()->sample(bRec, sampler->next2D());
                
                // We update the next ray
                inRay = Ray3f(its.p, its.shFrame.toWorld(bRec.wo));
                
                // Test if the indirect new ray hits a light source
                Intersection lightIts;
                if (scene->rayIntersect(inRay, lightIts) && lightIts.mesh->isEmitter()) {
                    float pBRDF = its.mesh->getBSDF()->pdf(bRec);
                    float pLight = pow(lightIts.t, 2) / (lightIts.mesh->surfaceArea() * m_emitters.size() * lightIts.shFrame.cosTheta(lightIts.shFrame.toLocal(-its.shFrame.toWorld(bRec.wo))));
                    Color3f fr = its.mesh->getBSDF()->eval(bRec);
                    Color3f Le = lightIts.mesh->getEmitter()->emittedRadiance(lightIts.shFrame.toLocal(-inRay.d));
                    radLd += fr * Le * abs(its.shFrame.n.dot((lightIts.p - its.p).normalized())) / (pLight + pBRDF);
                }
                
                // Account for emitted radiance only if the previous mesh was specular
                if (prevSpec)
                    add += mult * (radLe + radLd);
                else
                    add += mult * radLd;
                mult *= radLi / (1.f - q);
                
                prevSpec = false;
            }
            
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
    */
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    /// Compute direct illumination for a random light
    Color3f LdRandom(const Scene *scene, Sampler *sampler, const Ray3f &inRay, const Intersection &its) const {
        Point2f sample = sampler->next2D();
        Mesh *emitter = m_emitters[m_emitterDpdf->sampleReuse(sample.y())];
        return Ld(scene, sample, inRay, its, emitter);
    }
    
    /// Compute direct illumination from the given emitter (this will suppose the mesh is diffuse)
    Color3f Ld(const Scene *scene, Point2f sample, const Ray3f &inRay, const Intersection &its, Mesh *emitter) const {
        
        // This will store the probability for randomly choosing the point
        float lightPointPdf(1.f);
        
        // We sample a random point on its surface
        Intersection emitterIts = emitter->samplePosition(sample, lightPointPdf);
        
        // X->Y vector, and X->Y normed vector
        Vector3f xy(emitterIts.p - its.p);
        Vector3f xyN = xy.normalized();
        
        // Distance along the ray
        emitterIts.t = xy.norm();
        
        // Compute the ray to test visibility
        Ray3f emitterRay(its.p, xy, Epsilon * (1.f + its.p.norm()), 1.f - (Epsilon * (1.f + its.p.norm())));
        
        // Now we test visibility
        if(scene->rayIntersect(emitterRay))
            return Color3f(0.f);
        
        // Incident ray, in local coordinates
        Vector3f wi = its.shFrame.toLocal(-inRay.d);
        // Outgoing ray, in local coodinates
        Vector3f wo = its.shFrame.toLocal(xyN);
        // Now we create a bsdf query object
        BSDFQueryRecord bsdfQuery(wi, wo, EMeasure::ESolidAngle);
        bsdfQuery.uv = its.uv;
        
        // Compute fr value
        Color3f fr = its.mesh->getBSDF()->eval(bsdfQuery);
        
        if (!fr.isValid())
            throw NoriException("invalid fr value\nfr=%s\nemitter=%s\nits=%s", fr.toString(), emitter->toString(), its.toString());
        
        // Compute emmited light in the direction of the surface
        Color3f Le = emitterIts.mesh->getEmitter()->emittedRadiance(emitterIts.shFrame.toLocal(-xy));
        
        if (!Le.isValid())
            throw NoriException("invalid Le value\nfr=%s\nemitter=%s\nits=%s", Le.toString(), emitter->toString(), its.toString());
        
        // Now it is time to cumpute the probability with which this direction would have been sampled (using the brdf)
        float pBRDF = its.mesh->getBSDF()->pdf(bsdfQuery);
        
        // Convert from surface to solid angles pdf, and don't forget to take into account our uniform distribution when choosing the light
        lightPointPdf /= m_emitters.size();
        float pLight = lightPointPdf * xy.squaredNorm() / abs(emitterIts.shFrame.n.dot(-xyN));
        
        // Now we can finally compute the value
        return abs(its.shFrame.n.dot(xyN)) * fr * Le / (pLight + pBRDF);
    }
    
    /// Return a human-readable description for debugging purposes
    std::string toString() const {
        return "Path_mis_integrator[]";
    }
    
private:
    float m_probStopQuery;
    std::vector<Mesh*> m_emitters;
    DiscretePDF *m_emitterDpdf;
    DiscretePDF *m_triangleDpdf;
    float m_totalSurface;
    uint32_t m_triangleCount;
};

NORI_REGISTER_CLASS(Path_mis, "path_mis");
NORI_NAMESPACE_END
