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

#include <nori/bsdf.h>
#include <nori/frame.h>
#include <nori/warp.h>
#include <nori/common.h>

NORI_NAMESPACE_BEGIN

class Microfacet : public BSDF {
public:
    Microfacet(const PropertyList &propList) {
        /* RMS surface roughness */
        m_alpha = propList.getFloat("alpha", 0.1f);

        /* Interior IOR (default: BK7 borosilicate optical glass) */
        m_intIOR = propList.getFloat("intIOR", 1.5046f);

        /* Exterior IOR (default: air) */
        m_extIOR = propList.getFloat("extIOR", 1.000277f);

        /* Albedo of the diffuse base material (a.k.a "kd") */
        m_kd = propList.getColor("kd", Color3f(0.5f));

        /* To ensure energy conservation, we must scale the 
           specular component by 1-kd. 

           While that is not a particularly realistic model of what 
           happens in reality, this will greatly simplify the 
           implementation. Please see the course staff if you're 
           interested in implementing a more realistic version 
           of this BRDF. */
        m_ks = 1 - m_kd.maxCoeff();
    }

    /// Evaluate the BRDF for the given pair of directions
    Color3f eval(const BSDFQueryRecord &bRec) const {
        
        // If the incident ray comes from the back side of the triangle, we eval to 0
        if (Frame::cosTheta(bRec.wi) <= 0.f)
            return Color3f(0.f);
        
        // Sum vector (normalized)
        Vector3f wh = (bRec.wi + bRec.wo).normalized();
        
        // Beckmann term
        float D = Warp::squareToBeckmannPdf(wh, m_alpha);
        
        // Fresnel term
        float F = fresnel(wh.dot(bRec.wi), m_extIOR, m_intIOR);
        
        // Shadowing term
        float G = G1(bRec.wi, wh) * G1(bRec.wo, wh);
        
        // Final brdf
        Color3f brdf = m_kd * INV_PI + m_ks * D * F * G / (4.f*Frame::cosTheta(bRec.wi)*Frame::cosTheta(bRec.wo)*Frame::cosTheta(wh));
        
        if (!brdf.isValid())
            throw NoriException("invalid brdf value\nbrdf=%s\nwi=%s\nwo=%s\nwh=%s\nD=%s\nF=%s\nG=%s\nalpha=%s", brdf.toString(), bRec.wi.toString(), bRec.wo.toString(), wh.toString(), std::to_string(D), std::to_string(F), std::to_string(G), std::to_string(m_alpha));
        
        return brdf;
        
    }

    /// Evaluate the sampling density of \ref sample() wrt. solid angles
    float pdf(const BSDFQueryRecord &bRec) const {
        
        if (Frame::cosTheta(bRec.wo) <= 0.f)
            return 0.f;
     
        // Sum vector (normalized)
        Vector3f wh = (bRec.wi + bRec.wo).normalized();
        
        // Beckmann term
        float D = Warp::squareToBeckmannPdf(wh, m_alpha);
        
        // Jacobian term
        float J = 1.f / (4*wh.dot(bRec.wo));
        
        // Final density
        return m_ks * D * J + (1.f - m_ks) * Frame::cosTheta(bRec.wo) / M_PI;
    }

    /// Sample the BRDF
    Color3f sample(BSDFQueryRecord &bRec, const Point2f &_sample) const {

        // Note: Once you have implemented the part that computes the scattered
        // direction, the last part of this function should simply return the
        // BRDF value divided by the solid angle density and multiplied by the
        // cosine factor from the reflection equation, i.e.
        // return eval(bRec) * Frame::cosTheta(bRec.wo) / pdf(bRec);
        
        // First test if the ray comes from the correct side of the surface
        if (Frame::cosTheta(bRec.wi) < 0)
            return Color3f(0.0f);
        
        // Set the measure type
        bRec.measure = ESolidAngle;
        
        // Relative index of refraction: no change
        bRec.eta = 1.0f;
        
        // Sample value we will have to update
        Point2f s = _sample;
        
        // Diffuse reflection
        if (s.x() > m_ks) {
            // In order to reuse the sample
            s.x() = (s.x() - m_ks) / (1.f - m_ks);
            
            /* Warp a uniformly distributed sample on [0,1]^2
             to a direction on a cosine-weighted hemisphere */
            bRec.wo = Warp::squareToCosineHemisphere(s);
        }
        
        else {
            // In order to reuse the sample
            s.x() /= m_ks;
            
            // New normal frame
            Frame normalFrame(Warp::squareToBeckmann(s, m_alpha));
            
            // Reflection using the new normal
            Vector3f reflection = normalFrame.toLocal(bRec.wi);
            reflection.x() = -reflection.x();
            reflection.y() = -reflection.y();
            bRec.wo = normalFrame.toWorld(reflection);
        }
        
        if (Frame::cosTheta(bRec.wo) <= 0.f)
            return 0.f;
        
        float prob = pdf(bRec);
        if (prob == 0.f)
            throw NoriException("Math error in microfacet.cpp, pdf(bRec) == 0.f");
        
        return eval(bRec) * Frame::cosTheta(bRec.wo) / prob;
    }

    bool isDiffuse() const {
        /* While microfacet BRDFs are not perfectly diffuse, they can be
           handled by sampling techniques for diffuse/non-specular materials,
           hence we return true here */
        return true;
    }

    std::string toString() const {
        return tfm::format(
            "Microfacet[\n"
            "  alpha = %f,\n"
            "  intIOR = %f,\n"
            "  extIOR = %f,\n"
            "  kd = %s,\n"
            "  ks = %f\n"
            "]",
            m_alpha,
            m_intIOR,
            m_extIOR,
            m_kd.toString(),
            m_ks
        );
    }
private:
    
    // A simple method used to evaluate the shadowing term, inputs are in local coordinates
    float G1(Vector3f wv, Vector3f wh) const {
        float a = wv.dot(wh)/wv.z();
        if (a <= 0.f)
            return 0.f;
        float b = 1.f/(m_alpha*Frame::tanTheta(wv));
        if (b >= 1.6f)
            return 1.f;
        return (3.535f*b + 2.181f*b*b)/(1.f+2.276f*b+2.577f*b*b);
    }
    
    float m_alpha;
    float m_intIOR, m_extIOR;
    float m_ks;
    Color3f m_kd;
};

NORI_REGISTER_CLASS(Microfacet, "microfacet");
NORI_NAMESPACE_END
