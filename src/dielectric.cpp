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
#include <nori/common.h>

NORI_NAMESPACE_BEGIN

/// Ideal dielectric BSDF
class Dielectric : public BSDF {
public:
    Dielectric(const PropertyList &propList) {
        /* Interior IOR (default: BK7 borosilicate optical glass) */
        m_intIOR = propList.getFloat("intIOR", 1.5046f);

        /* Exterior IOR (default: air) */
        m_extIOR = propList.getFloat("extIOR", 1.000277f);
    }

    Color3f eval(const BSDFQueryRecord &) const {
        /* Discrete BRDFs always evaluate to zero in Nori */
        return Color3f(0.0f);
    }

    float pdf(const BSDFQueryRecord &) const {
        /* Discrete BRDFs always evaluate to zero in Nori */
        return 0.0f;
    }

    Color3f sample(BSDFQueryRecord &bRec, const Point2f &sample) const {
        
        float cosTheta = Frame::cosTheta(bRec.wi);
        float n_out = m_extIOR;
        float n_in = m_intIOR;
        
        // If the ray comes from inside
        if (Frame::cosTheta(bRec.wi) <= 0.f) {
            cosTheta = - Frame::cosTheta(bRec.wi);
            n_out = m_intIOR;
            n_in = m_extIOR;
        }
                
        // Reflected amount
        float Fr = fresnel(cosTheta, n_out, n_in);
        
        // Reflect the ray
        if (sample.x() <= Fr) {
            // Reflection in local coordinates
            bRec.wo = Vector3f(
                               -bRec.wi.x(),
                               -bRec.wi.y(),
                               bRec.wi.z()
                               );
            bRec.measure = EDiscrete;
        }
        
        // Refract the ray
        else {
            float n1n2 = n_out/n_in;
            
            float z;
            if (Frame::cosTheta(bRec.wi) <= 0.f)
                z = sqrt(1.f - pow(n1n2, 2)*(1.f - pow(bRec.wi.z(), 2)));
            else
                z = -sqrt(1.f - pow(n1n2, 2)*(1.f - pow(bRec.wi.z(), 2)));
            
            bRec.wo = Vector3f(
                               -n1n2*bRec.wi.x(),
                               -n1n2*bRec.wi.y(),
                               z
                               );
        }
        
        return Color3f(1.f);
    }

    std::string toString() const {
        return tfm::format(
            "Dielectric[\n"
            "  intIOR = %f,\n"
            "  extIOR = %f\n"
            "]",
            m_intIOR, m_extIOR);
    }
private:
    float m_intIOR, m_extIOR;
};

NORI_REGISTER_CLASS(Dielectric, "dielectric");
NORI_NAMESPACE_END
