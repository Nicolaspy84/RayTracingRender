//
//  merl.cpp
//  nori
//
//  Created by Nicolas Buisson on 03/05/2018.
//

#include <nori/bsdf.h>
#include <nori/warp.h>
#include <nori/bitmap.h>
#include <filesystem/resolver.h>
#include "merlLoader.cpp"

NORI_NAMESPACE_BEGIN

class Merl : public BSDF {
public:
    Merl(const PropertyList &propList) {
        /* File name */
        m_fileName = getFileResolver()->resolve(propList.getString("filename", "")).str();
        m_fileName = getFileResolver()->resolve(propList.getString("name", "")).str();
        if (m_fileName != getFileResolver()->resolve("").str())
            read_brdf(m_fileName.c_str(), m_binary);
    }
    
    /// Evaluate the BRDF for the given pair of directions
    Color3f eval(const BSDFQueryRecord &bRec) const {
        if (bRec.measure != ESolidAngle
            || Frame::cosTheta(bRec.wi) <= 0
            || Frame::cosTheta(bRec.wo) <= 0)
            return Color3f(0.0f);
        
        double red, green, blue;
        lookup_brdf_val(m_binary, bRec.wi, bRec.wo, red, green, blue);
        
        return Color3f(red, green, blue);
    }
    
    /// Evaluate the sampling density of \ref sample() wrt. solid angles
    float pdf(const BSDFQueryRecord &bRec) const {
        
        /* This is a smooth BRDF -- return zero if the measure
         is wrong, or when queried for illumination on the backside */
        if (bRec.measure != ESolidAngle
            || Frame::cosTheta(bRec.wi) <= 0
            || Frame::cosTheta(bRec.wo) <= 0)
            return 0.0f;
        
        return INV_PI * Frame::cosTheta(bRec.wi);
    }
    
    /// Sample the BRDF
    Color3f sample(BSDFQueryRecord &bRec, const Point2f &sample) const {
        if (Frame::cosTheta(bRec.wi) <= 0)
            return Color3f(0.0f);
        
        bRec.measure = ESolidAngle;
        
        /* Warp a uniformly distributed sample on [0,1]^2
         to a direction on a cosine-weighted hemisphere */
        bRec.wo = Warp::squareToCosineHemisphere(sample);
        
        /* Relative index of refraction: no change */
        bRec.eta = 1.0f;
        
        /* eval() / pdf() * cos(theta) = albedo. */
        return eval(bRec) * Frame::cosTheta(bRec.wi) / pdf(bRec);
    }
    
    bool isDiffuse() const {
        return true;
    }
    
    std::string toString() const {
        return tfm::format(
                           "Merl[\n"
                           "  fileName = %s,\n"
                           "]",
                           m_fileName
                           );
    }
    
private:
    
    std::string m_fileName;
    double *m_binary = nullptr;
};

NORI_REGISTER_CLASS(Merl, "merl");
NORI_NAMESPACE_END
