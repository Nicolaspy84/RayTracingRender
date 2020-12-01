//
//  area.cpp
//  nori
//
//  Created by Nicolas Buisson on 24/03/2018.
//

#include <nori/emitter.h>

NORI_NAMESPACE_BEGIN

class AreaLight: public Emitter {
public:
    AreaLight(const PropertyList &props) {
        m_radiance = props.getColor("radiance");
    }
    
    Color3f emittedRadiance(Vector3f direction) const {
        // Check if the direction is positive
        if (direction.z() <= 0)
            return Color3f(0.f);
        
        // Return a constant (same radiance in any direction over the hemisphere)
        return m_radiance;
    }
    
    /// Return a human-readable description for debugging purposes
    std::string toString() const {
        return tfm::format(
                           "AreaLight[\n"
                           "  radiance = %s,\n"
                           "]",
                           m_radiance.toString()
                           );
    }
    
private:
    Color3f m_radiance;
};

NORI_REGISTER_CLASS(AreaLight, "area");
NORI_NAMESPACE_END
