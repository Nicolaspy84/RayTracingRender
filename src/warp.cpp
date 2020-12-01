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

#include <nori/warp.h>
#include <nori/vector.h>
#include <nori/frame.h>
#include <math.h>

NORI_NAMESPACE_BEGIN

Point2f Warp::squareToUniformSquare(const Point2f &sample) {
    return sample;
}

float Warp::squareToUniformSquarePdf(const Point2f &sample) {
    return ((sample.array() >= 0).all() && (sample.array() <= 1).all()) ? 1.0f : 0.0f;
}

Point2f Warp::squareToTent(const Point2f &sample) {
    Point2f dist;
    
    // Different cases are possible (see the mathematical justification in report.html 3)
    
    if (sample.x() < 0.5f)
        dist.x() = sqrt(2*sample.x()) - 1.f;
    else if (sample.x() >= 0.5f)
        dist.x() = 1.f - sqrt(2.f*(1.f-sample.x()));
    
    if (sample.y() < 0.5f)
        dist.y() = sqrt(2.f*sample.y()) - 1.f;
    else if (sample.y() >= 0.5f)
        dist.y() = 1.f - sqrt(2.f*(1.f-sample.y()));
    
    return dist;
}

float Warp::squareToTentPdf(const Point2f &p) {
    float px = 0.f;
    float py = 0.f;
    
    if (p.x() >= -1.f && p.x() <= 1.f)
        px = 1.f - abs(p.x());
    if (p.y() >= -1.f && p.y() <= 1.f)
        py = 1.f - abs(p.y());
    
    return px*py;
}

Point2f Warp::squareToUniformDisk(const Point2f &sample) {
    Point2f dist;
    float sqrtX = sqrt(sample.x());
    float piY = 2.f * M_PI * sample.y();
    
    dist.x() = sqrtX * cos(piY);
    dist.y() = sqrtX * sin(piY);
    
    return dist;
}

float Warp::squareToUniformDiskPdf(const Point2f &p) {
    if (pow(p.x(), 2)+pow(p.y(), 2) <= 1.f)
        return 1.f/(M_PI);
    return 0.f;
}

Vector3f Warp::squareToUniformSphere(const Point2f &sample) {
    Vector3f dist;
    
    float cosTheta = 1.f - 2.f*sample.x();
    float sinTheta = 2.f * sqrt(sample.x() - pow(sample.x(), 2));
    float phi = sample.y() * 2 * M_PI;
    float sinPhi = sin(phi);
    float cosPhi = cos(phi);
    dist.x() = sinTheta * cosPhi;
    dist.y() = sinTheta * sinPhi;
    dist.z() = cosTheta;
    return dist;
}

float Warp::squareToUniformSpherePdf(const Vector3f &v) {
    float eps = 0.000001f;
    float val = pow(v.x(), 2) + pow(v.y(), 2) + pow(v.z(), 2);
    if (val < 1.f + eps && val > 1.f - eps)
        return 1.f/(4.f*M_PI);
    return 0.f;
}

Vector3f Warp::squareToUniformHemisphere(const Point2f &sample) {
    Vector3f dist;
    dist = squareToUniformSphere(sample);
    // We simply do the same as for a sphere, but only sample positive z points)
    if (dist.z() < 0.f)
        dist.z() = -dist.z();
    return dist;
}

float Warp::squareToUniformHemispherePdf(const Vector3f &v) {
    float eps = 0.000001f;
    float val = pow(v.x(), 2) + pow(v.y(), 2) + pow(v.z(), 2);
    if (val < 1.f + eps && val > 1.f - eps && v.z() >= 0.f)
        return 1.f/(2.f*M_PI);
    return 0.f;
}

Vector3f Warp::squareToCosineHemisphere(const Point2f &sample) {
    Vector3f dist;
    Point2f circleDist = squareToUniformDisk(sample);
    // We directly sample points uniformly on the disk, and then project this on the sphere
    dist.x() = circleDist.x();
    dist.y() = circleDist.y();
    dist.z() = sqrt(1.f - pow(dist.x(), 2) - pow(dist.y(), 2));
    return dist;
}

float Warp::squareToCosineHemispherePdf(const Vector3f &v) {
    float eps = 0.000001f;
    float cosTheta = v.z();
    float val = pow(v.x(), 2) + pow(v.y(), 2) + pow(v.z(), 2);
    if (val < 1.f + eps && val > 1.f - eps && v.z() >= 0.f)
        return cosTheta/(M_PI);
    return 0.f;
}

Vector3f Warp::squareToBeckmann(const Point2f &sample, float alpha) {
    Vector3f dist;
    float cosTheta;
    if (sample.x() == 1.f)
        cosTheta = 0.f;
    cosTheta = 1.f/sqrt(1.f-pow(alpha, 2)*log(1.f-sample.x()));
    float sinTheta = sqrt(1.f - pow(cosTheta, 2));
    float phi = sample.y()*2*M_PI;
    dist.x() = sinTheta * cos(phi);
    dist.y() = sinTheta * sin(phi);
    dist.z() = cosTheta;
    return dist;
}

float Warp::squareToBeckmannPdf(const Vector3f &m, float alpha) {
    float eps = 0.000001f;
    float val = pow(m.x(), 2) + pow(m.y(), 2) + pow(m.z(), 2);
    if (val < 1.f + eps && val > 1.f - eps && m.z() > 0.f)
        return exp((1.f - 1.f/pow(m.z(), 2)) / pow(alpha, 2)) / (M_PI * pow(alpha, 2) * pow(m.z(), 3));
    return 0.f;
}

NORI_NAMESPACE_END
