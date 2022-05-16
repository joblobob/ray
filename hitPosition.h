#ifndef HITPOSITION_H
#define HITPOSITION_H

#include <Ray.h>

struct hitPosition {
    hitPosition(const QVector3D& point, const QVector3D& normal, double t)
        : point(point)
        , normal(normal)
        , colorAttenuation()
        , scatteredRay()
        , t(t) {};
    QVector3D point;
    QVector3D normal;

    QVector3D colorAttenuation;
    Ray scatteredRay;

    double t;
    bool frontFace;
    inline void setFaceNormal(const Ray& ray, const QVector3D& outwardNormal)
    {
        frontFace = QVector3D::dotProduct(ray.direction, outwardNormal) < 0.0f;
        normal = frontFace ? outwardNormal : -outwardNormal;
    }
};

#endif // HITPOSITION_H
