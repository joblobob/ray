#ifndef RAY_H
#define RAY_H

#include <QVector3D>

struct Ray {
    Ray() {}
	Ray(const QVector3D& origin, const QVector3D& direction) : pointOfOrigin(origin), direction(direction), DirlengthSquared(direction.lengthSquared())
    {
    }
    QVector3D at(double t) const { return pointOfOrigin + t * direction; }
    QVector3D pointOfOrigin;
    QVector3D direction;

    float DirlengthSquared;
};

#endif // RAY_H
