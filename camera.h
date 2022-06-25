#ifndef CAMERA_H
#define CAMERA_H

#include <Ray.h>
#include <calc.h>

class camera {

public:
    camera(
        QVector3D lookfrom,
        QVector3D lookat,
        QVector3D vup,
        double vfov = 90.0, // vertical field-of-view in degrees
        double aspect_ratio = img::aspect_ratio)
    {
        resetCam(lookfrom, lookat, vup, vfov, aspect_ratio);
    }

    Ray get_ray(double s, double t) const
    {
        return Ray(origin, lower_left_corner + s * horizontal + t * vertical - origin);
    }

    void resetCam(QVector3D lookfrom,
        QVector3D lookat,
        QVector3D vup,
        double vfov = 90.0,
        double aspect_ratio = img::aspect_ratio)
    {
        auto theta = calc::degrees_to_radians(vfov);
        auto h = tan(theta / 2.0);
        auto viewport_height = 2.0 * h;
        auto viewport_width = aspect_ratio * viewport_height;
        auto w = QVector3D(lookfrom - lookat).normalized();
        auto u = QVector3D(QVector3D::crossProduct(vup, w)).normalized();
        auto v = QVector3D::crossProduct(w, u);

        origin = lookfrom;
        horizontal = viewport_width * u;
        vertical = viewport_height * v;
        lower_left_corner = origin - horizontal / 2.0 - vertical / 2.0 - w;
    }

private:
    QVector3D origin;
    QVector3D lower_left_corner;
    QVector3D horizontal;
    QVector3D vertical;
};

#endif // CAMERA_H
