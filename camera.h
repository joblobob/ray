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
        double focus_dist,
        double vfov = 40.0, // vertical field-of-view in degrees
        double aspect_ratio = img::aspect_ratio, double aperture = 0.1)
    {
        resetCam(lookfrom, lookat, vup, focus_dist, vfov, aspect_ratio, aperture);
    }

    Ray get_ray(double s, double t) const
    {
        QVector3D rd = lens_radius * calc::random_in_unit_disk();
        QVector3D offset = u * rd.x() + v * rd.y();

        return Ray(
            origin + offset,
            lower_left_corner + s * horizontal + t * vertical - origin - offset);
        return Ray(origin, lower_left_corner + s * horizontal + t * vertical - origin);
    }

    void resetCam(QVector3D lookfrom,
        QVector3D lookat,
        QVector3D vup,
        double focus_dist,
        double vfov = 40.0,
        double aspect_ratio = img::aspect_ratio, double aperture = 0.1)
    {
        auto theta = calc::degrees_to_radians(vfov);
        auto h = tan(theta / 2.0);
        auto viewport_height = 2.0 * h;
        auto viewport_width = aspect_ratio * viewport_height;

        w = QVector3D(lookfrom - lookat).normalized();
        u = QVector3D(QVector3D::crossProduct(vup, w)).normalized();
        v = QVector3D::crossProduct(w, u);

        origin = lookfrom;
        horizontal = focus_dist * viewport_width * u;
        vertical = focus_dist * viewport_height * v;
        lower_left_corner = origin - horizontal / 2.0f - vertical / 2.0f - focus_dist * w;

        lens_radius = aperture / 2.0f;
    }

private:
    QVector3D origin;
    QVector3D lower_left_corner;
    QVector3D horizontal;
    QVector3D vertical;

    QVector3D w, u, v;
    double lens_radius;
};

#endif // CAMERA_H
