#ifndef CAMERA_H
#define CAMERA_H

#include <Ray.h>

class camera {
public:
    camera()
    {
        auto aspect_ratio = 16.0 / 9.0;
        auto viewport_height = 2.0;
        auto viewport_width = aspect_ratio * viewport_height;
        float focal_length = 1.0;

        origin = QVector3D(0, 0, 0);
        horizontal = QVector3D(viewport_width, 0.0, 0.0);
        vertical = QVector3D(0.0, viewport_height, 0.0);
        lower_left_corner = origin - horizontal / 2 - vertical / 2 - QVector3D { 0, 0, focal_length };
    }

    Ray get_ray(double u, double v) const
    {
        return Ray(origin,
            lower_left_corner + u * horizontal + v * vertical - origin);
    }

private:
    QVector3D origin;
    QVector3D lower_left_corner;
    QVector3D horizontal;
    QVector3D vertical;
};

#endif // CAMERA_H
