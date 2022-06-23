#ifndef SHAPES_H
#define SHAPES_H

#include <optional>

#include <hitPosition.h>

struct shape {
    shape() {}
    shape(QVector3D center)
        : center(center)
    {
    }

    virtual std::optional<hit_record> hit(const Ray& r, double t_min, double t_max) const = 0;

    QVector3D center;
};

struct sphere : public shape {
    sphere(QVector3D center, float radius, std::shared_ptr<material> material)
        : shape(center)
        , radius(radius)
        , radiusSquared(radius * radius),
          material(material) {};

    // shpere hit code
    std::optional<hit_record> hit(const Ray& ray, double t_min,
        double t_max) const override
    {
        QVector3D originCenter = ray.pointOfOrigin - center;

        float a = ray.direction.lengthSquared();
        float b = QVector3D::dotProduct(originCenter, ray.direction);
        float c = originCenter.lengthSquared() - radiusSquared;

        float d = b * b - a * c; // sphere

        if (d < 0.0)
            return {};
        auto sqrtd = sqrt(d);
        auto root = (-b - sqrtd) / a;

        if (root < t_min || t_max < root) {
            root = (-b + sqrtd) / a;
            if (root < t_min || t_max < root)
                return {};
        }
        auto at = ray.at(root);
        QVector3D outwardNormal = (at - center) / radius;
        hit_record retVal;
        retVal.t = root;
        retVal.p = at;
        retVal.mat_ptr = material;
        //auto retVal = hitPosition { at, outwardNormal, root };
        retVal.setFaceNormal(ray, outwardNormal);
        return retVal;
    };

    float radius;
    float radiusSquared;
    std::shared_ptr<material> material;
};

struct rectangle : public shape {
    rectangle(QVector3D center, float w, float h, float d)
        : shape(center)
        , witdh(w)
        , height(h)
        , depth(d) {};

    //  hit code
    std::optional<hit_record> hit(const Ray& ray, double t_min,
        double t_max) const override
    {

        auto retVal = hit_record() ;

        return retVal;
    };

    float witdh;
    float height;
    float depth;
};

#endif // SHAPES_H
