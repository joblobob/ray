#ifndef HITPOSITION_H
#define HITPOSITION_H

#include <calc.h>
#include <Ray.h>

struct shape;

struct hit_record ;
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

class material {
    public:
        virtual bool scatter(
            const Ray& r_in, const hit_record& rec, QVector3D& attenuation, Ray& scattered
        ) const = 0;
};

struct hit_record {
    QVector3D p;
    QVector3D normal;
    std::shared_ptr<material> mat_ptr;
    double t;
    bool front_face;

    inline void setFaceNormal(const Ray& r, const QVector3D& outward_normal) {
        front_face = QVector3D::dotProduct(r.direction, outward_normal) < 0;
        normal = front_face ? outward_normal :-outward_normal;
    }
};


std::optional<hit_record>
hitFromList(const std::vector<std::unique_ptr<shape>>& sphereList,
    const Ray& ray,
    double t_min,
    double t_max);


class lambertian : public material {
    public:
        lambertian(const QVector3D& a) : albedo(a) {}

        virtual bool scatter(
            const Ray& r_in, const hit_record& rec, QVector3D& attenuation, Ray& scattered
        ) const override {
            auto scatter_direction = rec.normal + calc::random_unit_vector();
            if(calc::near_zero(scatter_direction))
                scatter_direction = rec.normal;
            scattered = Ray(rec.p, scatter_direction);
            attenuation = albedo;
            return true;
        }

    public:
        QVector3D albedo;
};

class metal : public material {
    public:
        metal(const QVector3D& a) : albedo(a) {}

        virtual bool scatter(
            const Ray& r_in, const hit_record& rec, QVector3D& attenuation, Ray& scattered
        ) const override {
            QVector3D reflected = calc::reflect(r_in.direction.normalized(), rec.normal);
            scattered = Ray(rec.p, reflected);
            attenuation = albedo;
            return (QVector3D::dotProduct(scattered.direction, rec.normal) > 0);
        }

    public:
        QVector3D albedo;
};


#endif // HITPOSITION_H
