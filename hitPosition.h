#ifndef HITPOSITION_H
#define HITPOSITION_H

#include <Ray.h>
#include <calc.h>
#include <optional>

struct shape;

struct hit_record;
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
        const Ray& r_in, const hit_record& rec, QVector3D& attenuation, Ray& scattered) const = 0;
};

struct hit_record {
    QVector3D p;
    QVector3D normal;
    std::shared_ptr<material> mat_ptr;
    double t;
    bool front_face;

    inline void setFaceNormal(const Ray& r, const QVector3D& outward_normal)
    {
        front_face = QVector3D::dotProduct(r.direction, outward_normal) < 0;
        normal = front_face ? outward_normal : -outward_normal;
    }
};

std::optional<hit_record>
hitFromList(const std::vector<std::unique_ptr<shape>>& sphereList,
    const Ray& ray,
    double t_min,
    double t_max);

class lambertian : public material {
public:
    lambertian(const QVector3D& a)
        : albedo(a)
    {
    }

    virtual bool scatter(
        const Ray& r_in, const hit_record& rec, QVector3D& attenuation, Ray& scattered) const override
    {
        auto scatter_direction = rec.normal + calc::random_unit_vector();
        if (calc::near_zero(scatter_direction))
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
    metal(const QVector3D& a, double f)
        : albedo(a)
        , fuzz(f < 1.0 ? f : 1.0)
    {
    }

    virtual bool scatter(
        const Ray& r_in, const hit_record& rec, QVector3D& attenuation, Ray& scattered) const override
    {
        QVector3D reflected = calc::reflect(r_in.direction.normalized(), rec.normal);
        scattered = Ray(rec.p, reflected + fuzz * calc::random_in_unit_sphere());
        attenuation = albedo;
        return (QVector3D::dotProduct(scattered.direction, rec.normal) > 0);
    }

public:
    QVector3D albedo;
    double fuzz;
};

class dielectric : public material {
public:
    dielectric(double index_of_refraction)
        : ir(index_of_refraction)
    {
    }

    virtual bool scatter(
        const Ray& r_in, const hit_record& rec, QVector3D& attenuation, Ray& scattered) const override
    {
        attenuation = img::bgColor;
        double refraction_ratio = rec.front_face ? (1.0f / ir) : ir;

        QVector3D unit_direction = r_in.direction.normalized();
        double cos_theta = fmin(QVector3D::dotProduct(-unit_direction, rec.normal), 1.0f);
        double sin_theta = sqrt(1.0f - cos_theta * cos_theta);

        bool cannot_refract = refraction_ratio * sin_theta > 1.0f;
        QVector3D direction;

        if (cannot_refract || reflectance(cos_theta, refraction_ratio) > calc::random_double01())
            direction = calc::reflect(unit_direction, rec.normal);
        else
            direction = calc::refract(unit_direction, rec.normal, refraction_ratio);

        scattered = Ray(rec.p, direction);
        return true;
    }

public:
    double ir; // Index of Refraction

private:
    static double reflectance(double cosine, double ref_idx)
    {
        // Use Schlick's approximation for reflectance.
        auto r0 = (1 - ref_idx) / (1 + ref_idx);
        r0 = r0 * r0;
        return r0 + (1 - r0) * pow((1 - cosine), 5);
    }
};

#endif // HITPOSITION_H
