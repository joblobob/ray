#ifndef HITPOSITION_H
#define HITPOSITION_H

#include <Ray.h>
#include <calc.h>
#include <optional>

struct shape;

class aabb;

struct hit_record;

class texture {
public:
    virtual QVector3D value(double u, double v, const QVector3D& p) const = 0;
};

class solid_color : public texture {
public:
    solid_color() {}
    solid_color(QVector3D c)
        : color_value(c)
    {
    }

    solid_color(double red, double green, double blue)
        : solid_color(QVector3D(red, green, blue))
    {
    }

    virtual QVector3D value(double u, double v, const QVector3D& p) const override
    {
        return color_value;
    }

private:
    QVector3D color_value;
};

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
    virtual QVector3D emitted(double u, double v, const QVector3D& p) const
    {
        return img::defaultVec;
    }
};

struct hit_record {
    QVector3D p;
    QVector3D normal;
    std::shared_ptr<material> mat_ptr;
    double t;
    double u;
    double v;
    bool front_face;

    inline void setFaceNormal(const Ray& r, const QVector3D& outward_normal)
    {
        front_face = QVector3D::dotProduct(r.direction, outward_normal) < 0.0f;
        normal = front_face ? outward_normal : -outward_normal;
    }
};

std::optional<hit_record>
hitFromList(const std::vector<std::shared_ptr<shape>>& sphereList,
    const Ray& ray,
    double t_min,
    double t_max);

bool boundingBoxFromList(const std::vector<std::shared_ptr<shape>>& list, aabb& output_box);

aabb surrounding_box(aabb box0, aabb box1);

class lambertian : public material {
public:
    lambertian(const QVector3D& a)
        : albedo(std::make_shared<solid_color>(a))
    {
    }
    lambertian(std::shared_ptr<texture> a)
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
        attenuation = albedo->value(rec.u, rec.v, rec.p);
        return true;
    }

public:
    std::shared_ptr<texture> albedo;
};

class metal : public material {
public:
    metal(const QVector3D& a, double f)
        : albedo(a)
        , fuzz(f < 1.0f ? f : 1.0f)
    {
    }

    virtual bool scatter(
        const Ray& r_in, const hit_record& rec, QVector3D& attenuation, Ray& scattered) const override
    {
        QVector3D reflected = calc::reflect(r_in.direction.normalized(), rec.normal);
        scattered = Ray(rec.p, reflected + fuzz * calc::random_in_unit_sphere());
        attenuation = albedo;
        return (QVector3D::dotProduct(scattered.direction, rec.normal) > 0.0f);
    }

public:
    QVector3D albedo;
    float fuzz;
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
        float cos_theta = fmin(QVector3D::dotProduct(-unit_direction, rec.normal), 1.0f);
        float sin_theta = sqrt(1.0f - cos_theta * cos_theta);

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
    static float reflectance(float cosine, float ref_idx)
    {
        // Use Schlick's approximation for reflectance.
        auto r0 = (1.0f - ref_idx) / (1.0f + ref_idx);
        r0 = r0 * r0;
        return r0 + (1.0f - r0) * pow((1.0f - cosine), 5);
    }
};

class diffuse_light : public material {
public:
    diffuse_light(std::shared_ptr<texture> a)
        : emiter(a)
    {
    }
    diffuse_light(QVector3D c)
        : emiter(std::make_shared<solid_color>(c))
    {
    }

    virtual bool scatter(
        const Ray& r_in, const hit_record& rec, QVector3D& attenuation, Ray& scattered) const override
    {
        return false;
    }

    virtual QVector3D emitted(double u, double v, const QVector3D& p) const override
    {
        return emiter->value(u, v, p);
    }

public:
    std::shared_ptr<texture> emiter;
};

class isotropic : public material {
public:
    isotropic(QVector3D c)
        : albedo(std::make_shared<solid_color>(c))
    {
    }
    isotropic(std::shared_ptr<texture> a)
        : albedo(a)
    {
    }

    virtual bool scatter(
        const Ray& r_in, const hit_record& rec, QVector3D& attenuation, Ray& scattered) const override
    {
        scattered = Ray(rec.p, calc::random_in_unit_sphere());
        attenuation = albedo->value(rec.u, rec.v, rec.p);
        return true;
    }

public:
    std::shared_ptr<texture> albedo;
};

#endif // HITPOSITION_H
