#ifndef SHAPES_H
#define SHAPES_H

#include <optional>

#include <hitPosition.h>

namespace RectangleNormal {
constexpr QVector3D xy { 0, 0, 1 };
constexpr QVector3D xz { 0, 1, 0 };
constexpr QVector3D yz { 1, 0, 0 };
}

class aabb {
public:
    aabb() {};
    aabb(const QVector3D& a, const QVector3D& b)
    {
        minimum = a;
        maximum = b;
    }

    QVector3D min() const { return minimum; }
    QVector3D max() const { return maximum; }

    bool hit(const Ray& r, double t_min, double t_max) const
    {
        for (int a = 0; a < 3; a++) {
            auto invD = 1.0f / r.direction[a];
            auto t0 = (min()[a] - r.pointOfOrigin[a]) * invD;
            auto t1 = (max()[a] - r.pointOfOrigin[a]) * invD;
            if (invD < 0.0f)
                std::swap(t0, t1);
            t_min = t0 > t_min ? t0 : t_min;
            t_max = t1 < t_max ? t1 : t_max;
            if (t_max <= t_min)
                return false;
        }
        return true;
    }

    QVector3D minimum;
    QVector3D maximum;
};

struct shape {
    shape() {};
    shape(QVector3D center)
        : center(center)
    {
    }

    virtual std::optional<hit_record> hit(const Ray& r, double t_min, double t_max) const = 0;
    virtual bool bounding_box(aabb& output_box) const = 0;

    virtual double pdf_value(const QVector3D& o, const QVector3D& v) const {
        return 0.0;
    }

    virtual QVector3D random(const QVector3D& o) const {
        return QVector3D(1, 0, 0);
    }

    QVector3D center;
};

struct sphere : public shape {
    sphere(QVector3D center, float radius, std::shared_ptr<material> material)
        : shape { center }
        , radius { radius }
        , radded { radius, radius, radius }
        , sphere_output_box { center - radded, center + radded }
        , radiusSquared { radius * radius }
        , material { material } {};

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

        hit_record rec;
        rec.t = root;
        rec.p = at;
        rec.mat_ptr = material;
        get_sphere_uv(outwardNormal, rec.u, rec.v);
        //auto retVal = hitPosition { at, outwardNormal, root };
        rec.setFaceNormal(ray, outwardNormal);
        return rec;
    };

    virtual bool bounding_box(aabb& output_box) const override
    {
        output_box = sphere_output_box;
        return true;
    }

    double pdf_value(const QVector3D& o, const QVector3D& v) const {
        hit_record rec;
        auto hitRec = this->hit(Ray(o, v), 0.001, calc::infinity);
        if (!hitRec.has_value())
            return 0.0;

        auto cos_theta_max = sqrt(1 - radius*radius/(center-o).lengthSquared());
        auto solid_angle = 2*calc::pi*(1-cos_theta_max);

        return  1 / solid_angle;
    }

    QVector3D random(const QVector3D& o) const {
        QVector3D direction = center - o;
        auto distance_squared = direction.lengthSquared();
        onb uvw;
        uvw.build_from_w(direction);
        return uvw.local(calc::random_to_sphere(radius, distance_squared));
    }

    const float radius;
    const QVector3D radded;
    const aabb sphere_output_box;
    const float radiusSquared;
    const std::shared_ptr<material> material;

private:
    static void get_sphere_uv(const QVector3D& p, double& u, double& v)
    {
        // p: a given point on the sphere of radius one, centered at the origin.
        // u: returned value [0,1] of angle around the Y axis from X=-1.
        // v: returned value [0,1] of angle from Y=-1 to Y=+1.
        //     <1 0 0> yields <0.50 0.50>       <-1  0  0> yields <0.00 0.50>
        //     <0 1 0> yields <0.50 1.00>       < 0 -1  0> yields <0.50 0.00>
        //     <0 0 1> yields <0.25 0.50>       < 0  0 -1> yields <0.75 0.50>

        auto theta = acos(-p.y());
        auto phi = atan2(-p.z(), p.x()) + calc::pi;

        u = phi / (2.0f * calc::pi);
        v = theta / calc::pi;
    }
};

struct rectangle : public shape {
    rectangle(double _x0, double _x1, double _y0, double _y1, double _k,
        std::shared_ptr<material> mat)
        : mp(mat)
        , x0(_x0)
        , x1(_x1)
        , y0(_y0)
        , y1(_y1)
        , k(_k)
        , xy_rect_output_box(QVector3D(x0, y0, k - 0.0001), QVector3D(x1, y1, k + 0.0001)) {};

    virtual std::optional<hit_record> hit(const Ray& r, double t_min, double t_max) const override
    {
        auto t = (k - r.pointOfOrigin.z()) / r.direction.z();
        if (t < t_min || t > t_max)
            return {};
        auto x = r.pointOfOrigin.x() + t * r.direction.x();
        auto y = r.pointOfOrigin.y() + t * r.direction.y();
        if (x < x0 || x > x1 || y < y0 || y > y1)
            return {};

        hit_record rec;
        rec.u = (x - x0) / (x1 - x0);
        rec.v = (y - y0) / (y1 - y0);
        rec.t = t;

        auto outward_normal = RectangleNormal::xy;
        rec.setFaceNormal(r, outward_normal);
        rec.mat_ptr = mp;
        rec.p = r.at(t);
        return rec;
    }

    virtual bool bounding_box(aabb& output_box) const override
    {
        // The bounding box must have non-zero width in each dimension, so pad the Z
        // dimension a small amount.
        output_box = xy_rect_output_box;
        return true;
    }

    virtual double pdf_value(const QVector3D& origin, const QVector3D& v) const override {
        hit_record rec;
        auto hitRec = this->hit(Ray(origin, v), 0.001, calc::infinity);
        if (!hitRec.has_value())
            return {};
        rec = hitRec.value();

        auto area = (x1-x0)*(y1-y0);
        auto distance_squared = rec.t * rec.t * v.lengthSquared();
        auto cosine = fabs(QVector3D::dotProduct(v, rec.normal) / v.length());

        return distance_squared / (cosine * area);
    }

    virtual QVector3D random(const QVector3D& origin) const override {
        auto random_point = QVector3D(calc::random_double(x0,x1), k, calc::random_double(y0,y1));
        return random_point - origin;
    }

    const std::shared_ptr<material> mp;
    const double x0, x1, y0, y1, k;
    const aabb xy_rect_output_box;
};

class xz_rect : public shape {
public:
    xz_rect(double _x0, double _x1, double _z0, double _z1, double _k,
        std::shared_ptr<material> mat)
        : mp(mat)
        , x0(_x0)
        , x1(_x1)
        , z0(_z0)
        , z1(_z1)
        , k(_k)
        , xz_rect_output_box(QVector3D(x0, k - 0.0001f, z0), QVector3D(x1, k + 0.0001f, z1)) {};

    virtual std::optional<hit_record> hit(const Ray& r, double t_min, double t_max) const override
    {
        auto t = (k - r.pointOfOrigin.y()) / r.direction.y();
        if (t < t_min || t > t_max)
            return {};

        auto x = r.pointOfOrigin.x() + t * r.direction.x();
        auto z = r.pointOfOrigin.z() + t * r.direction.z();
        if (x < x0 || x > x1 || z < z0 || z > z1)
            return {};

        hit_record rec;
        rec.u = (x - x0) / (x1 - x0);
        rec.v = (z - z0) / (z1 - z0);
        rec.t = t;
        auto outward_normal = RectangleNormal::xz;
        rec.setFaceNormal(r, outward_normal);
        rec.mat_ptr = mp;
        rec.p = r.at(t);
        return rec;
    };

    virtual bool bounding_box(aabb& output_box) const override
    {
        // The bounding box must have non-zero width in each dimension, so pad the Y
        // dimension a small amount.
        output_box = xz_rect_output_box;
        return true;
    }

    virtual double pdf_value(const QVector3D& origin, const QVector3D& v) const override {
        hit_record rec;
        auto hitRec = this->hit(Ray(origin, v), 0.001, calc::infinity);
        if (!hitRec.has_value())
            return {};
        rec = hitRec.value();
        auto area = (x1-x0)*(z1-z0);
        auto distance_squared = rec.t * rec.t * v.lengthSquared();
        auto cosine = fabs(QVector3D::dotProduct(v, rec.normal) / v.length());

        return distance_squared / (cosine * area);
    }

    virtual QVector3D random(const QVector3D& origin) const override {
        auto random_point = QVector3D(calc::random_double(x0,x1), k, calc::random_double(z0,z1));
        return random_point - origin;
    }

public:
    const std::shared_ptr<material> mp;
    const double x0, x1, z0, z1, k;
    const aabb xz_rect_output_box;
};

class yz_rect : public shape {
public:
    yz_rect(double _y0, double _y1, double _z0, double _z1, double _k,
        std::shared_ptr<material> mat)
        : mp(mat)
        , y0(_y0)
        , y1(_y1)
        , z0(_z0)
        , z1(_z1)
        , k(_k)
        , yz_rect_output_box(QVector3D(k - 0.0001, y0, z0), QVector3D(k + 0.0001, y1, z1)) {};

    virtual std::optional<hit_record> hit(const Ray& r, double t_min, double t_max) const override
    {
        auto t = (k - r.pointOfOrigin.x()) / r.direction.x();
        if (t < t_min || t > t_max)
            return {};
        auto y = r.pointOfOrigin.y() + t * r.direction.y();
        auto z = r.pointOfOrigin.z() + t * r.direction.z();
        if (y < y0 || y > y1 || z < z0 || z > z1)
            return {};

        hit_record rec;
        rec.u = (y - y0) / (y1 - y0);
        rec.v = (z - z0) / (z1 - z0);
        rec.t = t;
        auto outward_normal = RectangleNormal::yz;
        rec.setFaceNormal(r, outward_normal);
        rec.mat_ptr = mp;
        rec.p = r.at(t);
        return rec;
    };

    virtual bool bounding_box(aabb& output_box) const override
    {
        // The bounding box must have non-zero width in each dimension, so pad the X
        // dimension a small amount.
        output_box = yz_rect_output_box;
        return true;
    }

public:
    const std::shared_ptr<material> mp;
    double y0, y1, z0, z1, k;
    const aabb yz_rect_output_box;
};

class box : public shape {
public:
    box() {}
    box(const QVector3D& p0, const QVector3D& p1, std::shared_ptr<material> ptr)
    {
        box_min = p0;
        box_max = p1;

        sides.push_back(std::make_shared<rectangle>(p0.x(), p1.x(), p0.y(), p1.y(), p1.z(), ptr));
        sides.push_back(std::make_shared<rectangle>(p0.x(), p1.x(), p0.y(), p1.y(), p0.z(), ptr));

        sides.push_back(std::make_shared<xz_rect>(p0.x(), p1.x(), p0.z(), p1.z(), p1.y(), ptr));
        sides.push_back(std::make_shared<xz_rect>(p0.x(), p1.x(), p0.z(), p1.z(), p0.y(), ptr));

        sides.push_back(std::make_shared<yz_rect>(p0.y(), p1.y(), p0.z(), p1.z(), p1.x(), ptr));
        sides.push_back(std::make_shared<yz_rect>(p0.y(), p1.y(), p0.z(), p1.z(), p0.x(), ptr));
    }

    virtual std::optional<hit_record> hit(const Ray& r, double t_min, double t_max) const override
    {
        return hitFromList(sides, r, t_min, t_max);
    }

    virtual bool bounding_box(aabb& output_box) const override
    {
        output_box = aabb(box_min, box_max);
        return boundingBoxFromList(sides, output_box);
    }

public:
    QVector3D box_min;
    QVector3D box_max;
    std::vector<std::shared_ptr<shape>> sides;
};

class constant_medium : public shape {
public:
    constant_medium(std::shared_ptr<shape> b, double d, std::shared_ptr<texture> a)
        : boundary(b)
        , neg_inv_density(-1 / d)
        , phase_function(std::make_shared<isotropic>(a))
    {
    }

    constant_medium(std::shared_ptr<shape> b, double d, QVector3D c)
        : boundary(b)
        , neg_inv_density(-1 / d)
        , phase_function(std::make_shared<isotropic>(c))
    {
    }

    virtual std::optional<hit_record> hit(
        const Ray& r, double t_min, double t_max) const override
    {
        // Print occasional samples when debugging. To enable, set enableDebug true.
        const bool enableDebug = false;
        const bool debugging = enableDebug && calc::random_double01() < 0.00001;

        hit_record rec, rec1, rec2;

        auto didHit1 = boundary->hit(r, -calc::infinity, calc::infinity);
        if (!didHit1.has_value())
            return {};
        rec1 = didHit1.value();

        auto didHit2 = boundary->hit(r, rec1.t + 0.0001, calc::infinity);
        if (!didHit2.has_value())
            return {};
        rec2 = didHit2.value();

        // if (debugging)
        //     qCritical() << "\nt_min=" << rec1.t << ", t_max=" << rec2.t << '\n';

        if (rec1.t < t_min)
            rec1.t = t_min;
        if (rec2.t > t_max)
            rec2.t = t_max;

        if (rec1.t >= rec2.t)
            return {};

        if (rec1.t < 0)
            rec1.t = 0;

        const auto ray_length = r.direction.length();
        const auto distance_inside_boundary = (rec2.t - rec1.t) * ray_length;
        const auto hit_distance = neg_inv_density * log(calc::random_double01());

        if (hit_distance > distance_inside_boundary)
            return {};

        rec.t = rec1.t + hit_distance / ray_length;
        rec.p = r.at(rec.t);

        /* if (debugging) {
            std::cerr << "hit_distance = " << hit_distance << '\n'
                      << "rec.t = " << rec.t << '\n'
                      << "rec.p = " << rec.p << '\n';
        }*/

        rec.normal = QVector3D(1, 0, 0); // arbitrary
        rec.front_face = true; // also arbitrary
        rec.mat_ptr = phase_function;

        return rec;
    };

    virtual bool bounding_box(aabb& output_box) const override
    {
        return boundary->bounding_box(output_box);
    }

public:
    std::shared_ptr<shape> boundary;
    std::shared_ptr<material> phase_function;
    double neg_inv_density;
};


#endif // SHAPES_H
