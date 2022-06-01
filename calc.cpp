
#include <calc.h>

namespace calc {

constexpr double degrees_to_radians(double degrees)
{
    return degrees * pi / 180.0;
}

float random_double(float min, float max)
{
    std::random_device
        rd; // Will be used to obtain a seed for the random number engine
    /*std::mt19937 gen(rd());*/ // Standard mersenne_twister_engine seeded with
    // rd()
    // Changing mt19937 to minstd_rand makes the code run 88 times faster!
    std::minstd_rand gen(rd());
    std::uniform_real_distribution<> distribution(min, max);
    return distribution(gen);
}

QVector3D randomVec()
{
    return { random_double(), random_double(), random_double() };
}

QVector3D randomVec(float min, float max)
{
    return { random_double(min, max), random_double(min, max),
        random_double(min, max) };
}

QVector3D
random_in_unit_sphere()
{
    QVector3D point { -1.0f, -1.0f, -1.0f };
    while (point.lengthSquared() >= 1.0f) {
        point = calc::randomVec(-1.0f, 1.0f);
    }
    return point;
}

QVector3D
random_unit_vector()
{
    return random_in_unit_sphere().normalized();
}

QVector3D random_in_hemisphere(const QVector3D& normal)
{
    QVector3D in_unit_sphere = random_in_unit_sphere();
    if (QVector3D::dotProduct(in_unit_sphere, normal) > 0.0f) // In the same hemisphere as the normal
        return in_unit_sphere;
    else
        return -in_unit_sphere;
}

// couleur
QVector3D
ray_color(const Ray& r,
    const std::vector<std::unique_ptr<shape>>& worldObjects,
    int depth, bool drawOnlyColors)
{
    // If we've exceeded the ray bounce limit, no more light is gathered.
    if (depth <= 0.0f)
        return img::defaultVec;

    std::optional<hitPosition> hitRecord = calc::hitFromList(worldObjects, r, 0.001f, calc::infinity);
    if (!hitRecord.has_value()) {
        // pas de hit, background gradient
        QVector3D unit_direction = r.direction.normalized();
        float t = 0.5f * (unit_direction.y() + 1.0f);
        return (1.0f - t) * img::bgColor + t * img::gradientBgVec;
    }

    if (drawOnlyColors) {
        // juste couleur
        QVector3D N = QVector3D { r.at(hitRecord.value().t) - img::infiniteZ }.normalized();
        return 0.5f * QVector3D { N.x() + 1.0f, N.y() + 1.0f, N.z() + 1.0f };
    }

    // recursive bounce
    QVector3D target = hitRecord->point + calc::random_in_hemisphere(hitRecord->normal);
    return 0.5f * calc::ray_color({ hitRecord->point, target - hitRecord->point }, worldObjects, depth - 1, drawOnlyColors);
}

std::optional<hitPosition>
hitFromList(const std::vector<std::unique_ptr<shape>>& sphereList,
    const Ray& ray,
    double t_min,
    double t_max)
{
    bool hit_anything = false;
    auto closest_so_far = t_max;
    std::optional<hitPosition> retVal;

    for (const auto& object : sphereList) {
        std::optional<hitPosition> didHit = object->hit(ray, t_min, closest_so_far);
        if (didHit.has_value()) {
            hit_anything = true;
            closest_so_far = didHit->t;
            retVal = didHit.value();
        }
    }

    return retVal;
}

} // namespace calc
