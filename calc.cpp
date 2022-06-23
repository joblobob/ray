
#include <calc.h>

#include <shapes.h>
#include <hitPosition.h>

namespace calc {

constexpr double degrees_to_radians(double degrees)
{
    return degrees * pi / 180.0;
}

float random_double(float min, float max)
{
    std::uniform_real_distribution<float> distribution(min, max);
    return distribution(gen);
}

float random_double01()
{
    return simple_distribution(gen01);
}

QVector3D randomVec()
{
    return { random_double01(), random_double01(), random_double01() };
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
ray_color(const Ray& inboundRay,
    const std::vector<std::unique_ptr<shape>>& worldObjects,
    int depth, bool drawOnlyColors)
{
    // If we've exceeded the ray bounce limit, no more light is gathered.
    if (depth <= 0.0f)
        return img::defaultVec;

    std::optional<hit_record> hitRecord = hitFromList(worldObjects, inboundRay, 0.001f, calc::infinity);
    if (!hitRecord.has_value()) {
        // pas de hit, background gradient
        QVector3D unit_direction = inboundRay.direction.normalized();
        float t = 0.5f * (unit_direction.y() + 1.0f);
        return (1.0f - t) * img::bgColor + t * img::gradientBgVec;
    }

    if (drawOnlyColors) {
        // juste couleur
        QVector3D N = QVector3D { inboundRay.at(hitRecord.value().t) - img::infiniteZ }.normalized();
        return 0.5f * QVector3D { N.x() + 1.0f, N.y() + 1.0f, N.z() + 1.0f };
    }

    // recursive bounce
    Ray scattered;
    QVector3D attenuation;
    if (hitRecord.value().mat_ptr->scatter(inboundRay, hitRecord.value(), attenuation, scattered))
        return attenuation * ray_color(scattered, worldObjects, depth-1, drawOnlyColors);
    return img::defaultVec;
}

QVector3D rgbPerSamples(const QVector3D& pixel, int samples)
{
    // Divide the color by the number of samples. and gamma-correct for gamma=2.0.
    auto scale = 1.0f / samples;
    auto r = sqrt(scale * pixel.x());
    auto g = sqrt(scale * pixel.y());
    auto b = sqrt(scale * pixel.z());

    return { 255.0f * std::clamp(r, 0.0f, 0.999f), 255.0f * std::clamp(g, 0.0f, 0.999f), 255.0f * std::clamp(b, 0.0f, 0.999f) };
}

bool near_zero(const QVector3D& vec)  {
        // Return true if the vector is close to zero in all dimensions.
        const auto s = 1e-8;
        return (fabs(vec[0]) < s) && (fabs(vec[1]) < s) && (fabs(vec[2]) < s);
}

QVector3D reflect(const QVector3D &v, const QVector3D &n) {
    return v - 2*QVector3D::dotProduct(v,n)*n;
}

} // namespace calc
