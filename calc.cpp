
#include <calc.h>

#include <hitPosition.h>
#include <shapes.h>

namespace calc {

float random_double(float min, float max)
{
    std::uniform_real_distribution<float> distribution(min, max);
    return distribution(gen);
}

float random_double01()
{
    return simple_distribution(gen01);
}

float random_double11()
{
    return simple_distributionMinusOneToOne(gen11);
}

QVector3D randomVec()
{
    return { random_double01(), random_double01(), random_double01() };
}

QVector3D randomVec11()
{
    return { random_double11(), random_double11(), random_double11() };
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
        point = calc::randomVec11();
    }
    return point;
}

QVector3D
random_in_unit_disk()
{
    auto p = QVector3D(random_double11(), random_double11(), 0.0f);
    while (p.lengthSquared() >= 1.0f) {
        p = QVector3D(random_double11(), random_double11(), 0.0f);
    }
    return p;
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
ray_color(const Ray& inboundRay, const QVector3D& background,
          const std::vector<std::shared_ptr<shape>>& worldObjects, std::shared_ptr<shape>& lights,
    int depth, bool drawOnlyColors)
{
    // If we've exceeded the ray bounce limit, no more light is gathered.
    if (depth <= 0.0f)
        return img::defaultVec;

    std::optional<hit_record> hitRecord = hitFromList(worldObjects, inboundRay, calc::smallestVal, calc::infinity);
    if (!hitRecord.has_value()) {
        /*// pas de hit, background gradient
        QVector3D unit_direction = inboundRay.direction.normalized();
        float t = 0.5f * (unit_direction.y() + 1.0f);
        return (1.0f - t) * img::bgColor + t * img::gradientBgVec;*/
        // pas de hit, background
        return background;
    }

    if (drawOnlyColors) {
        // juste couleur
        QVector3D N = QVector3D { inboundRay.at(hitRecord.value().t) - img::infiniteZ }.normalized();
        return 0.5f * QVector3D { N.x() + 1.0f, N.y() + 1.0f, N.z() + 1.0f };
    }

    auto rec = hitRecord.value();

    scatter_record srec;
    QVector3D emitted = rec.mat_ptr->emitted(inboundRay, rec, rec.u, rec.v, rec.p);
    if (!rec.mat_ptr->scatter(inboundRay, rec, srec))
        return emitted;

    if (srec.is_specular) {
        return srec.attenuation
               * ray_color(srec.specular_ray, background, worldObjects, lights, depth-1);
    }

    auto light_ptr = std::make_shared<hittable_pdf>(lights, rec.p);
    mixture_pdf p(light_ptr, srec.pdf_ptr);

    Ray scattered = Ray(rec.p, p.generate());
    auto pdf_val = p.value(scattered.direction);

    return emitted
           + srec.attenuation * rec.mat_ptr->scattering_pdf(inboundRay, rec, scattered)
                 * ray_color(scattered, background, worldObjects, lights, depth-1) / pdf_val;
}

QVector3D rgbPerSamples(const QVector3D& pixel, int samples)
{
    auto r = pixel.x();
    auto g = pixel.y();
    auto b = pixel.z();

    // Replace NaN components with zero. See explanation in Ray Tracing: The Rest of Your Life.
    if (r != r) r = 0.0;
    if (g != g) g = 0.0;
    if (b != b) b = 0.0;

    // Divide the color by the number of samples. and gamma-correct for gamma=2.0.
    auto scale = 1.0f / samples;
     r = sqrt(scale * pixel.x());
     g = sqrt(scale * pixel.y());
     b = sqrt(scale * pixel.z());

    return { 255.0f * std::clamp(r, 0.0f, 0.999f), 255.0f * std::clamp(g, 0.0f, 0.999f), 255.0f * std::clamp(b, 0.0f, 0.999f) };
}

bool near_zero(const QVector3D& vec)
{
    // Return true if the vector is close to zero in all dimensions.
    const auto s = 1e-8;
    return (fabs(vec[0]) < s) && (fabs(vec[1]) < s) && (fabs(vec[2]) < s);
}

QVector3D reflect(const QVector3D& v, const QVector3D& n)
{
    return v - 2 * QVector3D::dotProduct(v, n) * n;
}

QVector3D refract(const QVector3D& uv, const QVector3D& n, double etai_over_etat)
{
    auto cos_theta = fmin(QVector3D::dotProduct(-uv, n), 1.0);
    QVector3D r_out_perp = etai_over_etat * (uv + cos_theta * n);
    QVector3D r_out_parallel = -sqrt(fabs(1.0 - r_out_perp.lengthSquared())) * n;
    return r_out_perp + r_out_parallel;
}

} // namespace calc
