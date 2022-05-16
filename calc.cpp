
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

} // namespace calc
