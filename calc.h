#ifndef CALC_H
#define CALC_H

#include <QVector3D>
#include <random>

#include <cmath>
#include <limits>

#include <Ray.h>

struct shape;
namespace img {
// img
constexpr float aspect_ratio = 1.0f;
constexpr int width = 500;
constexpr int height = static_cast<int>(width / aspect_ratio);
constexpr int totalPixels = width * height;
constexpr QVector3D defaultVec { 0.0f, 0.0f, 0.0f };
constexpr QVector3D infiniteZ { 0.0f, 0.0f, -1.0f };
constexpr QVector3D gradientBgVec { 0.7f, 0.70f, 0.70f };
constexpr QVector3D bgColor { 1.0f, 1.0f, 1.0f };
}

namespace calc {

// Constants

constexpr double infinity = std::numeric_limits<double>::infinity();
constexpr double pi = 3.1415926535897932385;
constexpr float smallestVal = 0.001f;

// Utility Functions

inline constexpr double degrees_to_radians(double degrees)
{
    return degrees * pi / 180.0;
}

float random_double(float min, float max);
float random_double01();
float random_double11();

QVector3D randomVec();
QVector3D randomVec11();
QVector3D randomVec(float min, float max);

QVector3D random_in_unit_sphere();
QVector3D random_in_unit_disk();
QVector3D random_unit_vector();
inline QVector3D random_in_hemisphere(const QVector3D& normal);

inline int random_int(int min, int max)
{
    // Returns a random integer in [min,max].
    return static_cast<int>(random_double(min, max + 1));
}

//ray calcs
QVector3D ray_color(const Ray& r, const QVector3D& background, const std::vector<std::shared_ptr<shape>>& worldObjects,
    int depth, bool drawOnlyColors = false);

QVector3D reflect(const QVector3D& v, const QVector3D& n);
QVector3D refract(const QVector3D& uv, const QVector3D& n, double etai_over_etat);

bool near_zero(const QVector3D& vec);

QVector3D rgbPerSamples(const QVector3D& pixel, int samples);

static std::random_device rd; // Will be used to obtain a seed for the random number engine
static std::minstd_rand0 gen(rd());
static std::random_device rd01; // Will be used to obtain a seed for the random number engine
static std::minstd_rand0 gen01(rd01());
static std::random_device rd11; // Will be used to obtain a seed for the random number engine
static std::minstd_rand0 gen11(rd11());
static std::uniform_real_distribution<float> simple_distribution(0.0f, 1.0f);

static std::uniform_real_distribution<float> simple_distributionMinusOneToOne(-1.0f, 1.0f);

} // namespace calc

#endif // CALC_H
