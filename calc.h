#ifndef CALC_H
#define CALC_H

#include <QVector3D>
#include <random>

#include <cmath>
#include <limits>

#include <Ray.h>
#include <hitPosition.h>
#include <shapes.h>

namespace img {
// img
constexpr float aspect_ratio = 16.0f / 9.0f;
constexpr int width = 600;
constexpr int height = static_cast<int>(width / aspect_ratio);
constexpr QVector3D defaultVec { 0.0f, 0.0f, 0.0f };
constexpr QVector3D infiniteZ { 0.0f, 0.0f, -1.0f };
constexpr QVector3D gradientBgVec { 0.5f, 0.7f, 1.0f };
constexpr QVector3D bgColor { 1.0f, 1.0f, 1.0f };
}

namespace calc {

// Constants

constexpr double infinity = std::numeric_limits<double>::infinity();
constexpr double pi = 3.1415926535897932385;

// Utility Functions

inline constexpr double degrees_to_radians(double degrees);

inline float random_double(float min = 0.0, float max = 1.0);

inline QVector3D randomVec();
inline QVector3D randomVec(float min, float max);

inline QVector3D random_in_unit_sphere();
inline QVector3D random_unit_vector();
inline QVector3D random_in_hemisphere(const QVector3D& normal);

//ray calcs
QVector3D ray_color(const Ray& r, const std::vector<std::unique_ptr<shape>>& worldObjects,
    int depth, bool drawOnlyColors = false);
inline std::optional<hitPosition> hitFromList(const std::vector<std::unique_ptr<shape>>& sphereList,
    const Ray& ray, double t_min,
    double t_max);

} // namespace calc

#endif // CALC_H
