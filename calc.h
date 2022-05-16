#ifndef CALC_H
#define CALC_H

#include <QVector3D>
#include <random>

#include <cmath>
#include <limits>

namespace calc {

// Constants

constexpr double infinity = std::numeric_limits<double>::infinity();
constexpr double pi = 3.1415926535897932385;

// Utility Functions

constexpr double degrees_to_radians(double degrees);

float random_double(float min = 0.0, float max = 1.0);

QVector3D randomVec();
QVector3D randomVec(float min, float max);

QVector3D random_in_unit_sphere();
QVector3D random_unit_vector();
QVector3D random_in_hemisphere(const QVector3D& normal);

} // namespace calc

#endif // CALC_H
