
#include <hitPosition.h>

#include <shapes.h>

std::optional<hit_record> hitFromList(const std::vector<std::unique_ptr<shape> > &sphereList, const Ray &ray, double t_min, double t_max)
{
    bool hit_anything = false;
    auto closest_so_far = t_max;
    std::optional<hit_record> retVal;

    for (const auto& object : sphereList) {
        std::optional<hit_record> didHit = object->hit(ray, t_min, closest_so_far);
        if (didHit.has_value()) {
            hit_anything = true;
            closest_so_far = didHit->t;
            retVal = didHit.value();
        }
    }

    return retVal;
}
