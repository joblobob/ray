
#include <hitPosition.h>

#include <shapes.h>

std::optional<hit_record> hitFromList(const std::vector<std::shared_ptr<shape>>& sphereList, const Ray& ray, double t_min, double t_max)
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

bool boundingBoxFromList(const std::vector<std::shared_ptr<shape>>& list, aabb& output_box)
{
    if (list.empty())
        return false;

    aabb temp_box;
    bool first_box = true;

    for (const auto& object : list) {
        if (!object->bounding_box(temp_box))
            return false;
        output_box = first_box ? temp_box : surrounding_box(output_box, temp_box);
        first_box = false;
    }

    return true;
}

aabb surrounding_box(aabb box0, aabb box1)
{
    QVector3D small(fmin(box0.min().x(), box1.min().x()),
        fmin(box0.min().y(), box1.min().y()),
        fmin(box0.min().z(), box1.min().z()));

    QVector3D big(fmax(box0.max().x(), box1.max().x()),
        fmax(box0.max().y(), box1.max().y()),
        fmax(box0.max().z(), box1.max().z()));

    return aabb(small, big);
}

double hittable_pdf::value(const QVector3D &direction) const {
    return ptr->pdf_value(o, direction);
}

QVector3D hittable_pdf::generate() const {
    return ptr->random(o);
}

double mixture_pdf::value(const QVector3D &direction) const {
    return 0.5 * p[0]->value(direction) + 0.5 *p[1]->value(direction);
}

QVector3D mixture_pdf::generate() const {
    if (calc::random_double01() < 0.5)
        return p[0]->generate();
    else
        return p[1]->generate();
}
