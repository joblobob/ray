#include "bvh_node.h"

#include <algorithm>

#include "qdebug.h"

inline bool box_compare(const std::shared_ptr<shape>& a, const std::shared_ptr<shape>& b, int axis)
{
    aabb box_a;
    aabb box_b;

    if (!a->bounding_box(box_a) || !b->bounding_box(box_b))
        qCritical() << "No bounding box in bvh_node constructor.\n";

    return box_a.min()[axis] < box_b.min()[axis];
}

bool box_x_compare(const std::shared_ptr<shape>& a, const std::shared_ptr<shape>& b)
{
    return box_compare(a, b, 0);
}

bool box_y_compare(const std::shared_ptr<shape>& a, const std::shared_ptr<shape>& b)
{
    return box_compare(a, b, 1);
}

bool box_z_compare(const std::shared_ptr<shape>& a, const std::shared_ptr<shape>& b)
{
    return box_compare(a, b, 2);
}

bvh_node::bvh_node(
    const std::vector<std::shared_ptr<shape>>& src_objects, size_t start, size_t end)
{
    auto objects = src_objects; // Create a modifiable array of the source scene objects

    int axis = calc::random_int(0, 2);
    auto comparator = (axis == 0) ? box_x_compare
                                  : (axis == 1) ? box_y_compare
                                                : box_z_compare;

    size_t object_span = end - start;

    if (object_span == 1) {
        left = right = objects[start];
    } else if (object_span == 2) {
        if (comparator(objects[start], objects[start + 1])) {
            left = objects[start];
            right = objects[start + 1];
        } else {
            left = objects[start + 1];
            right = objects[start];
        }
    } else {
        std::sort(objects.begin() + start, objects.begin() + end, comparator);

        auto mid = start + object_span / 2;
        left = std::make_shared<bvh_node>(objects, start, mid);
        right = std::make_shared<bvh_node>(objects, mid, end);
    }

    aabb box_left, box_right;

    if (!left->bounding_box(box_left)
        || !right->bounding_box(box_right))
        qCritical() << "No bounding box in bvh_node constructor.\n";

    box = surrounding_box(box_left, box_right);
}

std::optional<hit_record> bvh_node::hit(const Ray& r, double t_min, double t_max) const
{
    if (!box.hit(r, t_min, t_max))
        return {};

    auto hit_left = left->hit(r, t_min, t_max);
    hit_record rec;
    if (hit_left.has_value())
        rec = hit_left.value();

    auto hit_right = right->hit(r, t_min, hit_left.has_value() ? rec.t : t_max);

    if (hit_right.has_value())
        rec = hit_right.value();

    return rec;
}

bool bvh_node::bounding_box(aabb& output_box) const
{
    output_box = box;
    return true;
}
