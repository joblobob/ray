#ifndef BVH_H
#define BVH_H

#include "qdebug.h"
#include "shapes.h"

#include "calc.h"
#include "hitPosition.h"

class bvh_node : public shape {
public:
    bvh_node() {};

    bvh_node(const std::vector<std::shared_ptr<shape>>& list)
        : bvh_node(list, 0, list.size())
    {
    }

    bvh_node(const std::vector<std::shared_ptr<shape>>& src_objects, size_t start, size_t end);

    virtual std::optional<hit_record> hit(const Ray& r, double t_min, double t_max) const override;

    virtual bool bounding_box(aabb& output_box) const override;

public:
    std::shared_ptr<shape> left;
    std::shared_ptr<shape> right;
    aabb box;
};

#endif
