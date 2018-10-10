#ifndef CARIBOU_GEOMETRY_ENTITY_H
#define CARIBOU_GEOMETRY_ENTITY_H

#include <array>
#include <initializer_list>
#include <algorithm>

namespace caribou
{
namespace geometry
{

struct BaseData
{
    bool operator==(const BaseData &) const {
        return true;
    }
};


/**
 * An entity represent the base of every geometry classes.
 * @tparam Data A data attached to the entity. It can be used as a property or a container of properties for the entity.
 */
template<typename Data=BaseData>
class Entity
{
public:
    Entity() = default;
    Entity(const Data & d) : data(d) {}
    Data data;
};

} // namespace geometry

} // namespace caribou

#endif //CARIBOU_GEOMETRY_ENTITY_H
