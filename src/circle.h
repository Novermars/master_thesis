#ifndef NLO_WALBERLA_CIRCLE_H_
#define NLO_WALBERLA_CIRCLE_H_

#include "core/all.h"
#include <tuple>
#include <vector>
#include <iosfwd>

namespace NLO
{
class Circle
{
public:
    using Coords = std::tuple<walberla::cell_idx_t, walberla::cell_idx_t, walberla::cell_idx_t>;
private:
    std::vector<Coords> coordinates_;
    Coords middleCell_;
    walberla::uint_t radius_;
    walberla::uint_t radiusY_;
    walberla::uint_t radiusZ_;
public:
    Circle();
    Circle(std::vector<Coords> const& coordinates);
    Coords middle() const;
    walberla::uint_t radius() const;
    walberla::uint_t radiusY() const;
    walberla::uint_t radiusZ() const;
private:
    Coords findMiddle();
};
} // end namespace NLO
#endif