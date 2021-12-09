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
    Circle();
    Circle(std::vector<Coords> const& coordinates);
    Coords middle() const;
    walberla::real_t radius() const;
    //friend std::ostream& operator<< (std::ostream& stream, const Circle& matrix);
private:
    std::vector<Coords> coordinates_;
    Coords middleCell_;
    walberla::real_t radius_;

    Coords findMiddle();
};
} // end namespace NLO
#endif