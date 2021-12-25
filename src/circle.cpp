#include "circle.h"
#include <limits>

using namespace NLO;

Circle::Circle(std::vector<Coords> const& coordinates)
:
    coordinates_{coordinates}
{
    middleCell_ = findMiddle();
}

Circle::Coords Circle::findMiddle()
{
    int highest_x = 0;
    int lowest_x = INT_MAX;
    int highest_y = 0;
    int lowest_y = INT_MAX;
    int highest_z = 0;
    int lowest_z = INT_MAX;
    for (auto const& coord: coordinates_)
    {
        if (std::get<0>(coord) > highest_x)
            highest_x = std::get<0>(coord);
        if (std::get<0>(coord) < lowest_x)
            lowest_x = std::get<0>(coord);
        if (std::get<1>(coord) > highest_y)
            highest_y = std::get<1>(coord);
        if (std::get<1>(coord) < lowest_y)
            lowest_y = std::get<1>(coord);
        if (std::get<2>(coord) > highest_z)
            highest_z = std::get<2>(coord);
        if (std::get<2>(coord) < lowest_z)
            lowest_z = std::get<2>(coord);
    }
    int middle_x = 0;
    int middle_y = 0;
    int middle_z = 0;

    if ((highest_x - lowest_x) % 2 == 0)
        middle_x = (highest_x - lowest_x) / 2;
    else
        middle_x = (highest_x - lowest_x) / 2 + 1;

    if ((highest_y - lowest_y) % 2 == 0)
        middle_y = (highest_y - lowest_y) / 2;
    else
        middle_y = (highest_y - lowest_y) / 2 + 1;

    if ((highest_z - lowest_z) % 2 == 0)
        middle_z = (highest_z - lowest_z) / 2;
    else
        middle_z = (highest_z - lowest_z) / 2 + 1;

    //std::cout << middle_x << ' ' << middle_y << ' ' << middle_z << '\n';
    radiusY_ = middle_y;
    radiusZ_ = middle_z;
    // Assumes that the x coordinates are constant!!
    radius_ = std::max(middle_y, middle_z);
    std::cout << "Radius inside circle: " << radius_ << '\n';
    return {middle_x + lowest_x, middle_y + lowest_y, middle_z + lowest_z};
}

Circle::Coords Circle::middle() const
{
    return middleCell_;
}

walberla::uint_t Circle::radius() const
{
    return radius_;
}

walberla::uint_t Circle::radiusY() const
{
    return radiusY_;
}

walberla::uint_t Circle::radiusZ() const
{
    return radiusZ_;
}