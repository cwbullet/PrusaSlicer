#ifndef BOUNDINGCIRCLE_HPP
#define BOUNDINGCIRCLE_HPP

#include <array>

#include "Point.hpp"
#include "Eigen/Dense"

namespace Slic3r {

template<int N, class T> struct BoundingCircle
{
    Vec<N, T> center;
    T R = 0.f;

    BoundingCircle(const Vec<N, T> &p1, const Vec<N, T> &p2, const Vec<N, T> &p3)
    {
        static constexpr T ALMOST_ZERO = 1.;

        auto a = p1 - p3, b = p2 - p3;
        auto aXb = a.cross(b);
        T  sq_aXb = aXb.squaredNorm();

        if (sq_aXb < ALMOST_ZERO) {
            // Degenerate triangle, lets fall back to just a bounding box query

            Eigen::AlignedBox<T, N> box(p1);
            box.extend(p2); box.extend(p3);

            center = box.center();
        } else {
            // https://en.wikipedia.org/wiki/Circumscribed_circle
            // section "Higher dimensions":
            center = p3 +
                     (a.squaredNorm() * b - b.squaredNorm() * a).cross(aXb) /
                     (2 * sq_aXb);
        }

        R = std::abs((p1 - center).norm());
    }

    explicit BoundingCircle(const std::array<Vec<N, T>, 3> &pts)
        : BoundingCircle(pts[0], pts[1], pts[2])
    {}
};

using BoundingCircle3f = BoundingCircle<3, float>;
using BoundingCircle3d = BoundingCircle<3, double>;
using BoundingCircle2f = BoundingCircle<2, float>;
using BoundingCircle2d = BoundingCircle<2, double>;

}

#endif // BOUNDINGCIRCLE_HPP
