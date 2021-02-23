#include <catch2/catch.hpp>
#include "test_utils.hpp"

#include "libslic3r/SLA/Hollowing.hpp"
//#include <openvdb/tools/Filter.h>


//#include <libnest2d/tools/benchmark.h>

//#include <libslic3r/SimplifyMesh.hpp>
#include <libslic3r/BoundingCircle.hpp>

TEST_CASE("Bounding circle for right angle triangle, Z: 0", "[BoundingCircle]") {
    using namespace Slic3r;

    std::array<Vec3f, 3> pts = {{
        { .0f, .0f, .0f},
        {10.f,  0.f, 0.f},
        { 0.f, 5.f, 0.f}
    }};

    BoundingCircle bc{pts};

    Vec3f refc = (pts[1] + pts[2]) / 2.f;

    REQUIRE(bc.R == Approx((pts[1] + pts[2]).norm() / 2.));
    REQUIRE(bc.center.x() == Approx(refc.x()));
    REQUIRE(bc.center.y() == Approx(refc.y()));
    REQUIRE(bc.center.z() == Approx(refc.z()));
}

TEST_CASE("Bounding circle for degenerate triangle", "[BoundingCircle]") {
    using namespace Slic3r;

    std::array<Vec3f, 3> pts = {{
        { .0f, .0f, .0f},
        {10.f,  0.f, 0.f},
        { 0.f, 0.f, 0.f}
    }};

    BoundingCircle bc{pts};

    Vec3f refc = (pts[1] + pts[2]) / 2.f;

    REQUIRE(bc.R == Approx((pts[1] + pts[2]).norm() / 2.));
    REQUIRE(bc.center.x() == Approx(refc.x()));
    REQUIRE(bc.center.y() == Approx(refc.y()));
    REQUIRE(bc.center.z() == Approx(refc.z()));
}

TEST_CASE("Hollow two overlapping spheres", "[Hollowing]") {
    using namespace Slic3r;

    TriangleMesh sphere1 = make_sphere(10., 2 * PI / 20.), sphere2 = sphere1;

    sphere1.translate(-5.f, 0.f, 0.f);
    sphere2.translate( 5.f, 0.f, 0.f);

    sphere1.merge(sphere2);
    sphere1.require_shared_vertices();

    sla::hollow_mesh(sphere1, sla::HollowingConfig{}, sla::HollowingFlags::hfRemoveInsideTriangles);

    sphere1.WriteOBJFile("twospheres.obj");
}

