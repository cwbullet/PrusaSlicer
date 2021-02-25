#include <catch2/catch.hpp>
//#include "test_utils.hpp"

#include "libslic3r/SLA/Hollowing.hpp"
//#include <openvdb/tools/Filter.h>


//#include <libnest2d/tools/benchmark.h>

//#include <libslic3r/SimplifyMesh.hpp>

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

