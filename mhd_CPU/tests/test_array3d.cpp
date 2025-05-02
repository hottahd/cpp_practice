#include <catch2/catch_test_macros.hpp>
#include "array3d.hpp"

TEST_CASE("Array3D constructor and access") {
    Array3D<double> arr(3, 4, 5);

    // Check dimensions
    REQUIRE(arr.size_x() == 3);
    REQUIRE(arr.size_y() == 4);
    REQUIRE(arr.size_z() == 5);

    // Check access
    arr(1, 2, 3) = 42.0;
    REQUIRE(arr(1, 2, 3) == 42.0);
}