// #include <catch2/catch_test_macros.hpp>
// #include "model.hpp"
// #include "boundary_condition.hpp"

// TEST_CASE("Test BoundaryCondition class") {
//     Model<double> model = Model<double>::from_config_file("./config.json");

//     // // Create a BoundaryCondition object
//     // BoundaryCondition<double> bc(model);

//     // // Test the bnd_range function
//     // int i0, i1, j0, j1, k0, k1;
//     // bc.bnd_range(i0, i1, j0, j1, k0, k1, "x");
//     // REQUIRE(i0 == 0);
//     // REQUIRE(i1 == model.grid.i_margin);
//     // REQUIRE(j0 == 0);
//     // REQUIRE(j1 == model.grid.j_total);
//     // REQUIRE(k0 == 0);
//     // REQUIRE(k1 == model.grid.k_total);

//     // // Test the bnd_symmetric function
//     // Array3D<double> arr(model.grid.i_total, model.grid.j_total, model.grid.k_total);
//     // Array3D<double>* fac = nullptr;
//     // bc.bnd_symmetric(arr, fac, -1.0, "x", "inner");

// }
