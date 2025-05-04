
#include <memory>
#include <catch2/catch_test_macros.hpp>

#include "array3d.hpp"
#include "config.hpp"
#include "model.hpp"
#include "boundary_condition.hpp"


TEST_CASE("Test BoundaryCondition class") {
    Model<double> model = Model<double>::from_config_file("../tests/config.json");
    Grid<Real>& grid = model.grid;

    // Create a BoundaryCondition object
    BoundaryCondition<double> bc(model);

    // Test the bnd_range function
    int i0_, i1_, j0_, j1_, k0_, k1_;
    
    bc.bnd_range(i0_, i1_, j0_, j1_, k0_, k1_, "x");
    REQUIRE(i0_ == 0);
    REQUIRE(i1_ == grid.i_margin);
    REQUIRE(j0_ == 0);
    REQUIRE(j1_ == grid.j_total);
    REQUIRE(k0_ == 0);
    REQUIRE(k1_ == grid.k_total);

    bc.bnd_range(i0_, i1_, j0_, j1_, k0_, k1_, "y");
    REQUIRE(i0_ == 0);
    REQUIRE(i1_ == grid.i_total);
    REQUIRE(j0_ == 0);
    REQUIRE(j1_ == grid.j_margin);
    REQUIRE(k0_ == 0);
    REQUIRE(k1_ == grid.k_total);

    bc.bnd_range(i0_, i1_, j0_, j1_, k0_, k1_, "z");
    REQUIRE(i0_ == 0);
    REQUIRE(i1_ == grid.i_total);
    REQUIRE(j0_ == 0);
    REQUIRE(j1_ == grid.j_total);
    REQUIRE(k0_ == 0);
    REQUIRE(k1_ == grid.k_margin);


    // test for a margin = 2 case
    Array3D<double> arr(grid.i_total, grid.j_total, grid.k_total);

    // x boundary test
    for (int i = 0; i < grid.i_total; ++i) {
        for (int j = 0; j < grid.j_total; ++j) {
            for (int k = 0; k < grid.k_total; ++k) {
                arr(i, j, k) = i;
            }
        }
    }

    for (int j = 0; j < grid.j_total; ++j) {
        for (int k = 0; k < grid.k_total; ++k) {
            REQUIRE(arr(0, j, k) != arr(3, j, k));
            REQUIRE(arr(1, j, k) != arr(2, j, k));            
            REQUIRE(arr(grid.i_total - 1, j, k) != arr(grid.i_total-4, j, k));
            REQUIRE(arr(grid.i_total - 2, j, k) != arr(grid.i_total-3, j, k));
        }
    }

    bc.bnd_symmetric(arr, nullptr, 1.0, "x", "inner");
    for (int j = 0; j < grid.j_total; ++j) {
        for (int k = 0; k < grid.k_total; ++k) {
            REQUIRE(arr(0, j, k) == arr(3, j, k));
            REQUIRE(arr(1, j, k) == arr(2, j, k));            
        }
    }

    bc.bnd_symmetric(arr, nullptr, -1.0, "x", "inner");
    for (int j = 0; j < grid.j_total; ++j) {
        for (int k = 0; k < grid.k_total; ++k) {
            REQUIRE(arr(0, j, k) == -arr(3, j, k));
            REQUIRE(arr(1, j, k) == -arr(2, j, k));            
        }
    }

    bc.bnd_symmetric(arr, nullptr, 1.0, "x", "outer");
    for (int j = 0; j < grid.j_total; ++j) {
        for (int k = 0; k < grid.k_total; ++k) {
            REQUIRE(arr(grid.i_total - 1, j, k) == arr(grid.i_total-4, j, k));
            REQUIRE(arr(grid.i_total - 2, j, k) == arr(grid.i_total-3, j, k));
        }
    }

    bc.bnd_symmetric(arr, nullptr, -1.0, "x", "outer");
    for (int j = 0; j < grid.j_total; ++j) {
        for (int k = 0; k < grid.k_total; ++k) {
            REQUIRE(arr(grid.i_total - 1, j, k) == -arr(grid.i_total-4, j, k));
            REQUIRE(arr(grid.i_total - 2, j, k) == -arr(grid.i_total-3, j, k));
        }
    }

    // y boundary test
    for (int i = 0; i < grid.i_total; ++i) {
        for (int j = 0; j < grid.j_total; ++j) {
            for (int k = 0; k < grid.k_total; ++k) {
                arr(i, j, k) = j;
            }
        }
    }

    for (int i = 0; i < grid.i_total; ++i) {
        for (int k = 0; k < grid.k_total; ++k) {
            REQUIRE(arr(i, 0, k) != arr(i, 3, k));
            REQUIRE(arr(i, 1, k) != arr(i, 2, k));            
            REQUIRE(arr(i, grid.j_total - 1, k) != arr(i, grid.j_total-4, k));
            REQUIRE(arr(i, grid.j_total - 2, k) != arr(i, grid.j_total-3, k));
        }
    }

    bc.bnd_symmetric(arr, nullptr, 1.0, "y", "inner");
    for (int i = 0; i < grid.i_total; ++i) {
        for (int k = 0; k < grid.k_total; ++k) {
            REQUIRE(arr(i, 0, k) == arr(i, 3, k));
            REQUIRE(arr(i, 1, k) == arr(i, 2, k));            
        }
    }

    bc.bnd_symmetric(arr, nullptr, -1.0, "y", "inner");
    for (int i = 0; i < grid.i_total; ++i) {
        for (int k = 0; k < grid.k_total; ++k) {
            REQUIRE(arr(i, 0, k) == -arr(i, 3, k));
            REQUIRE(arr(i, 1, k) == -arr(i, 2, k));
        }
    }

    bc.bnd_symmetric(arr, nullptr, 1.0, "y", "outer");
    for (int i = 0; i < grid.i_total; ++i) {
        for (int k = 0; k < grid.k_total; ++k) {
            REQUIRE(arr(i, grid.j_total - 1, k) == arr(i, grid.j_total-4, k));
            REQUIRE(arr(i, grid.j_total - 2, k) == arr(i, grid.j_total-3, k));
        }
    }

    bc.bnd_symmetric(arr, nullptr, -1.0, "y", "outer");
    for (int i = 0; i < grid.i_total; ++i) {
        for (int k = 0; k < grid.k_total; ++k) {

            REQUIRE(arr(i, grid.j_total - 1, k) == -arr(i, grid.j_total - 4, k));
            REQUIRE(arr(i, grid.j_total - 2, k) == -arr(i, grid.j_total - 3, k));
        }
    }

    // z boundary test
    for (int i = 0; i < grid.i_total; ++i) {
        for (int j = 0; j < grid.j_total; ++j) {
            for (int k = 0; k < grid.k_total; ++k) {
                arr(i, j, k) = k;
            }
        }
    }

    for (int i = 0; i < grid.i_total; ++i) {
        for (int j = 0; j < grid.j_total; ++j) {
            REQUIRE(arr(i, j, 0) != arr(i, j, 3));
            REQUIRE(arr(i, j, 1) != arr(i, j, 2));
            REQUIRE(arr(i, j, grid.k_total - 1) != arr(i, j, grid.k_total-4));
            REQUIRE(arr(i, j, grid.k_total - 2) != arr(i, j, grid.k_total-3));
        }
    }

    bc.bnd_symmetric(arr, nullptr, 1.0, "z", "inner");
    for (int i = 0; i < grid.i_total; ++i) {
        for (int j = 0; j < grid.j_total; ++j) {
            REQUIRE(arr(i, j, 0) == arr(i, j, 3));
            REQUIRE(arr(i, j, 1) == arr(i, j, 2));            
        }
    }

    bc.bnd_symmetric(arr, nullptr, -1.0, "z", "inner");
    for (int i = 0; i < grid.i_total; ++i) {
        for (int j = 0; j < grid.j_total; ++j) {
            REQUIRE(arr(i, j, 0) == -arr(i, j, 3));
            REQUIRE(arr(i, j, 1) == -arr(i, j, 2));
        }
    }

    bc.bnd_symmetric(arr, nullptr, 1.0, "z", "outer");
    for (int i = 0; i < grid.i_total; ++i) {
        for (int j = 0; j < grid.j_total; ++j) {
            REQUIRE(arr(i, j, grid.k_total - 1) == arr(i, j, grid.k_total-4));
            REQUIRE(arr(i, j, grid.k_total - 2) == arr(i, j, grid.k_total-3));
        }
    }

    bc.bnd_symmetric(arr, nullptr, -1.0, "z", "outer");
    for (int i = 0; i < grid.i_total; ++i) {
        for (int j = 0; j < grid.j_total; ++j) {

            REQUIRE(arr(i, j, grid.k_total - 1) == -arr(i, j, grid.k_total - 4));
            REQUIRE(arr(i, j, grid.k_total - 2) == -arr(i, j, grid.k_total - 3));
        }
    }
    
    int i;
    int i_total = 10;
    int i_margin = 2;
    int i_ghst, i_trgt;

    i = 0;
    bc.bnd_symmetric_index(i, i_total, i_margin, i_ghst, i_trgt, "inner");
    REQUIRE(i_ghst == 0);
    REQUIRE(i_trgt == 3);

    bc.bnd_symmetric_index(i, i_total, i_margin, i_ghst, i_trgt, "outer");
    REQUIRE(i_ghst == 8);
    REQUIRE(i_trgt == 7);

    i = 1;
    bc.bnd_symmetric_index(i, i_total, i_margin, i_ghst, i_trgt, "inner");
    REQUIRE(i_ghst == 1);
    REQUIRE(i_trgt == 2);

    bc.bnd_symmetric_index(i, i_total, i_margin, i_ghst, i_trgt, "outer");
    REQUIRE(i_ghst == 9);
    REQUIRE(i_trgt == 6);

}

//     // // Test the bnd_symmetric function
//     // Array3D<double> arr(model.grid.i_total, model.grid.j_total, model.grid.k_total);
//     // Array3D<double>* fac = nullptr;
//     // bc.bnd_symmetric(arr, fac, -1.0, "x", "inner");

// }
