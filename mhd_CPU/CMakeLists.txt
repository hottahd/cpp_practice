cmake_minimum_required(VERSION 3.10)
project(mhd_CPU)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

find_package(nlohmann_json REQUIRED)
find_package(OpenMP REQUIRED)

file(GLOB SOURCES "src/*.cpp")

include_directories(include)
add_executable(a.out ${SOURCES})
target_link_libraries(a.out
    PRIVATE
    OpenMP::OpenMP_CXX
    nlohmann_json::nlohmann_json
)

# Testing
include(FetchContent)
FetchContent_Declare(
    Catch2 
    GIT_REPOSITORY https://github.com/catchorg/Catch2.git
    GIT_TAG v3.5.0
)
FetchContent_MakeAvailable(Catch2)

file(GLOB TEST_SOURCES "tests/*.cpp")

add_executable(tests ${TEST_SOURCES})
target_link_libraries(tests 
    PRIVATE
    Catch2::Catch2WithMain
    OpenMP::OpenMP_CXX
    nlohmann_json::nlohmann_json
)

target_include_directories(tests PRIVATE
    ${CMAKE_SOURCE_DIR}/include
    ${Catch2_SOURCE_DIR}/src # for Catch2
    )
include(CTest)
include(Catch)
catch_discover_tests(tests)