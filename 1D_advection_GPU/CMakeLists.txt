cmake_minimum_required(VERSION 3.10)
project(Advection LANGUAGES CXX CUDA)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)
set(CMAKE_CUDA_ARCHITECTURES 89)

find_package(nlohmann_json REQUIRED)

file(GLOB SOURCES "src/*.cpp" "src/*.cu")

include_directories(include)
add_executable(a.out ${SOURCES})
target_link_libraries(a.out)