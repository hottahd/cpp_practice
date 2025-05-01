cmake_minimum_required(VERSION 3.10)
project(mhd_CPU)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

find_package(nlohmann_json REQUIRED)

file(GLOB SOURCES "src/*.cpp")

include_directories(include)
add_executable(a.out ${SOURCES})
target_link_libraries(a.out nlohmann_json::nlohmann_json)