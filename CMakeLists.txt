cmake_minimum_required(VERSION 3.14)
project(CGAL_TDA)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_BUILD_TYPE Release)

# Tell CMake where to find CGAL on M1 Mac
set(CMAKE_PREFIX_PATH "/opt/homebrew")
set(CGAL_DIR "/opt/homebrew/lib/cmake/CGAL")

find_package(CGAL REQUIRED)
find_package(Boost REQUIRED)

add_executable(cgal_tda prototype.cpp)

target_link_libraries(cgal_tda CGAL::CGAL)

# Gudhi headers location found on your Mac
target_include_directories(cgal_tda PRIVATE
    /opt/homebrew/include
    /Users/souhardyatola/Desktop/gudhi-devel/src/Simplex_tree/include
    /Users/souhardyatola/Desktop/gudhi-devel/src/common/include
    /Users/souhardyatola/Desktop/gudhi-devel/src/Persistent_cohomology/include
)