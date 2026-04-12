# Topological Data Analysis (TDA) Prototype

This repository contains a from-scratch C++ implementation of a Topological Data Analysis pipeline, developed as a proof-of-concept for Google Summer of Code (GSoC) 2026 with CGAL.

## Architecture & Integration
Unlike standard baseline forks, this prototype was built from the ground up to ensure deep understanding of the toolchain:
* **Custom CMake:** A handcrafted `CMakeLists.txt` designed to dynamically link complex mathematical dependencies.
* **CGAL 6.0.3:** Utilized for core geometric data structures and mesh processing.
* **Gudhi:** Integrated for advanced topological computations (Simplex Trees, Persistence Homology).

## Benchmarking Models
The repository includes custom `.off` (Object File Format) files created specifically to benchmark the algorithm's performance on edge cases:
* `angle_test.off`
* `noisy_plane_benchmark.off`

## Build Instructions
This project uses CMake. To compile natively on a Linux/macOS environment:

```bash
mkdir build
cd build
cmake ..
make -j4