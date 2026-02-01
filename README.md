# Topological Filtering Prototypes (GSoC 2026 Prep)

This repository contains research prototypes for the "Topological Filtering of Features" project proposed for CGAL.
The goal is to filter noise from 3D meshes using Persistent Homology (TDA) rather than local angle checks.

## The Files

### 1. benchmark_final.cpp
* **Purpose:** Demonstrates the failure case of the current PMP::detect_sharp_edges algorithm.
* **Findings:** The current algorithm detects **109 false positives** on a simple plane with random Z-noise (Amplitude 0.4).
* **Status:** Reported as a GitHub Issue on the main CGAL repository.
* **Methodology:** Generates a synthetic 20x20 grid, perturbs vertices with random noise, and counts sharp edges detected at a 45-degree threshold.

### 2. find_saddles.cpp (The 3D Solver)
* **Purpose:** A custom, lightweight implementation of the topological pipeline without external libraries.
* **Algorithm:**
    * **Flow:** Computes Steepest Descent for every vertex to determine flow direction.
    * **Basins:** Performs Watershed Segmentation using a Union-Find data structure to group vertices into catchment basins.
    * **Persistence:** Calculates the persistence value (Saddle_Height - Basin_Height) to distinguish significant features from topological noise.

## Goal
To integrate this persistence-based filtering into CGAL::Polygon_mesh_processing to robustly detect sharp features even on noisy input scans.
