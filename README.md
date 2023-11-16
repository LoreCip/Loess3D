# Fortran implementation for Loess algorithm in 3 dimensions

## Overview

This repository contains the Fortran program `loess3d`, a three-dimensional locally weighted polynomial smoothing (LOESS) implementation. The program utilizes OpenMP for parallelization and optionally supports MPI for distributed computing.

## Table of Contents

- [Fortran implementation for Loess algorithm in 3 dimensions](#fortran-implementation-for-loess-algorithm-in-3-dimensions)
  - [Overview](#overview)
  - [Table of Contents](#table-of-contents)
  - [Getting Started](#getting-started)
    - [Prerequisites](#prerequisites)
    - [Building](#building)

## Getting Started

### Prerequisites

To build and run the `loess3d` program, you need the following:

- Fortran compiler with OpenMP support (e.g., GNU Fortran)
- HDF5 library
- Optionally, MPI library for distributed computing

### Building

1. Clone the repository:

   ```bash
   git clone https://github.com/LoreCip/Loess3D.git
   cd loess3d

2. Edit the Makefile to set the appropriate compiler and compiler flags (especially for MPI).

3. Build the program:

    ```bash
    make
    ```

### Usage

#### Input

The program reads input data from an HDF5 file. The HDF5 file must contain the following datasets:

- n: Number of data points along the first dimension.
- m: Number of data points along the second dimension.
- l: Number of data points along the third dimension.
- Nth: Number of OpenMP threads to use.
- degree: Degree of the LOESS polynomial.
- frac: Fraction of total data points to use in the LOESS computation.
- xx: Data for the first dimension.
- yy: Data for the second dimension.
- zz: Data for the third dimension.
- O_data: The actual data for smoothing. This name can be chosen freely as it is also specified as an additional command-line argument.

#### Output

The program outputs the smoothed data to an HDF5 file. The output file will contain datasets with the same dimensions as the input datasets.

### Running with MPI

If MPI is enabled, the program can be run in a parallel, distributed fashion. To run with MPI, use the following command:

```bash
mpirun -n <num_processes> ./bin/run <input_file.h5> O_data
```
