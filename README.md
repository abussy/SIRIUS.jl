# SIRIUS

[![Build Status](https://github.com/abussy/SIRIUS.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/abussy/SIRIUS.jl/actions/workflows/CI.yml?query=branch%3Amain)

SIRIUS.jl is a Julia package providing wrappers for the C-style API of the 
[SIRIUS](https://github.com/electronic-structure/SIRIUS) library, a domain specific library for
electronic structure calculations.

The package relies on the SIRIUS binary distributed as [SIRIUS_jll.jl](SIRIUS_jll). The source code
is generated by the `gen/generator.jl` script in a two-stage process: first, Julia wrappers are 
automatically generated around the C-style API of SIRIUS in the `src/LibSirius.jl` file, using Clang.jl. 
Then, post-processing takes place for selected functions (as listed in the `gen/python/functions_to_parse.txt`) 
for ease of use and better compliance with Julia conventions, resulting in the `src/SIRIUS.jl` file.

The package can be tested with and without MPI support. Results for energies, forces and stresses are
compared to reference numbers obtained via the SIRIUS miniapp. 
To run the MPI tests: `import Pkg; Pkg.test("SIRIUS", test_args=["mpi"])`
