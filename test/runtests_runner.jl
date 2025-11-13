using TestItemRunner
using MPI
using SIRIUS

@run_package_tests verbose=true

# Finalize once and for all, for all tests
SIRIUS.finalize(call_mpi_fin=false)
MPI.Finalize()
