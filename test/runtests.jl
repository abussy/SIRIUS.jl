#TODO: run tests folder by folder, provide options for the test run:
#      with MPI or not, with GPU or not, set OMP_NUM_THREADS, etc
#      For now, only single rank, 2 OMP threads
#TODO: find a way to pass number of ranks and threads as test option
#TODO: try to initialize sirius and MPI here, so that we do not leave tests with SIRIUS not finalized?

using MKL
using MPI
using SIRIUS

ENV["OMP_NUM_THREADS"] = "2"
runfile = joinpath(@__DIR__, "runtests_runner.jl")
run(`$(mpiexec()) -n 2 $(Base.julia_cmd())  $runfile`)
