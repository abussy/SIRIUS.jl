#TODO: run tests folder by folder, provide options for the test run:
#      with MPI or not, with GPU or not, set OMP_NUM_THREADS, etc
#      For now, only single rank, 2 OMP threads
#TODO: find a way to pass the number of OMP_THREADS

using MKL
using MPI
using SIRIUS

ENV["OMP_NUM_THREADS"] = "2"
runfile = joinpath(@__DIR__, "runtests_runner.jl")
run(`$(mpiexec()) -n 2 $(Base.julia_cmd()) --check-bounds=yes --depwarn=yes --project=@ --color=yes --startup-file=no $runfile`)
