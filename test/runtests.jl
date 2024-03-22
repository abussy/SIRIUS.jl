using TestItemRunner
using MPI
using MKL
using SIRIUS

#TODO: run tests folder by folder, provide options for the test run:
#      with MPI or not, with GPU or not, set OMP_NUM_THREADS, etc
#      For now, only single rank, 2 OMP threads


#TODO: find a way to pass the number of OMP_THREADS
cd("test_scf")
   #include("./test_scf/test_scf.jl")
   @run_package_tests run(`$(Base.julia_cmd()) ./test_scf/test_scf.jl`)
cd("..")

cd("test_nlcg")
   #include("./test_nlcg/test_nlcg.jl")
   run(`$(Base.julia_cmd()) test_nlcg.jl`)
cd("..")
