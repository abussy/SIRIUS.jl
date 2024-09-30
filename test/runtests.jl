#TODO: run tests folder by folder, provide options for the test run:
#      with MPI or not, with GPU or not, set OMP_NUM_THREADS, etc
#      For now, only single rank, 2 OMP threads
#TODO: find a way to pass number of ranks and threads as test option
#TODO: try to initialize sirius and MPI here, so that we do not leave tests with SIRIUS not finalized?

using MPI

function parse_test_args()
    args = Symbol.(ARGS)
    if "SIRIUS_TEST_ARGS" in keys(ENV) && isempty(ARGS)
        args = Symbol.(split(ENV["SIRIUS_TEST_ARGS"], "-"))
    end

    tags = [:mpi, :serial]
    base_tag = filter(in(tags), args)

    if isempty(base_tag)
        base_tag = :serial
    elseif length(base_tag) > 1
        error("SIRIUS.jl testing takes only one argument.")
    else
        base_tag = base_tag[1]
    end

    return base_tag
end

base_tag = parse_test_args()
runfile = joinpath(@__DIR__, "runtests_runner.jl")

ENV["OMP_NUM_THREADS"] = "2"
if base_tag == :mpi
    run(`$(mpiexec()) -n 2 $(Base.julia_cmd())  $runfile`)
else
    run(`$(Base.julia_cmd())  $runfile`)
end
