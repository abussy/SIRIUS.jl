using Clang.Generators
using SIRIUS_jll
using PyCall

### Automatically generate the Julia wrappers for SIRIUS in LibSirius.jl
cd(@__DIR__)

include_dir = normpath(SIRIUS_jll.artifact_dir, "include/sirius/src/api")

# wrapper generator options
options = load_options(joinpath(@__DIR__, "generator.toml"))

# add compiler flags
args = get_default_args()
push!(args, "-I$include_dir")
push!(args, "-fparse-all-comments")

header_dir = include_dir
headers = [joinpath(header_dir, header) for header in readdir(header_dir) if endswith(header, ".h")]

# create context
ctx = create_context(headers, args, options)

# run generator
build!(ctx)

### Second layer of wrapping for ease of use
cd("./python")
@pyinclude("generate_module.py")
