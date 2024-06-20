using MPI
using MKL

export SIRIU
include("../python/SIRIU.jl")

MPI.Init()
comm = MPI.COMM_WORLD

if !SIRIU.is_initialized()
   SIRIU.initialize(false)
end
ctx = SIRIU.create_context_from_json(comm, "./sirius.json")
SIRIU.initialize_context(ctx)

# Initialize Kpoint Set
k_grid = Vector{Cint}(undef, 3)
k_shift = Vector{Cint}(undef, 3)
use_symmetry = SIRIU.get_kp_params_from_ctx!(ctx, k_grid, k_shift)
kps = SIRIU.create_kset_from_grid(ctx; k_grid, k_shift, use_symmetry)

# Initialize Ground State
gs = SIRIU.create_ground_state(kps)

# Solve and retrieve results
initial_guess = true
save_state = true
density_tol, energy_tol, iter_solver_tol, max_niter = SIRIU.get_scf_params_from_ctx(ctx)
SIRIU.find_ground_state(gs; density_tol, energy_tol, iter_solver_tol, initial_guess, max_niter, save_state) 

# Querry results
energy = SIRIU.get_energy(gs, "total")
@show energy

#TODO: make really sure we are consistent with SIRIUS in the ordering
natoms = SIRIU.get_num_atoms(gs)
forces = Matrix{Cdouble}(undef, 3, natoms)
SIRIU.get_forces!(gs, "total", forces)
@show forces

#TODO: make really sure we are consistent with SIRIUS in the ordering
stress = Matrix{Cdouble}(undef, 3, 3)
SIRIU.get_stress_tensor!(gs, "total", stress)
@show stress

# Finalize
SIRIU.free_ground_state_handler!(gs)
SIRIU.free_kpoint_set_handler!(kps)
SIRIU.free_context_handler!(ctx)
SIRIU.finalize(call_mpi_fin=false)
MPI.Finalize()
