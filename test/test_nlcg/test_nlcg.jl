using MPI
using MKL
using SIRIUS

MPI.Init()
comm = MPI.COMM_WORLD

# Initialize SIRIUS and context
SIRIUS.initialize(false)
ctx = SIRIUS.create_context_from_json(comm, "./sirius.json")
SIRIUS.initialize_context(ctx)

# Initialize Kpoint Set
k_grid, k_shift, use_symmetry = SIRIUS.get_kp_params_from_ctx(ctx)
kps = SIRIUS.create_kset_from_grid(ctx; k_grid, k_shift, use_symmetry)

# Initialize Ground State
gs = SIRIUS.create_ground_state(kps)

# Solve and retrieve results
density_tol, energy_tol, iter_solver_tol, max_niter = SIRIUS.get_scf_params_from_ctx(ctx)
SIRIUS.find_ground_state(gs, true, true; density_tol, energy_tol, iter_solver_tol, max_niter)

# Finish minimization with nlcg
@show temp, smearing, kappa, tau, tol, maxiter, restart, processing_unit = SIRIUS.get_nlcg_params_from_ctx(ctx)
@show converged = SIRIUS.nlcg(gs, kps; temp, smearing, kappa, tau, tol, maxiter, restart, processing_unit)

@show energy = SIRIUS.get_energy(gs, "total")
@show forces = SIRIUS.get_forces(gs, "total")
@show stress = SIRIUS.get_stress_tensor(gs, "total")

ref_energy = -164.63167507636578
ediff = abs(energy - ref_energy)
@show ediff, ediff < 1.0e-6

ref_forces = [0.0; 0.0; 0.0]

fdiff = maximum(abs.(forces - ref_forces))
@show fdiff, fdiff < 1.0e-8

stress_ref = [0.0002640500238659804 0.0 -1.807003620809175e-20; 0.0 0.0002640500238660082 3.61400724161835e-20; -1.807003620809175e-20 3.61400724161835e-20 0.00026405002386584164]
sdiff = maximum(stress - stress_ref)
@show sdiff, sdiff < 1.0e-8

# Finalize
SIRIUS.free_ground_state_handler(gs)
SIRIUS.free_kpoint_set_handler(kps)
SIRIUS.free_context_handler(ctx)
SIRIUS.finalize(false)

MPI.Finalize()
