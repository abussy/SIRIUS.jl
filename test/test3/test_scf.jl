@testitem "Testing SIRIUS SCF #3:" begin 
   using MPI
   using SIRIUS
   
   include("../test_utils.jl")

   MPI.Init()
   comm = MPI.COMM_WORLD
   
   if !SIRIUS.is_initialized()
      SIRIUS.initialize(false)
   end
   ctx = SIRIUS.create_context_from_json(comm, "./sirius.json")
   SIRIUS.initialize_context(ctx)
   
   # Initialize Kpoint Set
   k_grid = Vector{Cint}(undef, 3)
   k_shift = Vector{Cint}(undef, 3)
   use_symmetry = SIRIUS.get_kp_params_from_ctx!(ctx, k_grid, k_shift)
   kps = SIRIUS.create_kset_from_grid(ctx; k_grid, k_shift, use_symmetry)
   
   # Initialize Ground State
   gs = SIRIUS.create_ground_state(kps)
   
   # Solve and retrieve results
   initial_guess = true
   save_state = true
   density_tol, energy_tol, iter_solver_tol, max_niter = SIRIUS.get_scf_params_from_ctx(ctx)
   SIRIUS.find_ground_state(gs; density_tol, energy_tol, iter_solver_tol, initial_guess, max_niter, save_state) 
   
   @testset begin
      energy = SIRIUS.get_energy(gs, "total")
      ref_energy = get_ref_energy("output_ref.json", "total")
      ediff = abs(energy-ref_energy)
      @show ediff < 1.0e-8
      @test ediff < 1.0e-8
   end
   
   # Finalize
   SIRIUS.free_ground_state_handler!(gs)
   SIRIUS.free_kpoint_set_handler!(kps)
   SIRIUS.free_context_handler!(ctx)

   # Last test
   SIRIUS.finalize(call_mpi_fin=false)
   MPI.Finalize()
end
