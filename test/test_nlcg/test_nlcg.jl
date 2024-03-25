@testitem "Testing SIRIUS SCF:" begin
   using MPI
   using MKL
   using SIRIUS

   cd("./test_nlcg")
   
   MPI.Init()
   comm = MPI.COMM_WORLD
   
   # Initialize SIRIUS and context
   if !SIRIUS.is_initialized()
      SIRIUS.initialize(false)
   end
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
   temp, smearing, kappa, tau, tol, maxiter, restart, processing_unit = SIRIUS.get_nlcg_params_from_ctx(ctx)
   converged = SIRIUS.nlcg(gs, kps; temp, smearing, kappa, tau, tol, maxiter, restart, processing_unit)
   
   energy = SIRIUS.get_energy(gs, "total")
   forces = SIRIUS.get_forces(gs, "total")
   stress = SIRIUS.get_stress_tensor(gs, "total")
   
   @testset begin
      ref_energy = -156.3801929967038
      ediff = abs(energy - ref_energy)
      @test ediff < 1.0e-7
      
      ref_forces = [-0.005319386671138636 0.0037611865774594263 0.0015926647641683897 -5.2070785584471655e-5 2.2184428132771524e-5; -0.0016842404567648636 0.009553606795418727 -0.009796977930904309 0.011630693136752458 -0.009699366955623394; 0.0 0.0 0.0 0.0 0.0]
      fdiff = maximum(abs.(forces - ref_forces))
      @test fdiff < 1.0e-5
      
      stress_ref = [0.0003104844541582559 -8.878986858346005e-7 0.0; -8.878986858346005e-7 0.00028104998359693503 0.0; 0.0 0.0 0.000310236182434925]
      sdiff = maximum(stress - stress_ref)
      @test sdiff < 1.0e-7
   end
   
   # Finalize
   SIRIUS.free_ground_state_handler(gs)
   SIRIUS.free_kpoint_set_handler(kps)
   SIRIUS.free_context_handler(ctx)
   #SIRIUS.finalize(false)
   #MPI.Finalize()

   cd("..")
end
