@testitem "Testing SIRIUS SCF:" begin
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
   
   energy = SIRIUS.get_energy(gs, "total")
   forces = SIRIUS.get_forces(gs, "total")
   stress = SIRIUS.get_stress_tensor(gs, "total")
   
   @testset begin
      ref_energy = -156.38601370759218
      ediff = abs(energy - ref_energy)
      @test ediff < 1.0e-8
      
      ref_forces = [-0.005115449649280062 0.004278805994561478 0.000951668287306915 -5.6159597693654385e-5 -5.7713961009117934e-5; -0.0009363236595949704 0.008810093900525541 -0.009470905886296609 0.011037466381802311 -0.009439624040984902; 0.0 0.0 0.0 0.0 0.0]
      fdiff = maximum(abs.(forces - ref_forces))
      @test fdiff < 1.0e-8
      
      stress_ref = [0.00033618607681220464 -1.5106589957323968e-6 0.0; -1.5106589957323968e-6 0.00030408232425386574 0.0; 0.0 0.0 0.00033600383184123096]
      sdiff = maximum(stress - stress_ref)
      @test sdiff < 1.0e-8
   end
   
   # Finalize
   SIRIUS.free_ground_state_handler(gs)
   SIRIUS.free_kpoint_set_handler(kps)
   SIRIUS.free_context_handler(ctx)
   SIRIUS.finalize(false)
   
   MPI.Finalize()
end
