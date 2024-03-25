@testitem "Testing SIRIUS nlcg:" begin
   using MPI
   using MKL
   using SIRIUS
   
   cd("test_scf")
   
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
   
   energy = SIRIUS.get_energy(gs, "total")
   forces = SIRIUS.get_forces(gs, "total")
   stress = SIRIUS.get_stress_tensor(gs, "total")
   
   @testset begin
      ref_energy = -164.63200441770132
      ediff = abs(energy - ref_energy)
      @test ediff < 1.0e-7
      
      ref_forces = [0.0; 0.0; 0.0]
      fdiff = maximum(abs.(forces - ref_forces))
      @test fdiff < 1.0e-5
      
      stress_ref = [0.00018548995729250617 4.012354050806741e-35 4.517509051496563e-21; 4.012354050806741e-35 0.0001854899572924784 9.035018103783277e-21; 4.517509051496563e-21 9.035018103783277e-21 0.00018548995729275597]
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
