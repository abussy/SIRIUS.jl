module SIRIUS

using MPI
using SIRIUS_jll

# Note on documentation: The functions defined in this module are wrappers around the C-style API 
#                        of SIRIUS. Their precise documentation is available in the SIRIUS source
#                        code, in src/api/sirius_api.cpp

# Points to the shared library provided by SIRIUS_jll
# Library path can be changed externally by SIRIUS.libpath = "path"
# This can be useful for linking to a SIRIUS library that isoptimally built
libpath = libsirius

### Utility function that maps an integer to a C MPI_comm
comm2f(comm::MPI.Comm) = ccall((:MPI_Comm_c2f, MPI.libmpi), Cint, (MPI.MPI_Comm,), comm)

### SIRIUS library initialization and finalization
function initialize(call_mpi_init::Bool)
   call_mpi_init__ = Ref{Cuchar}(call_mpi_init)
   error_code__ = Ref{Cint}(0)
   @ccall libpath.sirius_initialize(call_mpi_init__::Ref{Cuchar}, error_code__::Ref{Cint})::Cvoid
   if error_code__[] != 0
      error("Sirius.initialize failed with error code", error_code__[])
   end
end

function finalize(call_mpi_fin::Bool, call_device_reset::Bool=true, call_fftw_fin::Bool=true)
   call_mpi_fin__ = Ref{Cuchar}(call_mpi_fin)
   call_device_reset__ = Ref{Cuchar}(call_device_reset)
   call_fftw_fin__ = Ref{Cuchar}(call_fftw_fin)
   error_code__ = Ref{Cint}(0)
   @ccall libpath.sirius_finalize(call_mpi_fin__::Ref{Cuchar}, call_device_reset__::Ref{Cuchar},
                              call_fftw_fin__::Ref{Cuchar}, error_code__::Ref{Cint})::Cvoid
   if error_code__[] != 0
      error("Sirius.finalize failed with error code", error_code__[])
   end
end

function is_initialized()
   status__ = Ref{Cuchar}(false)
   error_code__ = Ref{Cint}(0)
   @ccall libpath.sirius_is_initialized(status__::Ref{Cuchar}, error_code__::Ref{Cint})::Cvoid
   if error_code__[] != 0
      error("Sirius.is_initialized failed with error code", error_code__[])
   end
   return Bool(status__[])
end

### Defining types for the handler pointers of the C-API
mutable struct ContextHandler
   handler_ptr::Ref{Ptr{Cvoid}}
end

mutable struct GroundStateHandler
   handler_ptr::Ref{Ptr{Cvoid}}
end

mutable struct KpointSetHandler
   handler_ptr::Ref{Ptr{Cvoid}}
end

### Handler freeing function. Note: not added as finalizer as call order matters
function free_context_handler(ctx::ContextHandler)
   error_code__ = Ref{Cint}(0)
   @ccall libpath.sirius_free_object_handler(ctx.handler_ptr::Ref{Ptr{Cvoid}}, error_code__::Ref{Cint})::Cvoid
   if error_code__[] != 0
      error("Sirius.free_context_handler failed with error code", error_code__[])
   end
end

function free_ground_state_handler(gs::GroundStateHandler)
   error_code__ = Ref{Cint}(0)
   @ccall libpath.sirius_free_object_handler(gs.handler_ptr::Ref{Ptr{Cvoid}}, error_code__::Ref{Cint})::Cvoid
   if error_code__[] != 0
      error("Sirius.free_ground_state_handler failed with error code", error_code__[])
   end
end

function free_kpoint_set_handler(kps::KpointSetHandler)
   error_code__ = Ref{Cint}(0)
   @ccall libpath.sirius_free_object_handler(kps.handler_ptr::Ref{Ptr{Cvoid}}, error_code__::Ref{Cint})::Cvoid
   if error_code__[] != 0
      error("Sirius.free_kpoint_set_handler failed with error code", error_code__[])
   end
end

### Functions related to the SIRIUS simulation context
function create_context_from_json(comm::MPI.Comm, fname::String)
   ctx = ContextHandler(C_NULL)
   fcomm__::Int32 = comm2f(comm)
   error_code__ = Ref{Cint}(0)
   @ccall libpath.sirius_create_context_from_json(fcomm__::Int32, ctx.handler_ptr::Ref{Ptr{Cvoid}},
                                                  fname::Cstring, error_code__::Ref{Cint})::Cvoid
   if error_code__[] != 0
      error("Sirius.create_context_from_json failed with error code", error_code__[])
   end
   return ctx
end
 
function create_context_from_json(fname::String)
   ctx = ContextHandler(C_NULL)
   error_code__ = Ref{Cint}(0)
   @ccall libpath.sirius_create_context_from_json_commworld(ctx.handler_ptr::Ref{Ptr{Cvoid}},
                                                            fname::Cstring, error_code__::Ref{Cint})::Cvoid
   if error_code__[] != 0
      error("Sirius.create_context_from_json failed with error code", error_code__[])
   end
   return ctx
end

function initialize_context(ctx::ContextHandler)
   error_code__ = Ref{Cint}(0)
   @ccall libpath.sirius_initialize_context(ctx.handler_ptr::Ptr{Cvoid}, error_code__::Ref{Cint})::Cvoid
   if error_code__[] != 0
      error("Sirius.initialize_context failed with error code", error_code__[])
   end
end 

function get_kp_params_from_ctx(ctx::ContextHandler)
   k_grid__ = Vector{Cint}(undef, 3)
   k_shift__ = Vector{Cint}(undef, 3)
   use_symmetry__ = Ref{Cuchar}(false)
   error_code__ = Ref{Cint}(0)
   @ccall libpath.sirius_get_kp_params_from_ctx(ctx.handler_ptr::Ptr{Cvoid}, k_grid__::Ref{Cint},
                                                k_shift__::Ref{Cint}, use_symmetry__::Ref{Cuchar},
                                                error_code__::Ref{Cint})::Cvoid
   if error_code__[] != 0
      error("Sirius.get_kp_params_from_ctx failed with error code", error_code__[])
   end
   return k_grid__, k_shift__, Bool(use_symmetry__[])
end

function get_scf_params_from_ctx(ctx::ContextHandler)
   density_tol__ = Ref{Cdouble}(0.0)
   energy_tol__ = Ref{Cdouble}(0.0)
   iter_tol__ = Ref{Cdouble}(0.0)
   max_niter__ = Ref{Cint}(0)
   error_code__ = Ref{Cint}(0)
   @ccall libpath.sirius_get_scf_params_from_ctx(ctx.handler_ptr::Ptr{Cvoid}, density_tol__::Ref{Cdouble},
                                                 energy_tol__::Ref{Cdouble}, iter_tol__::Ref{Cdouble},
                                                 max_niter__::Ref{Cint}, error_code__::Ref{Cint})::Cvoid
   if error_code__[] != 0
      error("Sirius.get_scf_params_from_ctx failed with error code", error_code__[])
   end
   return density_tol__[], energy_tol__[], iter_tol__[], max_niter__[]
end

function get_nlcg_params_from_ctx(ctx::ContextHandler)
   temp__ = Ref{Cdouble}(0.0)
   smearing__ = Vector{UInt8}(undef, 32)
   kappa__ = Ref{Cdouble}(0.0)
   tau__ = Ref{Cdouble}(0.0)
   tol__ = Ref{Cdouble}(0.0)
   maxiter__ = Ref{Cint}(0.0)
   restart__ = Ref{Cint}(0)
   processing_unit__ = Vector{UInt8}(undef, 32)
   error_code__ = Ref{Cint}(0)
   @ccall libpath.sirius_get_nlcg_params_from_ctx(ctx.handler_ptr::Ptr{Cvoid}, temp__::Ref{Cdouble},
                                                  smearing__::Ptr{UInt8}, kappa__::Ref{Cdouble}, 
                                                  tau__::Ref{Cdouble}, tol__::Ref{Cdouble}, 
                                                  maxiter__::Ref{Cint}, restart__::Ref{Cint},
                                                  processing_unit__::Ptr{UInt8}, 
                                                  error_code__::Ref{Cint})::Cvoid
   if error_code__[] != 0
      error("Sirius.get_nlcg_params_from_ctx failed with error code", error_code__[])
   end
   smearing__[end] = 0
   smear = GC.@preserve smearing__ unsafe_string(pointer(smearing__))
   processing_unit__[end] = 0
   pu = GC.@preserve processing_unit__ unsafe_string(pointer(processing_unit__))
   if isempty(pu) pu = "none" end
   return temp__[], smear, kappa__[], tau__[], tol__[], maxiter__[], restart__[], pu
end

### Functions related to the Kpoint set
function create_kset_from_grid(ctx::ContextHandler; k_grid::Vector{Int32}, k_shift::Vector{Int32},
                               use_symmetry::Bool)
   kps = KpointSetHandler(C_NULL)
   error_code__ = Ref{Cint}(0)
   use_symmetry__::Ref{Cuchar} = use_symmetry
   @ccall libpath.sirius_create_kset_from_grid(ctx.handler_ptr::Ptr{Cvoid}, k_grid::Ptr{Cint},
                                               k_shift::Ptr{Cint}, use_symmetry__::Ref{Cuchar},
                                               kps.handler_ptr::Ptr{Cvoid}, error_code__::Ref{Cint})::Cvoid
   if error_code__[] != 0
      error("Sirius.create_kset_from_grid failed with error code", error_code__[])
   end
   return kps
end

function create_kset(ctx::ContextHandler; num_kp::Integer, k_coords::AbstractVector, k_weights::Vector{Float64})
   num_kpoints__ = Ref{Cint}(num_kp)
   kpoints__ = Vector{Cdouble}(undef, 3*num_kp)
   for (ikp, k_coord) in enumerate(k_coords)
      kpoints__[3*(ikp-1)+1:3*ikp] = k_coord[:]
   end
   init_kset__::Ref{Cuchar} = true
   kps = KpointSetHandler(C_NULL)
   error_code__ = Ref{Cint}(0)
   @ccall libpath.sirius_create_kset(ctx.handler_ptr::Ptr{Cvoid}, num_kpoints__::Ref{Cint}, kpoints__::Ptr{Cdouble},
                                     k_weights::Ptr{Cdouble}, init_kset__::Ref{Cuchar}, kps.handler_ptr::Ptr{Cvoid},
                                     error_code__::Ref{Cint})::Cvoid
   if error_code__[] != 0
      error("Sirius.create_kset failed with error code", error_code__[])
   end
   return kps
end

### Functions related to the SIRIUS ground state
function create_ground_state(kps::KpointSetHandler)
   gs = GroundStateHandler(C_NULL)
   error_code__ = Ref{Cint}(0)
   @ccall libpath.sirius_create_ground_state(kps.handler_ptr::Ptr{Cvoid}, gs.handler_ptr::Ptr{Cvoid},
                                         error_code__::Ref{Cint})::Cvoid
   if error_code__[] != 0
      error("Sirius.create_ground_state failed with error code", error_code__[])
   end
   return gs
end

function find_ground_state(gs::GroundStateHandler, initial_guess::Bool, save_state::Bool;
                           density_tol::Float64, energy_tol::Float64,
                           iter_solver_tol::Float64, max_niter::Integer)
   #input
   initial_guess__ = Ref{Cuchar}(initial_guess)
   save_state__ = Ref{Cuchar}(save_state)
   density_tol__ = Ref{Cdouble}(density_tol)
   energy_tol__ = Ref{Cdouble}(energy_tol)
   iter_solver_tol__ = Ref{Cdouble}(iter_solver_tol)
   max_niter__ = Ref{Cint}(max_niter)

   #output
   converged__ = Ref{Cuchar}(false)
   niter__ = Ref{Cint}(-1)
   rho_min__ = Ref{Cdouble}(-1.0)

   error_code__ = Ref{Cint}(0)
   @ccall libpath.sirius_find_ground_state(gs.handler_ptr::Ptr{Cvoid}, density_tol__::Ref{Cdouble},
                                       energy_tol__::Ref{Cdouble}, iter_solver_tol__::Ref{Cdouble},
                                       initial_guess__::Ref{Cuchar}, max_niter__::Ref{Cint},
                                       save_state__::Ref{Cuchar}, converged__::Ref{Cuchar},
                                       niter__::Ref{Cint}, rho_min__::Ref{Cdouble},
                                       error_code__::Ref{Cint})::Cvoid
   if error_code__[] != 0
      error("Sirius.find_ground_state failed with error code", error_code__[])
   end
   return Bool(converged__[]), niter__[], rho_min__[]
end

function nlcg(gs::GroundStateHandler, kps::KpointSetHandler; temp::Float64, smearing::String,
              kappa::Float64, tau::Float64, tol::Float64, maxiter::Integer, restart::Integer,
              processing_unit::String)

   temp__ = Ref{Cdouble}(temp)
   smearing__ = smearing
   kappa__ = Ref{Cdouble}(kappa)
   tau__ = Ref{Cdouble}(tau)
   tol__ = Ref{Cdouble}(tol)
   maxiter__ = Ref{Cint}(maxiter)
   restart__ = Ref{Cint}(restart)
   processing_unit__ = processing_unit

   converged__ = Ref{Cuchar}(false)
   error_code__ = Ref{Cint}(0)
   @ccall libpath.sirius_nlcg_params(gs.handler_ptr::Ptr{Cvoid}, kps.handler_ptr::Ptr{Cvoid},
                                     temp__::Ref{Cdouble}, smearing__::Cstring, kappa__::Ref{Cdouble},
                                     tau__::Ref{Cdouble}, tol__::Ref{Cdouble}, maxiter__::Ref{Cint},
                                     restart__::Ref{Cint}, processing_unit__::Cstring,
                                     converged__::Ref{Cuchar},error_code__::Ref{Cint})::Cvoid
   if error_code__[] != 0
      error("Sirius.nlcg failed with error code", error_code__[])
   end
   return Bool(converged__[])
end

### Functions for querrying the SIRIUS results
function get_num_atoms(gs::GroundStateHandler)
   num_atoms__ = Ref{Cint}(0)
   error_code__ = Ref{Cint}(0)
   @ccall libpath.sirius_get_num_atoms(gs.handler_ptr::Ptr{Cvoid}, num_atoms__::Ref{Cint},
                                       error_code__::Ref{Cint})::Cvoid
   if error_code__[] != 0
      error("Sirius.get_num_atoms failed with error code", error_code__[])
   end
   return num_atoms__[]
end

function get_energy(gs::GroundStateHandler, label::String)
   energy__ = Ref{Cdouble}(0.0)
   error_code__ = Ref{Cint}(0)
   @ccall libpath.sirius_get_energy(gs.handler_ptr::Ptr{Cvoid}, label::Cstring, energy__::Ref{Cdouble},
                                error_code__::Ref{Cint})::Cvoid
   if error_code__[] != 0
      error("Sirius.get_energy failed with error code", error_code__[])
   end
   return energy__[]
end

function get_forces(gs::GroundStateHandler, label::String)
   forces__ = Matrix{Cdouble}(undef, 3, get_num_atoms(gs))
   error_code__ = Ref{Cint}(0)
   @ccall libpath.sirius_get_forces(gs.handler_ptr::Ptr{Cvoid}, label::Cstring, forces__::Ref{Cdouble},
                                error_code__::Ref{Cint})::Cvoid
   if error_code__[] != 0
      error("Sirius.get_forces failed with error code", error_code__[])
   end
   return forces__
end

function get_stress_tensor(gs::GroundStateHandler,  label::String)
   stress__ = Matrix{Cdouble}(undef, 3, 3)
   error_code__ = Ref{Cint}(0)
   @ccall libpath.sirius_get_stress_tensor(gs.handler_ptr::Ptr{Cvoid}, label::Cstring, stress__::Ref{Cdouble},
                                       error_code__::Ref{Cint})::Cvoid
   if error_code__[] != 0
      error("Sirius.get_stress_tensor failed with error code", error_code__[])
   end
   return stress__
end

end
