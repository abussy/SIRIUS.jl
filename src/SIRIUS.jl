module SIRIUS

using MPI
using SIRIUS_jll

export LibSirius
include("LibSirius.jl")

### Hand written wrapper around the SIRIUS handlers (C void pointers)
mutable struct ContextHandler
   handler_ptr::Ref{Ptr{Cvoid}}
end

mutable struct GroundStateHandler
   handler_ptr::Ref{Ptr{Cvoid}}
end

mutable struct KpointSetHandler
   handler_ptr::Ref{Ptr{Cvoid}}
end

mutable struct HamiltonianHandler
   handler_ptr::Ref{Ptr{Cvoid}}
end

### Handler freeing function. Note: not added as finalizer as call order matters
function free_context_handler!(ctx::ContextHandler)
   error_code__ = Ref{Cint}(0)
   LibSirius.sirius_free_object_handler(ctx.handler_ptr, error_code__)
   if error_code__[] != 0
      error("SIRIUS.free_context_handler failed with error code", error_code__[])
   end
end

function free_ground_state_handler!(gs::GroundStateHandler)
   error_code__ = Ref{Cint}(0)
   LibSirius.sirius_free_object_handler(gs.handler_ptr, error_code__)
   if error_code__[] != 0
      error("SIRIUS.free_ground_state_handler failed with error code", error_code__[])
   end
end

function free_kpoint_set_handler!(kps::KpointSetHandler)
   error_code__ = Ref{Cint}(0)
   LibSirius.sirius_free_object_handler(kps.handler_ptr, error_code__)
   if error_code__[] != 0
      error("SIRIUS.free_kpoint_set_handler failed with error code", error_code__[])
   end
end

function free_hamiltonian_handler!(H0::HamiltonianHandler)
   error_code__ = Ref{Cint}(0)
   LibSirius.sirius_free_object_handler(H0.handler_ptr, error_code__)
   if error_code__[] != 0
      error("SIRIUS.free_hamiltonian_handler failed with error code", error_code__[])
   end
end

### Utility function that maps an integer to a C MPI_comm
comm2f(comm::MPI.Comm) = ccall((:MPI_Comm_c2f, MPI.libmpi), Cint, (MPI.MPI_Comm,), comm)

### Hard coded handler creation
function create_context_from_json(comm::MPI.Comm, fname)
   ctx = ContextHandler(C_NULL)
   fcomm__::Int32 = comm2f(comm)
   error_code__ = Ref{Cint}(0)
   LibSirius.sirius_create_context_from_json(fcomm__, ctx.handler_ptr, fname, error_code__)
   if error_code__[] != 0
      error("SIRIUS.create_context_from_json failed with error code", error_code__[])
   end
   return ctx
end

function create_kset_from_grid(ctx::ContextHandler; k_grid, k_shift, use_symmetry)
   kps = KpointSetHandler(C_NULL)
   error_code__ = Ref{Cint}(0)
   use_symmetry__::Ref{Bool} = use_symmetry
   LibSirius.sirius_create_kset_from_grid(ctx.handler_ptr, k_grid, k_shift, use_symmetry__,
                                          kps.handler_ptr, error_code__)
   if error_code__[] != 0
      error("SIRIUS.create_kset_from_grid failed with error code", error_code__[])
   end
   return kps
end

function create_ground_state(kps::KpointSetHandler)
   gs = GroundStateHandler(C_NULL)
   error_code__ = Ref{Cint}(0)
   LibSirius.sirius_create_ground_state(kps.handler_ptr, gs.handler_ptr, error_code__)
   if error_code__[] != 0
      error("SIRIUS.create_ground_state failed with error code", error_code__[])
   end
   return gs
end

### Generated wrapper code around LibSirius
function initialize(call_mpi_init)

   #input arguments (non-array)
   call_mpi_init__ = Ref{Bool}(call_mpi_init)

   #output arguments (non-array)

   error_code__ = Ref{Cint}(0)
   LibSirius.sirius_initialize(call_mpi_init__, error_code__)
   if error_code__[] != 0
      error("SIRIUS.initialize failed with error code", error_code__[])
   end

end

function finalize(; call_mpi_fin=nothing, call_device_reset=nothing, call_fftw_fin=nothing)

   #input arguments (non-array)
   if isnothing(call_mpi_fin)
      call_mpi_fin__ = Ptr{Bool}(C_NULL)
   else
      call_mpi_fin__ = Ref{Bool}(call_mpi_fin)
   end

   if isnothing(call_device_reset)
      call_device_reset__ = Ptr{Bool}(C_NULL)
   else
      call_device_reset__ = Ref{Bool}(call_device_reset)
   end

   if isnothing(call_fftw_fin)
      call_fftw_fin__ = Ptr{Bool}(C_NULL)
   else
      call_fftw_fin__ = Ref{Bool}(call_fftw_fin)
   end


   #output arguments (non-array)

   error_code__ = Ref{Cint}(0)
   LibSirius.sirius_finalize(call_mpi_fin__, call_device_reset__, call_fftw_fin__, error_code__)
   if error_code__[] != 0
      error("SIRIUS.finalize failed with error code", error_code__[])
   end

end

function is_initialized()

   #input arguments (non-array)

   #output arguments (non-array)
   status__ = Ref{Bool}(0)

   error_code__ = Ref{Cint}(0)
   LibSirius.sirius_is_initialized(status__, error_code__)
   if error_code__[] != 0
      error("SIRIUS.is_initialized failed with error code", error_code__[])
   end

   return status__[]
end

function initialize_context(ctx_handler)

   #input arguments (non-array)

   #output arguments (non-array)

   error_code__ = Ref{Cint}(0)
   LibSirius.sirius_initialize_context(ctx_handler.handler_ptr, error_code__)
   if error_code__[] != 0
      error("SIRIUS.initialize_context failed with error code", error_code__[])
   end

end

function get_kp_params_from_ctx!(ctx_handler, k_grid, k_shift)

   #input arguments (non-array)

   #output arguments (non-array)
   use_symmetry__ = Ref{Bool}(0)

   error_code__ = Ref{Cint}(0)
   LibSirius.sirius_get_kp_params_from_ctx(ctx_handler.handler_ptr, k_grid, k_shift, use_symmetry__, error_code__)
   if error_code__[] != 0
      error("SIRIUS.get_kp_params_from_ctx failed with error code", error_code__[])
   end

   return use_symmetry__[]
end

function get_scf_params_from_ctx(ctx_handler)

   #input arguments (non-array)

   #output arguments (non-array)
   density_tol____ = Ref{Cdouble}(0)
   energy_tol____ = Ref{Cdouble}(0)
   iter_solver_tol__ = Ref{Cdouble}(0)
   max_niter__ = Ref{Cint}(0)

   error_code__ = Ref{Cint}(0)
   LibSirius.sirius_get_scf_params_from_ctx(ctx_handler.handler_ptr, density_tol____, energy_tol____, iter_solver_tol__, max_niter__, error_code__)
   if error_code__[] != 0
      error("SIRIUS.get_scf_params_from_ctx failed with error code", error_code__[])
   end

   return density_tol____[], energy_tol____[], iter_solver_tol__[], max_niter__[]
end

function find_ground_state(gs_handler; density_tol=nothing, energy_tol=nothing, iter_solver_tol=nothing, initial_guess=nothing, max_niter=nothing, save_state=nothing)

   #input arguments (non-array)
   if isnothing(density_tol)
      density_tol__ = Ptr{Cdouble}(C_NULL)
   else
      density_tol__ = Ref{Cdouble}(density_tol)
   end

   if isnothing(energy_tol)
      energy_tol__ = Ptr{Cdouble}(C_NULL)
   else
      energy_tol__ = Ref{Cdouble}(energy_tol)
   end

   if isnothing(iter_solver_tol)
      iter_solver_tol__ = Ptr{Cdouble}(C_NULL)
   else
      iter_solver_tol__ = Ref{Cdouble}(iter_solver_tol)
   end

   if isnothing(initial_guess)
      initial_guess__ = Ptr{Bool}(C_NULL)
   else
      initial_guess__ = Ref{Bool}(initial_guess)
   end

   if isnothing(max_niter)
      max_niter__ = Ptr{Cint}(C_NULL)
   else
      max_niter__ = Ref{Cint}(max_niter)
   end

   if isnothing(save_state)
      save_state__ = Ptr{Bool}(C_NULL)
   else
      save_state__ = Ref{Bool}(save_state)
   end


   #output arguments (non-array)
   converged__ = Ref{Bool}(0)
   niter__ = Ref{Cint}(0)
   rho_min__ = Ref{Cdouble}(0)

   error_code__ = Ref{Cint}(0)
   LibSirius.sirius_find_ground_state(gs_handler.handler_ptr, density_tol__, energy_tol__, iter_solver_tol__, initial_guess__, max_niter__, save_state__, converged__, niter__, rho_min__, error_code__)
   if error_code__[] != 0
      error("SIRIUS.find_ground_state failed with error code", error_code__[])
   end

   return converged__[], niter__[], rho_min__[]
end

function get_num_atoms(gs_handler)

   #input arguments (non-array)

   #output arguments (non-array)
   num_atoms__ = Ref{Cint}(0)

   error_code__ = Ref{Cint}(0)
   LibSirius.sirius_get_num_atoms(gs_handler.handler_ptr, num_atoms__, error_code__)
   if error_code__[] != 0
      error("SIRIUS.get_num_atoms failed with error code", error_code__[])
   end

   return num_atoms__[]
end

function get_energy(gs_handler, label)

   #input arguments (non-array)

   #output arguments (non-array)
   energy__ = Ref{Cdouble}(0)

   error_code__ = Ref{Cint}(0)
   LibSirius.sirius_get_energy(gs_handler.handler_ptr, label, energy__, error_code__)
   if error_code__[] != 0
      error("SIRIUS.get_energy failed with error code", error_code__[])
   end

   return energy__[]
end

function get_forces!(gs_handler, label, forces)

   #input arguments (non-array)

   #output arguments (non-array)

   error_code__ = Ref{Cint}(0)
   LibSirius.sirius_get_forces(gs_handler.handler_ptr, label, forces, error_code__)
   if error_code__[] != 0
      error("SIRIUS.get_forces failed with error code", error_code__[])
   end

end

function get_stress_tensor!(gs_handler, label, stress_tensor)

   #input arguments (non-array)

   #output arguments (non-array)

   error_code__ = Ref{Cint}(0)
   LibSirius.sirius_get_stress_tensor(gs_handler.handler_ptr, label, stress_tensor, error_code__)
   if error_code__[] != 0
      error("SIRIUS.get_stress_tensor failed with error code", error_code__[])
   end

end


end #SIRIUS module
