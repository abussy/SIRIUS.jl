module SIRIUS

using MPI
using OpenBLAS32
using SIRIUS_jll
using JSON3
using Libdl

export LibSirius
include("LibSirius.jl")

### Users can silence SIRIUS stdout output by calling SIRIUS.output_mode(true)
const silent = Ref(false)
function output_mode(;make_silent=true)
   silent[] = make_silent
end

function get_outstream()
   if silent[]
      return devnull
   else
      return nothing
   end
end

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
   redirect_stdio(;stdout=get_outstream()) do
      LibSirius.sirius_free_object_handler(ctx.handler_ptr, error_code__)
      Base.Libc.flush_cstdio()
   end
   if error_code__[] != 0
      error("SIRIUS.free_context_handler failed with error code", error_code__[])
   end
end

function free_ground_state_handler!(gs::GroundStateHandler)
   error_code__ = Ref{Cint}(0)
   redirect_stdio(;stdout=get_outstream()) do
      LibSirius.sirius_free_object_handler(gs.handler_ptr, error_code__)
      Base.Libc.flush_cstdio()
   end
   if error_code__[] != 0
      error("SIRIUS.free_ground_state_handler failed with error code", error_code__[])
   end
end

function free_kpoint_set_handler!(kps::KpointSetHandler)
   error_code__ = Ref{Cint}(0)
   redirect_stdio(;stdout=get_outstream()) do
      LibSirius.sirius_free_object_handler(kps.handler_ptr, error_code__)
      Base.Libc.flush_cstdio()
   end
   if error_code__[] != 0
      error("SIRIUS.free_kpoint_set_handler failed with error code", error_code__[])
   end
end

function free_hamiltonian_handler!(H0::HamiltonianHandler)
   error_code__ = Ref{Cint}(0)
   redirect_stdio(;stdout=get_outstream()) do
      LibSirius.sirius_free_object_handler(H0.handler_ptr, error_code__)
      Base.Libc.flush_cstdio()
   end
   if error_code__[] != 0
      error("SIRIUS.free_hamiltonian_handler failed with error code", error_code__[])
   end
end

### Utility function that maps an integer to a C MPI_comm
function comm2f(comm::MPI.Comm)
   # Some MPI libraries (cray-mpich) do not define MPI_Comm_c2f, in this case
   # the C MPI communicator is just an integer in the first place
   handle = Libdl.dlopen(MPI.libmpi)
   if isnothing(Libdl.dlsym(handle, :MPI_Comm_c2f; throw_error = false))
      comm.val
   else
      ccall((:MPI_Comm_c2f, MPI.libmpi), Cint, (MPI.MPI_Comm,), comm)
   end
end

### Hard coded handler creation
function create_context_from_json(comm::MPI.Comm, fname)
   ctx = ContextHandler(C_NULL)
   fcomm__::Int32 = comm2f(comm)
   error_code__ = Ref{Cint}(0)
   redirect_stdio(;stdout=get_outstream()) do
      LibSirius.sirius_create_context_from_json(fcomm__, ctx.handler_ptr, fname, error_code__)
      Base.Libc.flush_cstdio()
   end
   if error_code__[] != 0
      error("SIRIUS.create_context_from_json failed with error code", error_code__[])
   end
   return ctx
end

function create_context_from_dict(comm::MPI.Comm, dict)
   json_string = JSON3.write(dict)
   create_context_from_json(comm, json_string)
end

function create_kset_from_grid(ctx::ContextHandler; k_grid, k_shift, use_symmetry)
   kps = KpointSetHandler(C_NULL)
   error_code__ = Ref{Cint}(0)
   use_symmetry__::Ref{Bool} = use_symmetry
   redirect_stdio(;stdout=get_outstream()) do
      LibSirius.sirius_create_kset_from_grid(ctx.handler_ptr, k_grid, k_shift, use_symmetry__,
                                             kps.handler_ptr, error_code__)
      Base.Libc.flush_cstdio()
   end
   if error_code__[] != 0
      error("SIRIUS.create_kset_from_grid failed with error code", error_code__[])
   end
   return kps
end

function create_kset(ctx::ContextHandler; num_kp, k_coords, k_weights, init_kset=true)
   kps = KpointSetHandler(C_NULL)
   num_kpoints__ = Ref{Cint}(num_kp)
   kpoints__ = Vector{Cdouble}(undef, 3*num_kp)
   for (ikp, k_coord) in enumerate(k_coords)
      kpoints__[3*(ikp-1)+1:3*ikp] = k_coord[:]
   end
   init_kset__::Ref{Bool} = init_kset
   error_code__ = Ref{Cint}(0)      
   redirect_stdio(;stdout=get_outstream()) do
      LibSirius.sirius_create_kset(ctx.handler_ptr, num_kpoints__, kpoints__, k_weights, init_kset__,
                                   kps.handler_ptr, error_code__)
      Base.Libc.flush_cstdio()
   end
   if error_code__[] != 0
      error("SIRIUS.create_kset failed with error code", error_code__[])
   end
   return kps
end

function create_ground_state(kps::KpointSetHandler)
   gs = GroundStateHandler(C_NULL)
   error_code__ = Ref{Cint}(0)
   redirect_stdio(;stdout=get_outstream()) do
      LibSirius.sirius_create_ground_state(kps.handler_ptr, gs.handler_ptr, error_code__)
      Base.Libc.flush_cstdio()
   end
   if error_code__[] != 0
      error("SIRIUS.create_ground_state failed with error code", error_code__[])
   end
   return gs
end

function create_hamiltonian(gs::GroundStateHandler)
   H0 = HamiltonianHandler(C_NULL)
   error_code__ = Ref{Cint}(0)
   redirect_stdio(;stdout=get_outstream()) do
      LibSirius.sirius_create_hamiltonian(gs.handler_ptr, H0.handler_ptr, error_code__)
      Base.Libc.flush_cstdio()
   end
   if error_code__[] != 0
      error("SIRIUS.create_hamiltonian failed with error code", error_code__[])
   end
   return H0
end

### Generated wrapper code around LibSirius
function initialize(call_mpi_init)

   #input arguments (non-array)
   call_mpi_init__ = Ref{Bool}(call_mpi_init)

   #output arguments (non-array)

   error_code__ = Ref{Cint}(0)
   redirect_stdio(;stdout=get_outstream()) do
      LibSirius.sirius_initialize(call_mpi_init__, error_code__)
      Base.Libc.flush_cstdio()
   end
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
   redirect_stdio(;stdout=get_outstream()) do
      LibSirius.sirius_finalize(call_mpi_fin__, call_device_reset__, call_fftw_fin__, error_code__)
      Base.Libc.flush_cstdio()
   end
   if error_code__[] != 0
      error("SIRIUS.finalize failed with error code", error_code__[])
   end

end

function is_initialized()

   #input arguments (non-array)

   #output arguments (non-array)
   status__ = Ref{Bool}(0)

   error_code__ = Ref{Cint}(0)
   redirect_stdio(;stdout=get_outstream()) do
      LibSirius.sirius_is_initialized(status__, error_code__)
      Base.Libc.flush_cstdio()
   end
   if error_code__[] != 0
      error("SIRIUS.is_initialized failed with error code", error_code__[])
   end

   return status__[]
end

function initialize_context(ctx_handler)

   #input arguments (non-array)

   #output arguments (non-array)

   error_code__ = Ref{Cint}(0)
   redirect_stdio(;stdout=get_outstream()) do
      LibSirius.sirius_initialize_context(ctx_handler.handler_ptr, error_code__)
      Base.Libc.flush_cstdio()
   end
   if error_code__[] != 0
      error("SIRIUS.initialize_context failed with error code", error_code__[])
   end

end

function get_kp_params_from_ctx!(ctx_handler, k_grid, k_shift)

   #input arguments (non-array)

   #output arguments (non-array)
   use_symmetry__ = Ref{Bool}(0)

   error_code__ = Ref{Cint}(0)
   redirect_stdio(;stdout=get_outstream()) do
      LibSirius.sirius_get_kp_params_from_ctx(ctx_handler.handler_ptr, k_grid, k_shift, use_symmetry__, error_code__)
      Base.Libc.flush_cstdio()
   end
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
   redirect_stdio(;stdout=get_outstream()) do
      LibSirius.sirius_get_scf_params_from_ctx(ctx_handler.handler_ptr, density_tol____, energy_tol____, iter_solver_tol__, max_niter__, error_code__)
      Base.Libc.flush_cstdio()
   end
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
   redirect_stdio(;stdout=get_outstream()) do
      LibSirius.sirius_find_ground_state(gs_handler.handler_ptr, density_tol__, energy_tol__, iter_solver_tol__, initial_guess__, max_niter__, save_state__, converged__, niter__, rho_min__, error_code__)
      Base.Libc.flush_cstdio()
   end
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
   redirect_stdio(;stdout=get_outstream()) do
      LibSirius.sirius_get_num_atoms(gs_handler.handler_ptr, num_atoms__, error_code__)
      Base.Libc.flush_cstdio()
   end
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
   redirect_stdio(;stdout=get_outstream()) do
      LibSirius.sirius_get_energy(gs_handler.handler_ptr, label, energy__, error_code__)
      Base.Libc.flush_cstdio()
   end
   if error_code__[] != 0
      error("SIRIUS.get_energy failed with error code", error_code__[])
   end

   return energy__[]
end

function get_forces!(gs_handler, label, forces)

   #input arguments (non-array)

   #output arguments (non-array)

   error_code__ = Ref{Cint}(0)
   redirect_stdio(;stdout=get_outstream()) do
      LibSirius.sirius_get_forces(gs_handler.handler_ptr, label, forces, error_code__)
      Base.Libc.flush_cstdio()
   end
   if error_code__[] != 0
      error("SIRIUS.get_forces failed with error code", error_code__[])
   end

end

function get_stress_tensor!(gs_handler, label, stress_tensor)

   #input arguments (non-array)

   #output arguments (non-array)

   error_code__ = Ref{Cint}(0)
   redirect_stdio(;stdout=get_outstream()) do
      LibSirius.sirius_get_stress_tensor(gs_handler.handler_ptr, label, stress_tensor, error_code__)
      Base.Libc.flush_cstdio()
   end
   if error_code__[] != 0
      error("SIRIUS.get_stress_tensor failed with error code", error_code__[])
   end

end

function initialize_kset(ks_handler; count=C_NULL)

   #input arguments (non-array)

   #output arguments (non-array)

   error_code__ = Ref{Cint}(0)
   redirect_stdio(;stdout=get_outstream()) do
      LibSirius.sirius_initialize_kset(ks_handler.handler_ptr, count, error_code__)
      Base.Libc.flush_cstdio()
   end
   if error_code__[] != 0
      error("SIRIUS.initialize_kset failed with error code", error_code__[])
   end

end

function set_periodic_function(gs_handler, label; f_mt=C_NULL, lmmax=nothing, nrmtmax=nothing, num_atoms=nothing, f_rg=C_NULL, size_x=nothing, size_y=nothing, size_z=nothing, offset_z=nothing)

   #input arguments (non-array)
   if isnothing(lmmax)
      lmmax__ = Ptr{Cint}(C_NULL)
   else
      lmmax__ = Ref{Cint}(lmmax)
   end

   if isnothing(nrmtmax)
      nrmtmax__ = Ptr{Cint}(C_NULL)
   else
      nrmtmax__ = Ref{Cint}(nrmtmax)
   end

   if isnothing(num_atoms)
      num_atoms__ = Ptr{Cint}(C_NULL)
   else
      num_atoms__ = Ref{Cint}(num_atoms)
   end

   if isnothing(size_x)
      size_x__ = Ptr{Cint}(C_NULL)
   else
      size_x__ = Ref{Cint}(size_x)
   end

   if isnothing(size_y)
      size_y__ = Ptr{Cint}(C_NULL)
   else
      size_y__ = Ref{Cint}(size_y)
   end

   if isnothing(size_z)
      size_z__ = Ptr{Cint}(C_NULL)
   else
      size_z__ = Ref{Cint}(size_z)
   end

   if isnothing(offset_z)
      offset_z__ = Ptr{Cint}(C_NULL)
   else
      offset_z__ = Ref{Cint}(offset_z)
   end


   #output arguments (non-array)

   error_code__ = Ref{Cint}(0)
   redirect_stdio(;stdout=get_outstream()) do
      LibSirius.sirius_set_periodic_function(gs_handler.handler_ptr, label, f_mt, lmmax__, nrmtmax__, num_atoms__, f_rg, size_x__, size_y__, size_z__, offset_z__, error_code__)
      Base.Libc.flush_cstdio()
   end
   if error_code__[] != 0
      error("SIRIUS.set_periodic_function failed with error code", error_code__[])
   end

end

function get_periodic_function!(gs_handler, label; f_mt=C_NULL, lmmax=nothing, nrmtmax=nothing, num_atoms=nothing, f_rg=C_NULL, size_x=nothing, size_y=nothing, size_z=nothing, offset_z=nothing)

   #input arguments (non-array)
   if isnothing(lmmax)
      lmmax__ = Ptr{Cint}(C_NULL)
   else
      lmmax__ = Ref{Cint}(lmmax)
   end

   if isnothing(nrmtmax)
      nrmtmax__ = Ptr{Cint}(C_NULL)
   else
      nrmtmax__ = Ref{Cint}(nrmtmax)
   end

   if isnothing(num_atoms)
      num_atoms__ = Ptr{Cint}(C_NULL)
   else
      num_atoms__ = Ref{Cint}(num_atoms)
   end

   if isnothing(size_x)
      size_x__ = Ptr{Cint}(C_NULL)
   else
      size_x__ = Ref{Cint}(size_x)
   end

   if isnothing(size_y)
      size_y__ = Ptr{Cint}(C_NULL)
   else
      size_y__ = Ref{Cint}(size_y)
   end

   if isnothing(size_z)
      size_z__ = Ptr{Cint}(C_NULL)
   else
      size_z__ = Ref{Cint}(size_z)
   end

   if isnothing(offset_z)
      offset_z__ = Ptr{Cint}(C_NULL)
   else
      offset_z__ = Ref{Cint}(offset_z)
   end


   #output arguments (non-array)

   error_code__ = Ref{Cint}(0)
   redirect_stdio(;stdout=get_outstream()) do
      LibSirius.sirius_get_periodic_function(gs_handler.handler_ptr, label, f_mt, lmmax__, nrmtmax__, num_atoms__, f_rg, size_x__, size_y__, size_z__, offset_z__, error_code__)
      Base.Libc.flush_cstdio()
   end
   if error_code__[] != 0
      error("SIRIUS.get_periodic_function failed with error code", error_code__[])
   end

end

function fft_transform(gs_handler, label, direction)

   #input arguments (non-array)
   direction__ = Ref{Cint}(direction)

   #output arguments (non-array)

   error_code__ = Ref{Cint}(0)
   redirect_stdio(;stdout=get_outstream()) do
      LibSirius.sirius_fft_transform(gs_handler.handler_ptr, label, direction__, error_code__)
      Base.Libc.flush_cstdio()
   end
   if error_code__[] != 0
      error("SIRIUS.fft_transform failed with error code", error_code__[])
   end

end

function diagonalize_hamiltonian(ctx_handler, gs_handler, H0_handler, iter_solver_tol, max_steps; converge_by_energy=nothing, exact_diagonalization=nothing)

   #input arguments (non-array)
   iter_solver_tol__ = Ref{Cdouble}(iter_solver_tol)
   max_steps__ = Ref{Cint}(max_steps)
   if isnothing(converge_by_energy)
      converge_by_energy__ = Ptr{Cint}(C_NULL)
   else
      converge_by_energy__ = Ref{Cint}(converge_by_energy)
   end

   if isnothing(exact_diagonalization)
      exact_diagonalization__ = Ptr{Bool}(C_NULL)
   else
      exact_diagonalization__ = Ref{Bool}(exact_diagonalization)
   end


   #output arguments (non-array)
   converged__ = Ref{Bool}(0)
   niter__ = Ref{Cint}(0)

   error_code__ = Ref{Cint}(0)
   redirect_stdio(;stdout=get_outstream()) do
      LibSirius.sirius_diagonalize_hamiltonian(ctx_handler.handler_ptr, gs_handler.handler_ptr, H0_handler.handler_ptr, iter_solver_tol__, max_steps__, converge_by_energy__, exact_diagonalization__, converged__, niter__, error_code__)
      Base.Libc.flush_cstdio()
   end
   if error_code__[] != 0
      error("SIRIUS.diagonalize_hamiltonian failed with error code", error_code__[])
   end

   return converged__[], niter__[]
end

function get_band_energies!(ks_handler, ik, ispn, band_energies)

   #input arguments (non-array)
   ik__ = Ref{Cint}(ik)
   ispn__ = Ref{Cint}(ispn)

   #output arguments (non-array)

   error_code__ = Ref{Cint}(0)
   redirect_stdio(;stdout=get_outstream()) do
      LibSirius.sirius_get_band_energies(ks_handler.handler_ptr, ik__, ispn__, band_energies, error_code__)
      Base.Libc.flush_cstdio()
   end
   if error_code__[] != 0
      error("SIRIUS.get_band_energies failed with error code", error_code__[])
   end

end

function get_psi!(ks_handler, ik, ispin, psi)

   #input arguments (non-array)
   ik__ = Ref{Cint}(ik)
   ispin__ = Ref{Cint}(ispin)

   #output arguments (non-array)

   error_code__ = Ref{Cint}(0)
   redirect_stdio(;stdout=get_outstream()) do
      LibSirius.sirius_get_psi(ks_handler.handler_ptr, ik__, ispin__, psi, error_code__)
      Base.Libc.flush_cstdio()
   end
   if error_code__[] != 0
      error("SIRIUS.get_psi failed with error code", error_code__[])
   end

end

function set_band_energies(ks_handler, ik, ispn, band_energies)

   #input arguments (non-array)
   ik__ = Ref{Cint}(ik)
   ispn__ = Ref{Cint}(ispn)

   #output arguments (non-array)

   error_code__ = Ref{Cint}(0)
   redirect_stdio(;stdout=get_outstream()) do
      LibSirius.sirius_set_band_energies(ks_handler.handler_ptr, ik__, ispn__, band_energies, error_code__)
      Base.Libc.flush_cstdio()
   end
   if error_code__[] != 0
      error("SIRIUS.set_band_energies failed with error code", error_code__[])
   end

end

function set_band_occupancies(ks_handler, ik, ispn, band_occupancies)

   #input arguments (non-array)
   ik__ = Ref{Cint}(ik)
   ispn__ = Ref{Cint}(ispn)

   #output arguments (non-array)

   error_code__ = Ref{Cint}(0)
   redirect_stdio(;stdout=get_outstream()) do
      LibSirius.sirius_set_band_occupancies(ks_handler.handler_ptr, ik__, ispn__, band_occupancies, error_code__)
      Base.Libc.flush_cstdio()
   end
   if error_code__[] != 0
      error("SIRIUS.set_band_occupancies failed with error code", error_code__[])
   end

end

function get_band_occupancies!(ks_handler, ik, ispn, band_occupancies)

   #input arguments (non-array)
   ik__ = Ref{Cint}(ik)
   ispn__ = Ref{Cint}(ispn)

   #output arguments (non-array)

   error_code__ = Ref{Cint}(0)
   redirect_stdio(;stdout=get_outstream()) do
      LibSirius.sirius_get_band_occupancies(ks_handler.handler_ptr, ik__, ispn__, band_occupancies, error_code__)
      Base.Libc.flush_cstdio()
   end
   if error_code__[] != 0
      error("SIRIUS.get_band_occupancies failed with error code", error_code__[])
   end

end

function generate_density(gs_handler; add_core=nothing, transform_to_rg=nothing, paw_only=nothing)

   #input arguments (non-array)
   if isnothing(add_core)
      add_core__ = Ptr{Bool}(C_NULL)
   else
      add_core__ = Ref{Bool}(add_core)
   end

   if isnothing(transform_to_rg)
      transform_to_rg__ = Ptr{Bool}(C_NULL)
   else
      transform_to_rg__ = Ref{Bool}(transform_to_rg)
   end

   if isnothing(paw_only)
      paw_only__ = Ptr{Bool}(C_NULL)
   else
      paw_only__ = Ref{Bool}(paw_only)
   end


   #output arguments (non-array)

   error_code__ = Ref{Cint}(0)
   redirect_stdio(;stdout=get_outstream()) do
      LibSirius.sirius_generate_density(gs_handler.handler_ptr, add_core__, transform_to_rg__, paw_only__, error_code__)
      Base.Libc.flush_cstdio()
   end
   if error_code__[] != 0
      error("SIRIUS.generate_density failed with error code", error_code__[])
   end

end

function find_band_occupancies(ks_handler)

   #input arguments (non-array)

   #output arguments (non-array)

   error_code__ = Ref{Cint}(0)
   redirect_stdio(;stdout=get_outstream()) do
      LibSirius.sirius_find_band_occupancies(ks_handler.handler_ptr, error_code__)
      Base.Libc.flush_cstdio()
   end
   if error_code__[] != 0
      error("SIRIUS.find_band_occupancies failed with error code", error_code__[])
   end

end

function set_num_bands(ctx_handler, num_bands)

   #input arguments (non-array)
   num_bands__ = Ref{Cint}(num_bands)

   #output arguments (non-array)

   error_code__ = Ref{Cint}(0)
   redirect_stdio(;stdout=get_outstream()) do
      LibSirius.sirius_set_num_bands(ctx_handler.handler_ptr, num_bands__, error_code__)
      Base.Libc.flush_cstdio()
   end
   if error_code__[] != 0
      error("SIRIUS.set_num_bands failed with error code", error_code__[])
   end

end

function generate_initial_density(gs_handler)

   #input arguments (non-array)

   #output arguments (non-array)

   error_code__ = Ref{Cint}(0)
   redirect_stdio(;stdout=get_outstream()) do
      LibSirius.sirius_generate_initial_density(gs_handler.handler_ptr, error_code__)
      Base.Libc.flush_cstdio()
   end
   if error_code__[] != 0
      error("SIRIUS.generate_initial_density failed with error code", error_code__[])
   end

end

function get_gkvec!(ks_handler, ik, gvec)

   #input arguments (non-array)
   ik__ = Ref{Cint}(ik)

   #output arguments (non-array)

   error_code__ = Ref{Cint}(0)
   redirect_stdio(;stdout=get_outstream()) do
      LibSirius.sirius_get_gkvec(ks_handler.handler_ptr, ik__, gvec, error_code__)
      Base.Libc.flush_cstdio()
   end
   if error_code__[] != 0
      error("SIRIUS.get_gkvec failed with error code", error_code__[])
   end

end

function set_energy_fermi(ks_handler, energy_fermi)

   #input arguments (non-array)
   energy_fermi__ = Ref{Cdouble}(energy_fermi)

   #output arguments (non-array)

   error_code__ = Ref{Cint}(0)
   redirect_stdio(;stdout=get_outstream()) do
      LibSirius.sirius_set_energy_fermi(ks_handler.handler_ptr, energy_fermi__, error_code__)
      Base.Libc.flush_cstdio()
   end
   if error_code__[] != 0
      error("SIRIUS.set_energy_fermi failed with error code", error_code__[])
   end

end

function initialize_subspace(gs_handler, ks_handler)

   #input arguments (non-array)

   #output arguments (non-array)

   error_code__ = Ref{Cint}(0)
   redirect_stdio(;stdout=get_outstream()) do
      LibSirius.sirius_initialize_subspace(gs_handler.handler_ptr, ks_handler.handler_ptr, error_code__)
      Base.Libc.flush_cstdio()
   end
   if error_code__[] != 0
      error("SIRIUS.initialize_subspace failed with error code", error_code__[])
   end

end

function set_atom_vector_field(ctx_handler, ia, vector_field)

   #input arguments (non-array)
   ia__ = Ref{Cint}(ia)

   #output arguments (non-array)

   error_code__ = Ref{Cint}(0)
   redirect_stdio(;stdout=get_outstream()) do
      LibSirius.sirius_set_atom_vector_field(ctx_handler.handler_ptr, ia__, vector_field, error_code__)
      Base.Libc.flush_cstdio()
   end
   if error_code__[] != 0
      error("SIRIUS.set_atom_vector_field failed with error code", error_code__[])
   end

end

function nlcg_params(gs_handler, ks_handler, temp, smearing, kappa, tau, tol, maxiter, restart, processing_unit)

   #input arguments (non-array)
   temp__ = Ref{Cdouble}(temp)
   kappa__ = Ref{Cdouble}(kappa)
   tau__ = Ref{Cdouble}(tau)
   tol__ = Ref{Cdouble}(tol)
   maxiter__ = Ref{Cint}(maxiter)
   restart__ = Ref{Cint}(restart)

   #output arguments (non-array)
   converged__ = Ref{Bool}(0)

   error_code__ = Ref{Cint}(0)
   redirect_stdio(;stdout=get_outstream()) do
      LibSirius.sirius_nlcg_params(gs_handler.handler_ptr, ks_handler.handler_ptr, temp__, smearing, kappa__, tau__, tol__, maxiter__, restart__, processing_unit, converged__, error_code__)
      Base.Libc.flush_cstdio()
   end
   if error_code__[] != 0
      error("SIRIUS.nlcg_params failed with error code", error_code__[])
   end

   return converged__[]
end

function apply_h(ks_handler, H0_handler, ik, nbands, phi, hpsi)

   #input arguments (non-array)
   ik__ = Ref{Cint}(ik)
   nbands__ = Ref{Cint}(nbands)

   #output arguments (non-array)

   error_code__ = Ref{Cint}(0)
   redirect_stdio(;stdout=get_outstream()) do
      LibSirius.sirius_apply_h(ks_handler.handler_ptr, H0_handler.handler_ptr, ik__, nbands__, phi, hpsi, error_code__)
      Base.Libc.flush_cstdio()
   end
   if error_code__[] != 0
      error("SIRIUS.apply_h failed with error code", error_code__[])
   end

end


end #SIRIUS module
