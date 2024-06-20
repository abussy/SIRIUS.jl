module SIRIU

using MPI
using SIRIUS_jll

export LibSirius
include("../LibSirius.jl")

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
function create_context_from_json(comm::MPI.Comm, fname::AbstractString)
   ctx = ContextHandler(C_NULL)
   fcomm__::Int32 = comm2f(comm)
   error_code__ = Ref{Cint}(0)
   LibSirius.sirius_create_context_from_json(fcomm__, ctx.handler_ptr, fname, error_code__)
   if error_code__[] != 0
      error("SIRIUS.create_context_from_json failed with error code", error_code__[])
   end
   return ctx
end

function create_kset_from_grid(ctx::ContextHandler; k_grid::Vector{Int32}, k_shift::Vector{Int32},
                               use_symmetry::Bool)
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
__insert_generated_code_here__
end #SIRIUS module
