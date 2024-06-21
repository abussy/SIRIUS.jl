### Some utilites to parse reference SIRIUS output
using JSON3

function get_ref_energy(fname, label)
   json_string = read(fname, String)
   ref_output = JSON3.read(json_string)
   return ref_output["ground_state"]["energy"][label]
end

function get_ref_forces(fname; label="total")
   if label != "total"
      error("Only total forces availble in reference output")
   end
   json_string = read(fname, String)
   ref_output = JSON3.read(json_string)
   json_forces = ref_output["ground_state"]["forces"]

   #Convert to Julia Matrix and transpose
   natoms = length(json_forces)
   forces = Matrix{Float64}(undef, 3, natoms)
   for iat = 1:natoms
      for dir = 1:3
         forces[dir, iat] = json_forces[iat][dir]
      end
   end
   return forces
end

function get_ref_stress(fname; label="total")
   if label != "total"
      error("Only total stress availble in reference output")
   end
   json_string = read(fname, String)
   ref_output = JSON3.read(json_string)
   json_stress = ref_output["ground_state"]["stress"]

   #Convert to Julia Matrix
   stress = Matrix{Float64}(undef, 3, 3)
   for i = 1:3
      for j = 1:3
         stress[i, j] = json_stress[j][i]
      end
   end
   return stress
end
