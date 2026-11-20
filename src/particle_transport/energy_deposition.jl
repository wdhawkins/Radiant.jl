"""
    energy_deposition(cross_sections::Cross_Sections,geometry::Geometry,solvers::Solvers,
    sources::Fixed_Sources,flux::Flux,particles::Vector{<:Particle})

Calculate and extract the energy deposition. 

# Input Argument(s)
- `cross_sections::Cross_Sections`: cross section informations.
- `geometry::Geometry`: geometry informations.
- `solvers::Solvers`: solvers informations.
- `sources::Sources`: sources informations.
- `flux::Flux`: flux informations.
- `particles::Vector{<:Particle}`: particle list.

# Output Argument(s)
- `Dtot::Array{Float64}`: energy deposition per voxel [MeV/g × cmⁿ per particle], where
   n is the geometry dimension. 

# Reference(s)
- Morel et al. (1996), A Hybrid Multigroup/Continuous-Energy Monte Carlo Method for
  Solving the Boltzmann-Fokker-Planck Equation.

"""
function energy_deposition(cross_sections::Cross_Sections,geometry::Geometry,solvers::Solvers,sources::Fixed_Sources,flux::Flux,particles::Vector{<:Particle})

if cross_sections.get_opensn_csda() error("The cross sections exclude the continuous slowing-down contributions to energy deposition (OpenSn CSDA library). Radiant cannot compute energy deposition with them.") end

#----
# Extract geometry data
#----
Ndims = geometry.get_dimension()
Ns = geometry.get_number_of_voxels()
mat = geometry.get_material_per_voxel()

#----
# Calculate the total energy deposition
#----
Dtot = zeros(Ns[1],Ns[2],Ns[3])
for part in particles

    # Extract discrete_ordinates data
    discrete_ordinates = solvers.get_method(part)
    _,isCSD = discrete_ordinates.get_solver_type()

    norm = sources.get_normalization_factor()
    D = zeros(Ns[1],Ns[2],Ns[3])

    #----
    # Extract flux data
    #----
    𝚽l = flux.get_flux(part)
    if isCSD 𝚽cutoff = flux.get_flux_cutoff(part) end

    #----
    # Extract cross sections data
    #----
    Nmat = cross_sections.get_number_of_materials()
    Ng = cross_sections.get_number_of_groups(part)
    ρ = cross_sections.get_densities()

    # Extract energy deposition cross sections
    Σe = cross_sections.get_energy_deposition(part)

    #----
    # Energy deposition calculations
    #----
    for ix in range(1,Ns[1]), iy in range(1,Ns[2]), iz in range(1,Ns[3])

        # In-group energy deposition
        for ig in range(1,Ng)
            D[ix,iy,iz] += Σe[ig,mat[ix,iy,iz]] * 𝚽l[ig,1,1,ix,iy,iz]
        end
        if isCSD D[ix,iy,iz] += Σe[Ng+1,mat[ix,iy,iz]] * 𝚽cutoff[1,1,ix,iy,iz] end

        # Normalization
        D[ix,iy,iz] /= ρ[mat[ix,iy,iz]] * norm

    end
    Dtot += D
end

if Ndims == 1
    return Dtot[:,1,1]
elseif Ndims == 2
    return Dtot[:,:,1]
elseif Ndims == 3
    return Dtot
end

end