using Radiant
using Test
import Radiant

# Builds a small library with 8 logarithmic groups down to 10 keV and Legendre order 3: aluminum
# electrons, copper electrons and photons (`coupled`), or copper electrons, photons, and positrons
# (`positrons`). The positron library stays below the pair-production threshold.
function small_library(is_opensn_csda::Bool; coupled::Bool=false, positrons::Bool=false, inelastic=Inelastic_Collision())
    coupled = coupled || positrons
    mat = coupled ? Material("Cu") : Material("Al")
    mat.add_element(coupled ? "Cu" : "Al")
    particles = coupled ? Radiant.Particle[Electron(),Photon()] : Radiant.Particle[Electron()]
    if positrons push!(particles,Positron()) end
    interactions = Radiant.Interaction[inelastic,Elastic_Collision(),Bremsstrahlung()]
    if coupled append!(interactions,[Compton(),Photoelectric()]) end
    if positrons push!(interactions,Annihilation()) end
    cs = Cross_Sections()
    cs.set_source("physics-models")
    cs.set_materials([mat])
    cs.set_particles(particles)
    cs.set_group_structure("log",8,positrons ? 0.5 : 1.0,0.01)
    cs.set_interactions(interactions)
    cs.set_legendre_order(3)
    cs.set_opensn_csda(is_opensn_csda)
    cs.build()
    return cs,particles,mat
end

# Charge deposited per particle removed, in electrons: +1 for an electron, -1 for a positron.
deposited_charge(p::Radiant.Particle) = -Radiant.get_charge(p)

# Returns the largest |c_g - (q_g sigma_t,g - sum_g' q_g' sigma_0(g->g'))| / sigma_t,g over the
# groups of particle `pin`. Charge is conserved only if each charge-deposition coefficient equals
# the charge removed by the collisions in its group.
function charge_identity_error(cs,pin,particles)
    Σt = cs.get_total(pin)[:,1]
    expected = deposited_charge(pin) .* Σt
    for pout in particles
        Σs0 = cs.get_scattering(pin,pout,cs.get_legendre_order())[1,:,:,1]
        expected .-= deposited_charge(pout) .* vec(sum(Σs0,dims=2))
    end
    c = cs.get_charge_deposition(pin)[1:end-1,1]
    return maximum(abs.(c .- expected) ./ Σt)
end

@testset "Radiant.jl" begin

@testset "OpenSn CSDA mode" begin
    std,pstd,_ = small_library(false)
    os,pos,_ = small_library(true)
    e_std,e_os = pstd[1],pos[1]

    # Only the momentum transfer and the deposition data change.
    S = std.get_stopping_powers(e_std)
    @test all(S .> 0)
    @test os.get_stopping_powers(e_os) == S
    @test any(std.get_momentum_transfer(e_std) .!= 0)
    @test all(os.get_momentum_transfer(e_os) .== 0)

    # OpenSn adds the continuous slowing-down deposition itself.
    Ee_std = std.get_energy_deposition(e_std)
    Ee_os = os.get_energy_deposition(e_os)
    @test Ee_std[1:end-1,:] .- Ee_os[1:end-1,:] ≈ S rtol=1e-12
    @test Ee_std[end,1] != 0
    @test Ee_os[end,1] == 0
    Ec_std = std.get_charge_deposition(e_std)
    Ec_os = os.get_charge_deposition(e_os)
    @test Ec_os[1:end-1,:] == Ec_std[1:end-1,:]
    @test Ec_std[end,1] != 0
    @test Ec_os[end,1] == 0

    # Unsupported combinations are rejected.
    @test_throws ErrorException os.write(joinpath(mktempdir(),"library.matxs"),"matxs")
    fp = Inelastic_Collision()
    fp.set_scattering_model("FP")
    @test_throws ErrorException small_library(true;inelastic=fp)
end

@testset "Charge deposition conserves charge" begin
    for (coupled,positrons) in ((false,false),(true,false),(true,true))
        cs,particles,_ = small_library(true;coupled=coupled,positrons=positrons)
        for p in particles
            @test charge_identity_error(cs,p,particles) < 1e-12
        end
    end
end

@testset "FMAC-M round trip" begin
    for is_opensn_csda in (true,false)
        cs,particles,mat = small_library(is_opensn_csda;coupled=true)
        file = joinpath(mktempdir(),"library.txt")
        cs.write(file)
        rd = Cross_Sections()
        rd.set_source("fmac-m")
        rd.set_file(file)
        rd.set_materials([mat])
        rd.set_particles(particles)
        rd.build()
        @test rd.get_opensn_csda() == is_opensn_csda
        # OpenSn CSDA libraries are written with full precision, standard ones with six digits.
        rtol = is_opensn_csda ? 0.0 : 1e-5
        for p in particles
            @test isapprox(rd.get_total(p),cs.get_total(p);rtol=rtol)
            @test isapprox(rd.get_energy_deposition(p)[1:end-1,:],cs.get_energy_deposition(p)[1:end-1,:];rtol=rtol,atol=is_opensn_csda ? 0.0 : 1e-12)
            @test isapprox(rd.get_charge_deposition(p)[1:end-1,:],cs.get_charge_deposition(p)[1:end-1,:];rtol=rtol,atol=is_opensn_csda ? 0.0 : 1e-12)
            for q in particles
                L = cs.get_legendre_order()
                @test isapprox(rd.get_scattering(p,q,L),cs.get_scattering(p,q,L);rtol=rtol,atol=is_opensn_csda ? 0.0 : 1e-8)
            end
            if is_opensn_csda
                @test rd.get_stopping_powers(p) == cs.get_stopping_powers(p)
            end
        end
    end
end

end
nothing
