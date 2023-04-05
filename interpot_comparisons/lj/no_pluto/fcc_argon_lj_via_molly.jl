using Molly
using InteratomicPotentials
using AtomsIO # .CondaPkg/ comes from this
using Unitful
using UnitfulAtomic
using AtomsBase
 
# Wrapper struct enabling potentials implemented in InteratomicPotentials.jl to be used in Molly
struct InteratomicPotentialInter{P<:AbstractPotential}
	potential::P
    eaf::Function
    # TODO: should also use the energy cache Jeremiah had, but ignoring that for now
    # TODO: Default constructur with energy_and_force probably
end

# Extending Molly.forces with above struct
function Molly.forces(inter::InteratomicPotentialInter, sys::AbstractSystem, neighbors=nothing)
    eandf = inter.eaf(sys, inter.potential)
    eandf.f
end

function molly_params(sys::AtomsBase.FlexibleSystem)
    coords = [atm.position for atm in sys.particles]

    atoms = [Molly.Atom(mass=atm.atomic_mass) for atm in sys.particles]
    atoms_data = [AtomData(element=string(atm.atomic_symbol)) for atm in sys.particles]
    velocities = [atm.velocity for atm in sys.particles]

    # TODO: Generalize this. Currently very fragile, assumes you read in a cubic box. Otherwise need to use TriclinicBoundary, but some limitations...
    # Also by default, assumes directions are periodic..., ignoring sys.boundary_conditions
    boundary=CubicBoundary(sys.box[1][1],sys.box[2][2],sys.box[3][3]) # FRAGILE! 

    return NamedTuple((atoms=atoms, atoms_data=atoms_data, coords=coords, velocities = velocities, boundary=boundary))
end 


############################################################

# load system with AtomsIO
sys = load_system("../ref_config/dump_final.extxyz")


# set up InteratomicPotentials LennardJones, based off of 10.1103/PhysRevB.54.340 
ϵ = 0.01032u"eV"
σ = 3.405u"Å"
rcut = 8.51u"Å"
species = [:Ar]
lj_p = InteratomicPotentials.LennardJones(ϵ,σ,rcut,species)

inter_lj = InteratomicPotentialInter(lj_p,energy_and_force)
general_inters = (inter_lj,)

# Check force with regular sys obtained with AtomsIO
f_p = Molly.forces(inter_lj, sys)
#f_p = [uconvert.(u"eV/Å", fi) for fi in f_p]

mp = molly_params(sys)

# Both options here cause the dimension error when pushing to force logger
m_sys = System(;mp...,
            general_inters=general_inters, 
            loggers=(force=ForceLogger(1),),
            #force_units=u"eV/Å",
            #force_units=NoUnits,
            )


#### run zero simulation
simulator = VelocityVerlet(
    dt=0.001u"ps",
    coupling=AndersenThermostat(80u"K", 1.0u"ps"),
)

simulate!(m_sys,simulator,0)
"""
if not using default forces, get the following error:
    ERROR: DimensionError: kJ nm⁻¹ mol⁻¹ and 3.85209493530676e-5 eV Å⁻¹ are not dimensionally compatible.
    Stacktrace:
      [1] convert(#unused#::Type{Quantity{Float64, 𝐋 𝐌 𝐍⁻¹ 𝐓⁻², Unitful.FreeUnits{(kJ, nm⁻¹, mol⁻¹), 𝐋 𝐌 𝐍⁻¹ 𝐓⁻², nothing}}}, x::Quantity{Float64, 𝐋 𝐌 𝐓⁻², Unitful.FreeUnits{(Å⁻¹, eV), 𝐋 𝐌 𝐓⁻², nothing}})
        @ Unitful ~/.julia/packages/Unitful/G8F13/src/conversion.jl:112
      [2] macro expansion
        @ ~/.julia/packages/StaticArraysCore/U2Z1K/src/StaticArraysCore.jl:81 [inlined]
      [3] convert_ntuple
        @ ~/.julia/packages/StaticArraysCore/U2Z1K/src/StaticArraysCore.jl:77 [inlined]
      [4] SVector{3, Quantity{Float64, 𝐋 𝐌 𝐍⁻¹ 𝐓⁻², Unitful.FreeUnits{(kJ, nm⁻¹, mol⁻¹), 𝐋 𝐌 𝐍⁻¹ 𝐓⁻², nothing}}}(x::Tuple{Quantity{Float64, 𝐋 𝐌 𝐓⁻², Unitful.FreeUnits{(Å⁻¹, eV), 𝐋 𝐌 𝐓⁻², nothing}}, Quantity{Float64, 𝐋 𝐌 𝐓⁻², Unitful.FreeUnits{(Å⁻¹, eV), 𝐋 𝐌 𝐓⁻², nothing}}, Quantity{Float64, 𝐋 𝐌 𝐓⁻², Unitful.FreeUnits{(Å⁻¹, eV), 𝐋 𝐌 𝐓⁻², nothing}}})
        @ StaticArraysCore ~/.julia/packages/StaticArraysCore/U2Z1K/src/StaticArraysCore.jl:113
      [5] convert
        @ ~/.julia/packages/StaticArrays/jA1zK/src/convert.jl:176 [inlined]
      [6] setindex!(A::Vector{SVector{3, Quantity{Float64, 𝐋 𝐌 𝐍⁻¹ 𝐓⁻², Unitful.FreeUnits{(kJ, nm⁻¹, mol⁻¹), 𝐋 𝐌 𝐍⁻¹ 𝐓⁻², nothing}}}}, x::SVector{3, Quantity{Float64, 𝐋 𝐌 𝐓⁻², Unitful.FreeUnits{(Å⁻¹, eV), 𝐋 𝐌 𝐓⁻², nothing}}}, i1::Int64)
        @ Base ./array.jl:966
      [7] _unsafe_copyto!(dest::Vector{SVector{3, Quantity{Float64, 𝐋 𝐌 𝐍⁻¹ 𝐓⁻², Unitful.FreeUnits{(kJ, nm⁻¹, mol⁻¹), 𝐋 𝐌 𝐍⁻¹ 𝐓⁻², nothing}}}}, doffs::Int64, src::Vector{SVector{3, Quantity{Float64, 𝐋 𝐌 𝐓⁻², Unitful.FreeUnits{(Å⁻¹, eV), 𝐋 𝐌 𝐓⁻², nothing}}}}, soffs::Int64, n::Int64)
        @ Base ./array.jl:253
      [8] unsafe_copyto!
        @ ./array.jl:307 [inlined]
      [9] _copyto_impl!
        @ ./array.jl:331 [inlined]
     [10] copyto!
        @ ./array.jl:317 [inlined]
     [11] copyto!
        @ ./array.jl:343 [inlined]
     [12] copyto_axcheck!
        @ ./abstractarray.jl:1127 [inlined]
     [13] Array
        @ ./array.jl:626 [inlined]
     [14] convert
        @ ./array.jl:617 [inlined]
     [15] push!(a::Vector{Vector{SVector{3, Quantity{Float64, 𝐋 𝐌 𝐍⁻¹ 𝐓⁻², Unitful.FreeUnits{(kJ, nm⁻¹, mol⁻¹), 𝐋 𝐌 𝐍⁻¹ 𝐓⁻², nothing}}}}}, item::Vector{SVector{3, Quantity{Float64, 𝐋 𝐌 𝐓⁻², Unitful.FreeUnits{(Å⁻¹, eV), 𝐋 𝐌 𝐓⁻², nothing}}}})
        @ Base ./array.jl:1057
     [16] log_property!(logger::GeneralObservableLogger{Vector{SVector{3, Quantity{Float64, 𝐋 𝐌 𝐍⁻¹ 𝐓⁻², Unitful.FreeUnits{(kJ, nm⁻¹, mol⁻¹), 𝐋 𝐌 𝐍⁻¹ 𝐓⁻², nothing}}}}, typeof(forces)}, s::System{3, false, Float64, false, Vector{Molly.Atom{Float64, Quantity{Float64, 𝐌, Unitful.FreeUnits{(u,), 𝐌, nothing}}, Quantity{Float64, 𝐋, Unitful.FreeUnits{(nm,), 𝐋, nothing}}, Quantity{Float64, 𝐋² 𝐌 𝐍⁻¹ 𝐓⁻², Unitful.FreeUnits{(kJ, mol⁻¹), 𝐋² 𝐌 𝐍⁻¹ 𝐓⁻², nothing}}}}, Vector{AtomData}, Tuple{}, Tuple{}, Tuple{InteratomicPotentialInter{InteratomicPotentials.LennardJones{Float64}}}, Tuple{}, Vector{SVector{3, Quantity{Float64, 𝐋, Unitful.FreeUnits{(Å,), 𝐋, nothing}}}}, Vector{SVector{3, Quantity{Float64, 𝐋 𝐓⁻¹, Unitful.FreeUnits{(a₀, s⁻¹), 𝐋 𝐓⁻¹, nothing}}}}, CubicBoundary{Quantity{Float64, 𝐋, Unitful.FreeUnits{(Å,), 𝐋, nothing}}}, NoNeighborFinder, NamedTuple{(:force,), Tuple{GeneralObservableLogger{Vector{SVector{3, Quantity{Float64, 𝐋 𝐌 𝐍⁻¹ 𝐓⁻², Unitful.FreeUnits{(kJ, nm⁻¹, mol⁻¹), 𝐋 𝐌 𝐍⁻¹ 𝐓⁻², nothing}}}}, typeof(forces)}}}, Unitful.FreeUnits{(Å⁻¹, eV), 𝐋 𝐌 𝐓⁻², nothing}, Unitful.FreeUnits{(kJ, mol⁻¹), 𝐋² 𝐌 𝐍⁻¹ 𝐓⁻², nothing}, Quantity{Float64, 𝐋² 𝐌 𝚯⁻¹ 𝐓⁻², Unitful.FreeUnits{(kJ, K⁻¹), 𝐋² 𝐌 𝚯⁻¹ 𝐓⁻², nothing}}}, neighbors::Nothing, step_n::Int64; n_threads::Int64, kwargs::Base.Pairs{Symbol, Union{}, Tuple{}, NamedTuple{(), Tuple{}}})
        @ Molly ~/.julia/packages/Molly/SIRC5/src/loggers.jl:72
     [17] #run_loggers!#195
        @ ~/.julia/packages/Molly/SIRC5/src/loggers.jl:31 [inlined]
     [18] simulate!(sys::System{3, false, Float64, false, Vector{Molly.Atom{Float64, Quantity{Float64, 𝐌, Unitful.FreeUnits{(u,), 𝐌, nothing}}, Quantity{Float64, 𝐋, Unitful.FreeUnits{(nm,), 𝐋, nothing}}, Quantity{Float64, 𝐋² 𝐌 𝐍⁻¹ 𝐓⁻², Unitful.FreeUnits{(kJ, mol⁻¹), 𝐋² 𝐌 𝐍⁻¹ 𝐓⁻², nothing}}}}, Vector{AtomData}, Tuple{}, Tuple{}, Tuple{InteratomicPotentialInter{InteratomicPotentials.LennardJones{Float64}}}, Tuple{}, Vector{SVector{3, Quantity{Float64, 𝐋, Unitful.FreeUnits{(Å,), 𝐋, nothing}}}}, Vector{SVector{3, Quantity{Float64, 𝐋 𝐓⁻¹, Unitful.FreeUnits{(a₀, s⁻¹), 𝐋 𝐓⁻¹, nothing}}}}, CubicBoundary{Quantity{Float64, 𝐋, Unitful.FreeUnits{(Å,), 𝐋, nothing}}}, NoNeighborFinder, NamedTuple{(:force,), Tuple{GeneralObservableLogger{Vector{SVector{3, Quantity{Float64, 𝐋 𝐌 𝐍⁻¹ 𝐓⁻², Unitful.FreeUnits{(kJ, nm⁻¹, mol⁻¹), 𝐋 𝐌 𝐍⁻¹ 𝐓⁻², nothing}}}}, typeof(forces)}}}, Unitful.FreeUnits{(Å⁻¹, eV), 𝐋 𝐌 𝐓⁻², nothing}, Unitful.FreeUnits{(kJ, mol⁻¹), 𝐋² 𝐌 𝐍⁻¹ 𝐓⁻², nothing}, Quantity{Float64, 𝐋² 𝐌 𝚯⁻¹ 𝐓⁻², Unitful.FreeUnits{(kJ, K⁻¹), 𝐋² 𝐌 𝚯⁻¹ 𝐓⁻², nothing}}}, sim::VelocityVerlet{Quantity{Float64, 𝐓, Unitful.FreeUnits{(ps,), 𝐓, nothing}}, AndersenThermostat{Quantity{Int64, 𝚯, Unitful.FreeUnits{(K,), 𝚯, nothing}}, Quantity{Float64, 𝐓, Unitful.FreeUnits{(ps,), 𝐓, nothing}}}}, n_steps::Int64; n_threads::Int64)
        @ Molly ~/.julia/packages/Molly/SIRC5/src/simulators.jl:143
     [19] simulate!(sys::System{3, false, Float64, false, Vector{Molly.Atom{Float64, Quantity{Float64, 𝐌, Unitful.FreeUnits{(u,), 𝐌, nothing}}, Quantity{Float64, 𝐋, Unitful.FreeUnits{(nm,), 𝐋, nothing}}, Quantity{Float64, 𝐋² 𝐌 𝐍⁻¹ 𝐓⁻², Unitful.FreeUnits{(kJ, mol⁻¹), 𝐋² 𝐌 𝐍⁻¹ 𝐓⁻², nothing}}}}, Vector{AtomData}, Tuple{}, Tuple{}, Tuple{InteratomicPotentialInter{InteratomicPotentials.LennardJones{Float64}}}, Tuple{}, Vector{SVector{3, Quantity{Float64, 𝐋, Unitful.FreeUnits{(Å,), 𝐋, nothing}}}}, Vector{SVector{3, Quantity{Float64, 𝐋 𝐓⁻¹, Unitful.FreeUnits{(a₀, s⁻¹), 𝐋 𝐓⁻¹, nothing}}}}, CubicBoundary{Quantity{Float64, 𝐋, Unitful.FreeUnits{(Å,), 𝐋, nothing}}}, NoNeighborFinder, NamedTuple{(:force,), Tuple{GeneralObservableLogger{Vector{SVector{3, Quantity{Float64, 𝐋 𝐌 𝐍⁻¹ 𝐓⁻², Unitful.FreeUnits{(kJ, nm⁻¹, mol⁻¹), 𝐋 𝐌 𝐍⁻¹ 𝐓⁻², nothing}}}}, typeof(forces)}}}, Unitful.FreeUnits{(Å⁻¹, eV), 𝐋 𝐌 𝐓⁻², nothing}, Unitful.FreeUnits{(kJ, mol⁻¹), 𝐋² 𝐌 𝐍⁻¹ 𝐓⁻², nothing}, Quantity{Float64, 𝐋² 𝐌 𝚯⁻¹ 𝐓⁻², Unitful.FreeUnits{(kJ, K⁻¹), 𝐋² 𝐌 𝚯⁻¹ 𝐓⁻², nothing}}}, sim::VelocityVerlet{Quantity{Float64, 𝐓, Unitful.FreeUnits{(ps,), 𝐓, nothing}}, AndersenThermostat{Quantity{Int64, 𝚯, Unitful.FreeUnits{(K,), 𝚯, nothing}}, Quantity{Float64, 𝐓, Unitful.FreeUnits{(ps,), 𝐓, nothing}}}}, n_steps::Int64)
        @ Molly ~/.julia/packages/Molly/SIRC5/src/simulators.jl:136
     [20] top-level scope
        @ ~/exploratory/public/interpot_comparisons/lj/no_pluto/fcc_argon_lj_via_molly.jl:730
"""

# letting the force_units be the default
# Hacky conversion. Need to strip it of the wrong units, append the correct units, then convert to eV/A 
f_check =  map(x->uconvert.(u"eV/Å", ustrip.(x)*u"Eh_au/a0_au"), m_sys.loggers.force.history[1])
"""
1536-element Vector{SVector{3, Quantity{Float64, 𝐋 𝐌 𝐓⁻², Unitful.FreeUnits{(Å⁻¹, eV), 𝐋 𝐌 𝐓⁻², nothing}}}}:
 [0.001980826856892926 eV Å⁻¹, 0.020532665135586135 eV Å⁻¹, -0.017177486149641766 eV Å⁻¹]
 [0.0159924490724765 eV Å⁻¹, 0.003581448193741032 eV Å⁻¹, -0.048376398005418646 eV Å⁻¹]
 [-0.02028532914817166 eV Å⁻¹, -0.005666241304493413 eV Å⁻¹, 0.020174725791425048 eV Å⁻¹]
 [0.04269416032056045 eV Å⁻¹, 0.01709052715183096 eV Å⁻¹, -0.07347766990200989 eV Å⁻¹]
 [-0.057566936261611516 eV Å⁻¹, 0.07178852927421989 eV Å⁻¹, -0.0682536204422086 eV Å⁻¹]
 [0.0277148269047553 eV Å⁻¹, -0.03688926738500493 eV Å⁻¹, 0.06014495409852314 eV Å⁻¹]
 [0.022093629592778875 eV Å⁻¹, -0.047387516149956586 eV Å⁻¹, -0.07037312491157584 eV Å⁻¹]
 [-0.026274454948404934 eV Å⁻¹, -0.012639241466099996 eV Å⁻¹, 0.036808828946277995 eV Å⁻¹]
 [-0.01913305345153546 eV Å⁻¹, 0.08943042296010875 eV Å⁻¹, -0.008002507976242357 eV Å⁻¹]
 [-0.06344539201413879 eV Å⁻¹, 0.0025539835735561337 eV Å⁻¹, 0.003882271799257547 eV Å⁻¹]
...
By eye, these match the forces in ../ref_config/raw_forces
"""