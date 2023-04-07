using Molly
using InteratomicPotentials
using AtomsIO # .CondaPkg/ comes from this, but shouldn't with newest version
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

inter_lj = InteratomicPotentialInter(lj_p,InteratomicPotentials.energy_and_force)
general_inters = (inter_lj,)

# Check force with regular sys obtained with AtomsIO
f_p = Molly.forces(inter_lj, sys)
#f_p = [uconvert.(u"eV/Å", fi) for fi in f_p]

mp = molly_params(sys)

m_sys = System(;mp...,
            general_inters=general_inters, 
            loggers=(force=ForceLogger(typeof(1.0u"eV/Å"), 1),),
            force_units=u"eV/Å",
            #loggers=(force=ForceLogger(Float32, 1),),
            #force_units=NoUnits,
            )

#### run zero simulation
simulator = VelocityVerlet(
    dt=0.001u"ps",
    coupling=AndersenThermostat(80u"K", 1.0u"ps"),
)

simulate!(m_sys,simulator,0)

f_p_molly = m_sys.loggers.force.history[1] # outputs units as eV/Å, but values weren't actually converted

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