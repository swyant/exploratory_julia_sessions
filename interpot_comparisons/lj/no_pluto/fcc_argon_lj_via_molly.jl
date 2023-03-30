using Molly
using InteratomicPotentials
using AtomsIO
using Unitful
using AtomsBase
 
# Wrapper struct enabling potentials implemented in InteratomicPotentials.jl to be used in Molly
struct InteratomicPotentialInter{P<:AbstractPotential}
	potential::P
    eaf::Function
    # TODO: should also use the energy cache Jeremiah had, but ignoring that for now
    # TODO: Default constructur with energy_and_force probably
end

# Extending Molly.forces with above struct
function Molly.forces(inter::InteratomicPotentialInter, sys, neighbors=nothing)
    eandf = inter.eaf(sys, inter.potential)
    eandf.f
end

function molly_params(sys::AtomsBase.FlexibleSystem)
    coords = [atm.position for atm in sys.particles]
    atoms  = [Atom(mass=atm.atomic_mass,atomic_symbol=atm.atomic_symbol) for atm in sys.particles] # clean Atoms, mass instead of atomic_mass
    velocities = [atm.velocity for atm in sys.particles]

    # TODO: Generalize this. Currently very fragile, assumes you read in a cubic box. Otherwise need to use TriclinicBoundary, but some limitations...
    # Also by default, assumes directions are periodic..., ignoring sys.boundary_conditions
    boundary=CubicBoundary(sys.box[1][1],sys.box[2][2],sy.box[3][3]) # FRAGILE! 

    return NamedTuple((atoms=atoms, coords=coords, velocities = velocities, boundary=boundary))
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
f_p = [uconvert.(u"eV/Å", fi) for fi in f_p]

mp = molly_params(sys)

m_sys = System(mp...,
            general_inters=general_inters, 
            loggers=(force=ForceLogger(1),)
            )


#### run zero simulation
simulator = VelocityVerlet(
    dt=0.001u"ps",
    coupling=AndersenThermostat(80u"K", 1.0u"ps"),
)

simulate!(sys,simulator,0)