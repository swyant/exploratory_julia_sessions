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
function Molly.forces(inter::InteratomicPotentialInter, msys::Molly.System, neighbors=nothing)
    sys = convertMollySys(msys)
    eandf = inter.eaf(sys, inter.potential)
    eandf.f
end

function molly_params(sys::AtomsBase.FlexibleSystem)
    coords = [atm.position for atm in sys.particles]

    # Molly has its own Atom, but want to use AtomsBase version to carry around element info (for converting back to FlexibleSystem)
    # this also means that it has to carry around a position (AtomsBase.Atom constructor requires it), which is NOT UPDATED during simulate!. basically ignore its
    atoms  = [AtomsBase.Atom(atm.atomic_symbols, atm.position, mass=atm.atomic_mass) for atm in sys.particles] # not using Molly.Atom
    velocities = [atm.velocity for atm in sys.particles]

    # TODO: Generalize this. Currently very fragile, assumes you read in a cubic box. Otherwise need to use TriclinicBoundary, but some limitations...
    # Also by default, assumes directions are periodic..., ignoring sys.boundary_conditions
    boundary=CubicBoundary(sys.box[1][1],sys.box[2][2],sys.box[3][3]) # FRAGILE! 

    return NamedTuple((atoms=atoms, coords=coords, velocities = velocities, boundary=boundary))
end 


function convertMollySys(msys::Molly.System)
    atoms = []
    for i in 1:length(msys.atoms)
        coord = msys.coords[i]
        atm_species = msys.atoms[i].atomic_symbol 
        push!(atoms,AtomsBase.Atom(atm_species,coord))
    end

    # FRAGILE! Assuming CubicBoundary. And angstrom units...
    side_lengths = msys.boundary.side_lengths
    box = [[side_lengths[1],0.0,0.0],[0.0,side_lengths[2],0.0],[0.0,0.0,side_lengths[3]]]u"Å"

    # FRAGILE! assuming periodic boundary conditions 
    bcs = [Periodic(),Periodic(),Periodic()]
    print("hey")

    FlexibleSystem(atoms,box,bcs)
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
#f_p = Molly.forces(inter_lj, sys)
#f_p = [uconvert.(u"eV/Å", fi) for fi in f_p]

mp = molly_params(sys)


m_sys = System(;mp...,
            general_inters=general_inters, 
            loggers=(force=ForceLogger(1),)
            )


#### run zero simulation
simulator = VelocityVerlet(
    dt=0.001u"ps",
    coupling=AndersenThermostat(80u"K", 1.0u"ps"),
)

simulate!(m_sys,simulator,1)