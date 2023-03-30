using Molly
using InteratomicPotentials
using AtomsIO
using Unitful
 

struct InteratomicPotentialInter{P<:AbstractPotential}
	potential::P
    # should also use the energy cache Jeremiah had, but ignoring that for now
end

function Molly.forces(inter::InteratomicPotentialInter, sys, neighbors=nothing)
    eandf = energy_and_force(sys, inter.potential)
    eandf.f
end


sys = load_system("../ref_config/dump_final.extxyz")

ϵ = 0.01032u"eV"
σ = 3.405u"Å"
rcut = 8.51u"Å"
species = [:Ar]
lj_p = InteratomicPotentials.LennardJones(ϵ,σ,rcut,species)


inter_lj = InteratomicPotentialInter(lj_p)

general_inters = (inter_lj,)

f_p = Molly.forces(inter_lj, sys)

f_p = [uconvert.(u"eV/Å", fi) for fi in f_p]


#### run zero simulation
simulator = VelocityVerlet(
    dt=0.001u"ps",
    coupling=AndersenThermostat(80u"K", 1.0u"ps"),
)

simulate!(sys,simulator,0)