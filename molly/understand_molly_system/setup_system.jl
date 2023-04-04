using Molly 
using AtomsBase
using Unitful 


n_atoms = 200
mass   = 15.0u"u"
atoms  = [Molly.Atom(mass=mass, σ=0.3u"nm", ϵ=0.2u"kJ * mol^-1") for i in 1:n_atoms]

boundary = CubicBoundary(4.0u"nm", 4.0u"nm", 4.0u"nm")
coords   = place_atoms(n_atoms, boundary; min_dist=0.3u"nm")

temp       = 50.0u"K"
velocities = [velocity(mass,temp) for i in 1:n_atoms] # sample from M-B dist

pairwise_inters = (LennardJones(),)

sys = System(
	atoms     		= atoms,
	pairwise_inters = pairwise_inters,
	coords 		   	= coords,
	velocities     	= velocities,
	boundary 		= boundary,
	loggers 		= (temp=TemperatureLogger(10),
    					   coords=CoordinateLogger(5)),
)

simulator = VelocityVerlet(
			dt=0.002u"ps",
			coupling=AndersenThermostat(temp, 1.0u"ps")
)
simulate!(sys,simulator,5000)
