using ACE1 
using CSV, DataFrames
using JuLIP
using Unitful, UnitfulAtomic
using AtomsBase
using Plots

rpib = ACE1.rpi_basis(;
            species = [:Ti, :Al],
            N       = 2, 
            maxdeg  = 6,
            D       = ACE1.SparsePSHDegree(; wL = 1.5, csp = 1.0),
            r0      = 2.9,
            rin     = 0.65*2.9,
            rcut    = 5.5,
            pin     = 0,
       )

df = CSV.read("TiAl_example_N2_rinDefault_pin0_coeffs.csv", DataFrame, header=false)
coeffs = vec(Matrix(df))

pot_check = JuLIP.MLIPs.combine(rpib,coeffs)

box = [[50.0,0.0,0.0], [0.0, 50.0, 0.0], [0.0, 0.0, 50.0]]u"Å"
bcs = [AtomsBase.Periodic(),AtomsBase.Periodic(),AtomsBase.Periodic()]
test_sys_tmp = FlexibleSystem([ AtomsBase.Atom(:Ti, [10.,25.,25.]u"Å"),
                                AtomsBase.Atom(:Al, [13,25.0,25.]u"Å")],
                                box,bcs)

test_sys = JuLIP.Atoms(test_sys_tmp)

energy(pot_check,test_sys)


### Actual scan 
rijs = range(0.005,15.005,length=600)

pair_systems = Vector{JuLIP.Atoms}[]
energies = Vector{Float64}()
atom1_pos = [10.,25.,25.]
for rij in rijs
    atom2_pos = atom1_pos .+ [rij, 0.0, 0.0]
    sys_tmp = FlexibleSystem([ AtomsBase.Atom(:Ti, atom1_pos*u"Å"),
                                AtomsBase.Atom(:Al, atom2_pos*u"Å")],
                                box,bcs)
    sys = JuLIP.Atoms(sys_tmp)
    e = energy(pot_check,sys)
    push!(energies,e) 
end

plotly() #GR backend doesn't play well with the smaller limits?
plot(rijs,energies)
ylims!(-2.5,10)
xlims!(1.0,5.6)


### Nonsense potential, but checking what happens when pin=2
rpib2 = ACE1.rpi_basis(;
            species = [:Ti, :Al],
            N       = 2, 
            maxdeg  = 6,
            D       = ACE1.SparsePSHDegree(; wL = 1.5, csp = 1.0),
            r0      = 2.9,
            rin     = 0.65*2.9,
            rcut    = 5.5,
            pin     = 2,
       )
