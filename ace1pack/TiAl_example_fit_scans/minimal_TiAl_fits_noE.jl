using ACE1 
using ACE1pack
using ACEfit
using LinearAlgebra: I, Diagonal, pinv
using JuLIP

ds_path = "/Users/swyant/.julia/artifacts/b437d7d5fac4424b8203d0afc31732879d3da5b2/TiAl_tutorial.xyz" 
raw_data = read_extxyz(ds_path)

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

weights = Dict("default" => Dict("E" =>0.0,"F" => 1.0, "V"=>0.0))

data = [ AtomsData(at;  energy_key = "energy", force_key="force",
                        weights = weights) for at in raw_data[1:5:end]]

A, Y, W = ACEfit.assemble(data, rpib)
solver = ACEfit.QR()


Aw = Diagonal(W) * A
Yw = W .* Y
resultsw = ACEfit.solve(solver, Aw, Yw)

@show resultsw["C"]

potw = JuLIP.MLIPs.combine(rpib,resultsw["C"])

ACE1pack.linear_errors(data,potw)

### Manually making basis, haven't checked this 
dp = Float64[]
append!(dp, ACE1.scaling(rpib, 2, 1.0))
prior = Diagonal(1 .+ dp)

improper_coeffs = prior \ resultsw["C"]
pot_improper = JuLIP.MLIPs.combine(rpib,improper_coeffs)
ACE1pack.linear_errors(data,pot_improper)

Ap = Diagonal(W) * (A / prior)
Yp = W .* Y 
resultsp = ACEfit.solve(solver, Ap, Yp)

@show resultsp["C"]

coeffs_p = prior \ resultsp["C"]

pot_p = JuLIP.MLIPs.combine(rpib,coeffs_p)
ACE1pack.linear_errors(data, pot_p)

#### can I manually recreate the PotentialLearning.jl fit?

AtA = transpose(Aw)*Aw
Atb = transpose(Aw)*Yw
Q = pinv(AtA, 1e-6)
manual_coeffs = Q*Atb

pot_manual = JuLIP.MLIPs.combine(rpib,manual_coeffs)
ACE1pack.linear_errors(data,pot_manual)

### Need to check these errors 
#check_data = [ AtomsData(at;  energy_key = "energy", force_key="force",
#                weights = weights) for at in [raw_data[end]]]
#
#forces(pot_manual, raw_data[end])
#
#ACE1pack.linear_errors(check_data, pot_manual)

check_data = [ AtomsData(at;  energy_key = "energy", force_key="force",
                weights = weights) for at in raw_data[1:2]]

ACE1pack.linear_errors(check_data, pot_manual)

