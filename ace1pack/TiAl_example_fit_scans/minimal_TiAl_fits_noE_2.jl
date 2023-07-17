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

Aw = Diagonal(W) * A
Yw = W .* Y


### manual way to match PL.jl learn!(::UnivariateLinear)
AtA = transpose(Aw)*Aw
Atb = transpose(Aw)*Yw
Q = pinv(AtA, 1e-6)
manual_coeffs = Q*Atb
@show manual_coeffs

pot_manual = JuLIP.MLIPs.combine(rpib,manual_coeffs)
ACE1pack.linear_errors(data,pot_manual)
#=
  [ Info: RMSE Table
┌────────────┬────────────┬──────────┬─────────┐
│       Type │    E [meV] │ F [eV/A] │ V [meV] │
├────────────┼────────────┼──────────┼─────────┤
│   FLD_TiAl │ 850622.427 │    0.184 │   0.000 │
│ TiAl_T5000 │ 850471.352 │    0.399 │   0.000 │
├────────────┼────────────┼──────────┼─────────┤
│        set │ 850615.561 │    0.351 │   0.000 │
└────────────┴────────────┴──────────┴─────────┘
[ Info: MAE Table
┌────────────┬────────────┬──────────┬─────────┐
│       Type │    E [meV] │ F [eV/A] │ V [meV] │
├────────────┼────────────┼──────────┼─────────┤
│   FLD_TiAl │ 850622.402 │    0.124 │   0.000 │
│ TiAl_T5000 │ 850471.351 │    0.311 │   0.000 │
├────────────┼────────────┼──────────┼─────────┤
│        set │ 850615.536 │    0.257 │   0.000 │
└────────────┴────────────┴──────────┴─────────┘
=#


##### QR, no-prior
solver = ACEfit.QR()

resultsw = ACEfit.solve(solver, Aw, Yw)
@show resultsw["C"]

potw = JuLIP.MLIPs.combine(rpib,resultsw["C"])
ACE1pack.linear_errors(data,potw)
#=

  [ Info: RMSE Table
┌────────────┬────────────┬──────────┬─────────┐
│       Type │    E [meV] │ F [eV/A] │ V [meV] │
├────────────┼────────────┼──────────┼─────────┤
│   FLD_TiAl │ 854495.967 │    0.052 │   0.000 │
│ TiAl_T5000 │ 854364.199 │    0.266 │   0.000 │
├────────────┼────────────┼──────────┼─────────┤
│        set │ 854489.978 │    0.226 │   0.000 │
└────────────┴────────────┴──────────┴─────────┘
[ Info: MAE Table
┌────────────┬────────────┬──────────┬─────────┐
│       Type │    E [meV] │ F [eV/A] │ V [meV] │
├────────────┼────────────┼──────────┼─────────┤
│   FLD_TiAl │ 854495.933 │    0.035 │   0.000 │
│ TiAl_T5000 │ 854364.187 │    0.198 │   0.000 │
├────────────┼────────────┼──────────┼─────────┤
│        set │ 854489.944 │    0.151 │   0.000 │
└────────────┴────────────┴──────────┴─────────┘
=#

##### QR, w/ default prior (manually created, not checked!!!)
dp = Float64[]
append!(dp, ACE1.scaling(rpib, 2, 1.0))
prior = Diagonal(1 .+ dp)

Ap = Diagonal(W) * (A / prior)
Yp = W .* Y 
resultsp = ACEfit.solve(solver, Ap, Yp)

resultsp["C"]
@show coeffs_p = prior \ resultsp["C"]

pot_p = JuLIP.MLIPs.combine(rpib,coeffs_p)
ACE1pack.linear_errors(data, pot_p)
#=
  [ Info: RMSE Table
┌────────────┬────────────┬──────────┬─────────┐
│       Type │    E [meV] │ F [eV/A] │ V [meV] │
├────────────┼────────────┼──────────┼─────────┤
│   FLD_TiAl │ 854557.476 │    0.052 │   0.000 │
│ TiAl_T5000 │ 854446.331 │    0.266 │   0.000 │
├────────────┼────────────┼──────────┼─────────┤
│        set │ 854552.424 │    0.226 │   0.000 │
└────────────┴────────────┴──────────┴─────────┘
[ Info: MAE Table
┌────────────┬────────────┬──────────┬─────────┐
│       Type │    E [meV] │ F [eV/A] │ V [meV] │
├────────────┼────────────┼──────────┼─────────┤
│   FLD_TiAl │ 854557.441 │    0.034 │   0.000 │
│ TiAl_T5000 │ 854446.318 │    0.198 │   0.000 │
├────────────┼────────────┼──────────┼─────────┤
│        SET │ 854552.390 │    0.151 │   0.000 │
└────────────┴────────────┴──────────┴─────────┘
=#


#### Ridge-Regression (I think!!!) w/ QR solver
solver_reg = ACEfit.QR(;lambda=1e-2) # I varied lambda from 1e-2 to 1e-6, not much of a difference here

results_reg = ACEfit.solve(solver_reg, Aw, Yw)
@show results_reg["C"]

pot_reg = JuLIP.MLIPs.combine(rpib,results_reg["C"])
ACE1pack.linear_errors(data,pot_reg)
#=
 [ Info: RMSE Table
┌────────────┬────────────┬──────────┬─────────┐
│       Type │    E [meV] │ F [eV/A] │ V [meV] │
├────────────┼────────────┼──────────┼─────────┤
│   FLD_TiAl │ 854513.488 │    0.050 │   0.000 │
│ TiAl_T5000 │ 854385.359 │    0.266 │   0.000 │
├────────────┼────────────┼──────────┼─────────┤
│        set │ 854507.665 │    0.226 │   0.000 │
└────────────┴────────────┴──────────┴─────────┘
[ Info: MAE Table
┌────────────┬────────────┬──────────┬─────────┐
│       Type │    E [meV] │ F [eV/A] │ V [meV] │
├────────────┼────────────┼──────────┼─────────┤
│   FLD_TiAl │ 854513.455 │    0.034 │   0.000 │
│ TiAl_T5000 │ 854385.347 │    0.198 │   0.000 │
├────────────┼────────────┼──────────┼─────────┤
│        set │ 854507.632 │    0.150 │   0.000 │
=#


### Default RRQR 
solver_rrqr = ACEfit.RRQR()

results_rrqr = ACEfit.solve(solver_rrqr, Aw, Yw)
@show results_rrqr["C"]

pot_rrqr = JuLIP.MLIPs.combine(rpib,results_rrqr["C"])
ACE1pack.linear_errors(data,pot_rrqr)

#=
  [ Info: RMSE Table
┌────────────┬────────────┬──────────┬─────────┐
│       Type │    E [meV] │ F [eV/A] │ V [meV] │
├────────────┼────────────┼──────────┼─────────┤
│   FLD_TiAl │ 854525.077 │    0.051 │   0.000 │
│ TiAl_T5000 │ 854396.990 │    0.266 │   0.000 │
├────────────┼────────────┼──────────┼─────────┤
│        set │ 854519.255 │    0.226 │   0.000 │
└────────────┴────────────┴──────────┴─────────┘
[ Info: MAE Table
┌────────────┬────────────┬──────────┬─────────┐
│       Type │    E [meV] │ F [eV/A] │ V [meV] │
├────────────┼────────────┼──────────┼─────────┤
│   FLD_TiAl │ 854525.043 │    0.034 │   0.000 │
│ TiAl_T5000 │ 854396.978 │    0.198 │   0.000 │
├────────────┼────────────┼──────────┼─────────┤
│        set │ 854519.222 │    0.150 │   0.000 │
└────────────┴────────────┴──────────┴─────────┘
=#

### Default LSQR 
solver_lsqr = ACEfit.LSQR()

results_lsqr = ACEfit.solve(solver_lsqr, Aw, Yw)
@show results_lsqr["C"]

pot_lsqr = JuLIP.MLIPs.combine(rpib,results_lsqr["C"])
ACE1pack.linear_errors(data,pot_lsqr)
