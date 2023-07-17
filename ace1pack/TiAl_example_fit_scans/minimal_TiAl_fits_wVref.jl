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

weights = Dict("default" => Dict("E" =>1.0,"F" => 1.0, "V"=>0.0))

vref = JuLIP.OneBody([:Ti => -1586.0195, :Al => -105.5954]...)
data = [ AtomsData(at;  energy_key = "energy", force_key="force",
                        weights = weights, v_ref=vref) for at in raw_data[1:5:end]]

A, Y, W = ACEfit.assemble(data, rpib)

# Emmanuel approach
W_man = ones(Float64,length(W))
man_eman_coeffs = (A'* Diagonal(W_man) *A) \ (A'* Diagonal(W_man) * Y )
@show man_eman_coeffs

pot_emmanuel_temp = JuLIP.MLIPs.combine(rpib,man_eman_coeffs)
pot_emmanuel = JuLIP.MLIPs.SumIP(pot_emmanuel_temp, vref)
ACE1pack.linear_errors(data,pot_emmanuel)
#=
[ Info: RMSE Table
┌────────────┬─────────┬──────────┬─────────┐
│       Type │ E [meV] │ F [eV/A] │ V [meV] │
├────────────┼─────────┼──────────┼─────────┤
│   FLD_TiAl │  69.229 │    0.067 │   0.000 │
│ TiAl_T5000 │   1.725 │    0.273 │   0.000 │
├────────────┼─────────┼──────────┼─────────┤
│        set │  67.638 │    0.233 │   0.000 │
└────────────┴─────────┴──────────┴─────────┘
[ Info: MAE Table
┌────────────┬─────────┬──────────┬─────────┐
│       Type │ E [meV] │ F [eV/A] │ V [meV] │
├────────────┼─────────┼──────────┼─────────┤
│   FLD_TiAl │  49.294 │    0.042 │   0.000 │
│ TiAl_T5000 │   1.360 │    0.205 │   0.000 │
├────────────┼─────────┼──────────┼─────────┤
│        set │  47.116 │    0.158 │   0.000 │
└────────────┴─────────┴──────────┴─────────┘
=#

### Default QR
solver_qr = ACEfit.QR()

results_qr = ACEfit.solve(solver_qr,A,Y)

@show results_qr["C"]

pot_qr_tmp = JuLIP.MLIPs.combine(rpib,results_qr["C"])
pot_qr = JuLIP.MLIPs.SumIP(pot_qr_tmp,vref)
ACE1pack.linear_errors(data,pot_qr)
#=
[ Info: RMSE Table
┌────────────┬─────────┬──────────┬─────────┐
│       Type │ E [meV] │ F [eV/A] │ V [meV] │
├────────────┼─────────┼──────────┼─────────┤
│   FLD_TiAl │  69.801 │    0.066 │   0.000 │
│ TiAl_T5000 │   2.632 │    0.273 │   0.000 │
├────────────┼─────────┼──────────┼─────────┤
│        set │  68.199 │    0.233 │   0.000 │
└────────────┴─────────┴──────────┴─────────┘
[ Info: MAE Table
┌────────────┬─────────┬──────────┬─────────┐
│       Type │ E [meV] │ F [eV/A] │ V [meV] │
├────────────┼─────────┼──────────┼─────────┤
│   FLD_TiAl │  50.002 │    0.042 │   0.000 │
│ TiAl_T5000 │   2.462 │    0.205 │   0.000 │
├────────────┼─────────┼──────────┼─────────┤
│        set │  47.841 │    0.158 │   0.000 │
└────────────┴─────────┴──────────┴─────────┘
Dict{String, Dict{String, Any}} with 2 entries:
  "rmse" => Dict("set"=>Dict("V"=>0.0, "E"=>0.0681987, "F"=>0.232668), "FLD_TiAl"=>Dict("V"=>0.0, "E"=>0.0698013, "F"=>0.0660986), "TiAl_T500…
  "mae"  => Dict("set"=>Dict("V"=>0.0, "E"=>0.047841, "F"=>0.157662), "FLD_TiAl"=>Dict("V"=>0.0, "E"=>0.0500019, "F"=>0.0417805), "TiAl_T5000…
=# 

### RRQR
solver_rrqr = ACEfit.RRQR()

results_rrqr = ACEfit.solve(solver_rrqr,A,Y)

@show results_rrqr["C"]

pot_rrqr_tmp = JuLIP.MLIPs.combine(rpib,results_rrqr["C"])
pot_rrqr = JuLIP.MLIPs.SumIP(pot_rrqr_tmp,vref)
ACE1pack.linear_errors(data,pot_rrqr)

#=
[ Info: RMSE Table
┌────────────┬─────────┬──────────┬─────────┐
│       Type │ E [meV] │ F [eV/A] │ V [meV] │
├────────────┼─────────┼──────────┼─────────┤
│   FLD_TiAl │  69.229 │    0.067 │   0.000 │
│ TiAl_T5000 │   1.725 │    0.273 │   0.000 │
├────────────┼─────────┼──────────┼─────────┤
│        set │  67.638 │    0.233 │   0.000 │
└────────────┴─────────┴──────────┴─────────┘
[ Info: MAE Table
┌────────────┬─────────┬──────────┬─────────┐
│       Type │ E [meV] │ F [eV/A] │ V [meV] │
├────────────┼─────────┼──────────┼─────────┤
│   FLD_TiAl │  49.294 │    0.042 │   0.000 │
│ TiAl_T5000 │   1.360 │    0.205 │   0.000 │
├────────────┼─────────┼──────────┼─────────┤
│        set │  47.116 │    0.158 │   0.000 │
└────────────┴─────────┴──────────┴─────────┘
=#