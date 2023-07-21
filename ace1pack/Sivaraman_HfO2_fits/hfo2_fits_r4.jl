using ACE1 
using ACE1pack
using ACEfit
using LinearAlgebra: I, Diagonal, pinv
using JuLIP
using CSV, DataFrames
using Random

### functions
function manual_fit(A,Y,rpib,vref)
    man_coeffs = (A'*A) \ (A'* Y )
    @show man_coeffs 

    pot_man_temp = JuLIP.MLIPs.combine(rpib,man_coeffs)
    pot_man = JuLIP.MLIPs.SumIP(pot_man_temp, vref)
   
    pot_man
end

function reg_qr_fit(A,Y,rpib,vref; lambda=1e-3)
    solver_qr = ACEfit.QR(;lambda=lambda)

    results_qr = ACEfit.solve(solver_qr,A,Y)
   
    @show results_qr["C"]

    pot_qr_tmp = JuLIP.MLIPs.combine(rpib,results_qr["C"])
    pot_qr = JuLIP.MLIPs.SumIP(pot_qr_tmp,vref)
 
    pot_qr
end

#####
ds_path = "../../../../datasets/HfO2_Sivaraman/prl_2021/raw/train.xyz" 
raw_data = read_extxyz(ds_path)

train_idxs = vec(Matrix(CSV.read("Siviraman_HfO2_my_train_idxs.csv",DataFrame,header=false)))
val_idxs   = vec(Matrix(CSV.read("Siviraman_HfO2_my_val_idxs.csv",DataFrame,header=false)))

rpib1 = ACE1.rpi_basis(;
            species = [:Hf, :O],
            N       = 3, 
            maxdeg  = 8,
            D       = ACE1.SparsePSHDegree(; wL = 1.5, csp = 1.0),
            r0      = 2.15,
            rin     = 0.65*2.15,  # Does this meaningfully get used when pin=0?
            rcut    = 4.0, # the only thing that has changed
            pin     = 0,
       )

@show length(rpib1) #498

weights = Dict("default" => Dict("E" =>1.0,"F" => 1.0, "V"=>0.0))

vref = JuLIP.OneBody([:Hf => -2.70516846, :O => -0.01277342]...)
train_data = [ AtomsData(at;  energy_key = "energy", force_key="forces",
                        weights = weights, v_ref=vref) for at in raw_data[train_idxs]]

val_data = [ AtomsData(at;  energy_key = "energy", force_key="forces",
                        weights = weights, v_ref=vref) for at in raw_data[val_idxs]]

A1, Y1, W1 = ACEfit.assemble(train_data, rpib1)

#### Manual approach
pot_man1 = manual_fit(A1,Y1,rpib1,vref)
ACE1pack.linear_errors(train_data,pot_man1)
#=
[ Info: RMSE Table
┌─────────┬─────────┬──────────┬─────────┐
│    Type │ E [meV] │ F [eV/A] │ V [meV] │
├─────────┼─────────┼──────────┼─────────┤
│  amorph │   8.358 │    0.357 │   0.000 │
│ crystal │   5.643 │    0.265 │   0.000 │
├─────────┼─────────┼──────────┼─────────┤
│     set │   6.584 │    0.295 │   0.000 │
└─────────┴─────────┴──────────┴─────────┘
[ Info: MAE Table
┌─────────┬─────────┬──────────┬─────────┐
│    Type │ E [meV] │ F [eV/A] │ V [meV] │
├─────────┼─────────┼──────────┼─────────┤
│  amorph │   6.467 │    0.268 │   0.000 │
│ crystal │   4.445 │    0.191 │   0.000 │
├─────────┼─────────┼──────────┼─────────┤
│     set │   5.057 │    0.214 │   0.000 │
└─────────┴─────────┴──────────┴─────────┘
=#
ACE1pack.linear_errors(val_data,pot_man1)
#=
[ Info: RMSE Table
┌─────────┬─────────┬──────────┬─────────┐
│    Type │ E [meV] │ F [eV/A] │ V [meV] │
├─────────┼─────────┼──────────┼─────────┤
│  amorph │   8.581 │    0.367 │   0.000 │
│ crystal │   6.052 │    0.263 │   0.000 │
├─────────┼─────────┼──────────┼─────────┤
│     set │   6.976 │    0.300 │   0.000 │
└─────────┴─────────┴──────────┴─────────┘
[ Info: MAE Table
┌─────────┬─────────┬──────────┬─────────┐
│    Type │ E [meV] │ F [eV/A] │ V [meV] │
├─────────┼─────────┼──────────┼─────────┤
│  amorph │   6.620 │    0.277 │   0.000 │
│ crystal │   4.900 │    0.186 │   0.000 │
├─────────┼─────────┼──────────┼─────────┤
│     set │   5.460 │    0.215 │   0.000 │
└─────────┴─────────┴──────────┴─────────┘
=#


###### Reg. QR

pot_qr_reg1 = reg_qr_fit(A1,Y1,rpib1,vref)

ACE1pack.linear_errors(train_data,pot_qr_reg1)
#=
[ Info: RMSE Table
┌─────────┬─────────┬──────────┬─────────┐
│    Type │ E [meV] │ F [eV/A] │ V [meV] │
├─────────┼─────────┼──────────┼─────────┤
│  amorph │   8.422 │    0.358 │   0.000 │
│ crystal │   5.582 │    0.267 │   0.000 │
├─────────┼─────────┼──────────┼─────────┤
│     set │   6.573 │    0.297 │   0.000 │
└─────────┴─────────┴──────────┴─────────┘
[ Info: MAE Table
┌─────────┬─────────┬──────────┬─────────┐
│    Type │ E [meV] │ F [eV/A] │ V [meV] │
├─────────┼─────────┼──────────┼─────────┤
│  amorph │   6.538 │    0.270 │   0.000 │
│ crystal │   4.399 │    0.192 │   0.000 │
├─────────┼─────────┼──────────┼─────────┤
│     set │   5.047 │    0.215 │   0.000 │
└─────────┴─────────┴──────────┴─────────┘
=#


rpib3 = ACE1.rpi_basis(;
            species = [:Hf, :O],
            N       = 3, 
            maxdeg  = 10,
            D       = ACE1.SparsePSHDegree(; wL = 1.5, csp = 1.0),
            r0      = 2.15,
            rin     = 0.65*2.15,  
            rcut    = 4.0,
            pin     = 0,
)

A3, Y3, W3 = ACEfit.assemble(train_data, rpib3)


#### Manual approach
pot_man3 = manual_fit(A3,Y3,rpib3,vref)
ACE1pack.linear_errors(train_data,pot_man3)
#=
[ Info: RMSE Table
┌─────────┬─────────┬──────────┬─────────┐
│    Type │ E [meV] │ F [eV/A] │ V [meV] │
├─────────┼─────────┼──────────┼─────────┤
│  amorph │   7.292 │    0.337 │   0.000 │
│ crystal │   5.372 │    0.243 │   0.000 │
├─────────┼─────────┼──────────┼─────────┤
│     set │   6.018 │    0.274 │   0.000 │
└─────────┴─────────┴──────────┴─────────┘
[ Info: MAE Table
┌─────────┬─────────┬──────────┬─────────┐
│    Type │ E [meV] │ F [eV/A] │ V [meV] │
├─────────┼─────────┼──────────┼─────────┤
│  amorph │   5.739 │    0.254 │   0.000 │
│ crystal │   4.226 │    0.176 │   0.000 │
├─────────┼─────────┼──────────┼─────────┤
│     set │   4.684 │    0.199 │   0.000 │
└─────────┴─────────┴──────────┴─────────┘
=#

ACE1pack.linear_errors(val_data,pot_man3)
#=
[ Info: RMSE Table
┌─────────┬─────────┬──────────┬─────────┐
│    Type │ E [meV] │ F [eV/A] │ V [meV] │
├─────────┼─────────┼──────────┼─────────┤
│  amorph │   7.527 │    0.346 │   0.000 │
│ crystal │   5.727 │    0.240 │   0.000 │
├─────────┼─────────┼──────────┼─────────┤
│     set │   6.369 │    0.278 │   0.000 │
└─────────┴─────────┴──────────┴─────────┘
[ Info: MAE Table
┌─────────┬─────────┬──────────┬─────────┐
│    Type │ E [meV] │ F [eV/A] │ V [meV] │
├─────────┼─────────┼──────────┼─────────┤
│  amorph │   5.755 │    0.261 │   0.000 │
│ crystal │   4.550 │    0.172 │   0.000 │
├─────────┼─────────┼──────────┼─────────┤
│     set │   4.942 │    0.200 │   0.000 │
└─────────┴─────────┴──────────┴─────────┘
=#


rpib4 = ACE1.rpi_basis(;
            species = [:Hf, :O],
            N       = 4, 
            maxdeg  = 9,
            D       = ACE1.SparsePSHDegree(; wL = 1.5, csp = 1.0),
            r0      = 2.15,
            rin     = 0.65*2.15,  
            rcut    = 4.0,
            pin     = 0,
)

A4, Y4, W4 = ACEfit.assemble(train_data, rpib4)

pot_man4 = manual_fit(A4,Y4,rpib4,vref)
ACE1pack.linear_errors(train_data,pot_man4)
#=
[ Info: RMSE Table
┌─────────┬─────────┬──────────┬─────────┐
│    Type │ E [meV] │ F [eV/A] │ V [meV] │
├─────────┼─────────┼──────────┼─────────┤
│  amorph │   7.393 │    0.345 │   0.000 │
│ crystal │   5.451 │    0.245 │   0.000 │
├─────────┼─────────┼──────────┼─────────┤
│     set │   6.104 │    0.278 │   0.000 │
└─────────┴─────────┴──────────┴─────────┘
[ Info: MAE Table
┌─────────┬─────────┬──────────┬─────────┐
│    Type │ E [meV] │ F [eV/A] │ V [meV] │
├─────────┼─────────┼──────────┼─────────┤
│  amorph │   5.660 │    0.259 │   0.000 │
│ crystal │   4.356 │    0.178 │   0.000 │
├─────────┼─────────┼──────────┼─────────┤
│     set │   4.751 │    0.201 │   0.000 │
=#


#### sanity check against just forces
rpib_orig = ACE1.rpi_basis(;
            species = [:Hf, :O],
            N       = 3, 
            maxdeg  = 10,
            D       = ACE1.SparsePSHDegree(; wL = 1.5, csp = 1.0),
            r0      = 2.15,
            rin     = 0.65*2.15,  # Does this meaningfully get used when pin=0?
            rcut    = 5.0, # the only thing that has changed
            pin     = 0,
       )


weights_noE = weights = Dict("default" => Dict("E" =>0.0,"F" => 1.0, "V"=>0.0))

train_noE_data = [ AtomsData(at;  energy_key = "eenergy", force_key="forces",  # purposefully!!! misspelling energy_key
                        weights = weights_noE, v_ref=vref) for at in raw_data[train_idxs]]

val_noE_data = [ AtomsData(at;  energy_key = "eenergy", force_key="forces",
                        weights = weights_noE, v_ref=vref) for at in raw_data[val_idxs]]

A_noE, Y_noE, W_noE = ACEfit.assemble(train_noE_data, rpib_orig)

man_coeffs = (A_noE'*A_noE) \ (A_noE'* Y_noE)
pot_man_noE = manual_fit(A_noE,Y_noE,rpib_orig,vref)

ACE1pack.linear_errors(train_noE_data, pot_man_noE)
#=
[ Info: RMSE Table
┌─────────┬─────────┬──────────┬─────────┐
│    Type │ E [meV] │ F [eV/A] │ V [meV] │
├─────────┼─────────┼──────────┼─────────┤
│  amorph │   0.000 │    0.267 │   0.000 │
│ crystal │   0.000 │    0.193 │   0.000 │
├─────────┼─────────┼──────────┼─────────┤
│     set │   0.000 │    0.218 │   0.000 │
└─────────┴─────────┴──────────┴─────────┘
[ Info: MAE Table
┌─────────┬─────────┬──────────┬─────────┐
│    Type │ E [meV] │ F [eV/A] │ V [meV] │
├─────────┼─────────┼──────────┼─────────┤
│  amorph │   0.000 │    0.201 │   0.000 │
│ crystal │   0.000 │    0.142 │   0.000 │
├─────────┼─────────┼──────────┼─────────┤
│     set │   0.000 │    0.159 │   0.000 │
└─────────┴─────────┴──────────┴─────────┘
=#


train_noE2_data = [ AtomsData(at;  energy_key = "energy", force_key="forces",  # purposefully!!! misspelling energy_key
                        weights = weights_noE, v_ref=vref) for at in raw_data[train_idxs]]

val_noE2_data = [ AtomsData(at;  energy_key = "energy", force_key="forces",
                        weights = weights_noE, v_ref=vref) for at in raw_data[val_idxs]]

A_noE2, Y_noE2, W_noE2 = ACEfit.assemble(train_noE2_data, rpib_orig)

man_coeffs2 = (A_noE2'* Diagonal(W_noE2) * A_noE2) \ (A_noE2'* Diagonal(W_noE2) * Y_noE2)
pot_man_2 = JuLIP.MLIPs.combine(rpib_orig,man_coeffs2)

ACE1pack.linear_errors(train_noE_data, pot_man_2)
#=
[ Info: RMSE Table
┌─────────┬─────────┬──────────┬─────────┐
│    Type │ E [meV] │ F [eV/A] │ V [meV] │
├─────────┼─────────┼──────────┼─────────┤
│  amorph │   0.000 │    0.267 │   0.000 │
│ crystal │   0.000 │    0.193 │   0.000 │
├─────────┼─────────┼──────────┼─────────┤
│     set │   0.000 │    0.218 │   0.000 │
└─────────┴─────────┴──────────┴─────────┘
[ Info: MAE Table
┌─────────┬─────────┬──────────┬─────────┐
│    Type │ E [meV] │ F [eV/A] │ V [meV] │
├─────────┼─────────┼──────────┼─────────┤
│  amorph │   0.000 │    0.201 │   0.000 │
│ crystal │   0.000 │    0.142 │   0.000 │
├─────────┼─────────┼──────────┼─────────┤
│     set │   0.000 │    0.159 │   0.000 │
└─────────┴─────────┴──────────┴─────────┘
=#
