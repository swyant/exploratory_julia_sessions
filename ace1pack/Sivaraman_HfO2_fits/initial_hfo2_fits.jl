using ACE1 
using ACE1pack
using ACEfit
using LinearAlgebra: I, Diagonal, pinv
using JuLIP
using CSV, DataFrames
using Random

ds_path = "../../../../datasets/HfO2_Sivaraman/prl_2021/raw/train.xyz" 
raw_data = read_extxyz(ds_path)

### Just did this once
#total_num_configs = length(raw_data)
#num_train = Int(floor(0.9*total_num_configs))
#num_val = total_num_configs-num_train
#
#perm_indxs = Random.randperm(total_num_configs)
#
#train_indxs = perm_indxs[begin:1:num_train]
#val_indxs = perm_indxs[num_train+1:1:end]
#
#CSV.write("Siviraman_HfO2_my_train_idxs.csv", DataFrame(Tables.table(train_indxs)),header=false)
#CSV.write("Siviraman_HfO2_my_val_idxs.csv", DataFrame(Tables.table(val_indxs)),header=false)

train_idxs = vec(Matrix(CSV.read("Siviraman_HfO2_my_train_idxs.csv",DataFrame,header=false)))
val_idxs   = vec(Matrix(CSV.read("Siviraman_HfO2_my_val_idxs.csv",DataFrame,header=false)))

rpib1 = ACE1.rpi_basis(;
            species = [:Hf, :O],
            N       = 3, 
            maxdeg  = 8,
            D       = ACE1.SparsePSHDegree(; wL = 1.5, csp = 1.0),
            r0      = 2.15,
            rin     = 0.65*2.15,  # Does this meaningfully get used in pin=0?
            rcut    = 5.0,
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

# Emmanuel approach
W_man = ones(Float64,length(W1))
man_eman_coeffs = (A1'* Diagonal(W_man) *A1) \ (A1'* Diagonal(W_man) * Y1 )
@show man_eman_coeffs

pot_emmanuel_temp = JuLIP.MLIPs.combine(rpib1,man_eman_coeffs)
pot_emmanuel = JuLIP.MLIPs.SumIP(pot_emmanuel_temp, vref)
ACE1pack.linear_errors(train_data,pot_emmanuel)
#=
[ Info: RMSE Table
┌─────────┬─────────┬──────────┬─────────┐
│    Type │ E [meV] │ F [eV/A] │ V [meV] │
├─────────┼─────────┼──────────┼─────────┤
│  amorph │   6.702 │    0.294 │   0.000 │
│ crystal │   5.417 │    0.220 │   0.000 │
├─────────┼─────────┼──────────┼─────────┤
│     set │   5.836 │    0.244 │   0.000 │
└─────────┴─────────┴──────────┴─────────┘
[ Info: MAE Table
┌─────────┬─────────┬──────────┬─────────┐
│    Type │ E [meV] │ F [eV/A] │ V [meV] │
├─────────┼─────────┼──────────┼─────────┤
│  amorph │   5.401 │    0.221 │   0.000 │
│ crystal │   4.418 │    0.162 │   0.000 │
├─────────┼─────────┼──────────┼─────────┤
│     set │   4.715 │    0.179 │   0.000 │
└─────────┴─────────┴──────────┴─────────┘
=#
ACE1pack.linear_errors(val_data,pot_emmanuel)
#=
[ Info: RMSE Table
┌─────────┬─────────┬──────────┬─────────┐
│    Type │ E [meV] │ F [eV/A] │ V [meV] │
├─────────┼─────────┼──────────┼─────────┤
│  amorph │   6.516 │    0.303 │   0.000 │
│ crystal │   5.867 │    0.216 │   0.000 │
├─────────┼─────────┼──────────┼─────────┤
│     set │   6.086 │    0.247 │   0.000 │
└─────────┴─────────┴──────────┴─────────┘
[ Info: MAE Table
┌─────────┬─────────┬──────────┬─────────┐
│    Type │ E [meV] │ F [eV/A] │ V [meV] │
├─────────┼─────────┼──────────┼─────────┤
│  amorph │   5.244 │    0.228 │   0.000 │
│ crystal │   4.789 │    0.157 │   0.000 │
├─────────┼─────────┼──────────┼─────────┤
│     set │   4.937 │    0.179 │   0.000 │
└─────────┴─────────┴──────────┴─────────┘
=#

### Default QR
solver_qr = ACEfit.QR()

results_qr = ACEfit.solve(solver_qr,A1,Y1)

@show results_qr["C"]

pot_qr_tmp = JuLIP.MLIPs.combine(rpib1,results_qr["C"])
pot_qr = JuLIP.MLIPs.SumIP(pot_qr_tmp,vref)
ACE1pack.linear_errors(train_data,pot_qr)
#=
[ Info: RMSE Table
┌─────────┬─────────┬──────────┬─────────┐
│    Type │ E [meV] │ F [eV/A] │ V [meV] │
├─────────┼─────────┼──────────┼─────────┤
│  amorph │   6.707 │    0.293 │   0.000 │
│ crystal │   5.439 │    0.220 │   0.000 │
├─────────┼─────────┼──────────┼─────────┤
│     set │   5.852 │    0.244 │   0.000 │
└─────────┴─────────┴──────────┴─────────┘
[ Info: MAE Table
┌─────────┬─────────┬──────────┬─────────┐
│    Type │ E [meV] │ F [eV/A] │ V [meV] │
├─────────┼─────────┼──────────┼─────────┤
│  amorph │   5.403 │    0.220 │   0.000 │
│ crystal │   4.443 │    0.162 │   0.000 │
├─────────┼─────────┼──────────┼─────────┤
│     set │   4.734 │    0.179 │   0.000 │
└─────────┴─────────┴──────────┴─────────┘
=#
ACE1pack.linear_errors(val_data,pot_qr)
#=
[ Info: RMSE Table
┌─────────┬─────────┬──────────┬─────────┐
│    Type │ E [meV] │ F [eV/A] │ V [meV] │
├─────────┼─────────┼──────────┼─────────┤
│  amorph │   6.493 │    0.302 │   0.000 │
│ crystal │   5.928 │    0.216 │   0.000 │
├─────────┼─────────┼──────────┼─────────┤
│     set │   6.117 │    0.246 │   0.000 │
└─────────┴─────────┴──────────┴─────────┘
[ Info: MAE Table
┌─────────┬─────────┬──────────┬─────────┐
│    Type │ E [meV] │ F [eV/A] │ V [meV] │
├─────────┼─────────┼──────────┼─────────┤
│  amorph │   5.250 │    0.227 │   0.000 │
│ crystal │   4.854 │    0.157 │   0.000 │
├─────────┼─────────┼──────────┼─────────┤
│     set │   4.983 │    0.179 │   0.000 │
└─────────┴─────────┴──────────┴─────────┘
=#


### QR w/ regularization
solver_reg = ACEfit.QR(; lambda=1e-2)

results_reg = ACEfit.solve(solver_reg,A1,Y1)

@show results_reg["C"]

pot_reg_tmp = JuLIP.MLIPs.combine(rpib1,results_reg["C"])
pot_reg = JuLIP.MLIPs.SumIP(pot_reg_tmp,vref)
ACE1pack.linear_errors(train_data,pot_reg)
#=
[ Info: RMSE Table
┌─────────┬─────────┬──────────┬─────────┐
│    Type │ E [meV] │ F [eV/A] │ V [meV] │
├─────────┼─────────┼──────────┼─────────┤
│  amorph │   6.702 │    0.294 │   0.000 │
│ crystal │   5.456 │    0.221 │   0.000 │
├─────────┼─────────┼──────────┼─────────┤
│     set │   5.861 │    0.245 │   0.000 │
└─────────┴─────────┴──────────┴─────────┘
[ Info: MAE Table
┌─────────┬─────────┬──────────┬─────────┐
│    Type │ E [meV] │ F [eV/A] │ V [meV] │
├─────────┼─────────┼──────────┼─────────┤
│  amorph │   5.388 │    0.221 │   0.000 │
│ crystal │   4.471 │    0.162 │   0.000 │
├─────────┼─────────┼──────────┼─────────┤
│     set │   4.748 │    0.180 │   0.000 │
└─────────┴─────────┴──────────┴─────────┘
=#

ACE1pack.linear_errors(val_data,pot_reg)
#=
[ Info: RMSE Table
┌─────────┬─────────┬──────────┬─────────┐
│    Type │ E [meV] │ F [eV/A] │ V [meV] │
├─────────┼─────────┼──────────┼─────────┤
│  amorph │   6.472 │    0.303 │   0.000 │
│ crystal │   5.916 │    0.217 │   0.000 │
├─────────┼─────────┼──────────┼─────────┤
│     set │   6.103 │    0.247 │   0.000 │
└─────────┴─────────┴──────────┴─────────┘
[ Info: MAE Table
┌─────────┬─────────┬──────────┬─────────┐
│    Type │ E [meV] │ F [eV/A] │ V [meV] │
├─────────┼─────────┼──────────┼─────────┤
│  amorph │   5.236 │    0.228 │   0.000 │
│ crystal │   4.865 │    0.157 │   0.000 │
├─────────┼─────────┼──────────┼─────────┤
│     set │   4.986 │    0.180 │   0.000 │
└─────────┴─────────┴──────────┴─────────┘
=#


####### RRQR 
solver_rrqr = ACEfit.RRQR()

results_rrqr = ACEfit.solve(solver_rrqr, A1, Y1)
@show results_rrqr["C"]

pot_rrqr_tmp = JuLIP.MLIPs.combine(rpib1,results_rrqr["C"])
pot_rrqr = JuLIP.MLIPs.SumIP(pot_rrqr_tmp,vref)

ACE1pack.linear_errors(train_data,pot_rrqr)
#=
[ Info: RMSE Table
┌─────────┬─────────┬──────────┬─────────┐
│    Type │ E [meV] │ F [eV/A] │ V [meV] │
├─────────┼─────────┼──────────┼─────────┤
│  amorph │   6.673 │    0.293 │   0.000 │
│ crystal │   5.423 │    0.220 │   0.000 │
├─────────┼─────────┼──────────┼─────────┤
│     set │   5.830 │    0.244 │   0.000 │
└─────────┴─────────┴──────────┴─────────┘
[ Info: MAE Table
┌─────────┬─────────┬──────────┬─────────┐
│    Type │ E [meV] │ F [eV/A] │ V [meV] │
├─────────┼─────────┼──────────┼─────────┤
│  amorph │   5.375 │    0.220 │   0.000 │
│ crystal │   4.421 │    0.162 │   0.000 │
├─────────┼─────────┼──────────┼─────────┤
│     set │   4.710 │    0.179 │   0.000 │
└─────────┴─────────┴──────────┴─────────┘
=#

ACE1pack.linear_errors(val_data,pot_rrqr)
#=
[ Info: RMSE Table
┌─────────┬─────────┬──────────┬─────────┐
│    Type │ E [meV] │ F [eV/A] │ V [meV] │
├─────────┼─────────┼──────────┼─────────┤
│  amorph │   6.491 │    0.302 │   0.000 │
│ crystal │   5.870 │    0.216 │   0.000 │
├─────────┼─────────┼──────────┼─────────┤
│     set │   6.079 │    0.246 │   0.000 │
└─────────┴─────────┴──────────┴─────────┘
[ Info: MAE Table
┌─────────┬─────────┬──────────┬─────────┐
│    Type │ E [meV] │ F [eV/A] │ V [meV] │
├─────────┼─────────┼──────────┼─────────┤
│  amorph │   5.228 │    0.227 │   0.000 │
│ crystal │   4.795 │    0.157 │   0.000 │
├─────────┼─────────┼──────────┼─────────┤
│     set │   4.936 │    0.179 │   0.000 │
=#


### Default LSQR 
solver_lsqr = ACEfit.LSQR()

results_lsqr = ACEfit.solve(solver_lsqr, A1, Y1)
@show results_lsqr["C"]

pot_lsqr_tmp = JuLIP.MLIPs.combine(rpib1,results_lsqr["C"])
pot_lsqr = JuLIP.MLIPs.SumIP(pot_lsqr_tmp,vref)

ACE1pack.linear_errors(train_data,pot_lsqr)
#=
[ Info: RMSE Table
┌─────────┬─────────┬──────────┬─────────┐
│    Type │ E [meV] │ F [eV/A] │ V [meV] │
├─────────┼─────────┼──────────┼─────────┤
│  amorph │   7.415 │    0.311 │   0.000 │
│ crystal │   5.206 │    0.243 │   0.000 │
├─────────┼─────────┼──────────┼─────────┤
│     set │   5.962 │    0.265 │   0.000 │
└─────────┴─────────┴──────────┴─────────┘
[ Info: MAE Table
┌─────────┬─────────┬──────────┬─────────┐
│    Type │ E [meV] │ F [eV/A] │ V [meV] │
├─────────┼─────────┼──────────┼─────────┤
│  amorph │   5.915 │    0.233 │   0.000 │
│ crystal │   4.205 │    0.175 │   0.000 │
├─────────┼─────────┼──────────┼─────────┤
│     set │   4.723 │    0.192 │   0.000 │
└─────────┴─────────┴──────────┴─────────┘
=#

ACE1pack.linear_errors(val_data,pot_lsqr)
#=
[ Info: RMSE Table
┌─────────┬─────────┬──────────┬─────────┐
│    Type │ E [meV] │ F [eV/A] │ V [meV] │
├─────────┼─────────┼──────────┼─────────┤
│  amorph │   6.942 │    0.317 │   0.000 │
│ crystal │   5.833 │    0.238 │   0.000 │
├─────────┼─────────┼──────────┼─────────┤
│     set │   6.215 │    0.266 │   0.000 │
└─────────┴─────────┴──────────┴─────────┘
[ Info: MAE Table
┌─────────┬─────────┬──────────┬─────────┐
│    Type │ E [meV] │ F [eV/A] │ V [meV] │
├─────────┼─────────┼──────────┼─────────┤
│  amorph │   5.594 │    0.239 │   0.000 │
│ crystal │   4.649 │    0.170 │   0.000 │
├─────────┼─────────┼──────────┼─────────┤
│     set │   4.956 │    0.192 │   0.000 │
└─────────┴─────────┴──────────┴─────────┘
=#



#### BLR 
solver_blr = ACEfit.BLR()

results_blr = ACEfit.solve(solver_blr, A1, Y1)
@show results_blr["C"]

pot_blr_tmp = JuLIP.MLIPs.combine(rpib1,results_blr["C"])
pot_blr = JuLIP.MLIPs.SumIP(pot_blr_tmp,vref)

ACE1pack.linear_errors(train_data,pot_blr)
#=
[ Info: RMSE Table
┌─────────┬─────────┬──────────┬─────────┐
│    Type │ E [meV] │ F [eV/A] │ V [meV] │
├─────────┼─────────┼──────────┼─────────┤
│  amorph │   6.675 │    0.293 │   0.000 │
│ crystal │   5.438 │    0.221 │   0.000 │
├─────────┼─────────┼──────────┼─────────┤
│     set │   5.840 │    0.244 │   0.000 │
└─────────┴─────────┴──────────┴─────────┘
[ Info: MAE Table
┌─────────┬─────────┬──────────┬─────────┐
│    Type │ E [meV] │ F [eV/A] │ V [meV] │
├─────────┼─────────┼──────────┼─────────┤
│  amorph │   5.366 │    0.221 │   0.000 │
│ crystal │   4.444 │    0.162 │   0.000 │
├─────────┼─────────┼──────────┼─────────┤
│     set │   4.723 │    0.179 │   0.000 │
└─────────┴─────────┴──────────┴─────────┘
=#

#################

rpib2 = ACE1.rpi_basis(;
            species = [:Hf, :O],
            N       = 3, 
            maxdeg  = 12,
            D       = ACE1.SparsePSHDegree(; wL = 1.5, csp = 1.0),
            r0      = 2.15,
            rin     = 0.65*2.15,  # Does this meaningfully get used in pin=0?
            rcut    = 5.0,
            pin     = 0,
)

@show length(rpib2) #2408


############## 
#A2, Y2, W2 = ACEfit.assemble(train_data, rpib2)

#CSV.write("rpib2_design_mat.csv", DataFrame(Tables.table(A2)),header=false) # This would've been like > 20 GBs

rpib3 = ACE1.rpi_basis(;
            species = [:Hf, :O],
            N       = 3, 
            maxdeg  = 10,
            D       = ACE1.SparsePSHDegree(; wL = 1.5, csp = 1.0),
            r0      = 2.15,
            rin     = 0.65*2.15,  # Does this meaningfully get used in pin=0?
            rcut    = 5.0,
            pin     = 0,
)

@show length(rpib3)

A3, Y3, W3 = ACEfit.assemble(train_data, rpib3)


#### Manual approach
W3_man = ones(Float64,length(W3))
man_eman_coeffs3 = (A3'* Diagonal(W3_man) *A3) \ (A3'* Diagonal(W3_man) * Y3 )
@show man_eman_coeffs3

pot_emmanuel_temp3 = JuLIP.MLIPs.combine(rpib3,man_eman_coeffs3)
pot_emmanuel3 = JuLIP.MLIPs.SumIP(pot_emmanuel_temp3, vref)
ACE1pack.linear_errors(train_data,pot_emmanuel3)
#=
[ Info: RMSE Table
┌─────────┬─────────┬──────────┬─────────┐
│    Type │ E [meV] │ F [eV/A] │ V [meV] │
├─────────┼─────────┼──────────┼─────────┤
│  amorph │   5.946 │    0.271 │   0.000 │
│ crystal │   5.021 │    0.195 │   0.000 │
├─────────┼─────────┼──────────┼─────────┤
│     set │   5.318 │    0.220 │   0.000 │
└─────────┴─────────┴──────────┴─────────┘
[ Info: MAE Table
┌─────────┬─────────┬──────────┬─────────┐
│    Type │ E [meV] │ F [eV/A] │ V [meV] │
├─────────┼─────────┼──────────┼─────────┤
│  amorph │   4.655 │    0.203 │   0.000 │
│ crystal │   3.920 │    0.143 │   0.000 │
├─────────┼─────────┼──────────┼─────────┤
│     set │   4.143 │    0.161 │   0.000 │
└─────────┴─────────┴──────────┴─────────┘
=# 
ACE1pack.linear_errors(val_data,pot_emmanuel3)
#=
[ Info: RMSE Table
┌─────────┬─────────┬──────────┬─────────┐
│    Type │ E [meV] │ F [eV/A] │ V [meV] │
├─────────┼─────────┼──────────┼─────────┤
│  amorph │   6.050 │    0.282 │   0.000 │
│ crystal │   5.480 │    0.192 │   0.000 │
├─────────┼─────────┼──────────┼─────────┤
│     set │   5.672 │    0.224 │   0.000 │
└─────────┴─────────┴──────────┴─────────┘
[ Info: MAE Table
┌─────────┬─────────┬──────────┬─────────┐
│    Type │ E [meV] │ F [eV/A] │ V [meV] │
├─────────┼─────────┼──────────┼─────────┤
│  amorph │   4.729 │    0.212 │   0.000 │
│ crystal │   4.298 │    0.139 │   0.000 │
├─────────┼─────────┼──────────┼─────────┤
│     set │   4.438 │    0.162 │   0.000 │
└─────────┴─────────┴──────────┴─────────┘
=#


#### regularized QR
solver_reg = ACEfit.QR(; lambda=1e-3)

results_reg3 = ACEfit.solve(solver_reg,A3,Y3)

@show results_reg3["C"]

pot_reg_tmp3 = JuLIP.MLIPs.combine(rpib3,results_reg3["C"])
pot_reg3 = JuLIP.MLIPs.SumIP(pot_reg_tmp3,vref)
ACE1pack.linear_errors(train_data,pot_reg3)
#=
[ Info: RMSE Table
┌─────────┬─────────┬──────────┬─────────┐
│    Type │ E [meV] │ F [eV/A] │ V [meV] │
├─────────┼─────────┼──────────┼─────────┤
│  amorph │   5.961 │    0.268 │   0.000 │
│ crystal │   4.987 │    0.194 │   0.000 │
├─────────┼─────────┼──────────┼─────────┤
│     set │   5.301 │    0.219 │   0.000 │
└─────────┴─────────┴──────────┴─────────┘
[ Info: MAE Table
┌─────────┬─────────┬──────────┬─────────┐
│    Type │ E [meV] │ F [eV/A] │ V [meV] │
├─────────┼─────────┼──────────┼─────────┤
│  amorph │   4.670 │    0.201 │   0.000 │
│ crystal │   3.900 │    0.142 │   0.000 │
├─────────┼─────────┼──────────┼─────────┤
│     set │   4.133 │    0.160 │   0.000 │
└─────────┴─────────┴──────────┴─────────┘
=#


rpib4 = ACE1.rpi_basis(;
            species = [:Hf, :O],
            N       = 4, 
            maxdeg  = 9,
            D       = ACE1.SparsePSHDegree(; wL = 1.5, csp = 1.0),
            r0      = 2.15,
            rin     = 0.65*2.15,  # Does this meaningfully get used in pin=0?
            rcut    = 5.0,
            pin     = 0,
)

@show length(rpib4)

A4, Y4, W4 = ACEfit.assemble(train_data, rpib4)

solver_reg = ACEfit.QR(; lambda=1e-3)

results_reg4 = ACEfit.solve(solver_reg,A4,Y4)

@show results_reg4["C"]

pot_reg_tmp4 = JuLIP.MLIPs.combine(rpib4,results_reg4["C"])
pot_reg4 = JuLIP.MLIPs.SumIP(pot_reg_tmp4,vref)
ACE1pack.linear_errors(train_data,pot_reg4)
#=
[ Info: RMSE Table
┌─────────┬─────────┬──────────┬─────────┐
│    Type │ E [meV] │ F [eV/A] │ V [meV] │
├─────────┼─────────┼──────────┼─────────┤
│  amorph │   5.689 │    0.267 │   0.000 │
│ crystal │   4.883 │    0.193 │   0.000 │
├─────────┼─────────┼──────────┼─────────┤
│     set │   5.140 │    0.218 │   0.000 │
└─────────┴─────────┴──────────┴─────────┘
[ Info: MAE Table
┌─────────┬─────────┬──────────┬─────────┐
│    Type │ E [meV] │ F [eV/A] │ V [meV] │
├─────────┼─────────┼──────────┼─────────┤
│  amorph │   4.513 │    0.200 │   0.000 │
│ crystal │   3.900 │    0.143 │   0.000 │
├─────────┼─────────┼──────────┼─────────┤
│     set │   4.085 │    0.160 │   0.000 │
└─────────┴─────────┴──────────┴─────────┘
=#
