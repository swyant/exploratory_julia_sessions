using ACE1pack
using CSV, DataFrames

ds_path = "../../../../datasets/HfO2_Sivaraman/prl_2021/raw/train.xyz" 
raw_data = read_extxyz(ds_path)

train_idxs = vec(Matrix(CSV.read("Siviraman_HfO2_my_train_idxs.csv",DataFrame,header=false)))
val_idxs   = vec(Matrix(CSV.read("Siviraman_HfO2_my_val_idxs.csv",DataFrame,header=false)))

data_train = raw_data[train_idxs]
data_val   = raw_data[val_idxs]

model = acemodel(elements = [:Hf, :O],
                 order = 3,
                 totaldegree = 10,
                 rcut = 5.0,
                 Eref = [:Hf => -2.70516846, :O => -0.01277342]
                )
#
solver_qr_reg = ACEfit.QR(;lambda=1e-4)
P = smoothness_prior(model; p=3)

weights = Dict("default" => Dict("E" =>1.0,"F" => 1.0, "V"=>0.0))

acefit!(model,data_train;solver=solver_qr_reg, prior=P, weights=weights, force_key ="forces")

@info("Training Error Table")
ACE1pack.linear_errors(data_train, model; weights=weights, force_key="forces")
#=
[ Info: RMSE Table
┌─────────┬─────────┬──────────┬─────────┐
│    Type │ E [meV] │ F [eV/A] │ V [meV] │
├─────────┼─────────┼──────────┼─────────┤
│  amorph │  12.819 │    0.417 │   0.000 │
│ crystal │   8.859 │    0.284 │   0.000 │
├─────────┼─────────┼──────────┼─────────┤
│     set │  10.221 │    0.329 │   0.000 │
└─────────┴─────────┴──────────┴─────────┘
[ Info: MAE Table
┌─────────┬─────────┬──────────┬─────────┐
│    Type │ E [meV] │ F [eV/A] │ V [meV] │
├─────────┼─────────┼──────────┼─────────┤
│  amorph │  10.349 │    0.311 │   0.000 │
│ crystal │   6.727 │    0.208 │   0.000 │
├─────────┼─────────┼──────────┼─────────┤
│     set │   7.823 │    0.238 │   0.000 │
└─────────┴─────────┴──────────┴─────────┘
=#

@info("Val Error Table")
ACE1pack.linear_errors(data_val, model; weights=weights, force_key="forces")
#= 
[ Info: RMSE Table
┌─────────┬─────────┬──────────┬─────────┐
│    Type │ E [meV] │ F [eV/A] │ V [meV] │
├─────────┼─────────┼──────────┼─────────┤
│  amorph │  14.682 │    0.425 │   0.000 │
│ crystal │   9.344 │    0.278 │   0.000 │
├─────────┼─────────┼──────────┼─────────┤
│     set │  11.358 │    0.332 │   0.000 │
└─────────┴─────────┴──────────┴─────────┘
[ Info: MAE Table
┌─────────┬─────────┬──────────┬─────────┐
│    Type │ E [meV] │ F [eV/A] │ V [meV] │
├─────────┼─────────┼──────────┼─────────┤
│  amorph │  11.679 │    0.319 │   0.000 │
│ crystal │   7.382 │    0.202 │   0.000 │
├─────────┼─────────┼──────────┼─────────┤
│     set │   8.780 │    0.240 │   0.000 │
└─────────┴─────────┴──────────┴─────────┘
=#