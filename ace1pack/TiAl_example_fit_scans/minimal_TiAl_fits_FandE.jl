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

data = [ AtomsData(at;  energy_key = "energy", force_key="force",
                        weights = weights) for at in raw_data[1:5:end]]


A, Y, W = ACEfit.assemble(data, rpib)

W_man = ones(Float64,length(W))
Aw_man = Diagonal(W_man) * A
Yw_man = W_man .* Y

#### Manual PL.jl approach
AtA = transpose(A)*Aw_man
Atb = transpose(A)*Yw_man
Q = pinv(AtA, 1e-6) . # So 1e-6 is probably cuts out too much, because I get a much worse fit here. 1e-10 almost matched though
manual_coeffs = Q*Atb
@show manual_coeffs

pot_manual = JuLIP.MLIPs.combine(rpib,manual_coeffs)
ACE1pack.linear_errors(data,pot_manual)
#=
pinv α=1e-6
[ Info: RMSE Table
┌────────────┬───────────┬──────────┬─────────┐
│       Type │   E [meV] │ F [eV/A] │ V [meV] │
├────────────┼───────────┼──────────┼─────────┤
│   FLD_TiAl │ 36227.875 │   53.727 │   0.000 │
│ TiAl_T5000 │  5805.594 │   88.854 │   0.000 │
├────────────┼───────────┼──────────┼─────────┤
│        set │ 35416.574 │   80.297 │   0.000 │
└────────────┴───────────┴──────────┴─────────┘
[ Info: MAE Table
┌────────────┬───────────┬──────────┬─────────┐
│       Type │   E [meV] │ F [eV/A] │ V [meV] │
├────────────┼───────────┼──────────┼─────────┤
│   FLD_TiAl │ 31777.996 │   36.428 │   0.000 │
│ TiAl_T5000 │  3924.128 │   70.158 │   0.000 │
├────────────┼───────────┼──────────┼─────────┤
│        set │ 30511.911 │   60.410 │   0.000 │
└────────────┴───────────┴──────────┴─────────┘

pinv α=1e-10
[ Info: RMSE Table
┌────────────┬───────────┬──────────┬─────────┐
│       Type │   E [meV] │ F [eV/A] │ V [meV] │
├────────────┼───────────┼──────────┼─────────┤
│   FLD_TiAl │ 10729.022 │    3.363 │   0.000 │
│ TiAl_T5000 │   160.986 │    8.948 │   0.000 │
├────────────┼───────────┼──────────┼─────────┤
│        set │ 10482.401 │    7.759 │   0.000 │
└────────────┴───────────┴──────────┴─────────┘
[ Info: MAE Table
┌────────────┬──────────┬──────────┬─────────┐
│       Type │  E [meV] │ F [eV/A] │ V [meV] │
├────────────┼──────────┼──────────┼─────────┤
│   FLD_TiAl │ 6171.808 │    1.934 │   0.000 │
│ TiAl_T5000 │  145.346 │    6.933 │   0.000 │
├────────────┼──────────┼──────────┼─────────┤
│        set │ 5897.878 │    5.488 │   0.000 │
└────────────┴──────────┴──────────┴─────────┘
=#



### Emmanuel approach 
man_eman_coeffs = (A'* Diagonal(W_man) *A) \ (A'* Diagonal(W_man) * Y )
@show man_eman_coeffs

pot_emmanuel = JuLIP.MLIPs.combine(rpib,man_eman_coeffs)
ACE1pack.linear_errors(data,pot_emmanuel)
#=
[ Info: RMSE Table
┌────────────┬───────────┬──────────┬─────────┐
│       Type │   E [meV] │ F [eV/A] │ V [meV] │
├────────────┼───────────┼──────────┼─────────┤
│   FLD_TiAl │ 10729.294 │    3.311 │   0.000 │
│ TiAl_T5000 │   159.834 │    8.944 │   0.000 │
├────────────┼───────────┼──────────┼─────────┤
│        set │ 10482.666 │    7.749 │   0.000 │
└────────────┴───────────┴──────────┴─────────┘
[ Info: MAE Table
┌────────────┬──────────┬──────────┬─────────┐
│       Type │  E [meV] │ F [eV/A] │ V [meV] │
├────────────┼──────────┼──────────┼─────────┤
│   FLD_TiAl │ 6182.621 │    1.894 │   0.000 │
│ TiAl_T5000 │  144.354 │    6.939 │   0.000 │
├────────────┼──────────┼──────────┼─────────┤
│        set │ 5908.155 │    5.481 │   0.000 │
└────────────┴──────────┴──────────┴─────────┘
=#
forces(pot_emmanuel, raw_data[2])
energy(pot_emmanuel, raw_data[2])
check_data = [ AtomsData(at;  energy_key = "energy", force_key="force",
                weights = weights) for at in [raw_data[2]]]
ACE1pack.linear_errors(check_data,pot_emmanuel)


# QR with no weights (i.e. W = Diag(ones)), no prior
solver_qr = ACEfit.QR()

results_qr = ACEfit.solve(solver_qr,A,Y)

@show results_qr["C"]

pot_qr = JuLIP.MLIPs.combine(rpib,results_qr["C"])
ACE1pack.linear_errors(data,pot_qr)
#=
[ Info: RMSE Table
┌────────────┬───────────┬──────────┬─────────┐
│       Type │   E [meV] │ F [eV/A] │ V [meV] │
├────────────┼───────────┼──────────┼─────────┤
│   FLD_TiAl │ 10652.180 │    3.389 │   0.000 │
│ TiAl_T5000 │   503.936 │    8.909 │   0.000 │
├────────────┼───────────┼──────────┼─────────┤
│        set │ 10407.824 │    7.730 │   0.000 │
└────────────┴───────────┴──────────┴─────────┘
[ Info: MAE Table
┌────────────┬──────────┬──────────┬─────────┐
│       Type │  E [meV] │ F [eV/A] │ V [meV] │
├────────────┼──────────┼──────────┼─────────┤
│   FLD_TiAl │ 6212.777 │    1.879 │   0.000 │
│ TiAl_T5000 │  367.370 │    6.915 │   0.000 │
├────────────┼──────────┼──────────┼─────────┤
│        set │ 5947.077 │    5.459 │   0.000 │
└────────────┴──────────┴──────────┴─────────┘
=#


# QR with regularization 
solver_reg = ACEfit.QR(; lambda=1e-3) # varying lambda didn't do too much

results_reg = ACEfit.solve(solver_reg,A,Y)

@show results_reg["C"]

pot_reg = JuLIP.MLIPs.combine(rpib,results_reg["C"])
ACE1pack.linear_errors(data,pot_reg)
#=
lambda=1e-3
[ Info: RMSE Table
┌────────────┬───────────┬──────────┬─────────┐
│       Type │   E [meV] │ F [eV/A] │ V [meV] │
├────────────┼───────────┼──────────┼─────────┤
│   FLD_TiAl │ 10729.291 │    3.311 │   0.000 │
│ TiAl_T5000 │   159.834 │    8.944 │   0.000 │
├────────────┼───────────┼──────────┼─────────┤
│        set │ 10482.663 │    7.749 │   0.000 │
└────────────┴───────────┴──────────┴─────────┘
[ Info: MAE Table
┌────────────┬──────────┬──────────┬─────────┐
│       Type │  E [meV] │ F [eV/A] │ V [meV] │
├────────────┼──────────┼──────────┼─────────┤
│   FLD_TiAl │ 6182.619 │    1.894 │   0.000 │
│ TiAl_T5000 │  144.354 │    6.939 │   0.000 │
├────────────┼──────────┼──────────┼─────────┤
│        set │ 5908.152 │    5.481 │   0.000 │
=#

# RRQR 
solver_rrqr = ACEfit.RRQR()

results_rrqr = ACEfit.solve(solver_rrqr, A, Y)
@show results_rrqr["C"]

pot_rrqr = JuLIP.MLIPs.combine(rpib,results_rrqr["C"])
ACE1pack.linear_errors(data,pot_rrqr)

#=
[ Info: RMSE Table
┌────────────┬───────────┬──────────┬─────────┐
│       Type │   E [meV] │ F [eV/A] │ V [meV] │
├────────────┼───────────┼──────────┼─────────┤
│   FLD_TiAl │ 10729.294 │    3.311 │   0.000 │
│ TiAl_T5000 │   159.834 │    8.944 │   0.000 │
├────────────┼───────────┼──────────┼─────────┤
│        set │ 10482.666 │    7.749 │   0.000 │
└────────────┴───────────┴──────────┴─────────┘
[ Info: MAE Table
┌────────────┬──────────┬──────────┬─────────┐
│       Type │  E [meV] │ F [eV/A] │ V [meV] │
├────────────┼──────────┼──────────┼─────────┤
│   FLD_TiAl │ 6182.621 │    1.894 │   0.000 │
│ TiAl_T5000 │  144.354 │    6.939 │   0.000 │
├────────────┼──────────┼──────────┼─────────┤
│        set │ 5908.155 │    5.481 │   0.000 │
└────────────┴──────────┴──────────┴─────────┘
=#

### Default LSQR 
solver_lsqr = ACEfit.LSQR()

results_lsqr = ACEfit.solve(solver_lsqr, A, Y)
@show results_lsqr["C"]

pot_lsqr = JuLIP.MLIPs.combine(rpib,results_lsqr["C"])
ACE1pack.linear_errors(data,pot_lsqr)

#=
[ Info: RMSE Table
┌────────────┬───────────┬──────────┬─────────┐
│       Type │   E [meV] │ F [eV/A] │ V [meV] │
├────────────┼───────────┼──────────┼─────────┤
│   FLD_TiAl │ 10735.034 │    3.378 │   0.000 │
│ TiAl_T5000 │   145.604 │    8.977 │   0.000 │
├────────────┼───────────┼──────────┼─────────┤
│        set │ 10488.265 │    7.784 │   0.000 │
└────────────┴───────────┴──────────┴─────────┘
[ Info: MAE Table
┌────────────┬──────────┬──────────┬─────────┐
│       Type │  E [meV] │ F [eV/A] │ V [meV] │
├────────────┼──────────┼──────────┼─────────┤
│   FLD_TiAl │ 6156.919 │    1.952 │   0.000 │
│ TiAl_T5000 │  136.136 │    6.977 │   0.000 │
├────────────┼──────────┼──────────┼─────────┤
│        set │ 5883.247 │    5.525 │   0.000 │
└────────────┴──────────┴──────────┴─────────┘
=#


# modifying weights 
weights2 = Dict(
        "FLD_TiAl" => Dict("E" => 60.0, "F" => 1.0 , "V" => 0.0 ),
        "TiAl_T5000" => Dict("E" => 60.0, "F" => 1.0 , "V" => 0.0 ));

data2 = [ AtomsData(at;  energy_key = "energy", force_key="force",
                    weights = weights2) for at in raw_data[1:5:end]]


A2, Y2, W2 = ACEfit.assemble(data2, rpib)

A_w2 = Diagonal(W2) * A2
Y_w2 = W2 .* Y2

results_w2 = ACEfit.solve(solver_qr,A_w2,Y_w2)

@show results_w2["C"]

pot_w2 = JuLIP.MLIPs.combine(rpib,results_w2["C"])
ACE1pack.linear_errors(data2,pot_w2) # personal sanity check: swapping data2 w/ data here changes nothing

#=

    "FLD_TiAl" => Dict("E" => 60.0, "F" => 1.0 , "V" => 0.0 ),
    "TiAl_T5000" => Dict("E" => 60.0, "F" => 1.0 , "V" => 0.0 )
[ Info: RMSE Table
┌────────────┬─────────┬──────────┬─────────┐
│       Type │ E [meV] │ F [eV/A] │ V [meV] │
├────────────┼─────────┼──────────┼─────────┤
│   FLD_TiAl │ 506.151 │    7.515 │   0.000 │
│ TiAl_T5000 │  55.992 │   23.803 │   0.000 │
├────────────┼─────────┼──────────┼─────────┤
│        set │ 494.658 │   20.474 │   0.000 │
└────────────┴─────────┴──────────┴─────────┘
[ Info: MAE Table
┌────────────┬─────────┬──────────┬─────────┐
│       Type │ E [meV] │ F [eV/A] │ V [meV] │
├────────────┼─────────┼──────────┼─────────┤
│   FLD_TiAl │ 342.348 │    4.360 │   0.000 │
│ TiAl_T5000 │  37.437 │   18.772 │   0.000 │
├────────────┼─────────┼──────────┼─────────┤
│        set │ 328.489 │   14.607 │   0.000 │
└────────────┴─────────┴──────────┴─────────┘

    "FLD_TiAl" => Dict("E" => 60.0, "F" => 1.0 , "V" => 0.0 ),
    "TiAl_T5000" => Dict("E" => 5.0, "F" => 1.0 , "V" => 0.0 )
[ Info: RMSE Table
┌────────────┬─────────┬──────────┬─────────┐
│       Type │ E [meV] │ F [eV/A] │ V [meV] │
├────────────┼─────────┼──────────┼─────────┤
│   FLD_TiAl │ 518.747 │    7.927 │   0.000 │
│ TiAl_T5000 │ 815.081 │   23.422 │   0.000 │
├────────────┼─────────┼──────────┼─────────┤
│        set │ 535.785 │   20.205 │   0.000 │
└────────────┴─────────┴──────────┴─────────┘
[ Info: MAE Table
┌────────────┬─────────┬──────────┬─────────┐
│       Type │ E [meV] │ F [eV/A] │ V [meV] │
├────────────┼─────────┼──────────┼─────────┤
│   FLD_TiAl │ 361.515 │    4.533 │   0.000 │
│ TiAl_T5000 │ 715.384 │   18.494 │   0.000 │
├────────────┼─────────┼──────────┼─────────┤
│        set │ 377.600 │   14.460 │   0.000 │
└────────────┴─────────┴──────────┴─────────┘
=#


### Finally, manually adding a prior, but using no weights
dp = Float64[]
append!(dp, ACE1.scaling(rpib, 3, 1.0)) # p=3 smoothness prior
prior = Diagonal(1 .+ dp)
   
Ap = A / prior

results_p =ACEfit.solve(solver_qr,Ap,Y)
coeffs_p = prior \ results_p["C"]

pot_p = JuLIP.MLIPs.combine(rpib,coeffs_p)
ACE1pack.linear_errors(data,pot_p) 

#=
[ Info: RMSE Table
┌────────────┬───────────┬──────────┬─────────┐
│       Type │   E [meV] │ F [eV/A] │ V [meV] │
├────────────┼───────────┼──────────┼─────────┤
│   FLD_TiAl │ 10560.882 │    3.196 │   0.000 │
│ TiAl_T5000 │   409.910 │    8.956 │   0.000 │
├────────────┼───────────┼──────────┼─────────┤
│        set │ 10318.440 │    7.745 │   0.000 │
└────────────┴───────────┴──────────┴─────────┘
[ Info: MAE Table
┌────────────┬──────────┬──────────┬─────────┐
│       Type │  E [meV] │ F [eV/A] │ V [meV] │
├────────────┼──────────┼──────────┼─────────┤
│   FLD_TiAl │ 6178.498 │    1.862 │   0.000 │
│ TiAl_T5000 │  318.111 │    6.924 │   0.000 │
├────────────┼──────────┼──────────┼─────────┤
│        set │ 5912.117 │    5.461 │   0.000 │
└────────────┴──────────┴──────────┴─────────┘
=#