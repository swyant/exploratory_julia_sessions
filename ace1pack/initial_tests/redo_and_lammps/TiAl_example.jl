using ACE1pack 

data_file = joinpath(ACE1pack.artifact("TiAl_tutorial"), "TiAl_tutorial.xyz")

data = read_extxyz(data_file)

model = acemodel(elements = [:Ti, :Al],
                order = 3,
                totaldegree = 6,
                rcut = 5.5,
                Eref = [:Ti => -1586.0195, :Al => -105.5954])

@show length(model.basis)

weights = Dict(
        "FLD_TiAl" => Dict("E" => 60.0, "F" => 1.0 , "V" => 1.0 ),
        "TiAl_T5000" => Dict("E" => 5.0, "F" => 1.0 , "V" => 1.0 ));

solver = ACEfit.LSQR(damp = 1e-2, atol = 1e-6);

P = smoothness_prior(model; p = 3)

data_train = data[1:5:end] # taking every 5 

#acefit!(model, data_train; solver=solver, prior = P);
acefit!(model, data_train; solver=solver, prior = P, weights=weights); #unsure if model "resets", so I re-ran the model command

@info("Training Error Table")
ACE1pack.linear_errors(data_train, model; weights=weights);

@info("Test Error Table")
test_data = data[2:10:end]
ACE1pack.linear_errors(test_data, model; weights=weights);

#export2lammps("./TiAl_tutorial_pot.yace", model)

@show data[end-1]
@show forces(model.potential,data[end-1])