using ACE1pack 

data_file = joinpath(ACE1pack.artifact("TiAl_tutorial"), "TiAl_tutorial.xyz")

data = read_extxyz(data_file)

model = acemodel(elements = [:Ti, :Al],
                order = 3,
                totaldegree = 6,
                rcut = 5.5,
                Eref = [:Ti => -1586.0195, :Al => -105.5954]);

@show model.basis.BB

model2 = acemodel(elements = [:Ti, :Al],
                order = 2,
                totaldegree = 6,
                rcut = 5.5,
                Eref = [:Ti => -1586.0195, :Al => -105.5954]);

@show model2.basis.BB

model3 = acemodel(elements = [:Ti, :Al],
                order = 1,
                totaldegree = 6,
                rcut = 5.5,
                Eref = [:Ti => -1586.0195, :Al => -105.5954]);

model4 = acemodel(elements = [:Ti, :Al],
                order = 2,
                totaldegree = 3,
                rcut = 5.5,
                Eref = [:Ti => -1586.0195, :Al => -105.5954]);



@show model4.basis.BB