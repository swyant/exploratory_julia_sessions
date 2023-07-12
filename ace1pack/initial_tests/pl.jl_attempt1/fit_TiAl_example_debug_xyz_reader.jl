using PotentialLearning
using InteratomicPotentials, InteratomicBasisPotentials
using Unitful, UnitfulAtomic
using AtomsBase 
using StaticArrays


ds_path = "/Users/swyant/.julia/artifacts/b437d7d5fac4424b8203d0afc31732879d3da5b2/TiAl_tutorial.xyz" 

#ds = load_data(ds_path, ExtXYZ(u"eV", u"Å"))

extxyz =ExtXYZ(u"eV", u"Å") 
configs = Configuration[]
T = Float64
open(ds_path, "r") do io
    count = 1
    while !eof(io)
        # Read info line
        line = readline(io)
        num_atoms = parse(Int, line)

        line = readline(io)
        lattice_line = match(r"Lattice=\"(.*?)\" ", line).captures[1]
        lattice = parse.(Float64, split(lattice_line)) * extxyz.distance_units
        box = [lattice[1:3], lattice[4:6], lattice[7:9]]
        energy = try
            energy_line = match(r"energy=(.*?) ", line).captures[1]
            energy = parse(Float64, energy_line)
            Energy(energy, extxyz.energy_units)
        catch
            Energy(NaN, extxyz.energy_units)
        end

        # try
        #     stress_line = match(r"stress=\"(.*?)\" ", line).captures[1]
        #     stress = parse.(Float64, split(stress_line))
        #     push!(stresses, SVector{6}(stress))
        # catch
        #     push!(stresses, SVector{6}( fill(NaN, (6,))))
        # end

        bc = []
        try
            bc_line = match(r"pbc=\"(.*?)\"", line).captures[1]
            bc = [t == "T" ? Periodic() : DirichletZero() for t in split(bc_line)]
        catch
            bc = [DirichletZero(), DirichletZero(), DirichletZero()]
        end

        properties = match(r"Properties=(.*?) ", line).captures[1]
        properties = split(properties, ":")
        properties = [properties[i:i+2] for i = 1:3:(length(properties)-1)] 
        atoms = Vector{AtomsBase.Atom}(undef, num_atoms)
        forces = Force[]
        for i = 1:num_atoms
            line = split(readline(io))
            line_count = 1
            position = 0.0
            element = 0.0
            data = Dict(())
            for prop in properties
                if prop[1] == "species"
                    element = Symbol(line[line_count])
                    line_count += 1
                elseif prop[1] == "pos"
                    position = SVector{3}(parse.(T, line[line_count:line_count+2]))
                    line_count += 3
                elseif prop[1] == "move_mask"
                    ft = Symbol(line[line_count])
                    line_count += 1
                elseif prop[1] == "tags"
                    ft = Symbol(line[line_count])
                    line_count += 1
#                elseif prop[1] == "forces" . # it's force, not forces in this xyz file
                elseif prop[1] == "forces" || prop[1] == "force"

                    @show check_parse = parse.(T, line[line_count:line_count+2])
                    @show check_force =Force(
                                        parse.(T, line[line_count:line_count+2]),
                                        extxyz.energy_units / extxyz.distance_units,
                                        ) 
                    push!(
                        forces,
                        Force(
                            parse.(T, line[line_count:line_count+2]),
                            extxyz.energy_units / extxyz.distance_units,
                        ),
                    )
                    line_count += 3
                else
                    length = parse(Int, prop[3])
                    if length == 1
                        data = merge(data, Dict((Symbol(prop[1]) => line[line_count])))
                    else
                        data = merge(
                            data,
                            Dict((
                                Symbol(prop[1]) =>
                                    line[line_count:line_count+length-1]
                            )),
                        )
                    end
                end
            end
            #@show data
            #@show position
            #@show box 
            #@show energy 
            #@show bc
            #@show properties
            #@show forces
            atoms[i] = AtomsBase.Atom(element,position .* extxyz.distance_units)
#            if isempty(data)
#                atoms[i] = AtomsBase.Atom(element, position .* extxyz.distance_units)
#            else
#                atoms[i] =
#                    AtomsBase.Atom(element, position .* extxyz.distance_units, data...) # This is definitely the issue
#            end
        end
    end
end

conf_train = ds[1:5:end]

ace = ACE(species           = [:Ti, :Al],
          body_order        = 3, 
          polynomial_degree = 6,
          wL                = 1.5,
          csp               = 1.0,
          r0                = 2.9,
          rcutoff           = 5.5 )

# just fitting to forces, so as to not trigger nonlinear MLE calculation
f_descr_train = compute_force_descriptors(conf_train, ace)

ds_train = DataSet(conf_train .+ f_descr_train)

lb = LBasisPotential(ace)
lb, Σ = learn!(lb, conf_train; α=1e-6)