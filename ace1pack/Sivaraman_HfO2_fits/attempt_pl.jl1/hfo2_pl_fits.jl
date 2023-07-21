using PotentialLearning
using InteratomicPotentials, InteratomicBasisPotentials
using Unitful, UnitfulAtomic
using AtomsBase
using StaticArrays
using LinearAlgebra: Diagonal
using CSV, DataFrames

ds_path = "../../../../../datasets/HfO2_Sivaraman/prl_2021/raw/train.xyz"
ds = fixed_load_data(ds_path, ExtXYZ(u"eV", u"Å"))

train_idxs = vec(Matrix(CSV.read("Siviraman_HfO2_my_train_idxs.csv",DataFrame,header=false)))
val_idxs   = vec(Matrix(CSV.read("Siviraman_HfO2_my_val_idxs.csv",DataFrame,header=false)))

data_train = ds[train_idxs]
data_val   = ds[val_idxs]

ace1 = ACE(species           = [:Hf, :O],
          body_order        = 4, 
          polynomial_degree = 10,
          wL                = 1.5,
          csp               = 1.0,
          r0                = 2.15,
          rcutoff           = 5.0 )

e_descr_train = compute_local_descriptors(data_train,ace1)
f_descr_train = compute_force_descriptors(data_train,ace1)
ds_train = DataSet(data_train .+ e_descr_train .+ f_descr_train)

lb_noE_wls = LBasisPotential(ace1)
learn_wls!(lb_noE_wls, ds_train; w_e=0.0)  # just forces, because no vref

f_train = get_all_forces(ds_train)
f_train_pred_noE = get_all_forces(ds_train,lb_noE_wls)
f_train_mae, f_train_rmse, f_train_rsq = calc_metrics(f_train_pred_noE, f_train)
#(0.15924797855448045, 0.21769682978711252, 0.984147636963468)

f_val = get_all_forces(data_val)
f_descr_val = compute_force_descriptors(data_val, ace1)
ds_val = DataSet(data_val .+ f_descr_val)
f_val_pred_noE = get_all_forces(ds_val, lb_noE_wls)
f_val_mae, f_val_rmse, f_val_rsq = calc_metrics(f_val_pred_noE, f_val)
# (0.16070235630031215, 0.2213594369588226, 0.9830181728879664)

# So I was going to try DPP, but I think it's only really set up for energies atm?
lb_dpp_500 = LBasisPotnetial(ace1)
dpp = kDPP(ds_train, GlobalMean(), DotProduct(); batch_size = 1000)
dpp = kDPP(ds_train, GlobalMean(), DotProduct(); batch_size = 500)

#### Custom functions 
function fixed_load_data(file, extxyz::ExtXYZ; T=Float64)
    configs = Configuration[]
    open(file, "r") do io
        count = 1
        while !eof(io)
            # Read info line
            line = readline(io)
            num_atoms = parse(Int, line)

            line = readline(io)
            lattice_line = match(r"Lattice=\"(.*?)\" ", line).captures[1]
            lattice = parse.(Float64, split(lattice_line)) * extxyz.distance_units
            box = SVector{3}.([lattice[1:3], lattice[4:6], lattice[7:9]]) # MODIFIED: added SVector here
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
                bc = [t == "T" ? AtomsBase.Periodic() : DirichletZero() for t in split(bc_line)]
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
                    elseif prop[1] == "forces" || prop[1] == "force" # MODIFIED: the property is force, not forces
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
                atoms[i] = AtomsBase.Atom(element,position .* extxyz.distance_units) # MODIFIED: not including data 
            end

            system = FlexibleSystem(atoms, box, bc)
            count += 1
            push!(configs, Configuration(system, energy, Forces(forces)))
        end
    end
    return DataSet(configs)
end


function learn_wls!(lb::LBasisPotential, ds::DataSet; w_e = 1.0, w_f = 1.0, intercept = false)
    lp = PotentialLearning.LinearProblem(ds)

    @views B_train = reduce(hcat, lp.B)'
    @views dB_train = reduce(hcat, lp.dB)'
    @views e_train = lp.e
    @views f_train = reduce(vcat, lp.f)

    # Calculate A and b.
    if intercept
        int_col = ones(size(B_train, 1)+size(dB_train, 1))
        @views A = hcat(int_col, [B_train; dB_train])
    else
        @views A = [B_train; dB_train]
    end
    @views b = [e_train; f_train]

    # Calculate coefficients β.
    Q = Diagonal([w_e * ones(length(e_train));
                  w_f * ones(length(f_train))])
    βs = (A'*Q*A) \ (A'*Q*b)

    if intercept
        copyto!(lb.β0,βs[1])
        copyto!(lb.β,βs[2:end])
    else
        copyto!(lb.β,βs)
    end
end
####