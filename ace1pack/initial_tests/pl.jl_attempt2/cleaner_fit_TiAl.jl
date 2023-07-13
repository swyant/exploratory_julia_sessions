using PotentialLearning
using InteratomicPotentials, InteratomicBasisPotentials
using Unitful, UnitfulAtomic
using AtomsBase 
using StaticArrays
using JuLIP
using ACE1pack, ACE1x
using ACE1
using Interpolations
using OrderedCollections
using YAML

#### Custom functions #####
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

# TODO: update this to auto-detect elements from rpib
# TODO: auto-delete .table files since they're junk
function mock_export2lammps(fname, lb::LinearBasisPotential)
    @assert typeof(lb.basis) == ACE

    ### This is all so that I can generated a pair basis, which export2lammps needed
    ### should be refactored once it's actually added to our codebases
    mock_kwargs = Dict( :pair_basis     => :legendre,
                        :elements       => [:Ti,:Al],
                        :pair_degree    => 6,
                        :pair_rcut      => 5.5,
                        :pair_transform => (:agnesi,1,3),
                        :pair_envelope  => (:r,2),
                        :r0             => :bondlen,
                        :rcut           => 5.5)
    
    pairb = ACE1x._pair_basis(mock_kwargs)
    pairb_length = length(pairb)
    
    rpib  = lb.basis.rpib
    super_basis = JuLIP.MLIPs.IPSuperBasis([pairb,rpib])

    rpib_params = copy(lb.β)
    params = prepend!(rpib_params,randn(pairb_length))
    pot_novref = JuLIP.MLIPs.combine(super_basis,params)

    Eref = [:Ti => 0.0, :Al => 0.0]
    vref = JuLIP.OneBody(Eref...)

    pot = JuLIP.MLIPs.SumIP([pot_novref.components...,vref])

    mock_ace1model = ACE1x.ACE1Model(super_basis,params,vref,pot,Dict())
    export2lammps(fname, mock_ace1model.potential, rpib) 
    
end

function export2lammps(fname, IP,rpib::RPIBasis)

    if length(IP.components) != 3
        throw("IP must have three components which are OneBody, pair potential, and ace")
    end

    ordered_components = []

    for target_type in [OneBody, PolyPairPot, PIPotential]
        did_not_find = true
        for i = 1:3
            if typeof(IP.components[i]) <: target_type
                push!(ordered_components, IP.components[i])
                did_not_find = false
            end
        end

        if did_not_find
            throw("IP must have three components which are OneBody, pair potential, and ace")
        end
    end

    V1 = ordered_components[1]
    V2 = ordered_components[2]
    V3 = ordered_components[3]

    species = collect(string.(chemical_symbol.(V3.pibasis.zlist.list)))
    species_dict = Dict(zip(collect(0:length(species)-1), species))
    reversed_species_dict = Dict(zip(species, collect(0:length(species)-1)))


    elements = Vector(undef, length(species))
    E0 = zeros(length(elements))

    for (index, element) in species_dict
        E0[index+1] = V1(Symbol(element))
        elements[index+1] = element
    end

    # V1 and V3  (V2 handled below)

    # Begin assembling data structure for YAML
    data = OrderedDict()
    data["elements"] = elements

    data["E0"] = E0

    # embeddings
    data["embeddings"] = Dict()
    for species_ind1 in sort(collect(keys(species_dict)))
        data["embeddings"][species_ind1] = Dict(
            "ndensity" => 1,
            "FS_parameters" => [1.0, 1.0],
            "npoti" => "FinnisSinclairShiftedScaled",
            "drho_core_cutoff" => 1.000000000000000000,
            "rho_core_cutoff" => 100000.000000000000000000)
    end

    # bonds
    data["bonds"] = OrderedDict()
    basis1p = deepcopy(rpib.pibasis.basis1p)
    radialsplines = ACE1.Splines.RadialSplines(basis1p.J; zlist = basis1p.zlist, nnodes = 10000)
    ranges, nodalvals, zlist = ACE1.Splines.export_splines(radialsplines)
    # compute spline derivatives
    # TODO: move this elsewhere
    nodalderivs = similar(nodalvals)
    for iz1 in 1:size(nodalvals,2), iz2 in 1:size(nodalvals,3)
        for i in 1:size(nodalvals,1)
            range = ranges[i,iz1,iz2]
            spl = radialsplines.splines[i,iz1,iz2]
            deriv(r) = Interpolations.gradient(spl,r)[1]
            nodalderivs[i,iz1,iz2] = deriv.(range)
        end
    end
    # ----- end section to move
    for iz1 in 1:size(nodalvals,2), iz2 in 1:size(nodalvals,3)
        data["bonds"][[iz1-1,iz2-1]] = OrderedDict{Any,Any}(
            "radbasename" => "ACE.jl",
            "rcut" => ranges[1,iz1,iz2][end],         # note hardcoded 1
            "nradial" => length(V3.pibasis.basis1p.J.J.A),
            "nbins" => length(ranges[1,iz1,iz2])-1)   # note hardcoded 1
        nodalvals_map = OrderedDict([i-1 => nodalvals[i,iz1,iz2] for i in 1:size(nodalvals,1)])
        data["bonds"][[iz1-1,iz2-1]]["splinenodalvals"] = nodalvals_map
        nodalderivs_map = OrderedDict([i-1 => nodalderivs[i,iz1,iz2] for i in 1:size(nodalvals,1)])
        data["bonds"][[iz1-1,iz2-1]]["splinenodalderivs"] = nodalderivs_map
    end

    functions, lmax = ACE1pack.export_ACE_functions(V3, species, reversed_species_dict)
    data["functions"] = functions
    data["lmax"] = lmax
    YAML.write_file(fname, data)
end
#######

ds_path = "/Users/swyant/.julia/artifacts/b437d7d5fac4424b8203d0afc31732879d3da5b2/TiAl_tutorial.xyz" 

ds = fixed_load_data(ds_path, ExtXYZ(u"eV", u"Å"))
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
lb, Σ = learn!(lb, ds_train; α=1e-6)

mock_export2lammps("IBP_ACE_example_TiAl.yace", lb)

f_train = get_all_forces(ds_train)
f_train_pred = get_all_forces(ds_train,lb)
f_train_mae, f_train_rmse, f_train_rsq = calc_metrics(f_train_pred, f_train)

@show f_train_mae

test_ds = DataSet([ds[end-1]])
f_descr_test = compute_force_descriptors(test_ds, ace)
test_ds = DataSet(test_ds .+ f_descr_test)

test_forces = get_all_forces(test_ds, lb)
@show test_forces = transpose(reshape(test_forces,3,:))

check_potential = JuLIP.MLIPs.combine(lb.basis.rpib,lb.β)
data = JuLIP.read_extxyz(ds_path)
test_struct = data[end-1]
@show forces(check_potential,data[end-1])