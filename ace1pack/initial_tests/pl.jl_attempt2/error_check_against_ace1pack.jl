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

#### Custom functions (NEED TO MERGE THIS ASAP)
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
###############################


ds_path = "/Users/swyant/.julia/artifacts/b437d7d5fac4424b8203d0afc31732879d3da5b2/TiAl_tutorial.xyz" 
ds = fixed_load_data(ds_path, ExtXYZ(u"eV", u"Å"))

ace = ACE(species           = [:Ti, :Al],
          body_order        = 3, 
          polynomial_degree = 6,
          wL                = 1.5,
          csp               = 1.0,
          r0                = 2.9,
          rcutoff           = 5.5 )

coeffs = [0.025591381519026762, 0.03527125788791906, 0.030612348271342734, 0.06646319273427727, -0.007326117584533097, 0.0021476223879284516, 0.007005564677788114, 0.007769244005502122, 0.0033019143419533176, -0.024266397752728517, -0.03511058549188825, -0.004486106026460792, 0.07649832144216986, -0.017295005584346018, 0.0348519410800185, -0.0026935081045876344, -0.0036996351104090115, 0.04250332320742131, -0.06611126598243479, 0.07452744399669442, -0.08807382022058645, 0.006553101218837121, -0.02825330387435087, -0.005070437887508557, 0.017488241826946662, -0.041461388491636234, -0.050152804966179194, 0.014554551186620662, 0.005494466857846328, 0.03395869840669037, -0.12004390275966798, 0.07758118243125994, 0.024624168020804672, 0.0006581992277555695, -0.002196641935532242, 0.03231551745953874, 0.0005431753297032715, 0.009602374511533056, 0.028907266845791348, 0.03557855646347803, 0.000832998634326787, 0.019238505326450918, -0.007863457928993406, 0.03497657242548427, -0.058485491203844206, -0.025527625067137013, -0.003851837725125408, 0.019472633328804008, -0.04975455754968226, 0.008243807089528446, 0.020612783411412677, -0.07411984524326856, 0.007005564677788615, 0.0077692440055020795, 0.003301914341951156, -0.024266397752728874, -0.03511058549188732, -0.004486106026461139, 0.004050167320520888, 0.011275083723136878, -0.009533633282696359, -0.016652089366136488, 0.005947187981081792, -0.0086386178798077, 0.027556838876613178, -0.01794755394550558, 0.0328518497817209, -0.0444008944069233, -0.04142464521909121, 0.014939466653767677, -0.0013061815492572404, -0.008399904141687925, -0.013070571180237286, 0.07022679858374972, 0.03655463426663164, -0.02425878114877371, -0.013089322405632224, 0.007663504514768707, -0.0006932536563853398, -0.015392489165582057, -0.005333834581033578, 0.000966860983206308, -0.06259246382383571, 0.04896321372972445, -0.012976299346956766, 0.00575471543255263, -0.010710826925328487, 0.009130893987440367, 0.025455356200891677, 0.03737467186835743, -0.061410072131176816, 0.05891873535070835, 0.02179281899408886, 0.02640823532251172, 0.0002904232787473565, -0.028649695323579742, -0.018788426151163752, -0.004911376520223526, 0.034242726688142995, -0.008960425717451335, -0.04627434272845332, 0.05321984323617818, 0.013802856612787245, 0.023560957961937884]

lb = LBasisPotential(coeffs,ace)

test_ds = DataSet([ds[2]])
#test_ds = DataSet([config for config in ds[1:2]])
compute_force_descriptors(get_system.(test_ds)[1],ace)
f_descr_test = compute_force_descriptors(test_ds, ace)

test_ds = DataSet(test_ds .+ f_descr_test)

test_forces_pred = get_all_forces(test_ds, lb)
#@show fforces = transpose(reshape(test_forces_pred,3,:))

f_ref = get_all_forces(test_ds)
f_mae, f_rmse, f_rsq = calc_metrics(test_forces_pred, f_ref)

mock_export2lammps("IBP_ACE_example_TiAl_2.yace", lb)