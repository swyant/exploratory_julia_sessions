using ACE1 
using ACE1pack
using ACEfit
using LinearAlgebra: I, Diagonal, pinv
using JuLIP
using CSV, DataFrames
using Random
using Interpolations
using OrderedCollections
using YAML

ds_path = "../../../../../datasets/HfO2_Sivaraman/prl_2021/raw/train.xyz" 
raw_data = read_extxyz(ds_path)

train_idxs = vec(Matrix(CSV.read("../Siviraman_HfO2_my_train_idxs.csv",DataFrame,header=false)))
val_idxs   = vec(Matrix(CSV.read("../Siviraman_HfO2_my_val_idxs.csv",DataFrame,header=false)))


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

pair_basis_kwargs= Dict(:elements => [:Hf, :O],
                        :pair_rin => 0.65*2.15,
                        :pair_rcut => 5.0, 
                        :pair_degree => 10,
                        :pair_transform => (:agnesi, 1, 3),
                        :pair_basis => :legendre, 
                        :pair_envelope => (:r,2),
                        :r0 => 2.15,
                        :rcut => 5.0)

pairb = ACE1x._pair_basis(pair_basis_kwargs)

super_basis = JuLIP.MLIPs.IPSuperBasis([pairb,rpib3])

weights = Dict("default" => Dict("E" =>1.0,"F" => 1.0, "V"=>0.0))

vref = JuLIP.OneBody([:Hf => -2.70516846, :O => -0.01277342]...)
train_data = [ AtomsData(at;  energy_key = "energy", force_key="forces",
                        weights = weights, v_ref=vref) for at in raw_data[train_idxs]]

val_data = [ AtomsData(at;  energy_key = "energy", force_key="forces",
                        weights = weights, v_ref=vref) for at in raw_data[val_idxs]]

A3, Y3, W3 = ACEfit.assemble(train_data, super_basis)

solver_reg = ACEfit.QR(; lambda=1e-3)

results_reg3 = ACEfit.solve(solver_reg,A3,Y3)

@show results_reg3["C"]
CSV.write("N3_rcut5_maxdeg10_wpair_1e-3lambdaQR_fit_coeffs.csv", DataFrame(Tables.table(results_reg3["C"])),header=false)

pot_reg_tmp = JuLIP.MLIPs.combine(super_basis,results_reg3["C"])
pot_reg_check = JuLIP.MLIPs.SumIP(pot_reg_tmp,vref)
pot_reg = JuLIP.MLIPs.SumIP([pot_reg_tmp.components..., vref])
ACE1pack.linear_errors(train_data,pot_reg)
#=
[ Info: RMSE Table
┌─────────┬─────────┬──────────┬─────────┐
│    Type │ E [meV] │ F [eV/A] │ V [meV] │
├─────────┼─────────┼──────────┼─────────┤
│  amorph │   5.980 │    0.267 │   0.000 │
│ crystal │   4.924 │    0.194 │   0.000 │
├─────────┼─────────┼──────────┼─────────┤
│     set │   5.266 │    0.218 │   0.000 │
└─────────┴─────────┴──────────┴─────────┘
[ Info: MAE Table
┌─────────┬─────────┬──────────┬─────────┐
│    Type │ E [meV] │ F [eV/A] │ V [meV] │
├─────────┼─────────┼──────────┼─────────┤
│  amorph │   4.673 │    0.201 │   0.000 │
│ crystal │   3.859 │    0.142 │   0.000 │
├─────────┼─────────┼──────────┼─────────┤
│     set │   4.105 │    0.160 │   0.000 │
└─────────┴─────────┴──────────┴─────────┘
=#

@show forces(pot_reg,raw_data[395])
#=
108-element Vector{StaticArraysCore.SVector{3, Float64}}:
 [0.2969494827866761, 2.255410117613309, 3.577896436861522]
 [2.49443970460365, 0.31173802816798, -1.484831557535756]
 [0.98481132922835, 1.2632024631226813, -0.09026280094151673]
 [2.2786812602820703, 1.820217313009934, 2.9836212008800658]
 [-0.5284930874160665, 4.863032238242113, -1.7131618087960163]
 [0.24701731065108845, 1.5310641006049082, 0.6889844843288166]
 [-0.7892022433418191, -2.289563992545027, 0.3722118639678911]
 [-1.761117064258201, -0.13086628369248388, -0.554853493199273]
 ⋮
 [0.14830759652855363, 0.15143218011843373, -2.1417340397560167]
 [0.5171113756600221, 2.575954389908887, 2.766449926559403]
 [0.6205930248282812, 4.53788524223782, -3.7731398098862883]
 [2.023276844306821, -0.7078815672608698, -0.7340974880484623]
 [-2.466843787627681, 1.2670295666032985, 2.464812843720864]
 [-0.5421809682362984, 0.8508643035482351, -1.4293002993509205]
 [-0.1896024135990615, -1.772392177098567, -0.7859575681054025]
 [4.820024598136101, 0.38417497286799573, -1.3519309026613655]
=#

custom2_export2lammps("Siviraman_GAP_HfO2_N3_rcut5_maxdeg10_regQR_withpair.yace", pot_reg,rpib3)

function custom2_export2lammps(fname, IP,rpib::RPIBasis)

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

    # ----- 2body handled separately -----
    # writes a .table file, so for simplicity require that export fname is passed with
    # .yace extension, and we remove this and add the .table extension instead
    fname_stem = fname[1:end-5]
    write_pairpot_table(fname_stem, V2, species_dict)    
end

make_dimer(s1, s2, rr) = JuLIP.Atoms(
    [[0.0,0.0,0.0],[rr,0.0,0.0]],
    [[0.0,0.0,0.0],[0.0,0.0,0.0]],
    [atomic_mass(s1), atomic_mass(s2)],
    [AtomicNumber(s1), AtomicNumber(s2)],
    [100.0,100.0,100.0],
    [false, false, false])

function write_pairpot_table(fname, V2, species_dict)
    # fname is JUST THE STEM
    # write a pair_style table file for LAMMPS
    # the file has a seperate section for each species pair interaction
    # format of table pair_style is described at https://docs.lammps.org/pair_table.html

    # Create filename. Only the stem is specified
    fname = fname * "_pairpot.table"

    # enumerate sections
    species_pairs = []
    for i in 0:length(species_dict) - 1
        for j in i:length(species_dict) - 1
            push!(species_pairs, (species_dict[i], species_dict[j]))
        end
    end

    lines = Vector{String}()

    # make header. date is none since ACE1 current doesnt depend on time/dates package
    push!(lines, "# DATE: none UNITS: metal CONTRIBUTOR: ACE1.jl - https://github.com/ACEsuit/ACE1.jl")
    push!(lines, "# ACE1 pair potential")
    push!(lines, "")

    for spec_pair in species_pairs
        # make dimer
        dimer = make_dimer(Symbol(spec_pair[1]), Symbol(spec_pair[2]), 1.0)

        # get inner and outer cutoffs

        if typeof(V2.basis.J) <: JuLIP.SMatrix
            get_ru(jj) = jj.ru
            rus = get_ru.(V2.basis.J)
            rout = maximum(rus)
        else
            rout = V2.basis.J.ru
        end

        rin = 0.001
        spacing = 0.0005 # SPENCER EDIT
        rs = rin:spacing:rout

        # section header
        push!(lines, string(spec_pair[1], "_", spec_pair[2]))
        push!(lines, string("N ", length(rs)))
        push!(lines, "")

        # values
        for (index, R) in enumerate(rs)
            set_positions!(dimer, AbstractVector{JVec{Float64}}([[R,0.0,0.0], [0.0,0.0,0.0]]))
            E = energy(V2, dimer)
            F = forces(V2, dimer)[1][1]
            push!(lines, string(index, " ", R, " ", E, " ", F))
        end
        push!(lines, "")
    end

    # write
    open(fname, "w+") do io
        for line in lines
            write(io, line * "\n")
        end
    end

    return nothing
end        