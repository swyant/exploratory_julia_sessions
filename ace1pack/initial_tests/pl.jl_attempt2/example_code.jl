ace = ACE(species           = [:Hf, :O],
          body_order        = 3, 
          polynomial_degree = 6,
          r0                = 2.9,
          rcutoff           = 5.5 )

lb = LBasisPotential(ace)
lb, Σ = learn!(lb, ds_train; α=1e-6)

export2lammps("example_HfO2.yace", lb)


