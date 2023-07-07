[:] So two concerning things about this tutorials

- The weights are explicitly defined, but then they aren't passed to acefit!(), so that function ends up just using defaults weights? However, they are passed in at
ACE1pack.linear_errors(data_train, model; weights=weights);

- when setting up the smoothness prior, you have to explicitly set p (i.e., there is no autodetection based off of your model). seems potentially buggy
  P = smoothness_prior(model; p = 3)

  (Note that p is the polynomial size, and I think the point of the smoothness prior is to help dampen the higher order/frequency polynomials?


[:] running this fit
acefit!(model, data_train; solver=solver, prior = P);

┌────────────┬──────────┬───────┬────┬──────┬─────┐
│       Type │ #Configs │ #Envs │ #E │   #F │  #V │
├────────────┼──────────┼───────┼────┼──────┼─────┤
│   FLD_TiAl │       63 │   126 │ 63 │  378 │ 378 │
│ TiAl_T5000 │        3 │   310 │  3 │  930 │  18 │
├────────────┼──────────┼───────┼────┼──────┼─────┤
│      total │       66 │   436 │ 66 │ 1308 │ 396 │
│    missing │        0 │     0 │  0 │    0 │   0 │
└────────────┴──────────┴───────┴────┴──────┴─────┘
[ Info: Assembling linear problem.
[ Info:   - Creating feature matrix with size (1770, 190).
[ Info:   - Beginning assembly with processor count:  1.
Progress: 100%|███████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████| Time: 0:00:07
[ Info:   - Assembly completed.
┌ Warning: Need to apply preconditioner in LSQR.
└ @ ACEfit ~/.julia/packages/ACEfit/ID48n/src/solvers.jl:93
damp  0.01
atol  1.0e-6
maxiter  100000
Converged after 470 iterations.
relative RMS error  0.003664538640680456

~~~>
@info("Training Error Table")
ACE1pack.linear_errors(data_train, model; weights=weights);
[ Info: Training Error Table
[ Info: RMSE Table
┌────────────┬─────────┬──────────┬─────────┐
│       Type │ E [meV] │ F [eV/A] │ V [meV] │
├────────────┼─────────┼──────────┼─────────┤
│   FLD_TiAl │   9.765 │    0.079 │  97.640 │
│ TiAl_T5000 │   1.022 │    0.375 │  35.953 │
├────────────┼─────────┼──────────┼─────────┤
│        set │   9.543 │    0.319 │  95.703 │
└────────────┴─────────┴──────────┴─────────┘
[ Info: MAE Table
┌────────────┬─────────┬──────────┬─────────┐
│       Type │ E [meV] │ F [eV/A] │ V [meV] │
├────────────┼─────────┼──────────┼─────────┤
│   FLD_TiAl │   6.689 │    0.050 │  60.831 │
│ TiAl_T5000 │   1.008 │    0.291 │  25.885 │
├────────────┼─────────┼──────────┼─────────┤
│        set │   6.431 │    0.221 │  59.242 │
└────────────┴─────────┴──────────┴─────────┘


~~~> changing the evaluation to using default weights doesn't change anything
[ Info: RMSE Table
┌────────────┬─────────┬──────────┬─────────┐
│       Type │ E [meV] │ F [eV/A] │ V [meV] │
├────────────┼─────────┼──────────┼─────────┤
│   FLD_TiAl │   9.765 │    0.079 │  97.640 │
│ TiAl_T5000 │   1.022 │    0.375 │  35.953 │
├────────────┼─────────┼──────────┼─────────┤
│        set │   9.543 │    0.319 │  95.703 │
└────────────┴─────────┴──────────┴─────────┘
[ Info: MAE Table
┌────────────┬─────────┬──────────┬─────────┐
│       Type │ E [meV] │ F [eV/A] │ V [meV] │
├────────────┼─────────┼──────────┼─────────┤
│   FLD_TiAl │   6.689 │    0.050 │  60.831 │
│ TiAl_T5000 │   1.008 │    0.291 │  25.885 │
├────────────┼─────────┼──────────┼─────────┤
│        set │   6.431 │    0.221 │  59.242 │



---> vs. explicitly setting the weights when fitting the model, after resetting the models (by re-instantiating it from scratch)
acefit!(model, data_train; solver=solver, prior = P, weights=weights);
┌────────────┬──────────┬───────┬────┬──────┬─────┐
│       Type │ #Configs │ #Envs │ #E │   #F │  #V │
├────────────┼──────────┼───────┼────┼──────┼─────┤
│   FLD_TiAl │       63 │   126 │ 63 │  378 │ 378 │
│ TiAl_T5000 │        3 │   310 │  3 │  930 │  18 │
├────────────┼──────────┼───────┼────┼──────┼─────┤
│      total │       66 │   436 │ 66 │ 1308 │ 396 │
│    missing │        0 │     0 │  0 │    0 │   0 │
└────────────┴──────────┴───────┴────┴──────┴─────┘
[ Info: Assembling linear problem.
[ Info:   - Creating feature matrix with size (1770, 190).
[ Info:   - Beginning assembly with processor count:  1.
Progress: 100%|███████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████| Time: 0:00:00
[ Info:   - Assembly completed.
┌ Warning: Need to apply preconditioner in LSQR.
└ @ ACEfit ~/.julia/packages/ACEfit/ID48n/src/solvers.jl:93
damp  0.01
atol  1.0e-6
maxiter  100000
Converged after 438 iterations.
relative RMS error  0.0031575666515624534

ACE1pack.linear_errors(data_train, model; weights=weights);
[ Info: Training Error Table
[ Info: RMSE Table
┌────────────┬─────────┬──────────┬─────────┐
│       Type │ E [meV] │ F [eV/A] │ V [meV] │
├────────────┼─────────┼──────────┼─────────┤
│   FLD_TiAl │   6.824 │    0.062 │  88.961 │
│ TiAl_T5000 │  15.453 │    0.388 │  34.202 │
├────────────┼─────────┼──────────┼─────────┤
│        set │   7.437 │    0.329 │  87.221 │
└────────────┴─────────┴──────────┴─────────┘
[ Info: MAE Table
┌────────────┬─────────┬──────────┬─────────┐
│       Type │ E [meV] │ F [eV/A] │ V [meV] │
├────────────┼─────────┼──────────┼─────────┤
│   FLD_TiAl │   4.704 │    0.039 │  55.433 │
│ TiAl_T5000 │  14.712 │    0.302 │  26.172 │
├────────────┼─────────┼──────────┼─────────┤
│        set │   5.159 │    0.226 │  54.103 │
└────────────┴─────────┴──────────┴─────────┘


~~~> once again changing this evaluation to use default weights has no effect. Not sure what the role of the weights are here

ACE1pack.linear_errors(data_train, model; weights=ACE1pack.default_weights());
┌────────────┬─────────┬──────────┬─────────┐
│       Type │ E [meV] │ F [eV/A] │ V [meV] │
├────────────┼─────────┼──────────┼─────────┤
│   FLD_TiAl │   6.824 │    0.062 │  88.961 │
│ TiAl_T5000 │  15.453 │    0.388 │  34.202 │
├────────────┼─────────┼──────────┼─────────┤
│        set │   7.437 │    0.329 │  87.221 │
└────────────┴─────────┴──────────┴─────────┘
[ Info: MAE Table
┌────────────┬─────────┬──────────┬─────────┐
│       Type │ E [meV] │ F [eV/A] │ V [meV] │
├────────────┼─────────┼──────────┼─────────┤
│   FLD_TiAl │   4.704 │    0.039 │  55.433 │
│ TiAl_T5000 │  14.712 │    0.302 │  26.172 │
├────────────┼─────────┼──────────┼─────────┤
│        set │   5.159 │    0.226 │  54.103 │
└────────────┴─────────┴──────────┴─────────┘


~~~> also these fits aren't deterministic? Unless I'm accidentally changing somethign. But otherwise they seem to obtain slightly different results...

each time, I reset the model by re-evaluating the model instantiation line
e.g.

Attempt 1
[ Info: RMSE Table
┌────────────┬─────────┬──────────┬─────────┐
│       Type │ E [meV] │ F [eV/A] │ V [meV] │
├────────────┼─────────┼──────────┼─────────┤
│   FLD_TiAl │   6.685 │    0.060 │  90.583 │
│ TiAl_T5000 │  16.066 │    0.384 │  33.924 │
├────────────┼─────────┼──────────┼─────────┤
│        set │   7.375 │    0.326 │  88.796 │
└────────────┴─────────┴──────────┴─────────┘
[ Info: MAE Table
┌────────────┬─────────┬──────────┬─────────┐
│       Type │ E [meV] │ F [eV/A] │ V [meV] │
├────────────┼─────────┼──────────┼─────────┤
│   FLD_TiAl │   4.626 │    0.038 │  56.117 │
│ TiAl_T5000 │  15.836 │    0.298 │  26.707 │
├────────────┼─────────┼──────────┼─────────┤
│        set │   5.136 │    0.223 │  54.780 │
└────────────┴─────────┴──────────┴─────────┘


Attempt 2
[ Info: RMSE Table
┌────────────┬─────────┬──────────┬─────────┐
│       Type │ E [meV] │ F [eV/A] │ V [meV] │
├────────────┼─────────┼──────────┼─────────┤
│   FLD_TiAl │   6.787 │    0.061 │  91.752 │
│ TiAl_T5000 │  15.092 │    0.389 │  32.178 │
├────────────┼─────────┼──────────┼─────────┤
│        set │   7.370 │    0.330 │  89.904 │
└────────────┴─────────┴──────────┴─────────┘
[ Info: MAE Table
┌────────────┬─────────┬──────────┬─────────┐
│       Type │ E [meV] │ F [eV/A] │ V [meV] │
├────────────┼─────────┼──────────┼─────────┤
│   FLD_TiAl │   4.696 │    0.038 │  57.149 │
│ TiAl_T5000 │  14.666 │    0.302 │  24.569 │
├────────────┼─────────┼──────────┼─────────┤
│        set │   5.150 │    0.226 │  55.668 │
└────────────┴─────────┴──────────┴─────────┘

Attempt 3
[ Info: RMSE Table
┌────────────┬─────────┬──────────┬─────────┐
│       Type │ E [meV] │ F [eV/A] │ V [meV] │
├────────────┼─────────┼──────────┼─────────┤
│   FLD_TiAl │   6.867 │    0.064 │  89.057 │
│ TiAl_T5000 │  14.118 │    0.390 │  32.789 │
├────────────┼─────────┼──────────┼─────────┤
│        set │   7.354 │    0.330 │  87.289 │
└────────────┴─────────┴──────────┴─────────┘
[ Info: MAE Table
┌────────────┬─────────┬──────────┬─────────┐
│       Type │ E [meV] │ F [eV/A] │ V [meV] │
├────────────┼─────────┼──────────┼─────────┤
│   FLD_TiAl │   4.804 │    0.040 │  55.217 │
│ TiAl_T5000 │  13.786 │    0.303 │  24.711 │
├────────────┼─────────┼──────────┼─────────┤
│        set │   5.212 │    0.227 │  53.831 │
└────────────┴─────────┴──────────┴─────────┘


This seems weird because there's nothing particularly stochastic about this if I'm using LSQR right? the prior and the weights are fixed in each attempt, and I'm running each one on my mac (so same machine). Differences are small, but still a bit strange to me.



---> Anyways, final test fits (passing weights to acefit!(), so slightly different than example)
[ Info: Test Error Table
[ Info: RMSE Table
┌────────────┬─────────┬──────────┬─────────┐
│       Type │ E [meV] │ F [eV/A] │ V [meV] │
├────────────┼─────────┼──────────┼─────────┤
│   FLD_TiAl │  10.201 │    0.112 │  94.371 │
│ TiAl_T5000 │  13.567 │    0.428 │ 139.140 │
├────────────┼─────────┼──────────┼─────────┤
│        set │  10.319 │    0.301 │  96.034 │
└────────────┴─────────┴──────────┴─────────┘
[ Info: MAE Table
┌────────────┬─────────┬──────────┬─────────┐
│       Type │ E [meV] │ F [eV/A] │ V [meV] │
├────────────┼─────────┼──────────┼─────────┤
│   FLD_TiAl │   6.563 │    0.055 │  57.388 │
│ TiAl_T5000 │  13.567 │    0.341 │ 112.380 │
├────────────┼─────────┼──────────┼─────────┤
│        set │   6.775 │    0.186 │  59.054 │
└────────────┴─────────┴──────────┴─────────┘
