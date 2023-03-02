# Molly exploration with Pluto

For context, I'm simultaneously learning Molly, Pluto, the plotting ecosystem, and even some Julia basics. 

## 3-1-23

### Julia-focused

- Had to figure out how to "slice" a vector of StaticArrays. Naively did it the way I would do in python, didn't work. Have to use either getindex.() or array comprehension 

### Pluto issues
- The manifest resolution warnings that sometimes pop up (in the terminal) when using Pluto's package manager are a bit annoying. I think it might be specific to Julia 1.8. Can't quite tell what triggers it, because it didn't happen every time I loaded a package. Also a fresh Julia REPL helped. 

- My plots were originally not interactive, but then became interactive. Maybe because I called plotlyJs() at one point (then deleted, since I didn't have the relevant PlotlyJs package. Also was messing around with (already loaded) GLMakie commands (couldn't get the scatter to work because of a unit issue, even though it worked fine with Plots.plot3d), so maybe that changed things?

- I really expected that when i called simulate!() again, that the plots would change immediately, Perhaps they didn't because I wrapped key definitions in begin...end statements. But then even the length(sys.loggers.coords.history) didn't update. Not detecting the topology correctly? 

### Plot issues
- Aspect ratio for 3D scatter is weird, z-axis is much more spaced out even though axes limits are the same. Hopefully there's a better solution than what is provided in this 2021 [discourse discussion](https://discourse.julialang.org/t/plots-jl-aspect-ratio-1-0-for-3d-plot/58607)
