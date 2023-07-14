using ACE1 
using ACE1pack
using ACEfit
using JuLIP

ds_path = "/Users/swyant/.julia/artifacts/b437d7d5fac4424b8203d0afc31732879d3da5b2/TiAl_tutorial.xyz" 
raw_data = read_extxyz(ds_path)

rpib = ACE1.rpi_basis(;
            species = [:Ti, :Al],
            N       = 2, 
            maxdeg  = 6,
            D       = ACE1.SparsePSHDegree(; wL = 1.5, csp = 1.0),
            r0      = 2.9,
            rin     = 0.65*2.9,
            rcut    = 5.5,
            pin     = 0,
       )
manual_coeffs = [0.025591381519026762, 0.03527125788791906, 0.030612348271342734, 0.06646319273427727, -0.007326117584533097, 0.0021476223879284516, 0.007005564677788114, 0.007769244005502122, 0.0033019143419533176, -0.024266397752728517, -0.03511058549188825, -0.004486106026460792, 0.07649832144216986, -0.017295005584346018, 0.0348519410800185, -0.0026935081045876344, -0.0036996351104090115, 0.04250332320742131, -0.06611126598243479, 0.07452744399669442, -0.08807382022058645, 0.006553101218837121, -0.02825330387435087, -0.005070437887508557, 0.017488241826946662, -0.041461388491636234, -0.050152804966179194, 0.014554551186620662, 0.005494466857846328, 0.03395869840669037, -0.12004390275966798, 0.07758118243125994, 0.024624168020804672, 0.0006581992277555695, -0.002196641935532242, 0.03231551745953874, 0.0005431753297032715, 0.009602374511533056, 0.028907266845791348, 0.03557855646347803, 0.000832998634326787, 0.019238505326450918, -0.007863457928993406, 0.03497657242548427, -0.058485491203844206, -0.025527625067137013, -0.003851837725125408, 0.019472633328804008, -0.04975455754968226, 0.008243807089528446, 0.020612783411412677, -0.07411984524326856, 0.007005564677788615, 0.0077692440055020795, 0.003301914341951156, -0.024266397752728874, -0.03511058549188732, -0.004486106026461139, 0.004050167320520888, 0.011275083723136878, -0.009533633282696359, -0.016652089366136488, 0.005947187981081792, -0.0086386178798077, 0.027556838876613178, -0.01794755394550558, 0.0328518497817209, -0.0444008944069233, -0.04142464521909121, 0.014939466653767677, -0.0013061815492572404, -0.008399904141687925, -0.013070571180237286, 0.07022679858374972, 0.03655463426663164, -0.02425878114877371, -0.013089322405632224, 0.007663504514768707, -0.0006932536563853398, -0.015392489165582057, -0.005333834581033578, 0.000966860983206308, -0.06259246382383571, 0.04896321372972445, -0.012976299346956766, 0.00575471543255263, -0.010710826925328487, 0.009130893987440367, 0.025455356200891677, 0.03737467186835743, -0.061410072131176816, 0.05891873535070835, 0.02179281899408886, 0.02640823532251172, 0.0002904232787473565, -0.028649695323579742, -0.018788426151163752, -0.004911376520223526, 0.034242726688142995, -0.008960425717451335, -0.04627434272845332, 0.05321984323617818, 0.013802856612787245, 0.023560957961937884] 
pot_manual = JuLIP.MLIPs.combine(rpib, manual_coeffs)

weights = Dict("default" => Dict("E" =>0.0,"F" => 1.0, "V"=>0.0))

#check_data = [ AtomsData(at;  energy_key = "energy", force_key="force",
#                weights = weights) for at in raw_data[1:2]]

check_data = [ AtomsData(at;  energy_key = "energy", force_key="force",
                weights = weights) for at in [raw_data[end]]]
forces(pot_manual, raw_data[20])
ACE1pack.linear_errors(check_data, pot_manual)

