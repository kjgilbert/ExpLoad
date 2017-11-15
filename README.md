# Expansion Load Simulation Code

C++ code modified from [Peischl & Excoffier 2015](http://onlinelibrary.wiley.com/doi/10.1111/mec.13154/abstract). If you use this code, please cite [Gilbert *et al. In prep.*]().

This code simulates either a range expansion or a range shift over space, whereby a total of 1000 loci of deleterious or beneficial effects mutate to impact fitness. The code can be compiled without any dependencies. Simulations are run in the command line by invoking `./expload input_file.txt`, where input_file.txt is a plain text file containing the input parameters for a simulation.

Simulates either a 1-dimensional or 2-dimensional landscape of discrete demes in a stepping stone or lattice model, respectively. Landscape is of size `m1` by `m2` demes. An initial number of demes (`starting_demes`) is colonized with individuals from an ancestral, panmictic population of size `anc_pop_size`. In two dimensions, `starting_demes` should only be the number of demes along the second landscape dimension (`m2`, usually the longest axis), and the code automatically fills all demes in that width (first landscape dimension, `m1`). For a range shift, the width of the metapopulation must be defined by `niche_width`, which again in two dimensions automatically uses all demes along the width of this defined length.

A range shift, unlike an expansion, has extinction at the trailing edge of the moving metapopulation, and the rate of change for movement of the metapopulation over space is defined by parameter `theta`, which is the number of generations between each instance of the metapopulation moving forward by 1 deme and simultaneously receding by 1 deme. Additionally, an open-front shift is possible to simulate where only the trailing edge has this forced speed (see further parameter options below). This extinction/expansion in shifts is controlled by modifying the carrying capacity of a deme to either 0 (extinction) or *K*, the carrying capacity defined by the user with the capacity parameter.

All parameters are as follows, and must be provided in this exact order with no other preceding text in the input file. All text after the last parameter is ignored by the program.

* `m1` = width of the landscape (set to 1 for a linear stepping stone model)
* `m2` = length of the landscape
* `starting_demes` = the number of demes to be colonized from the ancestral population
* `niche_width` = the width of the metapopulation during a range shift, ignored during a standard expansion
* `capacity` = carrying capacity per deme
* `anc_pop_size` = the size of the ancestral population
* `burnin_time` = the number of generations to burn in the ancestral population
* `expansion_start` = the number of generations to burn in the initial colonized populations (starting_demes) on the landscape before the expansion/shift begins
* `theta` = the number of generations between each successive shift of the metapopulation forward during a range shift
* `generations` = the total number of generations to run the simulation from when the expansion/shift begins; if the landscape is fully crossed, the simulation will continue and shifting populations will exist in a niche_width by m1 sized grid of demes for the duration of the siimulation
* `snapshot` = output simulation results at every generation that is a multiple of this number
* `replicates` = the number of simulation replicates to perform
* `expansionMode` = set equal to 0 for a linear expansion across the length of the landscape, all other modes are deprecated in this version of the code
* `expansionModeKim` = set to 0 for a standard full expansion, 1 for an open-front shift, and 2 for a shift controlled at both the front and trailing edges
* `selectionMode` = 0 for soft selection, 1 for hard selection
* `mu` = genome-wide mutation rate
* `m` = migration rate
* `s` = mean effect size of deleterious mutations (effect sizes are fixed, but code can be modified to an exponential distribution with this set as mean)
* `phi` = the proportion of loci which are unconditionally deleterious mutations, all others are then given the opposite s as set above, i.e. a mean of -0.1 makes all deleterious mutations have s = -0.1 and all beneficial mutations s = +0.1

More details or description of the function of the code are found in [Gilbert *et al. In prep.*]().