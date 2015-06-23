The code simulates heat diffusion in a 3-dimensional world using
monte-carlo statistics. 

Main.c:

First, the user must specify how many iterations the simulation 
should be run for. This is the variable TIME, and for each time
step, the walkers are allowed to move.

Then, the Rates and Times arrays specify respectively the rate
at which we inject walkers into the system and the times at which 
we measure the walker distribution. 

We inject the walkers from the centers of the n_bin_side^2 bins
on the boundaries y=0 and y=1.

The entries in Rates are integers. For a rate x, if x

= 0 : All walkers are injected at time 0.
> 0 : x walkers are injected from each face every time step.
< 0 : 1 walker is injected from each face every x time steps.

We enforce the condition that Rates is an ordered list.

The world has a volume of 1, and each side is broken up into
n_bins_side intervals. The cubic world is also split up into
n_bins_side^3 cubic bins with a side length of 1/n_bins_side.

The processes are alternated by hot and cold walkers.