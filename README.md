#Goal
The code simulates heat diffusion in a 3-dimensional world using
monte-carlo statistics. It analyzes carbon nanotube composites,
where the nanotubes are represented by cylindrical tubes.

I wrote this code during a physics REU internship at the University of
Oklahoma. My advisor was Dr. Kieran Mullen. 

##Introduction

The general layout of the program is as follows:

- Create a cubical world and set apart rectangular segments as
*Carbon Nanotubes*.
- Simulate a constant heat flux by periodically injecting 
hot and cold *thermal walkers* from the rightmost and leftmost
faces respectively.
- Let the walkers perform their random walk and at various times,
plot the walker population vs the direction of heat flux.
- Calculate the slope of these plots, dT/dy, and thus calculate
the thermal conductivity k.

Under these assumptions, the problem is almost trivially parallelizable.
We simply generate enough random walks, which is easy because the 
walks are independent; then we collect the positions of the walkers
at various offsets to simulate the heat flux.


###Main.c:
This file spawns the MPI processes that generate the random walks and
also collects the positions of the walkers. 

Each side of the world is split up into *n_bin_side* intervals, and these
split up the world into *n_bin_side^3* cubical bins. In our choice of axes,
the walkers are injected on the faces y=0 and y=1, and we plot the number
of walkers along the bins across the y-direction. 

The variable *Rate* specifies the rate at which walkers are injected into the
system. Given that there are *n_bin_side^2* bins on the faces, *Rate* specifies
how many walkers are injected into each bin on the faces y=0 and y=1 every time step.
Specifically, the walkers are spawned from the centers of the sides of the bins. 

*Rate* is represented as an integer. For *Rate*

- = 0 : All walkers are injected at time 0.
- > 0 : *Rate* walkers are injected from each bin every time step.
- < 0 : 1 walker is injected from each bin every |*Rate*| time steps.

We record walker numbers along the cross sectional bins at the times
in the array *Times*, and print these numbers onto csv files.

*walks* represents how many random walk simulations a process will run. 
So for *walks*=*n*, each process will simulate _n * -2 * n_bin_side^2_
random walks. 

We require that 
- walks * n_processes >= Rate*Times    if Rate > 0 
- walks * n_processes >= Times/|Rates| if Rate < 0

Processes with higher rank will generate the thermal walkers later in 
time. 



