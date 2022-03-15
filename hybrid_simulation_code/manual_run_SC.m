%This script allows you to manually run the full metabolism model

clear;clc
par.N = 5; %Number of cells
par.gamma_tot = 100; %Total enzyme production rate
par.P = 1; %Permeability
par.c = 100; %Corresponds to kappa in manuscript
par.max_t = 2500; %Simulation time 
par.beta = 20; %Burst size
par.delta = 1; %Degradation rate

par.feedback = 1; %Whether growth feeds back onto enzyme dynamics
par.epsi = 1e-5; %Corresponds to g* in manuscript

par.maxchange = 0.03; %Maximum relative change in one ODE step
par.V = 10; %Corresponds to r_V in manuscript
par.overlay = 0; %Option to overlay m and E in plot 
par.n_store = 1e5; %Size of storage vector

%Set up auxotrophies
%We don't use this feature in the manuscript, but the code can support
%auxotrophies
par.aux_type1 = 0;
par.aux_type2 = 0;

%Set seed
rng(555)

%Simulate
sim_obj = hybrid_simulation_master(par);

%Plot
plot_hybrid_SC(sim_obj,par)

