clear;clc
par.N = 1;
par.gamma_tot = 100;
par.P = 1;
par.c = 100;
par.max_t = 2500;
par.beta = 2;
par.delta = 1;

par.feedback = 1;
par.epsi = 1e-5;
max_growth = par.epsi*par.gamma_tot*par.c/(1+par.delta);

par.maxchange = 0.03;
par.V = 10;
par.overlay = 0;
par.n_replicate = 5;
par.n_store = 1e5;
%Set up auxotrophies
par.aux_type1 = 0;
par.aux_type2 = 0;

rng(555)
sim_obj = hybrid_simulation_master(par);
sim_obj.growth_ints

plot_hybrid_SC(sim_obj,par)
pos = get(gcf,'Position');
pos(3) = 1.5*pos(3);
set(gcf,'Position',pos);
print(gcf, '-dpng','manual_plots/lab_meeting_example.png','-r400');

