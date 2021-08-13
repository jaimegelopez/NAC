function sim_obj = construct_hybrid_simulation(par)

%This script constructs the hybrid simulation object

%Initialize time
sim_obj.t = 0;

%Make vector of all variables, order goes 
%m1_ext m2_ext m1_rho m2r_rho E1_rho E2_rho 
ss_c = ceil((0.5*par.c*par.gamma*par.beta)/(1 + par.delta + par.P));
ss_E = ceil(par.gamma*par.beta/2);
var_vec = [ss_c; ss_c; ss_E; ss_E];
if par.P == 0
    sim_obj.var = [0;0;repmat(var_vec,par.N,1)];
else
    sim_obj.var = [ss_c;ss_c;repmat(var_vec,par.N,1)];
end

%Set up variables
%Type of compound, 1-metabolite, 2-enzyme
var_id1_vec = [1;1;2;2];
sim_obj.var_id1 = [1;1;repmat(var_id1_vec,par.N,1)];

%Indices of location, 0-external, 1-N-cells
sim_obj.var_id2 = [0;0;repelem((1:par.N)',length(var_vec),1)];

%Indices of compounds, i indices of compounds
var_id3_vec = [1;2;1;2];
sim_obj.var_id3 = [1;2;repmat(var_id3_vec,par.N,1)];

sim_obj.m1_ind = (sim_obj.var_id1 == 1) & (sim_obj.var_id3 == 1) & (sim_obj.var_id2 > 0);
sim_obj.m2_ind = (sim_obj.var_id1 == 1) & (sim_obj.var_id3 == 2) & (sim_obj.var_id2 > 0);
sim_obj.E1_ind = (sim_obj.var_id1 == 2) & (sim_obj.var_id3 == 1) & (sim_obj.var_id2 > 0);
sim_obj.E2_ind = (sim_obj.var_id1 == 2) & (sim_obj.var_id3 == 2) & (sim_obj.var_id2 > 0);

%Indices of reactions
%cell index 1-N
sim_obj.reac_id1 = repelem((1:par.N)',2);
%reaction class index 1 - make
sim_obj.reac_id2 = repmat([1;1],par.N,1);
%type indices
sim_obj.reac_id3 = repmat([1;2],par.N,1);

%Reaction matrix
sim_obj.reac_matrix = zeros(length(sim_obj.reac_id1),length(sim_obj.var));
for i = 1:length(sim_obj.reac_id1)
    ind1 = sim_obj.var_id2 == sim_obj.reac_id1(i);
    ind2 = sim_obj.var_id3 == sim_obj.reac_id3(i);
    ind3 = sim_obj.E1_ind | sim_obj.E2_ind;

    sim_obj.reac_matrix(i,ind1&ind2&ind3) = 1;

end

sim_obj.make_E1_ind = (sim_obj.reac_id2 == 1) & (sim_obj.reac_id3 == 1);
sim_obj.make_E2_ind = (sim_obj.reac_id2 == 1) & (sim_obj.reac_id3 == 2);

sim_obj.record_var = zeros(length(sim_obj.var),par.n_store);
sim_obj.record_t = zeros(1,par.n_store);
sim_obj.burst_t = zeros(1,par.n_store);

end
