function rates = hybrid_rates(sim_obj,par)

%This script computes the stochastic process rates for the full model

rates = zeros(size(sim_obj.reac_matrix,1),1);

m1 = sim_obj.var(sim_obj.m1_ind);
m2 = sim_obj.var(sim_obj.m2_ind);

con_rate = min([m1,m2],[],2);

E1_prod = m1 < m2;

%Randomize decision if nutrients are equal
equi_met = m2 == m1;
E1_prod(equi_met) = randi(2,sum(equi_met),1) -1;

%Set regulation to single value if auxotrophic
E1_prod(par.aux_var) = par.aux_id(par.aux_var);

if par.feedback
    make_E1 = con_rate.*E1_prod*par.epsi*par.gamma;
    make_E2 = con_rate.*(~E1_prod)*par.epsi*par.gamma;
else
    make_E1 = E1_prod*par.epsi*par.gamma;
    make_E2 = (~E1_prod)*par.epsi*par.gamma;
end

rates(sim_obj.make_E1_ind) = make_E1;
rates(sim_obj.make_E2_ind) = make_E2;

end

