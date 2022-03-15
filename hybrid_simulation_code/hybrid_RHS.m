function RHS = hybrid_RHS(var,sim_obj,par)

%This script computes the ODE RHS for the full metabolism model

RHS = zeros(size(var));

m1_ext = var(1);
m2_ext = var(2);

m1 = var(sim_obj.m1_ind);
m2 = var(sim_obj.m2_ind);

E1 = var(sim_obj.E1_ind);
E2 = var(sim_obj.E2_ind);

con_rate = min([m1,m2],[],2);

%If V = 0 simulate direct transfer between cells
if par.V == 0
    dm1 = E1*par.c - con_rate + par.P*(sum(m1) - par.N*m1) - par.delta*m1;
    dm2 = E2*par.c - con_rate + par.P*(sum(m2) - par.N*m2) - par.delta*m2;
else
    dm1 = E1*par.c - con_rate + par.P*(m1_ext/par.V - m1) - par.delta*m1;
    dm2 = E2*par.c - con_rate + par.P*(m2_ext/par.V - m2) - par.delta*m2;
    RHS(1) = sum(par.P*(m1 - m1_ext/par.V)) - par.delta*m1_ext;
    RHS(2) = sum(par.P*(m2 - m2_ext/par.V)) - par.delta*m2_ext;
end

if par.feedback
    dE1 = -con_rate.*par.epsi.*E1;
    dE2 = -con_rate.*par.epsi.*E2;
else
    dE1 = -par.epsi.*E1;
    dE2 = -par.epsi.*E2;
end

RHS(sim_obj.m1_ind) = dm1;
RHS(sim_obj.m2_ind) = dm2;

RHS(sim_obj.E1_ind) = dE1;
RHS(sim_obj.E2_ind) = dE2;

end