function sim_obj = hybrid_simulation_master(par)

%This is the master script that runs the full metabolism model

if ~isfield(par,'feedback')
    par.feedback = 1;
end

%Get true gamma
par.gamma = par.gamma_tot/par.beta;

%Construct aux vectors
if (par.aux_type1 == 0) && (par.aux_type2 == 0)
    par.aux_var =[];
    par.aux_id = [];
else
    par.aux_var = [zeros(par.N - par.aux_type1 - par.aux_type2,1); ...
        ones(par.aux_type1 + par.aux_type2,1)]; %1 means it is auxotrophic
    par.aux_var = logical(par.aux_var);
    
    par.aux_id = [zeros(par.N - par.aux_type1 - par.aux_type2,1); ...
        ones(par.aux_type1,1); zeros(par.aux_type2,1)];
    par.aux_id = logical(par.aux_id); %1 means it is auxotrophic for nut1, 0 means nut2
end

tic

%Construction simulation object
sim_obj = construct_hybrid_simulation(par);

u = 0;

tn = 1;
burst_n = 1;
sim_obj.initial_vars = sim_obj.var;
sim_obj.initial_t = sim_obj.t;
while u == 0
        
    %Set-up random variables and integrals for next Gillespie step
    rand_int = exprnd(1);
    rand_var = rand(1);
    total_int = 0;
    
    %Run continuous timestepping until next gillespie event
    o = 0;
    while o == 0
        
        %Get total rate
        total_rate = sum(hybrid_rates(sim_obj,par));
        
        %COMPUTE K1
        var_i = sim_obj.var;
        dvar_1 = hybrid_RHS(var_i,sim_obj,par);
        
        %Get first RHS and determine dt
        relchange = abs(dvar_1./var_i);
        relchange = max(relchange(~isinf(relchange)));
        dt = par.maxchange/relchange; %Set time step
        
        %TRUE RHS
        RHS = dvar_1;
        
        %Compute running left-side rectangle integral
        total_int = total_int + dt*total_rate;
        
        if total_int > rand_int
            %Get the right dt to make the condition exact
            old_total = total_int - dt*total_rate;
            dt = (rand_int - old_total)/total_rate;
            total_int = old_total + dt*total_rate;
        end
        
        sim_obj.var = sim_obj.var + dt*RHS;
        sim_obj.t = sim_obj.t + dt;
        
        %Check if time has been exceeded
        if sim_obj.t > par.max_t
            o = 1;
        end
                
        %Determine if integral condition is met
        if total_int == rand_int
            
            rates = hybrid_rates(sim_obj,par);
            total_rate = sum(rates);
            probs = rates/total_rate;
            cumulative_probs = cumsum(probs);
            reac_ind = find(cumulative_probs > rand_var,1,'first');
            reac_ind_id = sim_obj.reac_id2(reac_ind);
            
            %Perform gillespie step
            if reac_ind_id == 1
                sim_obj.var = sim_obj.var + poissrnd(par.beta).*sim_obj.reac_matrix(reac_ind,:)';
                sim_obj.burst_t(burst_n) = sim_obj.t;
                burst_n = burst_n + 1;
            else
                sim_obj.var = sim_obj.var + sim_obj.reac_matrix(reac_ind,:)';
            end
            
            o = 1;
            
        end
        
        %Record state variables
        sim_obj.record_var(:,tn) = sim_obj.var;
        sim_obj.record_t(tn) = sim_obj.t;
        
        %Advance timestep
        tn = tn + 1;
        
    end
    
    %End simulation if time limit has been reached
    if sim_obj.t > par.max_t
        u = 1;
        disp('FINAL TIME REACHED. SIMULATION TERMINATING.')
    end
    
end
toc

%Trim storage vectors
sim_obj.record_var = sim_obj.record_var(:,1:(tn-1));
sim_obj.record_t = sim_obj.record_t(1:(tn-1));
sim_obj.burst_t = sim_obj.burst_t(1:(burst_n-1));

%Get the growth integrals from the last three quarters
sim_obj.growth_ints = zeros(par.N,1);

%Get indices
type1_ind = sim_obj.var_id3 == 1;
type2_ind = sim_obj.var_id3 == 2;
met_ind = sim_obj.var_id1 == 1;

%Find time that best approximates 1/4
max_ind = length(sim_obj.record_t);
quarter_t = max(sim_obj.record_t)/4;
[~,quarter_ind] = min(abs(sim_obj.record_t - quarter_t));

%Loop through and compute statistics for each cell
for i = 1:par.N
    
    %Extract metabolites
    cell_ind = sim_obj.var_id2 == i;
    m1_vec = sim_obj.record_var(met_ind & cell_ind & type1_ind,:);
    m2_vec = sim_obj.record_var(met_ind & cell_ind & type2_ind,:);    
    
    %Compute mean growth
    growth = min([m1_vec;m2_vec]);
    sim_obj.growth_ints(i) = trapz(sim_obj.record_t(quarter_ind:max_ind),growth(quarter_ind:max_ind));
    sim_obj.growth_ints(i) = sim_obj.growth_ints(i)/(sim_obj.record_t(max_ind) - sim_obj.record_t(quarter_ind));
    
    %Get enzyme data
    E1_data = sim_obj.record_var(sim_obj.E1_ind,:);
    E2_data = sim_obj.record_var(sim_obj.E2_ind,:);
    s1 = E1_data(quarter_ind:max_ind);
    s2 = E2_data(quarter_ind:max_ind);
    
    %Compute enzyme variances and covariances
    [sim_obj.E_covs(i),E1_mean,E2_mean,E1_vari,E2_vari] = compute_continuous_variance(s1,s2,sim_obj.record_t(quarter_ind:max_ind));
    sim_obj.E_means(i,:) = [E1_mean, E2_mean];
    sim_obj.E_varis(i,:) = [E1_vari, E2_vari];
    
    %Compute metabolite variances and covariances
    m1_trunc = m1_vec(quarter_ind:max_ind);
    m2_trunc = m2_vec(quarter_ind:max_ind);
    [sim_obj.m_covs(i),m1_mean,m2_mean,m1_vari,m2_vari] = compute_continuous_variance(m1_trunc,m2_trunc,sim_obj.record_t(quarter_ind:max_ind));
    sim_obj.m_means(i,:) = [m1_mean, m2_mean];
    sim_obj.m_varis(i,:) = [m1_vari, m2_vari];
    sim_obj.m_CV(i,:) = sqrt(sim_obj.m_varis(i,:))./sim_obj.m_means(i,:);
    
end

end