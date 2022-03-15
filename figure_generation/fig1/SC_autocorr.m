function [corr_vec,mean_time_vec] = SC_autocorr(par,n_rep,window_size)

%This script computes autocorrelation (from initial condition) in
%the full stochastic metabolite model

%Simulate replicates
obj_store = cell(n_rep,1);
for i = 1:n_rep
    obj_store{i} = hybrid_simulation_master(par);
end

%Define left endpoints
left_endpoints = 0:window_size:(par.max_t-window_size);

%Housekeeping
E_ind = obj_store{1}.E1_ind | obj_store{1}.E2_ind;
corr_vec = zeros(size(left_endpoints));
mean_time_vec = zeros(size(left_endpoints));
store_len = 1e7;

%Compute autocorrelation
for i = 1:length(left_endpoints)
    
    %Get endpoint
    left_e = left_endpoints(i);
    right_e = left_e + window_size;
    E_vec = zeros(store_len,1);
    init_E_vec = zeros(store_len,1);
    t_vec = zeros(store_len,1);
    c_ind = 1;
    
    %Loop through replicates
    for j = 1:n_rep
        
        temp_obj = obj_store{j};
        
        %Get initial E value
        E_init = temp_obj.initial_vars(E_ind);
        
        %Identify times within window
        time_ind = (temp_obj.record_t >= left_e) & (temp_obj.record_t < right_e);
        
        %Get E values in time window
        E_mat = temp_obj.record_var(E_ind,time_ind);
        
        vec_len = sum(time_ind);
        
        %Add in only the first enzyme, as it is randomized
        E_vec(c_ind:(c_ind+vec_len-1)) = E_mat(1,:)';
        init_E_vec(c_ind:(c_ind+vec_len-1)) = repmat(E_init(1),vec_len,1);
        t_vec(c_ind:(c_ind+vec_len-1)) = temp_obj.record_t(time_ind);
        c_ind = c_ind + vec_len;
    end
    
    %Trim zeros
    E_vec = E_vec(1:(c_ind-1));
    init_E_vec = init_E_vec(1:(c_ind-1));
    t_vec = t_vec(1:(c_ind-1));
    
    %Compute mean time and correlation
    mean_time_vec(i) = mean(t_vec);
    corr_vec(i) = corr(E_vec,init_E_vec);
    
end

mean_time_vec = [0,mean_time_vec];
corr_vec = [1,corr_vec];

end

