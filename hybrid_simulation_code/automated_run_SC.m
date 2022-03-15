function automated_run_SC(table_file,id,results_folder)

%This script runs the full metabolism model in an automated fashion
%made to work with Slurm array jobs

%Set seed with time to ensure randomization
rng('shuffle')

%Read in default parameter structure
load('parameters/default_parameter_struct.mat');

%Load parameter table
par_table = readtable(table_file,'ReadRowNames',false,...
    'ReadVariableNames',true);

%Replace parameters
non_id_vars = par_table.Properties.VariableNames(1:(end-1));
matching_row = par_table.id == id;
for i = 1:length(non_id_vars)
    par.(non_id_vars{i}) = par_table{matching_row,non_id_vars{i}};
end

%Generate results file
results_file = [results_folder,'/','sim_id_',num2str(id),'.mat'];

%Run the simulation
growth_cell = cell(par.n_replicate,1);
E_mean_cell = cell(par.n_replicate,1);
E_vari_cell = cell(par.n_replicate,1);
E_cov_cell = cell(par.n_replicate,1);

m_mean_cell = cell(par.n_replicate,1);
m_vari_cell = cell(par.n_replicate,1);
m_cov_cell = cell(par.n_replicate,1);

%Loop through all replicates and run simulation
disp(['Running ',num2str(par.n_replicate), ' simulations for ',num2str(par.max_t),' units each.'])
for i = 1:par.n_replicate
    
    clear sim_obj
    sim_obj = hybrid_simulation_master(par);
    
    disp(['Simulation required ',num2str(length(sim_obj.record_t)),' time steps.'])
    
    growth_cell{i} = sim_obj.growth_ints;
    E_mean_cell{i} = sim_obj.E_means;
    E_vari_cell{i} = sim_obj.E_varis;
    E_cov_cell{i} = sim_obj.E_covs;

    m_mean_cell{i} = sim_obj.m_means;
    m_vari_cell{i} = sim_obj.m_varis;
    m_cov_cell{i} = sim_obj.m_covs;
    
end

save(results_file,'par','growth_cell','E_mean_cell','E_vari_cell',...
'E_cov_cell','m_mean_cell','m_vari_cell','m_cov_cell');

disp(['Raw results saved to ',results_file,'.'])


end
