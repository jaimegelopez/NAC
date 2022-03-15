function results_table = aggregate_SC_results(sim_name,parts,m_data)

%This script aggregrates results from automated runs of the full metabolite
%model

if ~exist('parts')
    parts = [];
end

if ~exist('m_data')
    m_data = 0;
end

%Generate file and folder names
results_dir = ['/home/jglopez/research/stochastic_crossfeeding/automated_runs/results/',sim_name,'/'];
par_file = ['/home/jglopez/research/stochastic_crossfeeding/automated_runs/parameters/',sim_name,'.csv'];

%Read results table
results_table = readtable(par_file,'ReadRowNames',false,...
    'ReadVariableNames',true);

%Make new columns
results_table.mean_growth = nan(size(results_table.id));
results_table.std_growth = nan(size(results_table.id));
results_table.m_mean = nan(size(results_table.id));
results_table.m_vari = nan(size(results_table.id));
results_table.m_cov = nan(size(results_table.id));
results_table.m_CV = nan(size(results_table.id));

%Loop through and get data

if isempty(parts)
    for i = 1:length(results_table.id)
        temp_id = results_table.id(i);
        filename = [results_dir,'sim_id_',num2str(temp_id),'.mat'];
        
        if isfile(filename)
            temp = load(filename);
            agg_growth = cellfun(@(x) mean(mean(x)),temp.growth_cell);
            
            results_table.mean_growth(i) = mean(agg_growth);
            results_table.std_growth(i) = std(agg_growth);
            
            if m_data
                results_table.m_mean(i) = mean(cellfun(@(x) mean(x(:)),temp.m_mean_cell));
                results_table.m_vari(i) = mean(cellfun(@(x) mean(x(:)),temp.m_vari_cell));
                results_table.m_cov(i) = mean(cellfun(@(x) mean(x(:)),temp.m_cov_cell));
                results_table.m_CV(i) = sqrt(results_table.m_vari(i))/results_table.m_mean(i);
            end
        end
        
    end
    
else
    
    for i = 1:length(results_table.id)
        temp_id = results_table.id(i);
        filename = [results_dir,'sim_id_',num2str(temp_id),'.mat'];
        temp = load(filename);
        agg_growth1 = cellfun(@(x) mean(mean(x(1:parts(1)))),temp.growth_cell);
        agg_growth2 = cellfun(@(x) mean(mean(x(parts(1)+1:parts(2)))),temp.growth_cell);
        
        results_table.mean_growth1(i) = mean(agg_growth1);
        results_table.std_growth1(i) = std(agg_growth1);
        
        results_table.mean_growth2(i) = mean(agg_growth2);
        results_table.std_growth2(i) = std(agg_growth2);
    end
end




end

