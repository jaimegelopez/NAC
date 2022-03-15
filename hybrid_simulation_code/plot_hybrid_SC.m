function plot_hybrid_SC(sim_obj,par)

%This script plots the output of the full metabolism model

figure
subplot(3,1,1);
hold on

if par.N == 1
    colors = [1 0 0];
else
    colors = distinguishable_colors(par.N);
end

linetypes = {'-','--'};

% Plot the enzymes and metabolites
enzyme_ind = sim_obj.var_id1 == 2;
type1_ind = sim_obj.var_id3 == 1;
type2_ind = sim_obj.var_id3 == 2;
met_ind = sim_obj.var_id1 == 1;
if ~par.overlay
    for i = 1:par.N
        cell_ind = sim_obj.var_id2 == i;
        plot(sim_obj.record_t,sim_obj.record_var(enzyme_ind & cell_ind & type1_ind,:),linetypes{1},'Color',colors(i,:));
        plot(sim_obj.record_t,sim_obj.record_var(enzyme_ind & cell_ind & type2_ind,:),linetypes{2},'Color',colors(i,:));
    end
    ylim([0,max(ylim)])
    title('Enzyme concentration')
    xlim([0,max(sim_obj.record_t)])
    
    subplot(3,1,2);
    hold on
    for i = 1:par.N
        cell_ind = sim_obj.var_id2 == i;
        plot(sim_obj.record_t,sim_obj.record_var(met_ind & cell_ind & type1_ind,:),linetypes{1},'Color',colors(i,:));
        plot(sim_obj.record_t,sim_obj.record_var(met_ind & cell_ind & type2_ind,:),linetypes{2},'Color',colors(i,:));
    end
    ylim([0,max(ylim)])
    title('Metabolite concentration')
    xlim([0,max(sim_obj.record_t)])
    
    
else
    
    for i = 1:par.N
        cell_ind = sim_obj.var_id2 == i;
        met = sim_obj.record_var(met_ind & cell_ind & type1_ind,:);
        met = met/max(met);
        plot(sim_obj.record_t,met,linetypes{1},'Color',colors(i,:));
       
        enzyme = sim_obj.record_var(enzyme_ind & cell_ind & type1_ind,:);
        enzyme = enzyme/max(enzyme);
        plot(sim_obj.record_t,enzyme,linetypes{2},'Color',colors(i,:));
    end
    ylim([0,max(ylim)])
    title('Enzyme and metabolite 1 concentration')
    xlim([0,par.max_t])
    legend('Metabolite','Enzyme')
    
    subplot(3,1,2);
    hold on
    for i = 1:par.N
        cell_ind = sim_obj.var_id2 == i;
        met = sim_obj.record_var(met_ind & cell_ind & type2_ind,:);
        met = met/max(met);
        plot(sim_obj.record_t,met,linetypes{1},'Color',colors(i,:));
       
        enzyme = sim_obj.record_var(enzyme_ind & cell_ind & type2_ind,:);
        enzyme = enzyme/max(enzyme);
        plot(sim_obj.record_t,enzyme,linetypes{2},'Color',colors(i,:));
    end
    ylim([0,max(ylim)])
    title('Enzyme and metabolite 2 concentration')
    xlim([0,max(sim_obj.record_t)])
end

subplot(3,1,3);
hold on
ext_ind = sim_obj.var_id2 == 0;
plot(sim_obj.record_t,sim_obj.record_var(ext_ind & type1_ind,:)/par.V,linetypes{1},'Color','k');
plot(sim_obj.record_t,sim_obj.record_var(ext_ind & type2_ind,:)/par.V,linetypes{2},'Color','k');
ylim([0,max(ylim)])
xlabel('Time')
title('External metabolite concentration')
xlim([0,max(sim_obj.record_t)])

end

