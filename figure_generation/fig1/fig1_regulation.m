function fig1_live_dead(Fig1ax,linewidth,fontsize,seed,cutoff)

%% Good regulation cell

rng(seed);

axes(Fig1ax(5))
par.N = 1;
par.gamma_tot = 50;
par.P = 1;
par.c = 100;
par.max_t = 2500;
par.beta = 2;
par.delta = 1;

par.feedback = 1;
par.epsi = 1e-5;
max_growth = par.epsi*par.gamma_tot*par.c/(1+par.delta);

par.maxchange = 0.005;
par.V = 10;
par.overlay = 0;
par.n_replicate = 5;
par.n_store = 1e5;

%Set up auxotrophies
par.aux_type1 = 0;
par.aux_type2 = 0;

color1 = [206 37 123]/255;
color2 = [15 104 194]/255;

sim_obj = hybrid_simulation_master(par);

plot(sim_obj.record_t(cutoff:end)*max_growth,sim_obj.record_var(sim_obj.m1_ind,cutoff:end),'-',...
    'Color',color1,'LineWidth',linewidth)
hold on
plot(sim_obj.record_t(cutoff:end)*max_growth,sim_obj.record_var(sim_obj.m2_ind,cutoff:end),'-',...
    'Color',color2,'LineWidth',linewidth)
xlim([0,20])
xticks([0,10,20])
set(gca,'TickLabelInterpreter','latex')
ylim([0,2350]);
set(gca,'FontSize',fontsize);
ylabel('Metabolite level','Interpreter','latex')

xlabel('Time','Interpreter','latex')
box off
text(-0.25,1.25,'\textbf{C}','Interpreter','latex','Units','normalized','FontSize',4)
text(2.5,2*970,'$\beta = 2$','Interpreter','latex','FontSize',4);


%% Bad regulation cell

axes(Fig1ax(6))

par.beta = 20;

sim_obj = hybrid_simulation_master(par);

plot(sim_obj.record_t(cutoff:end)*max_growth,sim_obj.record_var(sim_obj.m1_ind,cutoff:end),'-',...
    'Color',color1,'LineWidth',linewidth)
hold on
plot(sim_obj.record_t(cutoff:end)*max_growth,sim_obj.record_var(sim_obj.m2_ind,cutoff:end),'-',...
    'Color',color2,'LineWidth',linewidth)
xlim([0,20])
xticks([0,10,20])
set(gca,'TickLabelInterpreter','latex')
ylim([0,2350])
set(gca,'FontSize',fontsize);
xlabel('Time','Interpreter','latex')
box off

text(-0.1,1.25,'\textbf{D}','Interpreter','latex','Units','normalized','FontSize',4)
text(2.5,2*970,'$\beta = 20$','Interpreter','latex','FontSize',4);


end
