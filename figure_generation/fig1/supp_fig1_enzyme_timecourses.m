clear;clc
newfigure(3.42/2, (0.9/3)*3.42/3*(3.5/2));

gap = [0.12,0.11];
marg_h = [0.25,0.1];
marg_w = [0.1,0.05];
Fig1ax = tight_subplot(1,2,gap,marg_h,marg_w);
seed=666;
cutoff = 1;
fontsize = 4;
linewidth = 0.4;

%% Good regulation cell

rng(seed);

axes(Fig1ax(1))
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

plot(sim_obj.record_t*max_growth,sim_obj.record_var(sim_obj.E1_ind,:),'-',...
    'Color',color1,'LineWidth',linewidth)
hold on
plot(sim_obj.record_t*max_growth,sim_obj.record_var(sim_obj.E2_ind,:),'-',...
    'Color',color2,'LineWidth',linewidth)
xlim([0,20])
xticks([0,10,20])
set(gca,'TickLabelInterpreter','latex')
ylim([0,60])
yticks([0,20,40,60])
set(gca,'FontSize',fontsize);
ylabel('Enzyme level','Interpreter','latex')

xlabel('Time','Interpreter','latex')
box off
text(-0.2,1.07,'\textbf{A}','Interpreter','latex','Units','normalized','FontSize',4)
text(2.5,50,'$\beta = 2$','Interpreter','latex','FontSize',4);

good_reg_r = corr(sim_obj.record_var(sim_obj.E1_ind,:)',sim_obj.record_var(sim_obj.E2_ind,:)','type','Pearson');

%% Bad regulation cell

axes(Fig1ax(2))

par.beta = 20;

sim_obj = hybrid_simulation_master(par);

plot(sim_obj.record_t*max_growth,sim_obj.record_var(sim_obj.E1_ind,:),'-',...
    'Color',color1,'LineWidth',linewidth)
hold on
plot(sim_obj.record_t*max_growth,sim_obj.record_var(sim_obj.E2_ind,:),'-',...
    'Color',color2,'LineWidth',linewidth)
xlim([0,20])
xticks([0,10,20])
set(gca,'TickLabelInterpreter','latex')
ylim([0,60])
yticks([0,20,40,60])
yticklabels([])
set(gca,'FontSize',fontsize);
xlabel('Time','Interpreter','latex')
box off

text(-0.1,1.07,'\textbf{B}','Interpreter','latex','Units','normalized','FontSize',4)
text(2.5,50,'$\beta = 20$','Interpreter','latex','FontSize',4);

bad_reg_r = corr(sim_obj.record_var(sim_obj.E1_ind,:)',sim_obj.record_var(sim_obj.E2_ind,:)','type','Pearson');

print(gcf, '-dpng','supp_fig1_enzyme_timecourses.png','-r1500');
close all
