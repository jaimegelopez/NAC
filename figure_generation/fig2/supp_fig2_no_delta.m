%This script generates the no-delta supplementary figure

clear;clc

%% Set up for delta = 0 plotting

max_growth = 25*100;
colors = [0 0 0;...
    230 159 0;...
    86 180 233;...
    0 158 115; ...
    240 228 66;...
    213 94 0;...
    204 121 167]/255;

newfigure(3.42/2, 4*(1/3)*3.42/3*(3.5/2));

gap = [0.15,0.17]*3/4;
marg_h = [0.09,0.07]*3/4;
marg_w = [0.16,0.06]*3/4;

Figax = tight_subplot(4,2,gap,marg_h,marg_w);
fontsize = 4;
linewidth = 0.4;

%% Plot single cell results with variable beta

load('large_single_cell_sweep_no_delta.mat')

inds = [1,5];
n1 = 40;
n3 = [4,3];

axes(Figax(2))
hold on
for j = 1:n3(2)
    truei = inds(2)+j-1;
    curr_index = (truei-1)*n1 + 1:truei*n1;
    curr_mean = results_table.mean_growth(curr_index);
    curr_par = results_table{curr_index,'beta'};
    
    plot(curr_par,curr_mean./max_growth,'k-','LineWidth',linewidth,'Color',colors(j+4*(2-1),:));
end

xlabel('Burst size, $\beta$','Interpreter','latex')
set(gca,'XScale','log')
set(gca,'FontSize',fontsize);
xlim([1e-1,max(curr_par)])
ylim([-0.03,1])
yticks([0,0.5,1])

ylabel('Growth rate','Interpreter','latex')
yticklabels({'0','0.5','1'})

text(-0.25,1.15,'\textbf{B}','Interpreter','latex','Units','normalized','FontSize',4)


xticks([1e-1,1e1,1e3])
xticklabels({'$10^{-1}$','$10^1$','$10^3$'})
set(gca,'TickLabelInterpreter','latex')
set(gca,'XMinorTick','Off')

leg_top = 0.55;
leg_left = 30;
mylines = {'-','-','-','-'};

labels = {'$P = 0$','$P = 1$','$P = 10$'};
line_length = 30;
spacing_x = 0.3;
curr_y = leg_top;
for k = 1:length(labels)
    plot([leg_left,leg_left + line_length],[curr_y,curr_y],mylines{k},'Color',colors(k+4,:))
    text(leg_left + line_length + spacing_x,curr_y,labels{k},'Interpreter','latex','FontSize',4)
    curr_y = curr_y-0.12;
end

text(0.5,1.1,'Single cell','Interpreter','latex','FontSize',4,'Units','normalized','HorizontalAlignment','center')



%% Plot single cell results with variable P

load('no_delta_single_cell_P_only.mat')

inds = 1;
n1 = 80;
n3 = 4;

axes(Figax(1))
hold on
for j = 1:n3
    truei = j;
    curr_index = (j-1)*n1 + 1:j*n1;
    curr_mean = results_table.mean_growth(curr_index);
    curr_par = results_table{curr_index,'P'};
    
    plot(curr_par,curr_mean./max_growth,'k-','LineWidth',linewidth,'Color',colors(j+3,:));
end

xlabel('Cell permeability, $P$','Interpreter','latex')
set(gca,'XScale','log')
set(gca,'FontSize',fontsize);
xlim([1e-6,max(curr_par)])
ylim([-0.03,1])
yticks([0,0.5,1])

ylabel('Growth rate','Interpreter','latex')
yticklabels({'0','0.5','1'})

text(-0.25,1.15,'\textbf{A}','Interpreter','latex','Units','normalized','FontSize',4)


xticks([1e-6,1e-3,1e0,1e3])
xticklabels({'$10^{-6}$','$10^{-3}$','$10^0$','$10^3$'})
set(gca,'TickLabelInterpreter','latex')
set(gca,'XMinorTick','Off')

text(0.5,1.1,'Single cell','Interpreter','latex','FontSize',4,'Units','normalized','HorizontalAlignment','center')


%% Plot metabolites and enzymes during crash behavior
color1 = [206 37 123]/255;
color2 = [15 104 194]/255;
rng(666)
par.N = 1;
par.gamma_tot = 50;
par.P = 0;
par.c = 100;
par.max_t = 4000;
par.beta = 20;
par.delta = 0;
par.feedback = 1;
par.epsi = 1e-5;

par.maxchange = 0.01;
par.V = 10;
par.overlay = 0;
par.n_replicate = 5;
par.n_store = 1e6;
par.aux_type1 = 0;
par.aux_type2 = 0;
sim_obj = hybrid_simulation_master(par);

%Metabolites
cutoff = 1;
true_max_growth = par.epsi*par.gamma_tot*par.c/(2*(1+par.delta));
axes(Figax(4))
axis off
axes(Figax(3))
mult = 2.335;
pos = get(gca,'Position');
pos(3) = mult*pos(3);
set(gca,'Position',pos);
plot(sim_obj.record_t(cutoff:end)*true_max_growth,sim_obj.record_var(sim_obj.m1_ind,cutoff:end),'-',...
    'Color',color1,'LineWidth',linewidth)
hold on
plot(sim_obj.record_t(cutoff:end)*true_max_growth,sim_obj.record_var(sim_obj.m2_ind,cutoff:end),'-',...
    'Color',color2,'LineWidth',linewidth)
set(gca,'YScale','log')
xlim([0,75])
xticks([0,15,30,45,60,75])
set(gca,'TickLabelInterpreter','latex')
ylim([1e-4,1e8])
yticks([1e-4,1e0,1e4,1e8])
set(gca,'FontSize',fontsize);
box off
ylabel('Metabolite level','Interpreter','latex')
xlabel('Time','Interpreter','latex')
text(-0.25/mult,1.25,'\textbf{C}','Interpreter','latex','Units','normalized','FontSize',4)

%Enzymes
axes(Figax(6))
axis off
axes(Figax(5))
pos = get(gca,'Position');
pos(3) = mult*pos(3);
set(gca,'Position',pos);
plot(sim_obj.record_t(cutoff:end)*true_max_growth,sim_obj.record_var(sim_obj.E1_ind,cutoff:end),'-',...
    'Color',color1,'LineWidth',linewidth)
hold on
plot(sim_obj.record_t(cutoff:end)*true_max_growth,sim_obj.record_var(sim_obj.E2_ind,cutoff:end),'-',...
    'Color',color2,'LineWidth',linewidth)
set(gca,'YScale','log')
xlim([0,75])
xticks([0,15,30,45,60,75])
set(gca,'TickLabelInterpreter','latex')
ylim([1e-3,1e3])
yticks([1e-3,1e0,1e3])
set(gca,'FontSize',fontsize);
box off
ylabel('Enzyme level','Interpreter','latex')
xlabel('Time','Interpreter','latex')
text(-0.25/mult,1.25,'\textbf{D}','Interpreter','latex','Units','normalized','FontSize',4)


%% Plot multicell growth function
load('large_multicell_sweep_no_delta.mat')

par_cell = {'P'};
label_cell = {'Cell permeability, $P$'};
n3 = 4;
n1 = 40;

axes(Figax(7))
hold on
for j = 1:n3
    truei = j;
    curr_index = (truei-1)*n1 + 1:truei*n1;
    curr_mean = results_table.mean_growth(curr_index);
    curr_par = results_table{curr_index,par_cell{1}};
    
    plot(curr_par,curr_mean./max_growth,'k-','LineWidth',linewidth,'Color',colors(j,:));
end

xlabel(label_cell{1},'Interpreter','latex')
set(gca,'XScale','log')
set(gca,'FontSize',fontsize);
xlim([1e-1,max(curr_par)])
yticks([0,0.5,1])
yticklabels({'0','0.5','1'})
ylabel('Mean growth','Interpreter','latex')

text(-0.25,1.15,'\textbf{E}','Interpreter','latex','Units','normalized','FontSize',4)

xticks([1e-1,1e1,1e3])
xticklabels({'$10^{-1}$','$10^1$','$10^3$'})
set(gca,'TickLabelInterpreter','latex')
set(gca,'XMinorTick','Off')

ylim([-0.03,1])
yticks([0,0.5,1])
yticklabels({'0','0.5','1'})
ylabel({'Growth rate'},'Interpreter','latex')
yhandle = get(gca,'YLabel');
pos = yhandle.Position;
pos(1) = 1.5*pos(1);
set(yhandle,'Position',pos);

leg_top = 0.7;
leg_left = 10;
mylines = {'-','-','-','-'};
labels = {'$\beta = 2$','$\beta = 10$','$\beta = 20$','$\beta = 100$'};

line_length = 30;
spacing_x = 0.6;
curr_y = leg_top;
for k = 1:length(labels)
    plot([leg_left,leg_left + line_length],[curr_y,curr_y],mylines{k},'Color',colors(k,:))
    text(leg_left + line_length + spacing_x,curr_y,labels{k},'Interpreter','latex','FontSize',4)
    curr_y = curr_y-0.14;
end

text(0.5,1.1,'10 cells','Interpreter','latex','FontSize',4,'Units','normalized','HorizontalAlignment','center')


%% Plot multicell noise

par_cell = {'P'};
label_cell = {'Cell permeability, $P$'};
n3 = 3;
n1 = 40;

axes(Figax(8))
hold on
for j = 1:n3
    truei = j;
    curr_index = (truei-1)*n1 + 1:truei*n1;
    curr_mean = results_table.m_CV(curr_index);
    curr_par = results_table{curr_index,par_cell{1}};
    
    plot(curr_par,curr_mean,'k-','LineWidth',linewidth,'Color',colors(j,:));
end

xlabel(label_cell{1},'Interpreter','latex')
set(gca,'XScale','log')
set(gca,'FontSize',fontsize);
xlim([1e-1,max(curr_par)])
ylim([0,2.8])
yticks([0,1,2])
yticklabels({'0','1','2'})
ylabel('Metabolite CV','Interpreter','latex')
yhandle = get(gca,'YLabel');
pos = yhandle.Position;
pos(1) = 1.5*pos(1);
set(yhandle,'Position',pos);
text(-0.25,1.15,'\textbf{F}','Interpreter','latex','Units','normalized','FontSize',4)

xticks([1e-1,1e1,1e3])
xticklabels({'$10^{-1}$','$10^1$','$10^3$'})
set(gca,'TickLabelInterpreter','latex')
set(gca,'XMinorTick','Off')
text(0.5,1.1,'10 cells','Interpreter','latex','FontSize',4,'Units','normalized','HorizontalAlignment','center')


print(gcf, '-dpng','supp_fig2_no_delta.png','-r1500');

close all