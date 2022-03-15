%This script plots the metabolite CVs corresponding to figure 1

%% Load data

clear;clc

newfigure(0.5*3.42/2, (1.3/3)*3.42/3*(3.5/2));

load('large_single_cell_sweep.mat')

max_growth = 25*250/(1+0.5);
fontsize = 4;
linewidth = 0.4;
colors = [0 0 0;...
    230 159 0;...
    86 180 233;...
    0 158 115; ...
    240 228 66;...
    213 94 0;...
    204 121 167]/255;

par_cell = {'P'};
label_cell = {'Cell permeability, $P$'};
n3 = 3;
n1 = 40;

%% Plot CVs

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
ylim([0,1.2])
yticks([0,1,2])
yticklabels({'0','1','2'})
ylabel('Metabolite CV','Interpreter','latex')
yhandle = get(gca,'YLabel');
pos = yhandle.Position;
pos(1) = 2*pos(1);
set(yhandle,'Position',pos);
text(-0.25,1.15,'\textbf{B}','Interpreter','latex','Units','normalized','FontSize',4)

xticks([1e-1,1e1,1e3])
xticklabels({'$10^{-1}$','$10^1$','$10^3$'})
set(gca,'TickLabelInterpreter','latex')
set(gca,'XMinorTick','Off')


leg_top = 1.1;
leg_left = 10;
mylines = {'-','-','-','-'};
labels = {'$\beta = 2$','$\beta = 10$','$\beta = 20$'};

line_length = 30;
spacing_x = 0.6;
curr_y = leg_top;
for k = 1:length(labels)
    plot([leg_left,leg_left + line_length],[curr_y,curr_y],mylines{k},'Color',colors(k,:))
    text(leg_left + line_length + spacing_x,curr_y,labels{k},'Interpreter','latex','FontSize',4)
    curr_y = curr_y-0.12;
end

print(gcf, '-dpng','supp_fig1_metabolite_CV.png','-r1500');
close all
