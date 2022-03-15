%This script plots the mean growth as a function of volume

clear;clc

load('large_volume_sweep.mat')

max_growth = 25*100/(1+1);
max_perm_growth = 25*100/3;

colors = [0 0 0;...
    230 159 0;...
    86 180 233;...
    0 158 115; ...
    240 228 66;...
    213 94 0;...
    204 121 167]/255;

fontsize = 4;
linewidth = 0.4;

%% Plot

par_name = {'V'};
label = 'Extracellular volume ratio, $r_V$';
ind = 1;
n1 = 40;
n2 = 2;
n3 = 4;

i = 1;
newfigure(0.6*3.42/2, (1/3)*3.42/3*(3.5/2));
hold on
for j = 1:n3
    truei = ind+j-1;
    curr_index = (truei-1)*n1 + 1:truei*n1;
    curr_mean = results_table.mean_growth(curr_index);
    curr_par = results_table{curr_index,par_name};
    
    plot(curr_par,curr_mean./max_growth,'-','LineWidth',linewidth,'Color',colors(j,:));
end
plot([1e-1,1e3],[max_perm_growth,max_perm_growth]/max_growth,'k--','LineWidth',linewidth)

xlabel(label,'Interpreter','latex')
set(gca,'XScale','log')
set(gca,'FontSize',fontsize);
xlim([1e-1,1e3])
ylim([-0.03,1.1])
yticks([0,0.5,1])

    ylabel('Growth rate','Interpreter','latex')
    yticklabels({'0','0.5','1'})

xticks([1e-1,1e0,1e1,1e2,1e3,1e4])
xticklabels({'$10^{-1}$','$10^{0}$','$10^1$','$10^{2}$','$10^3$','$10^4$'})
set(gca,'TickLabelInterpreter','latex')
set(gca,'XMinorTick','Off')

leg_top = 1.1;
leg_left = 20;
mylines = {'-','-','-','-'};
    labels = {'$\beta = 2$','$\beta = 10$','$\beta = 20$','$\beta = 100$'};

line_length = 30;
spacing_x = 0.3;
curr_y = leg_top;
for k = 1:length(labels)
    plot([leg_left,leg_left + line_length],[curr_y,curr_y],mylines{k},'Color',colors(k+4*(i-1),:))
    text(leg_left + line_length + spacing_x,curr_y,labels{k},'Interpreter','latex','FontSize',4)
    curr_y = curr_y-0.11;
end
print(gcf, '-dpng','supp_fig1_volume.png','-r1500');
close all
