function fig2_par_sweep(Fig2ax,linewidth,fontsize)

%This script plots the mean growth and CVs for figure 2

load('large_multicell_sweep.mat')

max_growth = 25*100/(1+1);

colors = [0 0 0;...
    230 159 0;...
    86 180 233;...
    0 158 115; ...
    240 228 66;...
    213 94 0;...
    204 121 167]/255;

%% Plot growth function

par_cell = {'P'};
label_cell = {'Cell permeability, $P$'};
n3 = 4;
n1 = 40;

axes(Fig2ax(4))
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


text(-0.25,1.15,'\textbf{C}','Interpreter','latex','Units','normalized','FontSize',4)


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
text(0.5,0.9,'10 cells','Interpreter','latex','FontSize',4,'Units','normalized','HorizontalAlignment','center')

%% Plot noise

par_cell = {'P'};
label_cell = {'Cell permeability, $P$'};
n3 = 4;
n1 = 40;

axes(Fig2ax(3))
hold on
for j = 1:n3
    truei = j;
    curr_index = (truei-1)*n1 + 1:truei*n1;
    curr_mean = results_table.m_CV(curr_index);
    curr_par = results_table{curr_index,par_cell{1}};
    
    plot(curr_par,curr_mean,'k-','LineWidth',linewidth,'Color',colors(j,:));
    %plot(par_vec,m_CVs(:,j),'k-','LineWidth',linewidth,'Color',colors(j,:))
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
text(-0.25,1.15,'\textbf{B}','Interpreter','latex','Units','normalized','FontSize',4)

xticks([1e-1,1e1,1e3])
xticklabels({'$10^{-1}$','$10^1$','$10^3$'})
set(gca,'TickLabelInterpreter','latex')
set(gca,'XMinorTick','Off')
text(0.5,0.9,'10 cells','Interpreter','latex','FontSize',4,'Units','normalized','HorizontalAlignment','center')


leg_top = 1.7;
leg_left = 10;
mylines = {'-','-','-','-'};
labels = {'$\beta = 2$','$\beta = 10$','$\beta = 20$','$\beta = 100$'};

line_length = 30;
spacing_x = 0.6;
curr_y = leg_top;
for k = 1:length(labels)
    plot([leg_left,leg_left + line_length],[curr_y,curr_y],mylines{k},'Color',colors(k,:))
    text(leg_left + line_length + spacing_x,curr_y,labels{k},'Interpreter','latex','FontSize',4)
    curr_y = curr_y-0.3;
end

end

