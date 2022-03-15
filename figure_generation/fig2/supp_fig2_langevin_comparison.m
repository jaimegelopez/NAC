%This script compares the analytical Langevin prediction to the numerical
%results

clear;clc

%% Generate analytical prediction
load('langevin_comparison.mat');

load('default_parameter_struct.mat');

gmax = 0.5*par.epsi*par.gamma_tot*par.c/(1+par.delta);
Gamma = par.gamma_tot*gmax;
mu_E = gmax;
mu_m = par.delta+1;
kappa = par.c;

num = mu_E.*mu_m.*(kappa.*results_table.beta + 2.*mu_E + 2.*mu_m + kappa);
den = 2.*kappa.*Gamma.*results_table.N.*(mu_E + mu_m);
results_table.analytic_m_CV = sqrt((num./den));

%% Population size vs. noise

colors = [0 0 0;...
    230 159 0;...
    86 180 233;...
    0 158 115]./255;
newfigure(2*0.6*3.42/2, (1/3)*3.42/3*(3.5/2));
gap = [0.12,0.11];
marg_h = [0.23,0.07];
Figax = tight_subplot(1,2,gap,marg_h,[0.1,0.05]);
axes(Figax(1));
hold on
n = 20;
fontsize = 4;
linewidth = 0.6;
num_set = 4; 
for i = 1:(num_set-1)
    index = ((i-1)*n + 1):(i*n);
    N_vec = results_table.N(index);
    CV_vec = results_table.m_CV(index);
    analytic_vec = results_table.analytic_m_CV(index);
    include = CV_vec ~= 0;
    plot(1./sqrt(N_vec(include)),CV_vec(include),'-','Color',colors(i,:),'LineWidth',linewidth)
    plot(1./sqrt(N_vec(include)),CV_vec(include),'o',...
        'MarkerEdgeColor',colors(i,:),'MarkerFaceColor',colors(i,:),'MarkerSize',0.76)
    plot(1./sqrt(N_vec(include)),analytic_vec(include),'--','Color',colors(i,:),'LineWidth',linewidth)

end

xlabel('$1/\sqrt{\textrm{N}}$','Interpreter','latex','FontSize',fontsize)
ylabel('$\textrm{CV}_\textrm{m}$','Interpreter','latex','FontSize',fontsize)

xticks([0.2,0.6,1])
xticklabels({'0.2','0.6','1'});
xlim([0.2,1.02])
yticks([0,0.4,0.8])
yticklabels({'0','0.4','0.8'});
ylim([0,0.8])
set(gca,'TickLabelInterpreter','latex')
set(gca,'FontSize',fontsize);

leg_top = 0.73;
leg_left = 0.3;
mylines = {'-','-','-','-'};
labels = {'$\beta = 2$','$\beta = 10$','$\beta = 20$'};

line_length = 0.05;
spacing_x = 0.03;
curr_y = leg_top;
for k = 1:length(labels)
    plot([leg_left,leg_left + line_length],[curr_y,curr_y],mylines{k},'Color',colors(k,:))
    text(leg_left + line_length + spacing_x,curr_y,labels{k},'Interpreter','latex','FontSize',4)
    curr_y = curr_y-0.1;
end

text(-0.22,1.04,'\textbf{A}','Interpreter','latex','Units','normalized','FontSize',4)


%% Burst size vs. noise

axes(Figax(2))
hold on
beta_vec = results_table.beta(81:100);
CV_vec = results_table.m_CV(81:100);
analytic_vec = results_table.analytic_m_CV(81:100);
plot(sqrt(beta_vec),CV_vec,'r-','LineWidth',linewidth)
plot(sqrt(beta_vec),analytic_vec,'r--','LineWidth',linewidth)
xlabel('$\sqrt{\beta}$','Interpreter','latex','FontSize',fontsize)
ylabel('$\textrm{CV}_\textrm{m}$','Interpreter','latex','FontSize',fontsize)

xticks([0,5,10])
xticklabels({'0','5','10'});
xlim([0,10])
yticks([0,0.2,0.4])
yticklabels({'0','0.2','0.4'});
ylim([0,0.4])
set(gca,'TickLabelInterpreter','latex')
set(gca,'FontSize',fontsize);
text(-0.22,1.04,'\textbf{B}','Interpreter','latex','Units','normalized','FontSize',4)

print(gcf, '-dpng','supp_fig2_langevin_comparison.png','-r1500');

close all



