%This script generates the autocorrelation supp figure

%% Generate autocorrelation data

clear;clc

par.N = 1;
par.gamma_tot = 50;
par.P = 0;
par.c = 100;
par.max_t = 1000;
par.beta = 2;
par.delta = 1;
par.random_ICs = true;
par.feedback = 1;
par.epsi = 1e-5;
par.maxchange = 0.02;
par.V = 0;
par.overlay = 0;
par.n_replicate = 1;
par.n_store = 1e5;

%Set up auxotrophies
par.aux_type1 = 0;
par.aux_type2 = 0;

n_rep = 2500;
window_size = 1;

[corr_vec_low_beta,mean_time_vec_low_beta] = SC_autocorr(par,n_rep,window_size);

par.beta = 20;
[corr_vec_high_beta,mean_time_vec_high_beta] = SC_autocorr(par,n_rep,window_size);

save('supp_fig1_autocorrelation.mat','par','corr_vec_high_beta','corr_vec_low_beta',...
    'mean_time_vec_high_beta','mean_time_vec_low_beta','n_rep','window_size',...
    'max_growth')


%% Plot figure

clear;clc

load('supp_fig1_autocorrelation.mat')

colors = [0 0 0;...
    86 180 233]/255;

max_growth = 0.5*par.epsi*par.gamma_tot*par.c/(1+par.delta);

newfigure(2*0.5*3.42/2, 1.7*(1.3/3)*3.42/3*(3.5/2));

hold on
plot(mean_time_vec_low_beta.*max_growth,corr_vec_low_beta,'-','LineWidth',0.4,'color',colors(1,:))
plot(mean_time_vec_high_beta.*max_growth,corr_vec_high_beta,'-','LineWidth',0.4,'color',colors(2,:))

xlabel('Time lag','Interpreter','latex')
ylabel('Correlation to initial condition','Interpreter','latex')
set(gca,'TickLabelInterpreter','latex')
set(gca,'FontSize',6)
xlim([0,12])
xticks([0,6,12])

leg_top = 0.9;
leg_left = 6;
mylines = {'-','-'};
    labels = {'$\beta = 2$','$\beta = 20$'};

line_length = 1;
spacing_x = 0.2;
curr_y = leg_top;
for k = 1:length(labels)
    plot([leg_left,leg_left + line_length],[curr_y,curr_y],mylines{k},'Color',colors(k,:))
    text(leg_left + line_length + spacing_x,curr_y,labels{k},'Interpreter','latex','FontSize',6)
    curr_y = curr_y-0.14;
end


print(gcf, '-dpng','supp_fig1_autocorrelation.png','-r600');

